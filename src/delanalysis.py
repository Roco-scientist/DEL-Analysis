import os
import math
import numpy as np
import plotly.graph_objects as go
import re

from datetime import date
from io import StringIO
from pandas import read_csv, DataFrame, concat, merge, Series, notna

from plotly.subplots import make_subplots
from scipy.stats import zscore
from typing import Optional, Type, List


class DelData:
    """
    An object for working with data output from NGS-Barcode-Count
    """

    def __init__(self, data: DataFrame, data_type="Count", sample_name=None, unique_synthons: List[int] = None):
        self.data = data  # DataFrame of the data
        self.data_type = data_type  # The type of data which determines what actions can happen
        self.sample_name = sample_name  # The sample name which is only used with DelDataSample
        self.barcode_info = {}  # holds a dictionary of indexes and library sizes for each barcode pairing found
        # Unique synthons per barcode holds a dicionary of how many synthons per building block.  This is used to calculate library size for each condition
        self.unique_synthons_per_barcode = {}
        # If unique_synthons is supplied use those, otherwise infer by unique synthons whithin
        # barcode columns
        if unique_synthons is None:
            self._infer_number_unique_synthons()
        else:
            self.update_synthon_numbers(unique_synthons)

    def __repr__(self):
        """
        Return a string representation of the data.  Very similar to how pandas.DataFrame outputs
        the data
        """
        buf = StringIO("")
        self.data.to_string(
            buf=buf,
            max_rows=10,
            min_rows=10,
            max_cols=10,
            line_width=None,
            max_colwidth=None,
            show_dimensions=True,
        )

        return buf.getvalue()

    def __str__(self):
        """
        Return a string print of the data.  Very similar to how pandas.DataFrame outputs
        the data
        """
        buf = StringIO("")
        self.data.to_string(
            buf=buf,
            max_rows=10,
            min_rows=10,
            max_cols=10,
            line_width=None,
            max_colwidth=None,
            show_dimensions=True,
        )

        return buf.getvalue()

    def _zscore(self, counts: np.ndarray) -> np.ndarray:
        """
        Standard z-scoring of the counts data
        """
        if self.data_type == "zscore":
            raise Exception("Data is already zscored")
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        return zscore(counts)

    def _binomial_zscore(self, counts: np.ndarray, del_library_size: int) -> np.ndarray:
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as (observed_count - expected_count)/sqrt(total_counts * expected_probability * (1-expected_probability))
        Where the denominator is the binomial standard deviation
        Which equals sqrt(total_counts) * (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        """
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        expected_probability = 1/del_library_size  # get the probability as 1/library size
        total_counts = np.sum(counts, axis=0)  # Total counts is the sum of each column
        observed_probability = counts / total_counts
        binomial_sd = math.sqrt(expected_probability * (1 - expected_probability))
        return np.sqrt(total_counts) * (observed_probability - expected_probability) / binomial_sd

    def _binomial_zscore_sample_normalized(self, counts: np.ndarray, del_library_size: int) -> np.ndarray:
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as:
        (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        Where the denominator is the binomial standard deviation
        """
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        expected_probability = 1/del_library_size  # get the probability as 1/library size
        total_counts = np.sum(counts, axis=0)
        observed_probability = counts / total_counts
        binomial_sd = math.sqrt(expected_probability * (1 - expected_probability))
        return (observed_probability - expected_probability) / binomial_sd

    def _enrichment(self, counts: np.ndarray, del_library_size: int) -> np.ndarray:
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        total_counts = np.sum(counts, axis=0)
        return counts * del_library_size / total_counts

    def _infer_number_unique_synthons(self):
        """
        Infers the number of unique synthons by finding the amount of unique barcodes per barcode column
        """
        print("Total synthons inferred")
        for column_name in self.counted_barcode_columns():
            unique_synthons = [synthon for synthon in set(self.data[column_name]) if notna(synthon)]
            print(f"{column_name}: {len(unique_synthons)}")
            self.unique_synthons_per_barcode[column_name] = len(unique_synthons)
        print("If this is not correct, update with known number with update_synthon_numbers()")

    def update_synthon_numbers(self, unique_synthons_per_barcode: List[int]) -> None:
        """
        Allows user to update unique synthon numbers to fix what is inferred
        """
        barcode_columns = self.counted_barcode_columns()
        # Enforce the input as a list
        if not isinstance(unique_synthons_per_barcode, list):
            raise Exception(
                f"List of number of unique sythons required with the same length as counted barcodes, {len(barcode_columns)}")
        # Make sure there are enough for as many barcode columns exist
        if len(unique_synthons_per_barcode) != len(barcode_columns):
            raise Exception(
                f"Length of list, {len(unique_synthons_per_barcode)}, shorter than the number of counted barcodes {','.join(barcode_columns)}")
        # For each inputed unique synthon, replace the dictionary value already held
        for barcode_num, total_barcodes in zip(barcode_columns, unique_synthons_per_barcode):
            self.unique_synthons_per_barcode[barcode_num] = total_barcodes
        # Get the barcode groups and update the DEL library sizes based on the new wynthon numbers
        self._infer_barcode_groups()

    def _infer_barcode_groups(self):
        """
        Finds the grouping of barcode enrichments.  This could be all tri-synthons, bi-synthons
        within a tri-synthon design or mono-synthons within tri-synthon design etc.  Each group is
        normalized separately, therefor the indexes of where each group occurs is inferred here and
        stored within self.barcode_info.  The potential DEL library size for each group is also
        inferred by how many unique synthons exists with each building block that is used within the
        group and also stored in self.barcode_info.
        """
        barcode_columns = np.array(self.counted_barcode_columns())
        # Create an array of true/false of where a barcode is included. For di-synthon counts, only
        # two barcode are filled.  This true false is used to create a key for each row
        barcode_included = notna(self.data[barcode_columns]).values
        # For each row, create a key of comma separated barcode column names.  So for di-synthon
        # which is Barcode_1 and Barcode_3 count, this would create "Barcode_1,Barcode_3"
        barcode_keys = np.apply_along_axis(lambda row: ",".join(
            barcode_columns[row]), 1, barcode_included)
        # For each unique group counting, ie "Barcode_1,Barcode_3", get the indexes and the inferred
        # library size.  The index is used later for scoring the counts by the groups
        for barcode_key in set(barcode_keys):
            self.barcode_info[barcode_key] = {"indexes": np.where(barcode_keys == barcode_key)[0]}
            library_size = 1
            for barcode_num in barcode_key.split(","):
                library_size = library_size * self.unique_synthons_per_barcode[barcode_num]
            self.barcode_info[barcode_key]["library_size"] = library_size

    def data_columns(self):
        """
        Returns all column names that contain data that is not the barcodes
        """
        # Return any column that does not start with Barcode
        return [col for col in self.data.columns if not re.match(r"Barcode(_\d+){0,1}", col)]

    def counted_barcode_columns(self):
        """
        Returns all building block barcode column names
        """
        # Return any column that starts with Barcode and does or does not have a digit
        return [col for col in self.data if re.match(r"Barcode(_\d+){0,1}", col)]

    def number_synthons(self) -> List[int]:
        """
        Returns the bumber of barcodes for each row of the data.  This is used for the comparison
        graph
        """
        # First select only barcode data, then true false whether or not the barcode was used for
        # counts, the sum the rows
        barcode_data = self.data.loc[:, self.counted_barcode_columns()]
        return notna(barcode_data).sum(axis=1).values

    def to_csv(self, out_file: str):
        """
        Writes the current data to a csv file
        """
        self.data.to_csv(out_file, index=False)

    def data_descriptor(self):
        """
        Returns the data type without spaces for output files
        """
        return self.data_type.replace(" ", "_")


class DelDataMerged(DelData):
    """
    An object for working with merged data output from NGS-Barcode-Count
    """

    def quantile_normalize(self, inplace=False):
        """
        Qauntile normalizes the data.  Best if used after z-scoring the data
        """
        if self.data_type != "zscore":
            print("The data should probably be z-scored first")
        if "quantile_normalize" in self.data_type:
            raise Exception("Data is already quantile normalized")
        # Rank mean the data
        rank_mean = self.data.loc[:, self.data_columns()].stack().groupby(
            self.data.loc[:, self.data_columns()].rank(method='first').stack().astype(int)).mean()
        # Quantile normalize based on rank mean
        quantile_norm = self.data.loc[:, self.data_columns()].rank(
            method='min').stack().astype(int).map(rank_mean).unstack()
        # Create the dataframe with barcode columns and quantile normalized data
        quantile_norm_df = concat(
            [self.data.loc[:, self.counted_barcode_columns()], quantile_norm], ignore_index=True,
            sort=False, axis=1)
        # Add back the column names
        quantile_norm_df.columns = self.counted_barcode_columns() + self.data_columns()
        if inplace:
            self.data = quantile_norm_df
            self.data_type = f"quantile normalized {self.data_type}"
            return None
        else:
            return DelDataMerged(quantile_norm_df, f"quantile normalized {self.data_type}")

    def subtract_background(self, background_name: str, inplace=False):
        """
        Subtracts background
        """
        # Create a DataFrame without the background data
        del_data = self.data.loc[:, self.data_columns()].drop(columns=[background_name])
        # Subtract the background data from the created del_data DataFrame
        background_data = self.data[background_name]
        del_data_back_sub = del_data.sub(
            background_data, axis=0)
        # Create the dataframe with the barcode columns and the subtracted data
        del_data_back_sub_df = concat([self.data.loc[:, self.counted_barcode_columns()],
                                       del_data_back_sub], ignore_index=True, sort=False, axis=1)
        # Add back the column names
        del_data_back_sub_df.columns = self.counted_barcode_columns() + del_data_back_sub.columns.tolist()
        if inplace:
            self.data = del_data_back_sub_df
            self.data_type = f"{self.data_type} background subtracted"
            return None
        else:
            return DelDataMerged(del_data_back_sub_df, f"{self.data_type} background subtracted")

    def reduce(self, min_score: float, inplace=False):
        """
        Reduces the data to only include rows with at least one sample higher than min_score.  Will do so in
        place or return a new DelDataSample
        """
        # Check to make sure the cutoff score is not above the maximum of the data
        max_value = max(self.data[self.data_columns()].max())
        if min_score > max_value:
            raise Exception(f"Reduce value cutoff of {min_score} above the maximum score within the data of {max_value}\n\
                            Choose a lower cutoff")
        # Reduce by removing any row which does not contain at least one value greater than the
        # cutoff
        reduced_data = self.data[(self.data[self.data_columns()] >= min_score).any(1)]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelDataMerged(reduced_data, self.data_type)

    def merge(self, deldata, inplace=False):
        """
        Merges a DelDataMerged into the current DelDataMerged object to create more sample columns.  
        Does not work with DelDataSample at this point
        """
        if type(deldata) == DelDataSample:
            raise Exception("Only merged data is supported at this time")
        # Make sure all merged data is of the same type, e.g. zscored
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        # Merge the two together
        merged_data = merge(self.data, deldata.data,
                            on=["Barcode_1", "Barcode_2", "Barcode_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
            return None
        else:
            return DelDataMerged(merged_data, self.data_type)

    def concat(self, deldata, inplace=False):
        """
        Concats rows of a DelDataMerged into the current DelDataMerged object to create more rows of
        counted barcodes/synthons.  Especially used to concat di-synthon counts with tri and mono
        etc.
        """
        # If the column names don't make, raise an exception
        if set(deldata.data.columns) != set(self.data.columns):
            raise Exception("Data column mismatch")
        # If the data types are not the same, e.g. zscored, raise an exception
        if self.data_type != deldata.data_type:
            raise Exception(
                "Data types do not match.  Trying to concat {deldata.data_type} into {self.data_type}")
        # Add the new rows
        concat_data = concat([self.data, deldata.data], ignore_index=True, sort=False)
        if inplace:
            self.data = concat_data
            return None
        else:
            return DelDataMerged(concat_data, self.data_type)

    def sample_data(self, sample_name: str):
        """
        Outputs a DelDataSample object from the DelDataMerged object
        """
        sample_data = self.data.loc[:, self.counted_barcode_columns() + [sample_name]]
        # Rename the sample name column to the data type
        sample_data.rename({sample_name: self.data_type}, axis=1, inplace=True)
        return DelDataSample(sample_data, self.data_type, sample_name)

    def select_samples(self, sample_names: List[str], inplace=False):
        """
        Returns a subset of samples from a DelDataMerged object.  Especially used for the
        comparison_graph
        """
        if not isinstance(sample_names, list):
            raise Exception("sample_names needs to be a list of sample names")
        sample_data = self.data.loc[:, self.counted_barcode_columns() + sample_names]
        if inplace:
            self.data = sample_data
            return None
        else:
            return DelDataMerged(sample_data, self.data_type)

    def sample_enrichment(self, inplace=False):
        """
        (sample compound count / total sample counts) / (non-sample compound count / total non-sample counts)
        A cross sample enrichment metric.  It breaks though when there are only counts within the
        sample.  This is because you get a zero divsion error.  This is the metric used by dynabind
        for their fragment DEL
        """
        counts = self.data.loc[:, self.data_columns()].values
        sample_enrichment_df = DataFrame()
        for col_index, sample_id in enumerate(self.data_columns()):
            # Get row sumns of all columns except the sample column
            non_sample_counts = np.sum(np.delete(counts, col_index, 1), axis=1)
            # Get total non-sample counts
            total_non_sample_counts = np.sum(non_sample_counts)
            # Get the sample counts with the col_index
            sample_counts = counts[:, col_index]
            # Get the total sample counts
            total_sample_counts = np.sum(sample_counts)
            # compute the equation within the function description
            sample_data = (sample_counts / total_sample_counts) / \
                (non_sample_counts / total_non_sample_counts)
            sample_enrichment_df[sample_id] = sample_data
        if inplace:
            self.data = sample_enrichment_df
            self.data_type = "sample enrichment"
            return None
        else:
            return DelDataMerged(sample_enrichment_df, "sample enrichment")

    def zscore(self, inplace=False):
        """
        Z-scores the data.  Best used when there is good sequencing coverage and the DEL library is
        not too large.
        (count - mean_count) / stadard_deviation
        """
        zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes and normalize
        for barcode_group in self.barcode_info.keys():
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            zscore_df.loc[sel_indexes, self.data_columns()] = self._zscore(
                zscore_df.loc[sel_indexes, self.data_columns()].values)
        if inplace:
            self.data_type = "zscore"
            self.data = zscore_df
            return None
        else:
            return DelDataMerged(zscore_df, data_type="zscore")

    def binomial_zscore(self, inplace=False):
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as (observed_count - expected_count)/sqrt(total_counts * expected_probability * (1-expected_probability))
        Where the denominator is the binomial standard deviation
        Which equals sqrt(total_counts) * (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        """
        binomial_zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            binomial_zscore_df.loc[sel_indexes, self.data_columns()] = self._binomial_zscore(
                binomial_zscore_df.loc[sel_indexes, self.data_columns()].values, del_library_size)
        if inplace:
            self.data_type = "binomial_zscore"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataMerged(binomial_zscore_df, data_type="binomial_zscore")

    def binomial_zscore_sample_normalized(self, inplace=False):
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as:
        (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        Where the denominator is the binomial standard deviation
        """
        binomial_zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            binomial_zscore_df.loc[sel_indexes, self.data_columns()] = self._binomial_zscore_sample_normalized(
                binomial_zscore_df.loc[sel_indexes, self.data_columns()].values, del_library_size)
        if inplace:
            self.data_type = "binomial_zscore_sample_normalized"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataMerged(binomial_zscore_df, data_type="binomial_zscore_sample_normalized")

    def enrichment(self, inplace=False):
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        enrichment_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            enrichment_df.loc[sel_indexes, self.data_columns()] = self._enrichment(
                enrichment_df.loc[sel_indexes, self.data_columns()].values, del_library_size)
        if inplace:
            self.data_type = "enrichment"
            self.data = enrichment_df
            return None
        else:
            return DelDataMerged(enrichment_df, data_type="enrichment")

    def comparison_graph(self, x_sample: str, y_sample: str, out_dir, min_score=0):
        """
        Compares two samples on a single graph where x_sample is on the x axis with the enrichment
        or count and y_sample is on the y_axis
        """
        reduced_data = self.select_samples([x_sample, y_sample]).reduce(min_score)
        max_value = max(reduced_data.data[x_sample].tolist() + reduced_data.data[y_sample].tolist())
        fig = go.Figure()
        # Get the number of synthons for each row of data
        number_synthons = reduced_data.number_synthons()
        # Create a unique list to iterate over for graphing sepearately
        unique_number_synthons = list(set(number_synthons))
        unique_number_synthons.sort(reverse=True)
        # For tri, di, and mono synthon, graph seequentially in the colors above
        colors = [None, "orange", "green", "blue", "black", "yellow", "red"]
        for num_synthons in unique_number_synthons:
            # create a dataframe only containing the rows with the number of synthons
            reduced_data_synthon_num = reduced_data.data.iloc[np.where(
                number_synthons == num_synthons)[0], :]
            # create the hover text
            barcode_texts = ["<br>".join([f"<b>{barcode_name}:<b> {barcode_value}"for barcode_name, barcode_value in
                                         zip(self.counted_barcode_columns(), barcode_values)]) for barcode_values in
                             reduced_data_synthon_num[self.counted_barcode_columns()].values]
            hovertext = [f"<b>{x_sample}:<b> {round(x, 3)}<br><b>{y_sample}:<b> {round(y, 3)}<br>{barcode_text}"
                         for x, y, barcode_text in
                         zip(reduced_data_synthon_num[x_sample], reduced_data_synthon_num[y_sample], barcode_texts)]
            # Add the synthon point plot
            fig.add_trace(go.Scatter(
                name=f"{num_synthons} synthons",
                x=reduced_data_synthon_num[x_sample].round(3),
                y=reduced_data_synthon_num[y_sample].round(3),
                mode='markers',
                marker=dict(color=colors[num_synthons]),
                hovertemplate="%{text}",
                text=hovertext
            ))
        # Add the diagonal
        fig.add_shape(type='line',
                      x0=0,
                      y0=0,
                      x1=max_value,
                      y1=max_value,
                      line=dict(color='Red'))
        # Add the sample names to the axis
        fig.update_layout(
            xaxis_title=f"{x_sample} {reduced_data.data_type}",
            yaxis_title=f"{y_sample} {reduced_data.data_type}")
        # Output hte graph
        file_name = f"{date.today()}_{x_sample}_vs_{y_sample}.{self.data_descriptor()}.2d.html"
        fig.write_html(os.path.join(out_dir, file_name))


class DelDataSample(DelData):
    """
    An object for working with sample output from NGS-Barcode-Count
    """

    def merge(self, deldata):
        """
        Merges two DelDataSample objects and outputs a DelDataMerged object
        """
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        # Merge data and renmae the data column to the sample name for each
        merged_data = merge(self.data.rename({self.data_type: self.sample_name}, axis=1),
                            deldata.data.rename({deldata.data_type: deldata.sample_name}, axis=1),
                            on=["Barcode_1", "Barcode_2", "Barcode_3"],
                            how="outer").fillna(0)
        return DelDataMerged(merged_data, self.data_type)

    def reduce(self, min_score: float, inplace=False):
        """
        Reduces the data to only include data which is higher than the min_score.  Will do so in
        place or return a new DelDataSample
        """
        max_value = max(self.data[self.data_type])
        # Make sure hte data isn't being cocmpletely reduced because the min_score is too high
        if min_score > max_value:
            raise Exception(f"Reduce value cutoff of {min_score} above the maximum score within the data of {max_value}\n\
                            Choose a lower cutoff")
        # Reduce
        reduced_data = self.data[self.data[self.data_type] >= min_score]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelDataSample(reduced_data, self.data_type, self.sample_name)

    def max_score(self) -> float:
        """
        Returns the maximum score, whether that is a zscore or a count
        """
        return max(self.data[self.data_type])

    def data_column(self) -> str:
        """
        Returns the name of the data column from the DelDataSample object
        """
        return self.data_type

    def zscore(self, inplace=False):
        """
        Z-scores the data.  Best used when there is good sequencing coverage and the DEL library is
        not too large.
        (count - mean_count) / stadard_deviation
        """
        zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes then normalize
        for barcode_group in self.barcode_info.keys():
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            zscore_df.loc[sel_indexes, self.data_columns()] =\
                self._zscore(zscore_df.loc[sel_indexes, self.data_columns()].values)
        zscore_df.rename({self.data_type: "zscore"}, axis=1, inplace=True)
        if inplace:
            self.data_type = "zscore"
            self.data = zscore_df
            return None
        else:
            return DelDataSample(zscore_df, "zscore", self.sample_name)

    def binomial_zscore(self, inplace=False):
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as (observed_count - expected_count)/sqrt(total_counts * expected_probability * (1-expected_probability))
        Where the denominator is the binomial standard deviation
        Which equals sqrt(total_counts) * (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        """
        binomial_zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            binomial_zscore_df.loc[sel_indexes, self.data_columns()] =\
                self._binomial_zscore(
                    binomial_zscore_df.loc[sel_indexes, self.data_columns()].values, del_library_size)
        binomial_zscore_df.rename({self.data_type:  "binomial_zscore"}, axis=1, inplace=True)
        if inplace:
            self.data_type = "binomial_zscore"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataSample(binomial_zscore_df, "binomial_zscore", self.sample_name)

    def binomial_zscore_sample_normalized(self, inplace=False):
        """
        See: https://pubs.acs.org/doi/10.1021/acscombsci.8b00116#
        Calculated as:
        (observed_probability - expected_probability) / sqrt(expected_probability * (1 - expected_probability))
        Where the denominator is the binomial standard deviation
        """
        binomial_zscore_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            binomial_zscore_df.loc[sel_indexes, self.data_columns()] = self._binomial_zscore_sample_normalized(
                binomial_zscore_df.loc[sel_indexes, self.data_columns()].values, del_library_size)
        binomial_zscore_df.rename({self.data_type:  "binomial_zscore_sample_normalized"}, axis=1,
                                  inplace=True)
        if inplace:
            self.data_type = "binomial_zscore_sample_normalized"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataSample(binomial_zscore_df, "binomial_zscore_sample_normalized", self.sample_name)

    def enrichment(self, inplace=False):
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        enrichment_df = self.data.copy()
        if len(self.barcode_info.keys()) == 0:
            # get indexes for tri/di/mono synthon groups to normalize together along with library
            # sizes
            self._infer_barcode_groups()
        # For each tri/di/mono synthon group, get the indexes, library sizes then normalize
        for barcode_group in self.barcode_info.keys():
            del_library_size = self.barcode_info[barcode_group]["library_size"]
            sel_indexes = self.barcode_info[barcode_group]["indexes"]
            enrichment_df.loc[sel_indexes, self.data_columns()] =\
                self._enrichment(enrichment_df.loc[sel_indexes,
                                 self.data_columns()].values, del_library_size)
        enrichment_df.rename({self.data_type: "enrichment"}, axis=1, inplace=True)
        if inplace:
            self.data_type = "enrichment"
            self.data = enrichment_df
            return None
        else:
            return DelDataSample(enrichment_df, "enrichment", self.sample_name)

    def graph_2d(self, out_dir="./", min_score=0, barcodes: Optional[List[str]] = None):
        """
        Creates a 2d graph from DelDataSample object with x-axis being the combo building block and
        y-axis the single building block.  Currently only works for 3 barcode data
        """
        if barcodes is None:
            barcodes = self.counted_barcode_columns()
        # Reduce the data
        reduced_data = self.reduce(min_score)
        # Get the max to scale size and color
        max_score = reduced_data.max_score()
        max_point_size = 12
        # Scale the sizes to the max_point_size using the max value
        sizes = reduced_data.data[reduced_data.data_column()].apply(
            lambda score: max_point_size * (score - min_score + 1) / (max_score - min_score))
        # If there are more than 3 barcodes, create two 2d graphs
        if len(barcodes) >= 3:
            self._graph_2d_3_barcodes(reduced_data, out_dir, sizes, barcodes)
        # If there are only 2 barcodes, create a single 2d graph
        elif len(barcodes) == 2:
            self._graph_2d_2_barcodes(reduced_data, out_dir, sizes, barcodes)
        else:
            raise Exception(f"Only {len(barcodes)} counted barcodes not supported at this time")

    def _graph_2d_3_barcodes(self, reduced_data, out_dir: str, sizes: Series, barcodes: List[str]):
        """
        Creates a 2d graph from DelDataSample object with x-axis being the combo building block and
        y-axis the single building block when there are 3 barcodes. Creates two graphs to represent
        both possibilities

        :barcodes: list of 3 barcodes if the first three barcodes are not wanted
        """
        # Text data for the first two barcodes together
        ab = [f"{a},{b}" for a, b in zip(reduced_data.data[barcodes[0]],
                                         reduced_data.data[barcodes[1]])]
        # Text data for the second two barcodes together
        bc = [f"{b},{c}" for b, c in zip(reduced_data.data[barcodes[1]],
                                         reduced_data.data[barcodes[2]])]
        # Create a two plot figure
        fig = make_subplots(rows=1, cols=2)
        # Create the barcode_1 + barcode_2 vs barcode_3 graph
        fig.append_trace(go.Scatter(
            x=ab,
            y=reduced_data.data[barcodes[2]],
            mode='markers',
            hovertemplate="%{text}",
            marker=dict(
                size=sizes,
                color=sizes,
                colorscale='YlOrRd',
                showscale=True,
                cmax=max(sizes),
                cmin=min([0, min(sizes)])
            ),
            text=[f"<b>{barcodes[0]}, {barcodes[1]}:<b> {x}<br><b>{barcodes[2]}:<b> {y}<br>{reduced_data.data_column()}: {round(score, 3)}" for x, y, score in
                  zip(ab, reduced_data.data[barcodes[2]], reduced_data.data[reduced_data.data_column()])]

        ), row=1, col=1)
        # Update the layout to remove tickmarks and add axis lables
        fig["layout"]["xaxis"]["title"] = f"{barcodes[0]} and {barcodes[1]}"
        fig["layout"]["xaxis"]["showticklabels"] = False
        fig["layout"]["yaxis"]["title"] = barcodes[2]
        fig["layout"]["yaxis"]["showticklabels"] = False

        # Create the barcode_1 vs barcode_2 + barcode_3 graph
        fig.append_trace(go.Scatter(
            x=bc,
            y=reduced_data.data[barcodes[0]],
            mode='markers',
            hovertemplate="%{text}",
            marker=dict(
                size=sizes,
                color=sizes,
                colorscale='YlOrRd',
                showscale=True,
                cmax=max(sizes),
                cmin=min([0, min(sizes)])
            ),
            text=[f"<b>{barcodes[1]}, {barcodes[2]}:<b> {x}<br><b>{barcodes[0]}:<b> {y}<br>{reduced_data.data_column()}: {round(score, 3)}" for x, y, score in
                  zip(bc, reduced_data.data[barcodes[0]], reduced_data.data[reduced_data.data_column()])]

        ), row=1, col=2)
        # Update the layout to remove tickmarks and add axis lables
        fig["layout"]["xaxis2"]["title"] = f"{barcodes[1]} and {barcodes[2]}"
        fig["layout"]["xaxis2"]["showticklabels"] = False
        fig["layout"]["yaxis2"]["title"] = barcodes[0]
        fig["layout"]["yaxis2"]["showticklabels"] = False
        file_name = f"{date.today()}_{reduced_data.sample_name}.{reduced_data.data_descriptor()}.2d.html"
        # Create the interactive graph
        fig.write_html(os.path.join(out_dir, file_name))

    def _graph_2d_2_barcodes(self, reduced_data, out_dir: str, sizes: Series, barcodes: List[str]):
        """
        Creates a 2d graph from DelDataSample object with each axis being a single barcode/synthon. 
        """
        fig = go.Figure(data=go.Scatter(
            x=reduced_data.data[barcodes[0]],
            y=reduced_data.data[barcodes[1]],
            mode='markers',
            hovertemplate="%{text}",
            marker=dict(
                size=sizes,
                color=sizes,
                colorscale='YlOrRd',
                showscale=True,
                cmax=max(sizes),
                cmin=min([0, min(sizes)])
            ),
            text=[f"<b>{barcodes[0]}:<b> {x}<br><b>{barcodes[1]}:<b> {y}<br>{reduced_data.data_column()}: {round(score, 3)}" for x, y, score in
                  zip(reduced_data.data[barcodes[0]], reduced_data.data[barcodes[1]], reduced_data.data[reduced_data.data_column()])]

        ))
        # Update the axis titles
        fig.update_layout(
            xaxis_title=barcodes[0],
            yaxis_title=barcodes[1])
        # Remove the tickmarks
        fig["layout"]["xaxis"]["showticklabels"] = False
        fig["layout"]["yaxis"]["showticklabels"] = False
        file_name = f"{date.today()}_{reduced_data.sample_name}.{reduced_data.data_descriptor()}.2d.html"
        fig.write_html(os.path.join(out_dir, file_name))

    def graph_3d(self, out_dir="./", min_score=0, barcodes: Optional[List[str]] = None) -> None:
        """
        Creates a 3d graph from DelDataSample object with each axis being a building block.  Currently
        only works for 3 barcode data

        :barcodes: list of 3 barcodes if the first three barcodes are not wanted
        """
        if barcodes is None:
            barcodes = self.counted_barcode_columns()
        if not len(barcodes) >= 3:
            raise Exception("At least 3 counted barcoded needed for a 3d graph")
        # Reduce the data.  Otherwise it can become too many megabytes
        reduced_data = self.reduce(min_score)
        # get the max score to scale the size and colors
        max_score = reduced_data.max_score()
        max_point_size = 12
        # Scale the sizes to the max point size, max and min score
        sizes = reduced_data.data[reduced_data.data_column()].apply(
            lambda score: max_point_size * (score - min_score + 1) / (max_score - min_score))
        # Create the 3d graph
        fig = go.Figure(data=[go.Scatter3d(
            x=reduced_data.data[barcodes[0]],
            y=reduced_data.data[barcodes[1]],
            z=reduced_data.data[barcodes[2]],
            mode='markers',
            hovertemplate="%{text}",
            marker=dict(
                size=sizes,
                color=sizes,
                colorscale='YlOrRd',
                opacity=0.8,
                showscale=True,
                cmax=max(sizes),
                cmin=min([0, min(sizes)])
            ),
            text=[f"<b>{barcodes[0]}<b>: {x}<br><b>{barcodes[1]}<b>: {y}<br><b>{barcodes[2]}<b>: {z}<br><b>{reduced_data.data_column()}<b>: {round(score, 3)}" for x, y, z, score in
                  zip(reduced_data.data[barcodes[0]], reduced_data.data[barcodes[1]], reduced_data.data[barcodes[2]], reduced_data.data[reduced_data.data_column()])]
        )])
        # Remove tick labels
        fig.update_layout(
            scene=dict(
                xaxis=dict(showticklabels=False, title_text=barcodes[0]),
                yaxis=dict(showticklabels=False, title_text=barcodes[1]),
                zaxis=dict(showticklabels=False, title_text=barcodes[2]),
            )
        )
        file_name = f"{date.today()}_{self.sample_name}.{self.data_descriptor()}.3d.html"
        fig.write_html(os.path.join(out_dir, file_name))


def read_merged(file_path: str):
    """
    Reads in merged output data from NGS-Barcode-Count and creates a DelDataMerged object
    """
    data = read_csv(file_path)
    if "Count" in data.columns:
        raise Exception("Data type is sample. Use delanalysis.read_sample()")
    return DelDataMerged(data)


def read_sample(file_path: str, sample_name: str):
    """
    Reads in sample output data from NGS-Barcode-Count and creates a DelDataMerged object
    """
    data = read_csv(file_path)
    if "Count" not in data.columns:
        raise Exception(
            "Data type is not sample. Maybe use delanalysis.read_merge() if it is the merged output data")
    return DelDataSample(data, sample_name=sample_name)
