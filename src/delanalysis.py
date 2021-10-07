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
    An object for working with data output from DEL-Decode
    """

    def __init__(self, data, data_type="Count", sample_name=None):
        self.data = data
        self.data_type = data_type
        self.sample_name = sample_name

    def __repr__(self):
        """
        Return a string representation of the data.
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
        Return a string print of the data.
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

    def _zscore(self):
        if self.data_type == "zscore":
            raise Exception("Data is already zscored")
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        # Get the zscore of the count columns
        zscore_values = zscore(self.data.loc[:, self.data_columns()].values)
        # Create a dataframe to add the blocks back
        zscore_df = DataFrame(data=zscore_values, columns=self.data_columns())
        zscore_df_final = concat([self.data.loc[:, self.counted_barcode_columns()], zscore_df],
                                 ignore_index=True, sort=False, axis=1)
        columns = self.counted_barcode_columns() + self.data_columns()
        zscore_df_final.columns = columns
        zscore_df_final.rename({self.data_type: "zscore"}, axis=1, inplace=True)
        return zscore_df_final

    def _binomial_zscore(self, del_library_size: int):
        """
        Calculated as (observed_count - expected_count)/sqrt(total_counts * expected_probability * (1-expected_probability))
        Where the denominator is the binomial standard deviation
        """
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        expected_probability = 1/del_library_size  # get the probability as 1/library size
        total_counts = self.total_counts()
        expected_counts = total_counts * expected_probability
        binomial_sds = np.sqrt(total_counts * (expected_probability * (1 - expected_probability)))
        count_minus_expected = self.data[self.data_columns()].values - expected_counts
        zscore = count_minus_expected / binomial_sds
        new_df = concat([self.data[self.counted_barcode_columns()], DataFrame(data=zscore)],
                        ignore_index=True, sort=False, axis=1)
        new_df.columns = self.counted_barcode_columns() + self.data_columns()
        new_df.rename({self.data_type: "zscore"}, axis=1, inplace=True)
        return new_df

    def _enrichment(self, del_library_size: int, inplace=False):
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        if self.data_type != "Count":
            raise Exception("This calculation is meant for raw counts")
        enrichment_score = self.data[self.data_columns()].values * \
            del_library_size / self.total_counts()
        enrichment_df = concat([self.data[self.counted_barcode_columns()],
                                DataFrame(data=enrichment_score)],
                               ignore_index=True, sort=False, axis=1)
        enrichment_df.columns = self.counted_barcode_columns() + self.data_columns()
        enrichment_df.rename({"Count": "enrichment"}, axis=1, inplace=True)
        return enrichment_df

    def total_counts(self):
        return self.data[self.data_columns()].sum(axis=0).values

    def data_columns(self):
        """
        Returns all column names that contain data that is not the barcodes
        """
        return [col for col in self.data.columns if not re.search("^Barcode(_\d+){0,1}$", col)]

    def counted_barcode_columns(self):
        """
        Returns all building block barcode column names
        """
        return [col for col in self.data if re.search("^Barcode(_\d+){0,1}$", col)]

    def number_barcodes(self) -> List[int]:
        """
        Returns the bumber of barcodes for each row of the data
        """
        barcode_data = self.data.loc[:, self.counted_barcode_columns()]
        return notna(barcode_data).sum(axis=1).tolist()

    def to_csv(self, out_file: str):
        """
        Writes the current data to a csv file
        """
        self.data.to_csv(out_file, index=False)

    def data_descriptor(self):
        """
        Returns teh data type without spaces for output files
        """
        return self.data_type.replace(" ", "_")


class DelDataMerged(DelData):
    """
    An object for working with merged data output from DEL-Decode
    """

    def quantile_normalize(self, inplace=False):
        """
        Qauntile normalizes the data.  Best if used after z-scoring the data
        """
        if self.data_type != "zscore":
            print("The data should probably be z-scored first")
        if "quantile_normalize" in self.data_type:
            raise Exception("Data is already quantile normalized")
        rank_mean = self.data.loc[:, self.data_columns()].stack().groupby(
            self.data.loc[:, self.data_columns()].rank(method='first').stack().astype(int)).mean()
        quantile_norm = self.data.loc[:, self.data_columns()].rank(
            method='min').stack().astype(int).map(rank_mean).unstack()
        quantile_norm_df = concat(
            [self.data.loc[:, self.counted_barcode_columns()], quantile_norm], ignore_index=True,
            sort=False, axis=1)
        quantile_norm_df.columns = self.counted_barcode_columns() + self.data_columns()
        if inplace:
            self.data = quantile_norm_df
            self.data_type = f"quantile normalized {self.data_type}"
            return None
        else:
            return DelDataMerged(quantile_norm_df, f"quantile normalized {self.data_type}")

    def __subtract_within(self, sample_1: str, sample_2: str):
        """
        WORK IN PROGRESS
        """
        # TODO fill this in
        pass

    def __subtract_sample(self, sample_1: str, sample_2: str):
        """
        WORK IN PROGRESS
        """
        # TODO fill this in
        pass

    def subtract_background(self, background_name: str, inplace=False):
        """
        Subtracts background
        """
        del_data = self.data.loc[:, self.data_columns()].drop(columns=[background_name])
        background_data = self.data[background_name]
        del_data_back_sub = del_data.sub(
            background_data, axis=0)
        del_data_back_sub_df = concat([self.data.loc[:, self.counted_barcode_columns()],
                                       del_data_back_sub], ignore_index=True, sort=False, axis=1)
        del_data_back_sub_df.columns = self.counted_barcode_columns() + del_data_back_sub.columns.tolist()
        if inplace:
            self.data = del_data_back_sub_df
            self.data_type = f"{self.data_type} background subtracted"
            return None
        else:
            return DelDataMerged(del_data_back_sub_df, f"{self.data_type} background subtracted")

    def reduce(self, min_score: float, inplace=False):
        """
        Reduces the data to only include with at least on sample higher than min_score.  Will do so in
        place or return a new DelDataSample
        """
        reduced_data = self.data[(self.data[self.data_columns()] >= min_score).any(1)]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelDataMerged(reduced_data, self.data_type)

    def merge(self, deldata, inplace=False):
        if type(deldata) == DelDataSample:
            raise Exception("Only merged data is supported at this time")
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        merged_data = merge(self.data, deldata.data,
                            on=["Barcode_1", "Barcode_2", "Barcode_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
            return None
        else:
            return DelDataMerged(merged_data, self.data_type)

    def concat(self, deldata, inplace=False):
        if any([col not in deldata.data.columns for col in self.data.columns])\
                or any([col not in self.data.columns for col in deldata.data.columns]):
            raise Exception("Data column mismatch")
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
        sample_data = self.data.loc[:, ["Barcode_1", "Barcode_2", "Barcode_3", sample_name]]
        if self.data_type == "zscore":
            sample_data.rename({sample_name: "zscore"}, axis=1, inplace=True)
        else:
            sample_data.rename({sample_name: self.data_type}, axis=1, inplace=True)
        return DelDataSample(sample_data, self.data_type, sample_name)

    def sample_enrichment(self, inplace=False):
        """
        (sample compound count / total sample counts) / (non-sample compound count / total non-sample counts)
        """
        samples_data = self.data.loc[:, self.data_columns()]
        total_counts = {column: count for column, count in
                        zip(self.data_columns(), samples_data.sum(0))}
        sample_enrichment_df = DataFrame()
        for column in samples_data:
            sample_data = samples_data.loc[:, column].values / total_counts[column]
            total_other_reads = 0
            for other_column in total_counts.keys():
                if other_column != column:
                    total_other_reads += total_counts[other_column]
            all_other_data = (samples_data.drop(columns=[column]).sum(1) / total_other_reads).values
            sample_enrichment = np.divide(sample_data, all_other_data)
            sample_enrichment_df[column] = sample_enrichment
        if inplace:
            self.data = sample_enrichment_df
            self.data_type = "sample enrichment"
            return None
        else:
            return DelDataMerged(sample_enrichment_df, "sample enrichment")

    def select_samples(self, sample_names: List[str], inplace=False):
        """
        Returns a subset of samples from a DelDataMerged object
        """
        if not isinstance(sample_names, list):
            raise Exception("sample_names needs to be a list of sample names")
        sample_data = self.data.loc[:, ["Barcode_1", "Barcode_2", "Barcode_3"] + sample_names]
        if inplace:
            self.data = sample_data
            return None
        else:
            return DelDataMerged(sample_data, self.data_type)

    def zscore(self, inplace=False):
        zscore_df = self._zscore()
        if inplace:
            self.data_type = "zscore"
            self.data = zscore_df
            return None
        else:
            return DelDataMerged(zscore_df, data_type="zscore")

    def binomial_zscore(self, del_library_size: int, inplace=False):
        binomial_zscore_df = self._binomial_zscore(del_library_size)
        if inplace:
            self.data_type = "binomial_zscore"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataMerged(binomial_zscore_df, data_type="binomial_zscore")

    def enrichment(self, del_library_size: int, inplace=False):
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        enrichment_df = self._enrichment(del_library_size)
        if inplace:
            self.data_type = "enrichment"
            self.data = enrichment_df
            return None
        else:
            return DelDataMerged(enrichment_df, data_type="enrichment")


class DelDataSample(DelData):
    """
    An object for working with sample output from DEL-Decode
    """

    def merge(self, deldata):
        """
        Merges two DelDataSample objects and outputs a DelDataMerged object
        Still a work in progress
        """
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        merged_data = merge(self.data.rename({self.data_type: self.sample_name}, axis=1), deldata.data,
                            on=["Barcode_1", "Barcode_2", "Barcode_3"],
                            how="outer").fillna(0)
        return DelDataMerged(merged_data, self.data_type)

    def reduce(self, min_score: float, inplace=False):
        """
        Reduces the data to only include data which is higher than the min_score.  Will do so in
        place or return a new DelDataSample
        """
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
        zscore_df = self._zscore()
        if inplace:
            self.data_type = "zscore"
            self.data = zscore_df
            return None
        else:
            return DelDataSample(zscore_df, "zscore", self.sample_name)

    def binomial_zscore(self, del_library_size: int, inplace=False):
        binomial_zscore_df = self._binomial_zscore(del_library_size)
        if inplace:
            self.data_type = "binomial_zscore"
            self.data = binomial_zscore_df
            return None
        else:
            return DelDataSample(binomial_zscore_df, "binomial_zscore", self.sample_name)

    def enrichment(self, del_library_size: int, inplace=False):
        """
        From https://doi.org/10.1177%2F2472555218757718
        count * library diversity / total sample counts
        """
        enrichment_df = self._enrichment(del_library_size)
        if inplace:
            self.data_type = "enrichment"
            self.data = enrichment_df
            return None
        else:
            return DelDataSample(enrichment_df, "enrichment", self.sample_name)


def graph_3d(deldata, out_dir="./", min_score=0, barcodes: Optional[List[str]] = None) -> None:
    """
    Creates a 3d graph from DelDataSample object with each axis being a building block.  Currently
    only works for 3 barcode data
    """
    if not type(deldata) == DelDataSample:
        raise Exception(
            "Only sample data can be graphed.  Try merged_data.sample_data(<sample_name>)")
    if barcodes is None:
        barcodes = deldata.counted_barcode_columns()
    if not len(barcodes) >= 3:
        raise Exception("At least 3 counted barcoded needed for a 3d graph")
    reduced_data = deldata.reduce(min_score)
    max_score = reduced_data.max_score()
    max_point_size = 12
    sizes = reduced_data.data[reduced_data.data_column()].apply(
        lambda score: max_point_size * (score - min_score + 1) / (max_score - min_score))
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
    file_name = f"{date.today()}_{deldata.sample_name}.{deldata.data_descriptor()}.3d.html"
    fig.write_html(os.path.join(out_dir, file_name))


def graph_2d(deldata, out_dir="./", min_score=0, barcodes: Optional[List[str]] = None):
    """
    Creates a 2d graph from DelDataSample object with x-axis being the combo building block and
    y-axis the single building block.  Currently only works for 3 barcode data
    """
    if not type(deldata) == DelDataSample:
        raise Exception(
            "Only sample data can be graphed.  Try merged_data.sample_data(<sample_name>)")
    if barcodes is None:
        barcodes = deldata.counted_barcode_columns()
    reduced_data = deldata.reduce(min_score)
    max_score = reduced_data.max_score()
    max_point_size = 12
    sizes = reduced_data.data[reduced_data.data_column()].apply(
        lambda score: max_point_size * (score - min_score + 1) / (max_score - min_score))
    if len(barcodes) >= 3:
        _graph_2d_3_barcodes(reduced_data, out_dir, sizes, barcodes)
    elif len(barcodes) == 2:
        _graph_2d_2_barcodes(reduced_data, out_dir, sizes, barcodes)
    else:
        raise Exception(f"Only {len(barcodes)} counted barcodes not supported at this time")


def _graph_2d_3_barcodes(reduced_data, out_dir: str, sizes: Series, barcodes: List[str]):
    """
    Creates a 2d graph from DelDataSample object with x-axis being the combo building block and
    y-axis the single building block when there are 3 barcodes. 
    """
    ab = [f"{a},{b}" for a, b in zip(reduced_data.data[barcodes[0]],
                                     reduced_data.data[barcodes[1]])]
    bc = [f"{b},{c}" for b, c in zip(reduced_data.data[barcodes[1]],
                                     reduced_data.data[barcodes[2]])]
    fig = make_subplots(rows=1, cols=2)
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
    fig["layout"]["xaxis"]["title"] = f"{barcodes[0]} and {barcodes[1]}"
    fig["layout"]["xaxis"]["showticklabels"] = False
    fig["layout"]["yaxis"]["title"] = barcodes[2]
    fig["layout"]["yaxis"]["showticklabels"] = False
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
    fig["layout"]["xaxis2"]["title"] = f"{barcodes[1]} and {barcodes[2]}"
    fig["layout"]["xaxis2"]["showticklabels"] = False
    fig["layout"]["yaxis2"]["title"] = barcodes[0]
    fig["layout"]["yaxis2"]["showticklabels"] = False
    file_name = f"{date.today()}_{reduced_data.sample_name}.{reduced_data.data_descriptor()}.2d.html"
    fig.write_html(os.path.join(out_dir, file_name))


def _graph_2d_2_barcodes(reduced_data, out_dir: str, sizes: Series, barcodes: List[str]):
    """
    Creates a 2d graph from DelDataSample object with x-axis being the combo building block and
    y-axis the single building block when there are 2 barcodes. 
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
    fig.update_layout(
        xaxis_title=barcodes[0],
        yaxis_title=barcodes[1])
    fig["layout"]["xaxis"]["showticklabels"] = False
    fig["layout"]["yaxis"]["showticklabels"] = False
    file_name = f"{date.today()}_{reduced_data.sample_name}.{reduced_data.data_descriptor()}.2d.html"
    fig.write_html(os.path.join(out_dir, file_name))


def comparison_graph(deldatamerged, x_sample: str, y_sample: str, out_dir, min_score=0):
    if not type(deldatamerged) == DelDataMerged:
        raise Exception("Comparison graph only works for merged data")
    reduced_data = deldatamerged.select_samples([x_sample, y_sample]).reduce(min_score)
    max_value = max(reduced_data.data[x_sample].tolist() + reduced_data.data[y_sample].tolist())
    colors = [None, "yellow", "green", "blue", "black", "orange", "red"]
    fig = go.Figure(data=go.Scatter(
        x=reduced_data.data[x_sample].round(3),
        y=reduced_data.data[y_sample].round(3),
        mode='markers',
        marker=dict(color=list(
            map(lambda index: colors[index], reduced_data.number_barcodes()))),
        hovertemplate="%{text}",
        text=[f"<b>{x_sample}:<b> {round(x, 3)}<br><b>{y_sample}:<b> {round(y, 3)}<br><b>Barcode_1:<b> {bb_1}<br><b>Barcode_2:<b> {bb_2}<br><b>Barcode_3:<b> {bb_3}" for x, y, bb_1, bb_2, bb_3 in
              zip(reduced_data.data[x_sample], reduced_data.data[y_sample], reduced_data.data.Barcode_1, reduced_data.data.Barcode_2,
                  reduced_data.data.Barcode_3)]
    ))
    fig.add_shape(type='line',
                  x0=0,
                  y0=0,
                  x1=max_value,
                  y1=max_value,
                  line=dict(color='Red'))
    fig.update_layout(
        xaxis_title=f"{x_sample} {reduced_data.data_type}",
        yaxis_title=f"{y_sample} {reduced_data.data_type}")
    file_name = f"{date.today()}_{x_sample}_vs_{y_sample}.{deldatamerged.data_descriptor()}.2d.html"
    fig.write_html(os.path.join(out_dir, file_name))


def read_merged(file_path: str):
    """
    Reads in merged output data from DEL-Decode and creates a DelDataMerged object
    """
    data = read_csv(file_path)
    if "Count" in data.columns:
        raise Exception("Data type is sample. Use delanalysis.read_sample()")
    return DelDataMerged(data)


def read_sample(file_path: str, sample_name: str):
    """
    Reads in sample output data from DEL-Decode and creates a DelDataMerged object
    """
    data = read_csv(file_path)
    if "Count" not in data.columns:
        raise Exception(
            "Data type is not sample. Maybe use delanalysis.read_merge() if it is the merged output data")
    return DelDataSample(data, sample_name=sample_name)


def _test():
    "Setup for testing"

    data = read_sample("../../test_del/test_counts.csv", "test")
    data_merge = read_merged("../../test_del/test_counts.all.csv")
    print("Transforming data")
    print("zscore")
    # data_transformed = data_merge.zscore()
    data_transformed = data_merge.binomial_zscore(len(data_merge.data.index))
    breakpoint()
    print("quantile_normalize")
    data_transformed.quantile_normalize(inplace=True)
    print("Subtracting background")
    data_transformed.subtract_background("test_1", inplace=True)
    print("Graphing")
    comparison_graph(data_transformed, "test_2", "test_3", "../../test_del/", 4)
    sample_2 = data_transformed.sample_data("test_2")
    graph_2d(sample_2, "../../test_del/", 4)
    graph_2d(sample_2, "../../test_del/", 4, barcodes=["Barcode_1", "Barcode_2"])
    graph_3d(sample_2, "../../test_del/", 4)


def main():
    _test()


if __name__ == "__main__":
    main()
