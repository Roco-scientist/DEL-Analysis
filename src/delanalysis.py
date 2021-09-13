import os
import plotly.graph_objects as go
import re

from datetime import date
from io import StringIO
from pandas import read_csv, DataFrame, concat, merge
from scipy.stats import zscore
from typing import Optional, Type


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
        # Get the zscore of the count columns
        zscore_values = zscore(self.data.loc[:, self.data_columns()].values)
        # Create a dataframe to add the blocks back
        zscore_df = DataFrame(data=zscore_values, columns=self.data_columns())
        zscore_df_final = concat([self.data.loc[:, self.building_block_columns()], zscore_df],
                                 ignore_index=True, sort=False, axis=1)
        columns = self.building_block_columns() + self.data_columns()
        zscore_df_final.columns = columns
        zscore_df_final.rename({self.data_type: "zscore"}, axis=1, inplace=True)
        return zscore_df_final

    def data_columns(self):
        """
        Returns all column names that contain data that is not the barcodes
        """
        return [col for col in self.data.columns if not re.search("^BB_\d+$", col)]

    def building_block_columns(self):
        """
        Returns all building block barcode column names
        """
        return [col for col in self.data if re.search("^BB_\d+$", col)]


class DelDataMerged(DelData):
    """
    An object for working with merged data output from DEL-Decode
    """

    def quantile_normalize(self):
        """
        WORK IN PROGRESS
        """
        # TODO finish this
        pass

    def subtract_within(self, sample_1: str, sample_2: str):
        """
        WORK IN PROGRESS
        """
        # TODO fill this in
        pass

    def subtract_sample(self, sample_1: str, sample_2: str):
        """
        WORK IN PROGRESS
        """
        # TODO fill this in
        pass

    def merge(self, deldata, inplace=False):
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
            return None
        else:
            return DelDataMerged(merged_data, self.data_type)

    def sample_data(self, sample_name: str):
        """
        Outputs a DelDataSample object from the DelDataMerged object
        """
        sample_data = self.data.iloc[:, ["BB_1", "BB_2", "BB_3", sample_name]]
        if self.data_type == "zscore":
            sample_data.rename({sample_name: "zscore"}, axis=1, inplace=True)
        else:
            sample_data.rename({sample_name: self.data_type}, axis=1, inplace=True)
        return DelDataSample(sample_data, self.data_type, sample_name)

    def zscore(self, inplace=False):
        zscore_df = self._zscore()
        if inplace:
            self.data_type = "zscore"
            self.data = zscore_df
            return None
        else:
            return DelDataMerged(zscore_df, data_type="zscore")


class DelDataSample(DelData):
    """
    An object for working with sample output from DEL-Decode
    """

    def merge(self, deldata):
        """
        Merges two DelDataSample objects and outputs a DelDataMerged object
        """
        if self.data_type != deldata.data_type:
            raise Exception(
                f"Data types are not the same.  Trying to merge {self.data_type} into {deldata.data_type}")
        merged_data = merge(self.data.rename({self.data_type: self.sample_name}, axis=1), deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
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


def graph_3d(deldata, out_dir="./", min_score=0):
    """
    Creates a 3d graph from DelDataSample object with each axis being a building block.  Currently
    only works for 3 barcode data
    """
    if not type(deldata) == DelDataSample:
        raise Exception(
            "Only sample data can be graphed.  Try merged_data.sample_data(<sample_name>)")
    reduced_data = deldata.reduce(min_score)
    max_score = reduced_data.max_score()
    max_point_size = 12
    sizes = reduced_data.data[reduced_data.data_column()].apply(
        lambda score: max_point_size * (score - min_score + 1) / (max_score - min_score))
    fig = go.Figure(data=[go.Scatter3d(
        x=reduced_data.data.BB_1,
        y=reduced_data.data.BB_2,
        z=reduced_data.data.BB_3,
        mode='markers',
        hovertemplate="<b>BB_1<b>: %{x}<br><b>BB_2<b>: %{y}<br><b>BB_3<b>: %{z}<br>%{text}",
        marker=dict(
            size=sizes,
            color=sizes,
            colorscale='YlOrRd',
            opacity=0.8,
            showscale=True,
            cmax=max(sizes),
            cmin=min([0, min(sizes)])
        ),
        text=[f"{reduced_data.data_column()}: {score}" for score in
              reduced_data.data[reduced_data.data_column()]]
    )])
    # Remove tick labels
    fig.update_layout(
        scene=dict(
            xaxis=dict(showticklabels=False, title_text="BB_1"),
            yaxis=dict(showticklabels=False, title_text="BB_2"),
            zaxis=dict(showticklabels=False, title_text="BB_3"),
        )
    )
    file_name = f"{date.today()}_{deldata.sample_name}.html"
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


def main():
    "Setup for testing"
    data = read_sample("../../test_del/test_counts.csv", "test")
    graph_3d(data.zscore(), "../../", 4)
    breakpoint()
    pass


if __name__ == "__main__":
    main()
