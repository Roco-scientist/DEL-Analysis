import plotly.graph_objects as go
import re

from pandas import read_csv, DataFrame, concat, merge
from scipy.stats import zscore
from typing import Optional, Type


class DelData:
    """
    An object for working with data output from DEL-Decode
    """

    def __init__(self, data, zscored=False):
        self.data = data
        self.zscored = zscored

    def __repr__(self):
        # TODO fix below to prevent error
        return print(self.data)

    def __str__(self):
        # TODO fix below to prevent error
        return print(self.data)

    def calculate_zscore(self, inplace=False):
        if self.zscored:
            raise Exception("Data is already zscored")
        data_columns = [col for col in self.data if not re.search("^BB_\d+$", col)]
        building_block_columns = [col for col in self.data if re.search("^BB_\d+$", col)]
        # Get the zscore of the count columns
        zscore_values = zscore(self.data.loc[:, data_columns].values)
        # Create a dataframe to add the blocks back
        zscore_df = DataFrame(data=zscore_values, columns=data_columns)
        zscore_df_final = concat([self.data.iloc[:, building_block_columns], zscore_df], sort=False, axis=1)\
            .rename({"Count": "zscore"}, axis=1)
        if inplace:
            self.zscored = True
            self.data = zscore_df_final
            return None
        else:
            return DelData(zscore_df_final, zscored=True)


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
        if self.zscored and not deldata.zscored:
            raise Exception("Self is z-scored while merging data is not z-scored")
        if not self.zscored and deldata.zscored:
            raise Exception("Self is not z-scored while merging data is z-scored")
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
            return None
        else:
            return DelDataMerged(merged_data, self.zscored)

    def sample_data(self, sample_name: str):
        """
        Outputs a DelDataSample object from the DelDataMerged object
        """
        sample_data = self.data.iloc[:, ["BB_1", "BB_2", "BB_3", sample_name]]
        if self.zscored:
            sample_data.rename({sample_name: "zscore"}, axis=1, inplace=True)
        else:
            sample_data.rename({sample_name: "count"}, axis=1, inplace=True)
        return DelDataSample(sample_data, sample_name, self.zscored)


class DelDataSample(DelData):
    """
    An object for working with sample output from DEL-Decode
    """

    def __init__(self, data, sample_name: str, zscored=False):
        super().__init__(data, zscored)
        self.sample_name = sample_name

    def merge(self, deldata):
        """
        Merges two DelDataSample objects and outputs a DelDataMerged object
        """
        if self.zscored and not deldata.zscored:
            raise Exception("Self is z-scored while merging data is not z-scored")
        if not self.zscored and deldata.zscored:
            raise Exception("Self is not z-scored while merging data is z-scored")
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        return DelDataMerged(merged_data, self.zscored)

    def reduce(self, min_score: float, inplace=False):
        """
        Reduces the data to only include data which is higher than the min_score.  Will do so in
        place or return a new DelDataSample
        """
        if self.zscored:
            reduced_data = self.data[self.data.zscore >= min_score]
        else:
            reduced_data = self.data[self.data.Count >= min_score]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelDataSample(reduced_data, self.sample_name, self.zscored)

    def max_score(self) -> float:
        """
        Returns the maximum score, whether that is a zscore or a count
        """
        if self.zscored:
            return max(self.data.zscore)
        else:
            return max(self.data.Count)

    def data_column(self) -> str:
        """
        Returns the name of the data column from the DelDataSample object
        """
        if self.zscored:
            return "zscore"
        else:
            return "Count"


def graph_3d(deldata, out_prefix: str, min_score: float):
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
    fig.write_html(f"{out_prefix}.html")


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
    return DelDataSample(data, sample_name)


def main():
    "Setup for testing"
    data = read_sample("../../test_del/test_counts.csv", "test")
    breakpoint()
    pass


if __name__ == "__main__":
    main()
