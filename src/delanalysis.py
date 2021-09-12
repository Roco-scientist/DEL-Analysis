import re

from pandas import read_csv, DataFrame, concat, merge
from scipy.stats import zscore
from typing import Optional, Type


class DelData:
    def __init__(self, data, zscored=False):
        self.data = data
        self.zscored = zscored

    def __repr__(self):
        self.data

    def __str__(self):
        self.data

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
    def quantile_normalize(self):
        # TODO finish this
        pass

    def subtract_within(self, sample_1: str, sample_2: str):
        # TODO fill this in
        pass

    def subtract_sample(self, sample_1: str, sample_2: str):
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
        sample_data = self.data.iloc[:, ["BB_1", "BB_2", "BB_3", sample_name]]
        if self.zscored:
            sample_data.rename({sample_name: "zscore"}, axis=1, inplace=True)
        else:
            sample_data.rename({sample_name: "count"}, axis=1, inplace=True)
        return DelDataSample(sample_data, sample_name, self.zscored)


class DelDataSample(DelData):
    def __init__(self, data, sample_name: str, zscored=False):
        super().__init__(data, zscored)
        self.sample_name = sample_name

    def merge(self, deldata):
        if self.zscored and not deldata.zscored:
            raise Exception("Self is z-scored while merging data is not z-scored")
        if not self.zscored and deldata.zscored:
            raise Exception("Self is not z-scored while merging data is z-scored")
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        return DelDataMerged(merged_data, self.zscored)

    def reduce(self, min_score: float, inplace=False):
        if self.zscored:
            reduced_data = self.data[self.data.zscore >= min_score]
        else:
            reduced_data = self.data[self.data.Count >= min_score]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelDataSample(reduced_data, self.sample_name, self.zscored)


def graph_3d(deldata, out_path: str, min_score: float):
    reduced_data = deldata.reduce(min_score)
    pass


def read_merge(file_path: str):
    data = read_csv(file_path)
    if "Count" in data.columns:
        raise Exception("Data type is sample. Use delanalysis.read_sample()")
    return DelDataMerged(data)


def read_sample(file_path: str, sample_name: str):
    data = read_csv(file_path)
    if "Count" not in data.columns:
        raise Exception(
            "Data type is not sample. Maybe use delanalysis.read_merge() if it is the merged output data")
    return DelDataSample(data, sample_name)


def main():
    "Setup for testing"
    breakpoint()
    pass


if __name__ == "__main__":
    main()
