import re

from pandas import read_csv, DataFrame, concat, merge
from scipy.stats import zscore
from typing import Optional


class DelData:
    def __init__(self, data, zscored=False):
        self.data = data
        self.__determine_file_type()
        self.zscored = zscored

    def __repr__(self):
        self.data

    def __str__(self):
        self.data

    def __determine_file_type(self):
        if "Count" in self.data.columns:
            self.file_type = "sample"
        else:
            self.file_type = "merged"

    def calculate_zscore(self, inplace=False) -> Optional[DelData]:
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

    def reduce(self, min_score: float, inplace=False) -> Optional[DelData]:
        if self.file_type == "merge":
            raise Exception(
                "Can only reduce sample data.  Either call sample_data() or input sample data")
        if self.zscored:
            reduced_data = self.data[self.data.zscore >= min_score]
        else:
            reduced_data = self.data[self.data.Count >= min_score]
        if inplace:
            self.data = reduced_data
            return None
        else:
            return DelData(reduced_data)

    def sample_data(self, sample_name: str, inplace=False) -> Optional[DelData]:
        if self.file_type == "sample":
            raise Exception("Data is already sample data")
        sample_data = self.data.iloc[:, ["BB_1", "BB_2", "BB_3", sample_name]]
        if self.zscored:
            sample_data.rename({sample_name: "zscore"}, axis=1, inplace=True)
        else:
            sample_data.rename({sample_name: "count"}, axis=1, inplace=True)
        if inplace:
            self.data = sample_data
            self.file_type = "sample"
            return None
        else:
            return DelData(sample_data, self.zscored)

    def quantile_normalize(self):
        if self.file_type == "sample":
            raise Exception("Only z-scored merged counts can be  quantile normalized")
        if not self.zscored:
            raise Exception("Run calculate_zscore() before trying to quatile normalize")
        # TODO finish this
        pass

    def graph(self, output: str, cutoff: float, sample=None):
        # TODO fill this in
        pass

    def subtract(self, sample_1: str, sample_2: str):
        # TODO fill this in
        pass

    def merge(self, deldata: DelData, inplace=False) -> Optional[DelData]:
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
            return None
        else:
            return DelData(merged_data)


def graph_3d(deldata: DelData, out_path: str, min_score: float):
    reduced_data = deldata.reduce(min_score)
    pass


def read_merge(file_path: str) -> DelData:
    data = read_csv(file_path)
    return DelData(data)


def read_sample(file_path: str) -> DelData:
    data = read_csv(file_path)
    return DelData(data)
