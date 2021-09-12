import re

from pandas import read_csv, DataFrame, concat, merge
from scipy.stats import zscore


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
        else:
            return DelData(zscore_df_final, zscored=True)

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

    def merge(self, deldata: DelData, inplace=False):
        merged_data = merge(self.data, deldata.data,
                            on=["BB_1", "BB_2", "BB_3"],
                            how="outer").fillna(0)
        if inplace:
            self.data = merged_data
        else:
            return DelData(merged_data)


def read_merge(file_path: str) -> DelData:
    data = read_csv(file_path)
    return DelData(data)


def read_sample(file_path: str) -> DelData:
    data = read_csv(file_path)
    return DelData(data)
