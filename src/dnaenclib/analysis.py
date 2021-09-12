import re

from pandas import read_csv, DataFrame, concat
from scipy.stats import zscore


class DelData:
    def __init__(self, file_path: str):
        self.data = read_csv(file_path)
        self.__determine_file_type()
        self.zscore = None

    def __determine_file_type(self):
        if "Count" in self.data.columns:
            self.file_type = "sample"
        else:
            self.file_type = "merged"

    def calculate_zscore(self):
        data_columns = [col for col in self.data if not re.search("^BB_\d+$", col)]
        building_block_columns = [col for col in self.data if re.search("^BB_\d+$", col)]
        # Get the zscore of the count columns
        zscore_values = zscore(self.data.loc[:, data_columns].values)
        # Create a dataframe to add the blocks back
        zscore_df = DataFrame(data=zscore_values, columns=data_columns)
        self.zscore = concat([self.data.iloc[:, building_block_columns], zscore_df], sort=False, axis=1)\
            .rename({"Count": "zscore"}, axis=1)

    def quantile_normalize(self):
        if self.file_type == "sample":
            print("Can only quantile normalize merged count data")
            return
        if self.zscore == None:
            print("Run calculate_zscore() before trying to quatile normalize")
            return
        # TODO finish this
        pass

    def graph(self, output: str, cutoff: float, sample=None):
        # TODO fill this in
        pass

    def subtract(self, sample_1: str, sample_2: str):
        # TODO fill this in
        pass

    def merge(self, deldata: DelData):
        # TODO fill this in
        pass
