import pandas as pd
from .template import Processor


class SubsetSamples(Processor):

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame

    def main(
            self,
            count_df: pd.DataFrame,
            sample_info_df: pd.DataFrame) -> pd.DataFrame:

        self.count_df = count_df.copy()
        self.sample_info_df = sample_info_df.copy()

        for sample_id in self.sample_info_df.index:
            assert sample_id in self.count_df.columns, f'"{sample_id}" not in count_df.columns'

        self.count_df = self.count_df[self.sample_info_df.index]

        return self.count_df
