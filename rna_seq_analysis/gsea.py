import pandas as pd
from .template import Processor


class GSEA(Processor):

    count_df: pd.DataFrame
    gmt: str

    def main(
            self,
            count_df: pd.DataFrame,
            gmt: str):

        self.count_df = count_df
        self.gmt = gmt
