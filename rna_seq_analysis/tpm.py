import numpy as np
import pandas as pd
from .template import Processor


class TPM(Processor):

    df: pd.DataFrame
    gene_info_df: pd.DataFrame
    gene_length_column: str

    gene_kb_lengths: pd.Series

    def main(
            self,
            count_df: pd.DataFrame,
            gene_info_df: pd.DataFrame,
            gene_length_column: str) -> pd.DataFrame:

        self.df = count_df.copy()
        self.gene_info_df = gene_info_df
        self.gene_length_column = gene_length_column

        self.set_gene_lengths()
        self.normalize_by_gene_kb_lengths()
        self.normalize_by_million_reads_per_sample()

        return self.df

    def set_gene_lengths(self):
        df = self.df.merge(
            right=self.gene_info_df[[self.gene_length_column]],
            left_index=True,
            right_index=True,
            how='left'
        )
        self.gene_kb_lengths = df[self.gene_length_column] / 1000

    def normalize_by_gene_kb_lengths(self):
        for c in self.df.columns:
            self.df[c] = self.df[c] / self.gene_kb_lengths

    def normalize_by_million_reads_per_sample(self):
        million_reads_per_sample = np.sum(self.df, axis=0) / 1e+6
        self.df = self.df / million_reads_per_sample
