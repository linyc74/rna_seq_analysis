import pandas as pd
from .tpm import TPM
from .template import Processor


class RNASeqAnalysis(Processor):

    count_table: str
    sample_info_table: str
    gene_info_table: str
    gene_length_column: str
    sample_group_column: str

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_info_df: pd.DataFrame

    tpm_df: pd.DataFrame

    def main(
            self,
            count_table: str,
            sample_info_table: str,
            gene_info_table: str,
            gene_length_column: str,
            sample_group_column: str):

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.gene_info_table = gene_info_table
        self.gene_length_column = gene_length_column
        self.sample_group_column = sample_group_column

        self.read_tables()
        self.tpm()

    def read_tables(self):
        self.count_df = self.__read(self.count_table)
        self.sample_info_df = self.__read(self.sample_info_table)
        self.gene_info_df = self.__read(self.gene_info_table)

    def __read(self, file: str) -> pd.DataFrame:
        sep = ','
        for ext in ['.tsv', '.txt', '.tab']:
            if file.endswith(ext):
                sep = '\t'
                break
        return pd.read_csv(file, sep=sep, index_col=0)

    def tpm(self):
        self.tpm_df = TPM(self.settings).main(
            count_df=self.count_df,
            gene_info_df=self.gene_info_df,
            gene_length_column=self.gene_length_column)
