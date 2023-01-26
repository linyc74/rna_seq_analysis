import pandas as pd
from .tpm import TPM
from .pca import PCA
from .gsea import GSEA
from .deseq2 import DESeq2
from .heatmap import Heatmap
from .template import Processor


class RNASeqAnalysis(Processor):

    count_table: str
    sample_info_table: str
    gene_info_table: str
    gene_length_column: str
    heatmap_read_fraction: float
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_info_df: pd.DataFrame

    tpm_df: pd.DataFrame
    deseq2_normalized_count_df: pd.DataFrame
    deseq2_statistics_df: pd.DataFrame

    def main(
            self,
            count_table: str,
            sample_info_table: str,
            gene_info_table: str,
            gene_length_column: str,
            heatmap_read_fraction: float,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str):

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.gene_info_table = gene_info_table
        self.gene_length_column = gene_length_column
        self.heatmap_read_fraction = heatmap_read_fraction
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name

        self.read_tables()
        self.tpm()
        self.heatmap()
        self.deseq2()
        self.gsea()
        self.pca()

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

    def heatmap(self):
        Heatmap(self.settings).main(
            tpm_df=self.tpm_df,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def deseq2(self):
        self.deseq2_statistics_df, self.deseq2_normalized_count_df = DESeq2(self.settings).main(
            count_table=self.count_table,
            sample_info_table=self.sample_info_table,
            sample_group_column=self.sample_group_column,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name)

    def gsea(self):
        # self.deseq2_normalized_count_df
        GSEA(self.settings)

    def pca(self):
        # self.tpm_df
        # self.deseq2_normalized_count_df
        PCA(self.settings)
