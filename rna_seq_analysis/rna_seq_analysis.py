import os
import pandas as pd
from typing import Optional
from .tpm import TPM
from .pca import PCA
from .gsea import GSEA
from .deseq2 import DESeq2
from .tools import get_files
from .heatmap import Heatmap
from .template import Processor


class RNASeqAnalysis(Processor):

    count_table: str
    sample_info_table: str
    gene_info_table: str
    gene_sets_gmt: Optional[str]
    gene_length_column: str
    gene_name_column: str
    heatmap_read_fraction: float
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    skip_deseq2_gsea: bool
    gsea_input: str

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_info_df: pd.DataFrame

    tpm_df: pd.DataFrame
    deseq2_normalized_count_df: Optional[pd.DataFrame]

    def main(
            self,
            count_table: str,
            sample_info_table: str,
            gene_info_table: str,
            gene_sets_gmt: Optional[str],
            gene_length_column: str,
            gene_name_column: str,
            heatmap_read_fraction: float,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            skip_deseq2_gsea: bool,
            gsea_input: str):

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.gene_info_table = gene_info_table
        self.gene_sets_gmt = gene_sets_gmt
        self.gene_length_column = gene_length_column
        self.gene_name_column = gene_name_column
        self.heatmap_read_fraction = heatmap_read_fraction
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.skip_deseq2_gsea = skip_deseq2_gsea
        self.gsea_input = gsea_input

        self.read_tables()
        self.tpm()
        self.deseq2()
        self.heatmap()
        self.gsea()
        self.pca()
        self.clean_up()

    def read_tables(self):
        self.count_df = self.__read(self.count_table)
        self.sample_info_df = self.__read(self.sample_info_table)
        self.gene_info_df = self.__read(self.gene_info_table)

        for df in [self.count_df, self.sample_info_df, self.gene_info_df]:
            df.index.name = None  # make all final output files clean without index names

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

    def deseq2(self):
        if self.skip_deseq2_gsea:
            self.deseq2_normalized_count_df = None
        else:
            self.deseq2_normalized_count_df = DESeq2(self.settings).main(
                count_table=self.count_table,
                sample_info_table=self.sample_info_table,
                sample_group_column=self.sample_group_column,
                control_group_name=self.control_group_name,
                experimental_group_name=self.experimental_group_name,
                gene_info_df=self.gene_info_df,
                gene_name_column=self.gene_name_column)

    def heatmap(self):
        Heatmap(self.settings).main(
            tpm_df=self.tpm_df,
            deseq2_normalized_count_df=self.deseq2_normalized_count_df,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def gsea(self):
        if self.skip_deseq2_gsea or self.gene_sets_gmt is None:
            return
        GSEA(self.settings).main(
            count_df=self.tpm_df if self.gsea_input == 'tpm' else self.deseq2_normalized_count_df,
            gene_info_df=self.gene_info_df,
            sample_info_df=self.sample_info_df,
            gene_name_column=self.gene_name_column,
            sample_group_column=self.sample_group_column,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name,
            gene_sets_gmt=self.gene_sets_gmt)

    def pca(self):
        PCA(self.settings).main(
            tpm_df=self.tpm_df,
            deseq2_normalized_count_df=self.deseq2_normalized_count_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column)

    def clean_up(self):
        CleanUp(self.settings).main()


class CleanUp(Processor):

    def main(self):
        self.collect_log_files()
        self.remove_workdir()

    def collect_log_files(self):
        log_files = get_files(
            source=self.outdir,
            endswith='.log')

        if len(log_files) == 0:
            return

        os.makedirs(f'{self.outdir}/log', exist_ok=True)
        cmd = f'mv {self.outdir}/*.log {self.outdir}/log/'
        self.call(cmd)

    def remove_workdir(self):
        if not self.debug:
            self.call(f'rm -r {self.workdir}')
