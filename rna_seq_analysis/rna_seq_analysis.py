import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from typing import Optional, List, Tuple
from .tpm import TPM
from .pca import PCA
from .gsea import GSEA
from .deseq2 import DESeq2
from .tools import get_files
from .heatmap import Heatmap
from .template import Processor
from .subset_samples import SubsetSamples
from .batch_correction import BatchCorrection


class RNASeqAnalysis(Processor):

    count_table: str
    sample_info_table: str
    gene_info_table: str
    gene_sets_gmt: Optional[str]
    gene_length_column: str
    gene_name_column: str
    gene_description_column: Optional[str]
    heatmap_read_fraction: float
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    sample_batch_column: Optional[str]
    skip_deseq2_gsea: bool
    volcano_plot_label_genes: Optional[List[str]]
    gsea_input: str
    gsea_gene_name_keywords: Optional[List[str]]
    gsea_gene_set_name_keywords: Optional[List[str]]
    colormap: str
    invert_colors: bool

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_info_df: pd.DataFrame
    colors: List[Tuple[float, float, float, float]]

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
            gene_description_column: Optional[str],
            heatmap_read_fraction: float,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            sample_batch_column: Optional[str],
            skip_deseq2_gsea: bool,
            volcano_plot_label_genes: Optional[List[str]],
            gsea_input: str,
            gsea_gene_name_keywords: Optional[List[str]],
            gsea_gene_set_name_keywords: Optional[List[str]],
            colormap: str,
            invert_colors: bool):

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.gene_info_table = gene_info_table
        self.gene_sets_gmt = gene_sets_gmt
        self.gene_length_column = gene_length_column
        self.gene_name_column = gene_name_column
        self.gene_description_column = gene_description_column
        self.heatmap_read_fraction = heatmap_read_fraction
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.sample_batch_column = sample_batch_column
        self.skip_deseq2_gsea = skip_deseq2_gsea
        self.volcano_plot_label_genes = volcano_plot_label_genes
        self.gsea_input = gsea_input
        self.gsea_gene_name_keywords = gsea_gene_name_keywords
        self.gsea_gene_set_name_keywords = gsea_gene_set_name_keywords
        self.colormap = colormap
        self.invert_colors = invert_colors

        self.read_tables()
        self.subset_samples()
        self.set_colors()
        self.batch_correction()
        self.tpm()
        self.deseq2()
        self.heatmap()
        self.pca()
        self.gsea()
        self.clean_up()

    def read_tables(self):
        self.count_df = read(self.count_table)
        self.sample_info_df = read(self.sample_info_table)
        self.gene_info_df = read(self.gene_info_table)

        for df in [self.count_df, self.sample_info_df, self.gene_info_df]:
            df.index.name = None  # make all final output files clean without index names

    def subset_samples(self):
        self.count_df = SubsetSamples(self.settings).main(
            count_df=self.count_df,
            sample_info_df=self.sample_info_df)

    def set_colors(self):
        self.colors = GetColors(self.settings).main(
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            colormap=self.colormap,
            invert_colors=self.invert_colors)

    def batch_correction(self):
        if self.sample_batch_column is not None:
            self.count_df = BatchCorrection(self.settings).main(
                count_df=self.count_df,
                sample_info_df=self.sample_info_df,
                sample_batch_column=self.sample_batch_column)

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
                count_df=self.count_df,
                sample_info_df=self.sample_info_df,
                sample_group_column=self.sample_group_column,
                control_group_name=self.control_group_name,
                experimental_group_name=self.experimental_group_name,
                gene_info_df=self.gene_info_df,
                gene_name_column=self.gene_name_column,
                gene_description_column=self.gene_description_column,
                volcano_plot_label_genes=self.volcano_plot_label_genes,
                colors=self.colors)

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
            gene_sets_gmt=self.gene_sets_gmt,
            gene_name_keywords=self.gsea_gene_name_keywords,
            gene_set_name_keywords=self.gsea_gene_set_name_keywords)

    def pca(self):
        PCA(self.settings).main(
            tpm_df=self.tpm_df,
            deseq2_normalized_count_df=self.deseq2_normalized_count_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            colors=self.colors)

    def clean_up(self):
        CleanUp(self.settings).main()


class GetColors(Processor):

    sample_info_df: pd.DataFrame
    sample_group_column: str
    colormap: str
    invert_colors: bool

    def main(
            self,
            sample_info_df: pd.DataFrame,
            sample_group_column: str,
            colormap: str,
            invert_colors: bool) -> List[Tuple[float, float, float, float]]:

        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column
        self.colormap = colormap
        self.invert_colors = invert_colors

        n_groups = len(self.sample_info_df[self.sample_group_column].unique())

        if ',' in self.colormap:
            names = self.colormap.split(',')
            if len(names) != n_groups:
                self.logger.info(f'WARNING! Number of colors "{self.colormap}" does not match number of groups ({n_groups})')
            colors = [to_rgba(n) for n in names]
        else:
            cmap = plt.colormaps[self.colormap]
            colors = [cmap(i) for i in range(n_groups)]

        if self.invert_colors:
            colors = colors[::-1]

        return colors


class CleanUp(Processor):

    def main(self):
        self.collect_files(file_ext='log', dstdir_name='log')
        self.remove_workdir()

    def collect_files(self, file_ext: str, dstdir_name: str):
        files = get_files(
            source=self.outdir,
            endswith=file_ext)

        if len(files) == 0:
            return

        d = os.path.join(self.outdir, dstdir_name)
        os.makedirs(d, exist_ok=True)
        cmd = f'mv {self.outdir}/*.{file_ext} {d}/'
        self.call(cmd)

    def remove_workdir(self):
        if not self.debug:
            self.call(f'rm -r {self.workdir}')


def read(file: str) -> pd.DataFrame:
    sep = ','
    for ext in ['.tsv', '.txt', '.tab']:
        if file.endswith(ext):
            sep = '\t'
            break
    return pd.read_csv(file, sep=sep, index_col=0)
