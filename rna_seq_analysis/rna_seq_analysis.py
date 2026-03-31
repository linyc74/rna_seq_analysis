import os
import pandas as pd
import matplotlib.pyplot as plt
from copy import copy
from itertools import combinations
from matplotlib.colors import to_rgba
from typing import Optional, List, Tuple
from .tpm import TPM
from .gsea import GSEA
from .pca import PCA
from .deseq2 import DESeq2
from .tools import get_files
from .heatmap import Heatmap
from .template import Processor
from .batch_correction import BatchCorrection
from .cluster_profiler import ClusterProfiler


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
    control_group_name: Optional[str]
    experimental_group_name: Optional[str]
    sample_batch_column: Optional[str]
    skip_differential_analysis: bool
    volcano_plot_label_genes: Optional[List[str]]
    gsea_input: str
    gsea_gene_name_keywords: Optional[List[str]]
    gsea_gene_set_name_keywords: Optional[List[str]]
    gene_p_threshold: float
    gene_q_threshold: float
    pathway_p_threshold: float
    pathway_q_threshold: float
    organism: str
    show_n_pathways: int
    colormap: str
    invert_colors: bool

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_info_df: pd.DataFrame
    colors: List[Tuple[float, float, float, float]]

    tpm_df: pd.DataFrame
    deseq2_normalized_count_df: Optional[pd.DataFrame]
    deseq2_statistics_df: Optional[pd.DataFrame]

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
            control_group_name: Optional[str],
            experimental_group_name: Optional[str],
            sample_batch_column: Optional[str],
            skip_differential_analysis: bool,
            volcano_plot_label_genes: Optional[List[str]],
            gsea_input: str,
            gsea_gene_name_keywords: Optional[List[str]],
            gsea_gene_set_name_keywords: Optional[List[str]],
            gene_p_threshold: float,
            gene_q_threshold: float,
            pathway_p_threshold: float,
            pathway_q_threshold: float,
            organism: str,
            show_n_pathways: int,
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
        self.skip_differential_analysis = skip_differential_analysis
        self.volcano_plot_label_genes = volcano_plot_label_genes
        self.gsea_input = gsea_input
        self.gsea_gene_name_keywords = gsea_gene_name_keywords
        self.gsea_gene_set_name_keywords = gsea_gene_set_name_keywords
        self.gene_p_threshold = gene_p_threshold
        self.gene_q_threshold = gene_q_threshold
        self.pathway_p_threshold = pathway_p_threshold
        self.pathway_q_threshold = pathway_q_threshold
        self.organism = organism
        self.show_n_pathways = show_n_pathways
        self.colormap = colormap
        self.invert_colors = invert_colors

        self.preprocessing()
        self.differential_analysis()
        self.heatmap_and_pca_for_deseq2()
        CleanUp(self.settings).main()

    def preprocessing(self):
        self.count_df = read(self.count_table)
        self.sample_info_df = read(self.sample_info_table)
        self.gene_info_df = read(self.gene_info_table)

        for df in [self.count_df, self.sample_info_df, self.gene_info_df]:
            df.index.name = None  # make all final output files clean without index names

        self.count_df = SubsetSamples(self.settings).main(
            count_df=self.count_df,
            sample_info_df=self.sample_info_df)

        self.colors = GetColors(self.settings).main(
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            colormap=self.colormap,
            invert_colors=self.invert_colors)

        if self.sample_batch_column is not None:
            self.count_df = BatchCorrection(self.settings).main(
                count_df=self.count_df,
                sample_info_df=self.sample_info_df,
                sample_batch_column=self.sample_batch_column)

        self.tpm_df = TPM(self.settings).main(
            count_df=self.count_df,
            gene_info_df=self.gene_info_df,
            gene_length_column=self.gene_length_column)
        
        Heatmap(self.settings).main(
            feature_by_sample_df=self.tpm_df,
            heatmap_read_fraction=self.heatmap_read_fraction,
            fname='heatmap-tpm')

        PCA(self.settings).main(
            feature_by_sample_df=self.tpm_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            colors=self.colors,
            fname='pca-tpm')

    def differential_analysis(self):
        if self.skip_differential_analysis:
            self.deseq2_normalized_count_df = None
            self.deseq2_statistics_df = None
            return

        if self.control_group_name is None or self.experimental_group_name is None:
            comparisons = list(combinations(self.sample_info_df[self.sample_group_column].unique(), 2))
        else:
            comparisons = [(self.control_group_name, self.experimental_group_name)]
        
        msg = f'Running differential expression analysis for {len(comparisons)} comparisons:'
        for control, experimental in comparisons:
            msg += f'\n  "{control}" vs "{experimental}"'
        self.logger.info(msg)

        for control, experimental in comparisons:
            self.compare(control_group_name=control, experimental_group_name=experimental)

    def compare(self, control_group_name: str, experimental_group_name: str):
        self.logger.info(f'Running differential expression analysis for "{control_group_name}" vs "{experimental_group_name}"')

        c, e = control_group_name, experimental_group_name
        new_settings = copy(self.settings)
        new_settings.outdir = f'{self.settings.outdir}/{c}__vs__{e}'
        new_settings.workdir = f'{self.settings.workdir}/{c}__vs__{e}'
        for d in [new_settings.workdir, new_settings.outdir]:
            os.makedirs(d, exist_ok=True)

        self.deseq2_normalized_count_df, self.deseq2_statistics_df = DESeq2(new_settings).main(
            count_df=self.count_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            control_group_name=c,
            experimental_group_name=e,
            gene_info_df=self.gene_info_df,
            gene_name_column=self.gene_name_column,
            gene_description_column=self.gene_description_column,
            volcano_plot_label_genes=self.volcano_plot_label_genes,
            gene_p_threshold=self.gene_p_threshold,
            gene_q_threshold=self.gene_q_threshold,
            colors=self.colors)
    
        ClusterProfiler(new_settings).main(
            statistics_df=self.deseq2_statistics_df,
            organism=self.organism,
            control_group_name=c,
            experimental_group_name=e,
            gene_name_column=self.gene_name_column,
            gene_q_threshold=self.gene_q_threshold,
            pathway_p_threshold=self.pathway_p_threshold,
            pathway_q_threshold=self.pathway_q_threshold,
            show_n_pathways=self.show_n_pathways)
        
        if self.gene_sets_gmt is not None:
            GSEA(new_settings).main(
                count_df=self.tpm_df if self.gsea_input == 'tpm' else self.deseq2_normalized_count_df,
                gene_info_df=self.gene_info_df,
                sample_info_df=self.sample_info_df,
                gene_name_column=self.gene_name_column,
                sample_group_column=self.sample_group_column,
                control_group_name=c,
                experimental_group_name=e,
                gene_sets_gmt=self.gene_sets_gmt,
                gene_name_keywords=self.gsea_gene_name_keywords,
                gene_set_name_keywords=self.gsea_gene_set_name_keywords)

    def heatmap_and_pca_for_deseq2(self):
        if self.deseq2_normalized_count_df is None:
            return

        Heatmap(self.settings).main(
            feature_by_sample_df=self.deseq2_normalized_count_df,
            heatmap_read_fraction=self.heatmap_read_fraction,
            fname='heatmap-deseq2')

        PCA(self.settings).main(
            feature_by_sample_df=self.deseq2_normalized_count_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            colors=self.colors,
            fname='pca-deseq2')


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
