import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Tuple
from .template import Processor
from .tools import get_temp_path


class DESeq2(Processor):

    DSTDIR_NAME = 'deseq2'

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    gene_info_df: pd.DataFrame
    gene_name_column: str
    gene_description_column: Optional[str]
    volcano_plot_label_genes: Optional[List[str]]
    colors: List[Tuple[float, float, float, float]]

    count_csv: str
    sample_info_csv: str

    r_script: str
    statistics_csv: str
    normalized_count_csv: str
    statistics_df: pd.DataFrame
    normalized_count_df: pd.DataFrame

    def main(
            self,
            count_df: pd.DataFrame,
            sample_info_df: pd.DataFrame,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            gene_info_df: pd.DataFrame,
            gene_name_column: str,
            gene_description_column: Optional[str],
            volcano_plot_label_genes: Optional[List[str]],
            colors: List[Tuple[float, float, float, float]]) -> pd.DataFrame:

        self.count_df = count_df
        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_info_df = gene_info_df
        self.gene_name_column = gene_name_column
        self.gene_description_column = gene_description_column
        self.volcano_plot_label_genes = volcano_plot_label_genes
        self.colors = colors
        
        self.check_group_name()
        self.write_input_csvs()
        self.set_output_csvs()
        self.set_r_script()
        self.run_r_script()
        self.read_deseq2_output_csvs()
        self.add_gene_name_and_description_to_statistics_df()
        self.sort_statistics_df()
        self.rewrite_output_csvs()
        self.volcano_plot()

        return self.normalized_count_df

    def check_group_name(self):
        valid_group_names = set(self.sample_info_df[self.sample_group_column])
        for name in [self.control_group_name, self.experimental_group_name]:
            if name not in valid_group_names:
                msg = f'"{name}" does not exists in the "{self.sample_group_column}" column of the sample info table'
                raise AssertionError(msg)

    def write_input_csvs(self):
        self.count_csv = get_temp_path(
            prefix=f'{self.workdir}/deseq2-count-',
            suffix='.csv')
        self.sample_info_csv = get_temp_path(
            prefix=f'{self.workdir}/deseq2-sample-info-',
            suffix='.csv')

        self.count_df.to_csv(self.count_csv, index=True)
        self.sample_info_df.to_csv(self.sample_info_csv, index=True)

    def set_output_csvs(self):
        os.makedirs(f'{self.outdir}/{self.DSTDIR_NAME}', exist_ok=True)
        self.statistics_csv = f'{self.outdir}/{self.DSTDIR_NAME}/deseq2-statistics.csv'
        self.normalized_count_csv = f'{self.outdir}/{self.DSTDIR_NAME}/deseq2-normalized-count.csv'

    def set_r_script(self):
        self.r_script = f'''\
library(DESeq2)

count_df <- read.table(
    file='{self.count_csv}',
    header=TRUE,
    sep=',',
    row.names=1,
    check.names=FALSE
)

sample_sheet_df <- read.table(
    file='{self.sample_info_csv}',
    header=TRUE,
    sep=',',
    row.names=1,
    check.names=FALSE
)

# load data for deseq2
dataset <- DESeqDataSetFromMatrix(
    countData=count_df,
    colData=sample_sheet_df,
    design=~{self.sample_group_column}
)

# run deseq2
dataset <- DESeq(dataset)

# get deseq2 results
res <- results(
    dataset,
    contrast=c("{self.sample_group_column}", "{self.experimental_group_name}", "{self.control_group_name}")
)

# differential gene expression statistics
statistics_df <- data.frame(
    res,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

write.csv(
    statistics_df, 
    file = '{self.statistics_csv}'
)

# normalized count
count_df <- counts(dataset, normalized=TRUE)

write.csv(
    count_df,
    file = '{self.normalized_count_csv}'
)
'''

    def run_r_script(self):
        r_file = f'{self.outdir}/{self.DSTDIR_NAME}/deseq2.R'
        with open(r_file, 'w') as fh:
            fh.write(self.r_script)

        log = f'{self.outdir}/{self.DSTDIR_NAME}/deseq2.log'
        cmd = self.CMD_LINEBREAK.join([
            'Rscript',
            r_file,
            f'1> {log}',
            f'2> {log}'
        ])
        self.call(cmd)

    def read_deseq2_output_csvs(self):
        self.statistics_df = pd.read_csv(self.statistics_csv, index_col=0)
        self.normalized_count_df = pd.read_csv(self.normalized_count_csv, index_col=0)

    def add_gene_name_and_description_to_statistics_df(self):
        cols = [self.gene_name_column]
        if self.gene_description_column is not None:
            cols.append(self.gene_description_column)

        df = left_join(
            left=self.statistics_df,
            right=self.gene_info_df[cols]
        )

        columns = df.columns.tolist()
        if self.gene_description_column is not None:
            reordered = columns[-2:] + columns[:-2]
        else:
            reordered = columns[-1:] + columns[:-1]
        df = df[reordered]

        self.statistics_df = df

    def sort_statistics_df(self):
        self.statistics_df = self.statistics_df.sort_values(
            by=['padj', 'pvalue'],
            ascending=[True, True]
        )

    def rewrite_output_csvs(self):
        self.statistics_df.to_csv(self.statistics_csv, index=True)
        self.normalized_count_df.to_csv(self.normalized_count_csv, index=True)

    def volcano_plot(self):
        # the order of group names should be the same as the order of colors
        all_group_names = self.sample_info_df[self.sample_group_column].unique().tolist()
        up_color_index = all_group_names.index(self.experimental_group_name)
        down_color_index = all_group_names.index(self.control_group_name)

        for p_value_column in ['padj', 'pvalue']:
            volcano_plot(
            df=self.statistics_df,
            fold_change_column='log2FoldChange',
            p_value_column=p_value_column,
            gene_name_column=self.gene_name_column,
            png=f'{self.outdir}/{self.DSTDIR_NAME}/{p_value_column}-volcano-plot.png',
            genes_to_label=self.volcano_plot_label_genes,
            up_color=self.colors[up_color_index],
            down_color=self.colors[down_color_index],
        )


def left_join(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    n = len(left)
    merged = left.merge(
        right=right,
        right_index=True,
        left_index=True,
        how='left'
    )
    assert len(merged) == n
    return merged


FOLD_CHANGE_THRESHOLD = 1.0
P_THRESHOLD = 0.05
NON_SIGNIFICANT_COLOR = (0.8, 0.8, 0.8, 1.0)  # light grey
ALPHA = 0.5
SIGNIFICANT_POINT_SIZE = 12
NON_SIGNIFICANT_POINT_SIZE = 8
GENE_LABEL_HORIZONTAL_OFFSET = 6
GENE_LABEL_VERTICAL_OFFSET_SCALE = 1.0
FIGSIZE = (10 / 2.54, 8 / 2.54)
DPI = 600
FONT_SIZE = 8
GENE_LABEL_FONT_SIZE = 6


def volcano_plot(
        df: pd.DataFrame,
        fold_change_column: str,
        p_value_column: str,
        gene_name_column: str,
        png: str,
        genes_to_label: Optional[List[str]],
        up_color: Tuple[float, float, float, float],
        down_color: Tuple[float, float, float, float]):

    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[fold_change_column, p_value_column])
    df['neglog10p'] = -np.log10(df[p_value_column].clip(lower=np.finfo(float).tiny))

    x = df[fold_change_column].values
    y = df['neglog10p'].values

    significant = (np.abs(x) >= FOLD_CHANGE_THRESHOLD) & (df[p_value_column].values <= P_THRESHOLD)
    up = significant & (x > 0)
    down = significant & (x < 0)

    plt.rcParams.update({'font.size': FONT_SIZE})

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

    ax.scatter(x[~significant], y[~significant], color=NON_SIGNIFICANT_COLOR, s=NON_SIGNIFICANT_POINT_SIZE, alpha=ALPHA, edgecolor='none')
    ax.scatter(x[down], y[down], color=down_color, s=SIGNIFICANT_POINT_SIZE, alpha=ALPHA, edgecolor='none')
    ax.scatter(x[up],   y[up],   color=up_color,   s=SIGNIFICANT_POINT_SIZE, alpha=ALPHA, edgecolor='none')

    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel(f'-Log10({p_value_column})')

    # make the x axis symmetric
    right, left = ax.get_xlim()
    abs_max = max(abs(right), abs(left))
    ax.set_xlim(-abs_max, abs_max)

    if genes_to_label is not None:
        # small vertical offsets to reduce collisions
        v_offsets = np.linspace(1, 10, num=max(3, len(genes_to_label))) * GENE_LABEL_VERTICAL_OFFSET_SCALE
        i = 0
        gset = set(genes_to_label)
        subdf = df[df[gene_name_column].isin(gset)]

        for _, r in subdf.iterrows():
            gx, gy, name = float(r[fold_change_column]), float(r['neglog10p']), str(r[gene_name_column])
            dx = GENE_LABEL_HORIZONTAL_OFFSET if gx >= 0 else -GENE_LABEL_HORIZONTAL_OFFSET  # text horizontal offset (points)
            ha = 'left' if gx >= 0 else 'right'
            dy = v_offsets[i % len(v_offsets)] * (1 if (i % 2 == 0) else -1)  # alternate up/down
            i += 1

            ax.annotate(
                name,
                xy=(gx, gy),
                xytext=(dx, dy),
                textcoords='offset points',
                ha=ha,
                va='center',
                fontsize=GENE_LABEL_FONT_SIZE,
                arrowprops=dict(arrowstyle='-', lw=0.5, shrinkA=0, shrinkB=0)
            )

    plt.tight_layout()
    plt.savefig(png, dpi=DPI)
