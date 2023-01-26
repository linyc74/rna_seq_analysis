import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple
from .template import Processor


class Heatmap(Processor):

    LOG_PSEUDOCOUNT = True
    NORMALIZE_BY_SAMPLE_READS = False

    tpm_df: pd.DataFrame
    heatmap_read_fraction: float

    def main(
            self,
            tpm_df: pd.DataFrame,
            heatmap_read_fraction: float):
        self.tpm_df = tpm_df
        self.heatmap_read_fraction = heatmap_read_fraction

        self.filter_by_cumulative_reads()
        self.count_normalization()
        self.clustermap()

    def filter_by_cumulative_reads(self):
        self.tpm_df = FilterByCumulativeReads(self.settings).main(
            df=self.tpm_df,
            heatmap_read_fraction=self.heatmap_read_fraction)

    def count_normalization(self):
        self.tpm_df = CountNormalization(self.settings).main(
            df=self.tpm_df,
            log_pseudocount=self.LOG_PSEUDOCOUNT,
            by_sample_reads=self.NORMALIZE_BY_SAMPLE_READS)

    def clustermap(self):
        Clustermap(self.settings).main(
            data=self.tpm_df,
            output_prefix=f'{self.outdir}/heatmap'
        )


class FilterByCumulativeReads(Processor):

    ROW_SUM = 'Row Sum'
    CUMULATIVE_SUM = 'Cumulative Sum'

    df: pd.DataFrame
    heatmap_read_fraction: float

    def main(
            self,
            df: pd.DataFrame,
            heatmap_read_fraction: float) -> pd.DataFrame:

        self.df = df
        self.heatmap_read_fraction = heatmap_read_fraction

        assert 0 < self.heatmap_read_fraction < 1

        self.sum_each_row()
        self.sort_by_sum()
        self.cumulate_sum()
        self.divide_by_total()
        self.filter()
        self.clean_up()

        return self.df

    def sum_each_row(self):
        self.df[self.ROW_SUM] = np.sum(self.df.to_numpy(), axis=1)

    def sort_by_sum(self):
        self.df = self.df.sort_values(
            by=self.ROW_SUM,
            ascending=False)

    def cumulate_sum(self):
        l = []
        c = 0
        for s in self.df[self.ROW_SUM]:
            c += s
            l.append(c)
        self.df[self.CUMULATIVE_SUM] = l

    def divide_by_total(self):
        total = self.df[self.CUMULATIVE_SUM][-1]
        self.df[self.CUMULATIVE_SUM] = self.df[self.CUMULATIVE_SUM] / total

    def filter(self):
        total_rows = len(self.df)
        n_rows = total_rows
        for i, v in enumerate(self.df[self.CUMULATIVE_SUM]):
            if v > self.heatmap_read_fraction:
                n_rows = i
                break
        self.df = self.df.iloc[:n_rows, ]

        msg = f'''\
Total genes: {total_rows}
Keep the most abundant genes covering {self.heatmap_read_fraction * 100:.2f}% of TPM: {n_rows} genes'''
        self.logger.info(msg)

    def clean_up(self):
        self.df = self.df.drop(
            columns=[self.ROW_SUM, self.CUMULATIVE_SUM]
        )


class CountNormalization(Processor):

    PSEUDOCOUNT_FACTOR = 0.1  # the factor between min count and pseudocount

    df: pd.DataFrame
    log_pseudocount: bool
    by_sample_reads: bool
    sample_reads_unit: int

    def main(
            self,
            df: pd.DataFrame,
            log_pseudocount: bool,
            by_sample_reads: bool,
            sample_reads_unit: int = 1000000) -> pd.DataFrame:

        self.df = df
        self.log_pseudocount = log_pseudocount
        self.by_sample_reads = by_sample_reads
        self.sample_reads_unit = sample_reads_unit

        self.normalize_by_sample_reads()
        self.pseudocount_then_log10()

        return self.df

    def normalize_by_sample_reads(self):
        if self.by_sample_reads:
            sum_per_column = np.sum(self.df, axis=0) / self.sample_reads_unit
            self.df = np.divide(self.df, sum_per_column)

    def pseudocount_then_log10(self):
        non_zero_min = np.nanmin(self.df.replace({0: np.nan}))
        pseudocount = non_zero_min * self.PSEUDOCOUNT_FACTOR
        if self.log_pseudocount:
            self.df = np.log10(self.df + pseudocount)


class Clustermap(Processor):

    CLUSTER_COLUMNS = True
    COLORMAP = 'PuBu'
    Y_LABEL_CHAR_WIDTH = 0.08
    X_LABEL_CHAR_WIDTH = 0.08
    XTICKLABELS = True
    YTICKLABELS = False
    CELL_WIDTH = 0.3
    CELL_HEIGHT = 0.01
    LINEWIDTH = 0.
    DENDROGRAM_RATIO = (0.1, 0.05)
    COLORBAR_WIDTH = 0.01
    COLORBAR_HORIZONTAL_POSITION = 1.
    DPI = 600

    data: pd.DataFrame
    output_prefix: str

    x_label_padding: float
    y_label_padding: float
    figsize: Tuple[float, float]
    grid: sns.matrix.ClusterGrid

    def main(self, data: pd.DataFrame, output_prefix: str):
        self.data = data
        self.output_prefix = output_prefix

        self.set_figsize()
        self.clustermap()
        self.config_clustermap()
        self.save_fig()
        self.save_csv()

    def set_figsize(self):
        self.__set_x_y_label_padding()
        w = (len(self.data.columns) * self.CELL_WIDTH) + self.y_label_padding
        h = (len(self.data.index) * self.CELL_HEIGHT) + self.x_label_padding
        self.figsize = (w, h)

    def __set_x_y_label_padding(self):
        max_x_label_length = pd.Series(self.data.columns).apply(len).max()
        self.x_label_padding = max_x_label_length * self.X_LABEL_CHAR_WIDTH

        max_y_label_length = pd.Series(self.data.index).apply(len).max()
        self.y_label_padding = max_y_label_length * self.Y_LABEL_CHAR_WIDTH

    def clustermap(self):
        self.grid = sns.clustermap(
            data=self.data,
            cmap=self.COLORMAP,
            figsize=self.figsize,
            xticklabels=self.XTICKLABELS,
            yticklabels=self.YTICKLABELS,
            col_cluster=self.CLUSTER_COLUMNS,
            dendrogram_ratio=self.DENDROGRAM_RATIO,
            linewidth=self.LINEWIDTH)
        self.__set_plotted_data()

    def __set_plotted_data(self):
        self.data = self.grid.__dict__['data2d']

    def config_clustermap(self):
        self.__set_x_y_axes()
        self.__set_colorbar()

    def __set_x_y_axes(self):
        heatmap = self.grid.ax_heatmap

        n_rows = len(self.data)
        heatmap.set_ylim(n_rows, 0)  # first and last row not be chopped half

        heatmap.tick_params(
            axis='x',
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelrotation=90
        )
        heatmap.tick_params(
            axis='y',
            left=False,
            right=True,
            labelleft=False,
            labelright=True,
            labelrotation=0
        )

    def __set_colorbar(self):
        colorbar = self.grid.cax
        p = colorbar.get_position()
        colorbar.set_position([
            self.COLORBAR_HORIZONTAL_POSITION,  # left
            p.y0,  # bottom
            self.COLORBAR_WIDTH,  # width
            p.height  # height
        ])

    def save_fig(self):
        # must use grid.savefig(), but not plt.savefig()
        # plt.savefig() crops out the colorbar

        dpi = self.__downsize_dpi_if_too_large()
        for ext in ['pdf', 'png']:
            self.grid.savefig(f'{self.output_prefix}.{ext}', dpi=dpi)
        plt.close()

    def __downsize_dpi_if_too_large(self) -> int:
        longer_side = max(self.figsize)
        dpi = self.DPI
        while longer_side * dpi >= 2**16:  # 2^16 is the limit
            dpi = int(dpi/2)  # downsize
        return dpi

    def save_csv(self):
        self.data.to_csv(f'{self.output_prefix}.csv', index=True)
