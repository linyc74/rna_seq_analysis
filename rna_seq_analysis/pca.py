import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, Optional
from matplotlib.axes import Axes
from sklearn import decomposition
from .template import Processor


DSTDIR_NAME = 'pca'


class PCA(Processor):

    tpm_df: pd.DataFrame
    deseq2_normalized_count_df: Optional[pd.DataFrame]
    sample_info_df: pd.DataFrame
    sample_group_column: str

    def main(
            self,
            tpm_df: pd.DataFrame,
            deseq2_normalized_count_df: Optional[pd.DataFrame],
            sample_info_df: pd.DataFrame,
            sample_group_column: str):

        self.tpm_df = tpm_df
        self.deseq2_normalized_count_df = deseq2_normalized_count_df
        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column

        OnePCA(self.settings).main(
            feature_by_sample_df=self.tpm_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            fname='pca-tpm'
        )
        if self.deseq2_normalized_count_df is not None:
            OnePCA(self.settings).main(
                feature_by_sample_df=self.deseq2_normalized_count_df,
                sample_info_df=self.sample_info_df,
                sample_group_column=self.sample_group_column,
                fname='pca-deseq2'
            )


class OnePCA(Processor):

    XY_COLUMNS = ('PC 1', 'PC 2')

    feature_by_sample_df: pd.DataFrame  # row features x column samples
    sample_info_df: pd.DataFrame
    sample_group_column: str
    fname: str

    sample_coordinate_df: pd.DataFrame
    proportion_explained_series: pd.Series

    def main(
            self,
            feature_by_sample_df: pd.DataFrame,
            sample_info_df: pd.DataFrame,
            sample_group_column: str,
            fname: str):

        self.feature_by_sample_df = feature_by_sample_df
        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column
        self.fname = fname

        self.compute_pca()
        self.merge_group_info()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()
        self.write_proportion_explained()

    def compute_pca(self):
        self.sample_coordinate_df, self.proportion_explained_series = ComputePCA(self.settings).main(
            feature_by_sample_df=self.feature_by_sample_df)

    def merge_group_info(self):
        self.sample_coordinate_df = self.sample_coordinate_df.merge(
            right=self.sample_info_df,
            right_index=True,
            left_index=True,
            how='left'
        )

    def make_dstdir(self):
        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)

    def write_sample_coordinate(self):
        self.sample_coordinate_df.to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/{self.fname}-sample-coordinate.csv'
        )

    def plot_sample_coordinate(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.sample_group_column,
            x_label_suffix=f' ({self.proportion_explained_series[0]*100:.2f}%)',
            y_label_suffix=f' ({self.proportion_explained_series[1]*100:.2f}%)',
            fname=f'{self.fname}-sample-coordinate'
        )

    def write_proportion_explained(self):
        self.proportion_explained_series.to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/{self.fname}-proportion-explained.csv',
            header=['Proportion Explained']
        )


class ComputePCA(Processor):

    XY_COLUMNS = ('PC 1', 'PC 2')
    N_COMPONENTS = 2
    RANDOM_STATE = 1  # to ensure reproducible result

    feature_by_sample_df: pd.DataFrame

    sample_by_feature_df: pd.DataFrame
    embedding: decomposition.PCA
    sample_coordinate_df: pd.DataFrame
    proportion_explained_series: pd.Series

    def main(self, feature_by_sample_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:

        self.feature_by_sample_df = feature_by_sample_df

        self.transpose()
        self.set_embedding()
        self.fit_transform()
        self.set_proportion_explained_serise()

        return self.sample_coordinate_df, self.proportion_explained_series

    def transpose(self):
        self.sample_by_feature_df = self.feature_by_sample_df.transpose()

    def set_embedding(self):
        self.embedding = decomposition.PCA(
            n_components=self.N_COMPONENTS,
            copy=True,
            whiten=False,
            svd_solver='auto',
            tol=0.0,
            iterated_power='auto',
            random_state=self.RANDOM_STATE)

    def fit_transform(self):
        transformed = self.embedding.fit_transform(
            self.sample_by_feature_df.to_numpy()
        )

        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.XY_COLUMNS,
            index=self.sample_by_feature_df.index
        )

    def set_proportion_explained_serise(self):
        self.proportion_explained_series = pd.Series(self.embedding.explained_variance_ratio_)


class ScatterPlot(Processor):

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    fname: str

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            x_label_suffix: str,
            y_label_suffix: str,
            fname: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
        self.x_label_suffix = x_label_suffix
        self.y_label_suffix = y_label_suffix
        self.fname = fname

        self.set_figsize()
        self.set_parameters()
        self.init_figure()
        self.scatterplot()
        self.label_points()
        self.save_figure()

    def set_figsize(self):
        df, group = self.sample_coordinate_df, self.group_column
        max_legend_chrs = get_max_str_length(df[group])
        n_groups = len(df[group].unique())

        self.figsize = GetFigsize(self.settings).main(
            max_legend_chrs=max_legend_chrs, n_groups=n_groups)

    def set_parameters(self):
        if self.settings.for_publication:
            self.point_size = 20.
            self.marker_edge_color = 'white'
            self.line_width = 0.5
            self.fontsize = 7
            self.dpi = 600
        else:
            self.point_size = 30.
            self.marker_edge_color = 'white'
            self.line_width = 1.0
            self.fontsize = 10
            self.dpi = 300

    def init_figure(self):
        plt.figure(figsize=self.figsize, dpi=self.dpi)

    def scatterplot(self):
        self.ax = sns.scatterplot(
            data=self.sample_coordinate_df,
            x=self.x_column,
            y=self.y_column,
            hue=self.group_column)
        plt.gca().xaxis.set_tick_params(width=self.line_width)
        plt.gca().yaxis.set_tick_params(width=self.line_width)
        plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))  # scientific notation, single digit tick labels to avoid squeezing the rectangle
        plt.xlabel(f'{self.x_column}{self.x_label_suffix}')
        plt.ylabel(f'{self.y_column}{self.y_label_suffix}')
        legend = plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        legend.set_frame_on(False)

    def label_points(self):
        if self.settings.for_publication:
            return
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name
            )

    def save_figure(self):
        plt.tight_layout()
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.outdir}/{DSTDIR_NAME}/{self.fname}.{ext}', dpi=self.dpi)
        plt.close()


def get_max_str_length(series: pd.Series) -> int:
    return series.astype(str).apply(len).max()


class GetFigsize(Processor):

    def main(self, max_legend_chrs: int, n_groups: int) -> Tuple[float, float]:

        if self.settings.for_publication:
            base_width = 7.6 / 2.54
            chr_width = 0.15 / 2.54
            base_height = 6 / 2.54
            line_height = 0.4 / 2.54
        else:
            base_width = 14 / 2.54
            chr_width = 0.218 / 2.54
            base_height = 12 / 2.54
            line_height = 0.5 / 2.54

        w = base_width + (max_legend_chrs * chr_width)
        h = max(base_height, n_groups * line_height)

        return w, h
