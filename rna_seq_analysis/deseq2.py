import pandas as pd
from typing import Optional
from .template import Processor
from .tools import get_temp_path


class DESeq2(Processor):

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    gene_info_df: pd.DataFrame
    gene_name_column: str
    gene_description_column: Optional[str]

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
            gene_description_column: Optional[str]) -> pd.DataFrame:

        self.count_df = count_df
        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_info_df = gene_info_df
        self.gene_name_column = gene_name_column
        self.gene_description_column = gene_description_column

        self.check_group_name()
        self.write_input_csvs()
        self.set_output_csvs()
        self.set_r_script()
        self.run_r_script()
        self.read_deseq2_output_csvs()
        self.add_gene_name_and_description_to_statistics_df()
        self.rewrite_output_csvs()

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
        self.statistics_csv = f'{self.outdir}/deseq2-statistics.csv'
        self.normalized_count_csv = f'{self.outdir}/deseq2-normalized-count.csv'

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
        r_file = f'{self.outdir}/deseq2.R'
        with open(r_file, 'w') as fh:
            fh.write(self.r_script)

        log = f'{self.outdir}/deseq2.log'
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

    def rewrite_output_csvs(self):
        self.statistics_df.to_csv(self.statistics_csv, index=True)
        self.normalized_count_df.to_csv(self.normalized_count_csv, index=True)


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
