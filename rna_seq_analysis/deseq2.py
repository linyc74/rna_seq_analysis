import pandas as pd
from .template import Processor


class DESeq2(Processor):

    count_table: str
    sample_info_table: str
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    gene_info_df: pd.DataFrame
    gene_name_column: str

    r_script: str
    statistics_csv: str
    normalized_count_csv: str
    statistics_df: pd.DataFrame
    normalized_count_df: pd.DataFrame

    def main(
            self,
            count_table: str,
            sample_info_table: str,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            gene_info_df: pd.DataFrame,
            gene_name_column: str) -> pd.DataFrame:

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_info_df = gene_info_df
        self.gene_name_column = gene_name_column

        self.check_group_name()
        self.set_output_csvs()
        self.set_r_script()
        self.run_r_script()
        self.read_deseq2_output_csvs()
        self.add_gene_name_to_statistics_df()
        self.rewrite_output_csvs()

        return self.normalized_count_df

    def check_group_name(self):
        fpath = self.sample_info_table

        df = pd.read_csv(
            fpath,
            sep=',' if fpath.endswith('.csv') else '\\t',
            index_col=0
        )

        valid_group_names = set(df[self.sample_group_column])
        for name in [self.control_group_name, self.experimental_group_name]:
            if name not in valid_group_names:
                msg = f'"{name}" does not exists in the "{self.sample_group_column}" column of "{self.sample_info_table}"'
                raise AssertionError(msg)

    def set_output_csvs(self):
        self.statistics_csv = f'{self.outdir}/deseq2-statistics.csv'
        self.normalized_count_csv = f'{self.outdir}/deseq2-normalized-count.csv'

    def set_r_script(self):
        sep1 = ',' if self.count_table.endswith('.csv') else '\\t'
        sep2 = ',' if self.sample_info_table.endswith('.csv') else '\\t'
        self.r_script = f'''\
library(DESeq2)

count_df <- read.table(
    file='{self.count_table}',
    header=TRUE,
    sep='{sep1}',
    row.names=1,
    check.names=FALSE
)

sample_sheet_df <- read.table(
    file='{self.sample_info_table}',
    header=TRUE,
    sep='{sep2}',
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

    def add_gene_name_to_statistics_df(self):
        df = self.statistics_df

        n = len(df)
        df = df.merge(
            right=self.gene_info_df[[self.gene_name_column]],
            right_index=True,
            left_index=True,
            how='left'
        )
        assert len(df) == n

        columns = df.columns.tolist()
        reordered = columns[-1:] + columns[0:-1]
        df = df[reordered]

        self.statistics_df = df

    def rewrite_output_csvs(self):
        self.statistics_df.to_csv(self.statistics_csv, index=True)
        self.normalized_count_df.to_csv(self.normalized_count_csv, index=True)
