import pandas as pd
from typing import List, Any
from .template import Processor
from.tools import get_temp_path


class BatchCorrection(Processor):

    count_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    sample_batch_column: str

    batch_list: List[Any]
    corrected_csv: str

    def main(
            self,
            count_df: pd.DataFrame,
            sample_info_df: pd.DataFrame,
            sample_batch_column: str) -> pd.DataFrame:

        self.count_df = count_df
        self.sample_info_df = sample_info_df
        self.sample_batch_column = sample_batch_column

        self.set_batch_list()
        self.combat_seq()

        return pd.read_csv(self.corrected_csv, index_col=0)

    def set_batch_list(self):
        self.batch_list = []
        for sample_id in self.count_df.columns:
            batch = self.sample_info_df.loc[sample_id, self.sample_batch_column]
            self.batch_list.append(batch)

    def combat_seq(self):
        csv = get_temp_path(
            prefix=f'{self.workdir}/raw-count-',
            suffix='.csv')
        self.count_df.to_csv(csv, index=True)
        self.corrected_csv = ComBatSeq(self.settings).main(
            count_csv=csv,
            batch_list=self.batch_list)


class ComBatSeq(Processor):

    count_csv: str
    batch_list: List[Any]

    r_script: str

    corrected_csv: str

    def main(
            self,
            count_csv: str,
            batch_list: List[Any]) -> str:

        self.count_csv = count_csv
        self.batch_list = batch_list

        self.corrected_csv = f'{self.outdir}/batch-corrected-count.csv'
        self.write_r_script()
        self.run_r_script()

        return self.corrected_csv

    def write_r_script(self):
        b = str(self.batch_list)[1:-1]  # strip off open and close brackets '[' and ']'

        text = f'''\
library(sva)

count_df <- read.csv(
    file='{self.count_csv}',
    header=TRUE,
    row.names=1,
    check.names=FALSE
)

count_matrix = as.matrix(count_df)

batch <- c({b})

adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)

write.csv(
    adjusted, 
    file='{self.corrected_csv}'
)
'''
        self.r_script = f'{self.outdir}/combat-seq.R'
        with open(self.r_script, 'w') as fh:
            fh.write(text)

    def run_r_script(self):
        log = f'{self.outdir}/combat-seq.log'
        cmd = self.CMD_LINEBREAK.join([
            'Rscript',
            self.r_script,
            f'1> {log}',
            f'2> {log}'
        ])
        self.call(cmd)
