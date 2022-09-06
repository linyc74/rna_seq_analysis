from .template import Settings
from .tools import get_temp_path
from .rna_seq_analysis import RNASeqAnalysis


class Main:

    count_table: str
    sample_info_table: str
    gene_info_table: str
    gene_length_column: str
    sample_group_column: str

    settings: Settings

    def main(
            self,
            count_table: str,
            sample_info_table: str,
            gene_info_table: str,
            gene_length_column: str,
            sample_group_column: str,
            threads: int,
            debug: bool,
            outdir: str):

        self.count_table = count_table
        self.sample_info_table = sample_info_table
        self.gene_info_table = gene_info_table
        self.gene_length_column = gene_length_column
        self.sample_group_column = sample_group_column

        self.settings = Settings(
            workdir=get_temp_path(prefix='./rna_seq_analysis_workdir_'),
            outdir=outdir,
            threads=int(threads),
            debug=debug,
            mock=False)

        RNASeqAnalysis(self.settings).main(
            count_table=self.count_table,
            sample_info_table=self.sample_info_table,
            gene_info_table=self.gene_info_table,
            gene_length_column=self.gene_length_column,
            sample_group_column=self.sample_group_column)

