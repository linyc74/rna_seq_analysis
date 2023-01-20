from .template import Settings
from .tools import get_temp_path
from .rna_seq_analysis import RNASeqAnalysis


def main(
        count_table: str,
        sample_info_table: str,
        gene_info_table: str,
        gene_id_column: str,
        gene_length_column: str,
        sample_id_column: str,
        heatmap_read_fraction: float,
        sample_group_column: str,
        control_group_name: str,
        experimental_group_name: str,
        threads: int,
        debug: bool,
        outdir: str):

    settings = Settings(
        workdir=get_temp_path(prefix='./rna_seq_analysis_workdir_'),
        outdir=outdir,
        threads=int(threads),
        debug=debug,
        mock=False)

    RNASeqAnalysis(settings).main(
        count_table=count_table,
        sample_info_table=sample_info_table,
        gene_info_table=gene_info_table,
        gene_id_column=gene_id_column,
        gene_length_column=gene_length_column,
        sample_id_column=sample_id_column,
        heatmap_read_fraction=heatmap_read_fraction,
        sample_group_column=sample_group_column,
        control_group_name=control_group_name,
        experimental_group_name=experimental_group_name)
