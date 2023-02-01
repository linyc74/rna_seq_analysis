import os
import shutil

from .template import Settings
from .tools import get_temp_path
from .rna_seq_analysis import RNASeqAnalysis


def main(
        count_table: str,
        sample_info_table: str,
        gene_info_table: str,
        gene_sets_gmt: str,
        gene_length_column: str,
        gene_name_column: str,
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

    for d in [settings.workdir, settings.outdir]:
        os.makedirs(d, exist_ok=True)

    RNASeqAnalysis(settings).main(
        count_table=count_table,
        sample_info_table=sample_info_table,
        gene_info_table=gene_info_table,
        gene_sets_gmt=gene_sets_gmt,
        gene_length_column=gene_length_column,
        gene_name_column=gene_name_column,
        heatmap_read_fraction=heatmap_read_fraction,
        sample_group_column=sample_group_column,
        control_group_name=control_group_name,
        experimental_group_name=experimental_group_name)

    if not settings.debug:
        shutil.rmtree(settings.workdir)
