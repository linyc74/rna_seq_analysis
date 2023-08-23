import os
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
        gene_description_column: str,
        heatmap_read_fraction: float,
        sample_group_column: str,
        control_group_name: str,
        experimental_group_name: str,
        sample_batch_column: str,
        skip_deseq2_gsea: bool,
        gsea_input: str,
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
        gene_sets_gmt=None if gene_sets_gmt.lower() == 'none' else gene_sets_gmt,
        gene_length_column=gene_length_column,
        gene_name_column=gene_name_column,
        gene_description_column=None if gene_description_column.lower() == 'none' else gene_description_column,
        heatmap_read_fraction=heatmap_read_fraction,
        sample_group_column=sample_group_column,
        control_group_name=control_group_name,
        experimental_group_name=experimental_group_name,
        sample_batch_column=None if sample_batch_column.lower() == 'none' else sample_batch_column,
        skip_deseq2_gsea=skip_deseq2_gsea,
        gsea_input=gsea_input)
