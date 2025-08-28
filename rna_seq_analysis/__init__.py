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
        volcano_plot_label_genes: str,
        gsea_input: str,
        gsea_gene_name_keywords: str,
        gsea_gene_set_name_keywords: str,
        colormap: str,
        invert_colors: bool,
        publication_figure: bool,
        threads: int,
        debug: bool,
        outdir: str):

    settings = Settings(
        workdir=get_temp_path(prefix='./rna_seq_analysis_workdir_'),
        outdir=outdir,
        threads=int(threads),
        debug=debug,
        mock=False,
        for_publication=publication_figure)

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
        volcano_plot_label_genes=None if volcano_plot_label_genes.lower() == 'none' else volcano_plot_label_genes.split(','),
        gsea_input=gsea_input,
        gsea_gene_name_keywords=None if gsea_gene_name_keywords.lower() == 'none' else gsea_gene_name_keywords.split(','),
        gsea_gene_set_name_keywords=None if gsea_gene_set_name_keywords.lower() == 'none' else gsea_gene_set_name_keywords.split(','),
        colormap=colormap,
        invert_colors=invert_colors
    )
