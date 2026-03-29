import os
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from typing import List, Dict
from .template import Processor


# import R packages
r_cluster_profiler = importr('clusterProfiler')
r_enrichplot = importr('enrichplot')
r_ggplot2 = importr('ggplot2')


ORGANISM_TO_DB = {
    'human': 'org.Hs.eg.db',
    'mouse': 'org.Mm.eg.db',
    'rat':   'org.Rn.eg.db',
}
ORGANISM_TO_KEGG_CODE = {
    'human': 'hsa',
    'mouse': 'mmu',
    'rat':   'rno',
}


class ClusterProfiler(Processor):

    DSTDIR_NAME = 'clusterProfiler'

    statistics_df: pd.DataFrame
    organism: str
    control_group_name: str
    experimental_group_name: str
    gene_name_column: str
    gene_q_threshold: float
    pathway_p_threshold: float
    pathway_q_threshold: float
    show_n_pathways: int

    group_name_to_entrez_ids: Dict[str, List[str]]
    enrichment_name_to_result: Dict[str, ro.methods.RS4]  # enrichResult object from clusterProfiler
    
    def main(
            self,
            statistics_df: pd.DataFrame,
            organism: str,
            control_group_name: str,
            experimental_group_name: str,
            gene_name_column: str,
            gene_q_threshold: float,
            pathway_p_threshold: float,
            pathway_q_threshold: float,
            show_n_pathways: int):

        self.statistics_df = statistics_df
        self.organism = organism    
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_name_column = gene_name_column
        self.gene_q_threshold = gene_q_threshold
        self.pathway_p_threshold = pathway_p_threshold
        self.pathway_q_threshold = pathway_q_threshold
        self.show_n_pathways = show_n_pathways
        
        os.makedirs(f'{self.outdir}/{self.DSTDIR_NAME}', exist_ok=True)

        self.set_group_name_to_entrez_ids()

        self.enrichment_name_to_result = {}
        for group_name in self.group_name_to_entrez_ids.keys():
            self.go_enrichment(group_name)
            self.kegg_enrichment(group_name)

        for name, result in self.enrichment_name_to_result.items():
            self.bubble_plot(
                enrich_result=result,
                title=name,
                png=f'{self.outdir}/{self.DSTDIR_NAME}/{name}.png'
            )

    def set_group_name_to_entrez_ids(self):
        significant = self.statistics_df['padj'] < self.gene_q_threshold
        has_gene_symbol = self.statistics_df[self.gene_name_column].notna()
        upregulated = self.statistics_df['log2FoldChange'] >= 0
        downregulated = ~upregulated

        self.group_name_to_entrez_ids = {}

        control = significant & has_gene_symbol & downregulated
        gene_symbols = self.statistics_df.loc[control, self.gene_name_column].tolist()
        self.group_name_to_entrez_ids[self.control_group_name] = self.__to_entrez_ids(gene_symbols)

        experimental = significant & has_gene_symbol & upregulated
        gene_symbols = self.statistics_df.loc[experimental, self.gene_name_column].tolist()
        self.group_name_to_entrez_ids[self.experimental_group_name] = self.__to_entrez_ids(gene_symbols)

    def __to_entrez_ids(self, gene_symbols: List[str]) -> List[str]:
        result = r_cluster_profiler.bitr(
            ro.StrVector(gene_symbols),
            fromType = 'SYMBOL',
            toType   = 'ENTREZID',
            OrgDb    = ORGANISM_TO_DB[self.organism],
        )
        df = pandas2ri.rpy2py(result)
        return df['ENTREZID'].tolist()

    def go_enrichment(self, group_name: str):
        gene_vector = ro.StrVector(self.group_name_to_entrez_ids[group_name])
        shot_to_long = {
            'BP': 'GO Biological Process',
            'MF': 'GO Molecular Function',
            'CC': 'GO Cellular Component',
        }
        for name in ['BP', 'MF', 'CC']:
            result = r_cluster_profiler.enrichGO(
                gene          = gene_vector,
                OrgDb         = ORGANISM_TO_DB[self.organism],
                keyType       = 'ENTREZID',
                ont           = name,
                pAdjustMethod = 'BH',
                pvalueCutoff  = self.pathway_p_threshold,
                qvalueCutoff  = self.pathway_q_threshold,
            )
            enrichment_name = f'{group_name} - {shot_to_long[name]}'
            df = pandas2ri.rpy2py(result.slots['result'])
            df.to_csv(f'{self.outdir}/{self.DSTDIR_NAME}/{enrichment_name}.csv', index=True)
            self.enrichment_name_to_result[enrichment_name] = result

    def kegg_enrichment(self, group_name: str):
        gene_vector = ro.StrVector(self.group_name_to_entrez_ids[group_name])
        result = r_cluster_profiler.enrichKEGG(
            gene          = gene_vector,
            organism      = ORGANISM_TO_KEGG_CODE[self.organism],
            keyType       = 'ncbi-geneid',  # this is Entrez ID
            pAdjustMethod = 'BH',
            pvalueCutoff  = self.pathway_p_threshold,
            qvalueCutoff  = self.pathway_q_threshold,
        )
        df = pandas2ri.rpy2py(result.slots['result'])
        enrichment_name = f'{group_name} - KEGG'
        df.to_csv(f'{self.outdir}/{self.DSTDIR_NAME}/{enrichment_name}.csv', index=True)
        self.enrichment_name_to_result[enrichment_name] = result

    def bubble_plot(self, enrich_result: ro.methods.RS4, title: str, png: str):
        plot = r_enrichplot.dotplot(
            object       = enrich_result,
            x            = 'GeneRatio',
            color        = 'p.adjust',
            showCategory = self.show_n_pathways,
            title        = title,
            orderBy      = 'x',
        )
        r_ggplot2.ggsave(
            filename = png,
            plot    = plot,
            width   = 18 / 2.54,
            height  = get_height(self.show_n_pathways),
            dpi     = 600
        )


def get_height(n_pathways: int) -> float:
    padding = 2 / 2.54
    return padding + (n_pathways * 1 / 2.54)
