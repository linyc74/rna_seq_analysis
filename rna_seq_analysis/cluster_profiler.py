import os
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.sexp import NULLType
from typing import List, Dict, Optional
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
    enrichment_pathway_keywords: Optional[List[str]]
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
            enrichment_pathway_keywords: Optional[List[str]],
            show_n_pathways: int):

        self.statistics_df = statistics_df
        self.organism = organism    
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_name_column = gene_name_column
        self.gene_q_threshold = gene_q_threshold
        self.pathway_p_threshold = pathway_p_threshold
        self.pathway_q_threshold = pathway_q_threshold
        self.enrichment_pathway_keywords = enrichment_pathway_keywords
        self.show_n_pathways = show_n_pathways
        
        os.makedirs(f'{self.outdir}/{self.DSTDIR_NAME}', exist_ok=True)

        self.set_group_name_to_entrez_ids()

        self.enrichment_name_to_result = {}

        for group_name, entrez_ids in self.group_name_to_entrez_ids.items():
            if len(entrez_ids) == 0:
                self.logger.info(f'No significantly upregulated genes found for the group "{group_name}" with q-value ≤ {self.gene_q_threshold}. Skipping GO and KEGG analysis.')
                continue
            self.go_enrichment(group_name=group_name)
            self.kegg_enrichment(group_name=group_name)

        for name in self.enrichment_name_to_result.keys():
            self.filter_pathways_by_keywords(enrichment_name=name)
            self.save_csv(enrichment_name=name)
            self.bubble_plot(enrichment_name=name)

    def set_group_name_to_entrez_ids(self):
        significant = self.statistics_df['padj'] <= self.gene_q_threshold
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
        for ontology in ['BP', 'MF', 'CC']:
            result = r_cluster_profiler.enrichGO(
                gene          = gene_vector,
                OrgDb         = ORGANISM_TO_DB[self.organism],
                keyType       = 'ENTREZID',
                ont           = ontology,
                pAdjustMethod = 'BH',
                pvalueCutoff  = self.pathway_p_threshold,
                qvalueCutoff  = self.pathway_q_threshold,
            )
            enrichment_name = f'{group_name} - {shot_to_long[ontology]}'
            
            if result is None or isinstance(result, NULLType):
                self.logger.info(f'GO enrichment returned NULL for "{enrichment_name}"')
                continue
            
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

        if result is None or isinstance(result, NULLType):
            self.logger.info(f'KEGG enrichment returned NULL for "{group_name}"')
            return

        enrichment_name = f'{group_name} - KEGG'
        self.enrichment_name_to_result[enrichment_name] = result

    def filter_pathways_by_keywords(self, enrichment_name: str):
        if self.enrichment_pathway_keywords is None:
            return

        result = self.enrichment_name_to_result[enrichment_name]
        df = pandas2ri.rpy2py(result.slots['result'])

        before = len(df)
        df = df[df['Description'].apply(self.__contain_any_keyword)].copy()
        after = len(df)

        self.logger.info(f'Using keywords to filter "{enrichment_name}" pathways: {before} -> {after}')

        self.enrichment_name_to_result[enrichment_name].slots['result'] = pandas_df_to_r_df(df)

    def __contain_any_keyword(self, description: str) -> bool:
        for k in self.enrichment_pathway_keywords:
            if k.lower() in description.lower():  # should be case-insensitive
                return True
        return False

    def save_csv(self, enrichment_name: str):
        result = self.enrichment_name_to_result[enrichment_name]
        df = r_df_to_pandas_df(result.slots['result'])
        df.to_csv(f'{self.outdir}/{self.DSTDIR_NAME}/{enrichment_name}.csv', index=True)

    def bubble_plot(self, enrichment_name: str):
        enrich_result = self.enrichment_name_to_result[enrichment_name]
        plot = r_enrichplot.dotplot(
            object       = enrich_result,
            x            = 'GeneRatio',
            color        = 'p.adjust',
            showCategory = self.show_n_pathways,
            title        = enrichment_name,
            orderBy      = 'x',
            label_format = 1000,  # max number of characters of the pathway name before it is wrapped, 1000 is basically no wrapping
        )

        n_plotted = get_n_plotted_pathways(enrichplot=plot)
        if len(plotted_pathways) == 0:
            return  # no need to save the plot

        height = get_height(n_plotted=n_plotted)

        df = r_df_to_pandas_df(enrich_result.slots['result'])
        plotted_pathways = df['Description'].tolist()[0:n_plotted]
        longest_pathway_chars = max(len(p) for p in plotted_pathways)
        width = get_width(longest_pathway_chars=longest_pathway_chars)

        r_ggplot2.ggsave(
            filename = f'{self.outdir}/{self.DSTDIR_NAME}/{enrichment_name}.png',
            plot     = plot,
            width    = width,
            height   = height,
            dpi      = 600
        )


def get_n_plotted_pathways(enrichplot: ro.methods.RS4) -> int:
    build = r_ggplot2.ggplot_build(enrichplot)
    list_vector = build.slots['data']
    df = pandas2ri.rpy2py(list_vector[0])
    return df.shape[0]


def get_height(n_plotted: int) -> float:
    min_height = 7 / 2.54  # to accommodate p value bar and size circle labels
    padding = 2 / 2.54
    calculated = padding + (n_plotted * 1 / 2.54)
    return float(max(min_height, calculated))


def get_width(longest_pathway_chars: int) -> float:
    base_width = 12 / 2.54
    char_width = (5.5 / 30) / 2.54
    return float(base_width + (longest_pathway_chars * char_width))


def r_df_to_pandas_df(r_df: ro.vectors.DataFrame) -> pd.DataFrame:
    """
    Convert an R data.frame to a pandas DataFrame.
    """
    return pandas2ri.rpy2py(r_df)


def pandas_df_to_r_df(df: pd.DataFrame) -> ro.vectors.DataFrame:
    """
    Convert a pandas DataFrame to an R data.frame.

    Uses a local converter context so pandas Series columns are converted
    correctly under rpy2 >= 3.5.
    """
    with (ro.default_converter + pandas2ri.converter).context():
        return ro.conversion.get_conversion().py2rpy(df)
