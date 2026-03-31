import argparse
import warnings
import rna_seq_analysis
warnings.filterwarnings('ignore')


__VERSION__ = '1.2.0'


PROG = 'python rna_seq_analysis'
DESCRIPTION = f'Count-based RNA-seq analysis (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-c', '--count-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the count table (gene rows x sample columns)',
        }
    },
    {
        'keys': ['-s', '--sample-info-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the sample info table (sample rows), the first (index) column should be sample ID',
        }
    },
    {
        'keys': ['-g', '--gene-info-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the gene info table (gene rows), the first (index) column should be gene ID',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-m', '--gene-sets-gmt'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'path to the gene sets gmt file for GSEA, if None then skip GSEA (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-length-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'gene_length',
            'help': 'gene length column in the gene-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-name-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'gene_name',
            'help': 'gene name (aka symbol) column in the gene-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-description-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'gene description column in the gene-info-table, if None then no description will be used (default: %(default)s)',
        }
    },
    {
        'keys': ['--heatmap-read-fraction'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.8,
            'help': 'fraction of TPM reads to be included in the heatmap (default: %(default)s)',
        }
    },
    {
        'keys': ['--sample-group-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'group',
            'help': 'sample group column in the sample-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['--control-group-name'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'control group name in the "sample group column" (default: %(default)s)',
        }
    },
    {
        'keys': ['--experimental-group-name'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'experimental group name in the "sample group column" (default: %(default)s)',
        }
    },
    {
        'keys': ['--sample-batch-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'sample batch column in the sample-info-table, if None then skip batch correction (default: %(default)s)',
        }
    },
    {
        'keys': ['--skip-differential-analysis'],
        'properties': {
            'action': 'store_true',
            'help': 'completely skip differential analysis including DESeq2, GSEA, and pathway enrichment, i.e. no group contrast',
        }
    },
    {
        'keys': ['--volcano-plot-label-genes'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'comma-separated gene names to be labeled on the volcano plot (default: %(default)s)',
        }
    },
    {
        'keys': ['--gsea-input'],
        'properties': {
            'type': str,
            'required': False,
            'choices': ['deseq2', 'tpm'],
            'default': 'deseq2',
            'help': 'choice of GSEA input count table (default: %(default)s)',
        }
    },
    {
        'keys': ['--gsea-gene-name-keywords'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'comma-separated gene-name keywords for gene-set pre-filtering before GSEA (default: %(default)s)',
        }
    },
    {
        'keys': ['--gsea-gene-set-name-keywords'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'comma-separated gene-set-name keywords for gene-set pre-filtering before GSEA (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-p-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.05,
            'help': 'gene-level p-value threshold for differential expression analysis (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-q-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.1,
            'help': 'gene-level q-value threshold for differential expression analysis (default: %(default)s)',
        }
    },
    {
        'keys': ['--pathway-p-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.05,
            'help': 'pathway-level p-value threshold for enrichment analysis (default: %(default)s)',
        }
    },
    {
        'keys': ['--pathway-q-threshold'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.2,
            'help': 'pathway-level q-value threshold for enrichment analysis (default: %(default)s)',
        }
    },
    {
        'keys': ['--organism'],
        'properties': {
            'type': str,
            'required': False,
            'choices': ['human', 'mouse', 'rat'],
            'default': 'human',
            'help': 'organism for clusterProfiler (default: %(default)s)',
        }
    },
    {
        'keys': ['--show-n-pathways'],
        'properties': {
            'type': int,
            'required': False,
            'default': 20,
            'help': 'number of pathways to show in the enrichment analysis (default: %(default)s)',
        }
    },
    {
        'keys': ['--colormap'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Set1',
            'help': 'matplotlib colormap for plotting, or comma-separated color names, e.g. "darkred,lightgreen,skyblue" (default: %(default)s)',
        }
    },
    {
        'keys': ['--invert-colors'],
        'properties': {
            'action': 'store_true',
            'help': 'invert the order of colors',
        }
    },
    {
        'keys': ['--publication-figure'],
        'properties': {
            'action': 'store_true',
            'help': 'plot figures in the form and quality for paper publication',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'rna_seq_analysis_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running RNA-seq Analysis version {__VERSION__}\n', flush=True)
        rna_seq_analysis.main(
            count_table=args.count_table,
            sample_info_table=args.sample_info_table,
            gene_info_table=args.gene_info_table,
            gene_sets_gmt=args.gene_sets_gmt,
            gene_length_column=args.gene_length_column,
            gene_name_column=args.gene_name_column,
            gene_description_column=args.gene_description_column,
            heatmap_read_fraction=args.heatmap_read_fraction,
            sample_group_column=args.sample_group_column,
            control_group_name=args.control_group_name,
            experimental_group_name=args.experimental_group_name,
            sample_batch_column=args.sample_batch_column,
            skip_differential_analysis=args.skip_differential_analysis,
            volcano_plot_label_genes=args.volcano_plot_label_genes,
            gsea_input=args.gsea_input,
            gsea_gene_name_keywords=args.gsea_gene_name_keywords,
            gsea_gene_set_name_keywords=args.gsea_gene_set_name_keywords,
            gene_p_threshold=args.gene_p_threshold,
            gene_q_threshold=args.gene_q_threshold,
            pathway_p_threshold=args.pathway_p_threshold,
            pathway_q_threshold=args.pathway_q_threshold,
            organism=args.organism,
            show_n_pathways=args.show_n_pathways,
            colormap=args.colormap,
            invert_colors=args.invert_colors,
            publication_figure=args.publication_figure,
            threads=args.threads,
            debug=args.debug,
            outdir=args.outdir)


if __name__ == '__main__':
    EntryPoint().main()
