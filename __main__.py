import argparse
import rna_seq_analysis


__VERSION__ = '1.0.1'


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
    {
        'keys': ['-m', '--gene-sets-gmt'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the gene sets gmt file for GSEA',
        }
    },
]
OPTIONAL = [
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
            'default': 'normal',
            'help': 'control group name in the "sample group column" (default: %(default)s)',
        }
    },
    {
        'keys': ['--experimental-group-name'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'tumor',
            'help': 'experimental group name in the "sample group column" (default: %(default)s)',
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
            heatmap_read_fraction=args.heatmap_read_fraction,
            sample_group_column=args.sample_group_column,
            control_group_name=args.control_group_name,
            experimental_group_name=args.experimental_group_name,
            threads=args.threads,
            debug=args.debug,
            outdir=args.outdir)


if __name__ == '__main__':
    EntryPoint().main()
