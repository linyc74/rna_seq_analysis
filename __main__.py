import argparse


__VERSION__ = '1.0.0-beta'


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
            'help': 'path to the sample info table (sample rows)',
        }
    },
    {
        'keys': ['-g', '--gene-info-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the gene info table (gene rows)',
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
        'keys': ['--group-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'group',
            'help': 'group column in the sample-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'rna_seq_pipeline_outdir',
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


if __name__ == '__main__':
    EntryPoint().main()
