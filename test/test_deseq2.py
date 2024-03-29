import pandas as pd
from rna_seq_analysis.deseq2 import DESeq2
from .setup import TestCase


class TestDESeq2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        DESeq2(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
            sample_info_df=pd.read_csv(f'{self.indir}/sample_info_df.csv', index_col=0),
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_info_df=pd.read_csv(f'{self.indir}/gene_info_df.csv', index_col=0),
            gene_name_column='gene_name',
            gene_description_column='gene_description',
        )

    def test_invalid_group_name(self):
        invalid_group_name = 'X'
        with self.assertRaises(AssertionError):
            DESeq2(self.settings).main(
                count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
                sample_info_df=pd.read_csv(f'{self.indir}/sample_info_df.csv', index_col=0),
                sample_group_column='group',
                control_group_name='normal',
                experimental_group_name=invalid_group_name,
                gene_info_df=pd.read_csv(f'{self.indir}/gene_info_df.csv', index_col=0),
                gene_name_column='gene_name',
                gene_description_column=None,
            )
