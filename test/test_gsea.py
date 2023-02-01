import pandas as pd
from rna_seq_analysis.gsea import GSEA
from .setup import TestCase


class TestGSEA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):

        count_df = pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0)
        gene_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv', index_col=0)
        sample_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0)
        gene_sets_gmt = f'{self.indir}/h.all.v2022.1.Hs.symbols.gmt'

        GSEA(self.settings).main(
            count_df=count_df,
            gene_info_df=gene_info_df,
            sample_info_df=sample_info_df,
            gene_name_column='gene_name',
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_sets_gmt=gene_sets_gmt,
        )
