import pandas as pd
from rna_seq_analysis.gsea import GSEA, FilterGeneSets
from .setup import TestCase


class TestGSEA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_without_filtering_gene_sets(self):

        count_df = pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0)
        gene_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv', index_col=0)
        sample_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0)
        gene_sets_gmt = f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt'

        GSEA(self.settings).main(
            count_df=count_df,
            gene_info_df=gene_info_df,
            sample_info_df=sample_info_df,
            gene_name_column='gene_name',
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_sets_gmt=gene_sets_gmt,
            gene_name_keywords=None,
            gene_set_name_keywords=None,
        )

    def test_with_filtering_gene_sets(self):

        count_df = pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0)
        gene_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv', index_col=0)
        sample_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0)
        gene_sets_gmt = f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt'

        GSEA(self.settings).main(
            count_df=count_df,
            gene_info_df=gene_info_df,
            sample_info_df=sample_info_df,
            gene_name_column='gene_name',
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_sets_gmt=gene_sets_gmt,
            gene_name_keywords=['cdkn2a'],
            gene_set_name_keywords=['NFKB'],
        )

    def test_no_gene_set_passed_filter(self):

        count_df = pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0)
        gene_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv', index_col=0)
        sample_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0)
        gene_sets_gmt = f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt'

        GSEA(self.settings).main(
            count_df=count_df,
            gene_info_df=gene_info_df,
            sample_info_df=sample_info_df,
            gene_name_column='gene_name',
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_sets_gmt=gene_sets_gmt,
            gene_name_keywords=['XXXXX'],
            gene_set_name_keywords=None,
        )


class TestFilterGeneSets(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_no_filtering(self):
        gene_sets_gmt = f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt'
        actual = FilterGeneSets(self.settings).main(
            gene_sets_gmt=gene_sets_gmt,
            gene_name_keywords=None,
            gene_set_name_keywords=None,
        )
        self.assertEqual(actual, gene_sets_gmt)

    def test_filtering(self):
        actual = FilterGeneSets(self.settings).main(
            gene_sets_gmt=f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt',
            gene_name_keywords=['cdkn2a'],
            gene_set_name_keywords=['NFKB'],
        )
        self.assertFileEqual(actual, f'{self.indir}/pre-filtered-h.all.v2023.1.Hs.symbols.gmt')
