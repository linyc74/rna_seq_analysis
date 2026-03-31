import pandas as pd
from .setup import TestCase
from rna_seq_analysis.rna_seq_analysis import RNASeqAnalysis, GetColors, SubsetSamples


class TestRNASeqAnalysis(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        RNASeqAnalysis(self.settings).main(
            count_table=f'{self.indir}/22_1209_randomize_rna_seq_data_count.csv',
            sample_info_table=f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv',
            gene_info_table=f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv',
            gene_sets_gmt=f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt',
            gene_length_column='gene_length',
            gene_name_column='gene_name',
            gene_description_column='gene_description',
            heatmap_read_fraction=0.8,
            sample_group_column='group',
            control_group_name=None,
            experimental_group_name=None,
            sample_batch_column='batch',
            skip_differential_analysis=False,
            volcano_plot_label_genes=[
                'FAM238B',
                'RP1L1',
                'BOP1',
                'DRAIC',
                'EIF2B2',
            ],
            gsea_input='deseq2',
            gsea_gene_name_keywords=None,
            gsea_gene_set_name_keywords=None,
            gene_p_threshold=0.05,
            gene_q_threshold=0.5,
            pathway_p_threshold=0.05,
            pathway_q_threshold=0.5,
            organism='human',
            show_n_pathways=20,
            colormap='Set1',
            invert_colors=True
        )


class TestSubsetSamples(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = SubsetSamples(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count.csv', index_col=0),
            sample_info_df=pd.read_csv(f'{self.indir}/sample-info.csv', index_col=0)
        )
        expected = pd.read_csv(f'{self.indir}/count-subset.csv', index_col=0)
        self.assertDataFrameEqual(actual, expected)

    def test_wrong_sample_info(self):
        with self.assertRaises(AssertionError):
            SubsetSamples(self.settings).main(
                count_df=pd.read_csv(f'{self.indir}/count.csv', index_col=0),
                sample_info_df=pd.read_csv(f'{self.indir}/wrong-sample-info.csv', index_col=0)
            )


class TestGetColors(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GetColors(self.settings).main(
            sample_info_df=pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0),
            sample_group_column='group',
            colormap='Set1',
            invert_colors=True
        )
        expected = [  # two colors
            (0.21568627450980393, 0.49411764705882355, 0.7215686274509804, 1.0),
            (0.8941176470588236, 0.10196078431372549, 0.10980392156862745, 1.0),
        ]
        self.assertListEqual(expected, actual)

    def test_input_color_names(self):
        actual = GetColors(self.settings).main(
            sample_info_df=pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0),
            sample_group_column='group',
            colormap='red,green',
            invert_colors=False
        )
        expected = [
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 0.5019607843137255, 0.0, 1.0)
        ]
        self.assertListEqual(expected, actual)

    def test_input_hex_colors(self):
        actual = GetColors(self.settings).main(
            sample_info_df=pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0),
            sample_group_column='group',
            colormap='#59A257,#4A759D',
            invert_colors=False
        )
        expected = [
            (0.34901960784313724, 0.6352941176470588, 0.3411764705882353, 1.0),
            (0.2901960784313726, 0.4588235294117647, 0.615686274509804, 1.0)
        ]
        self.assertListEqual(expected, actual)

