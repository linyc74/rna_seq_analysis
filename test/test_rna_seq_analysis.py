from .setup import TestCase
from rna_seq_analysis.rna_seq_analysis import RNASeqAnalysis


class TestRNASeqAnalysis(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        RNASeqAnalysis(self.settings).main(
            count_table=f'{self.indir}/22_1209_randomize_rna_seq_data_count.csv',
            sample_info_table=f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv',
            gene_info_table=f'{self.indir}/22_1209_randomize_rna_seq_data_gene_info.csv',
            gene_length_column='gene_length',
            heatmap_read_fraction=0.8,
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer'
        )
