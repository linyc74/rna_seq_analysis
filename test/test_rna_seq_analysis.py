from .setup import TestCase
from rna_seq_analysis.rna_seq_analysis import RNASeqAnalysis


class TestRNASeqAnalysis(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        RNASeqAnalysis(self.settings).main(
            count_table=f'{self.indir}/count-table.csv',
            sample_info_table=f'{self.indir}/sample-info-table.csv',
            gene_info_table=f'{self.indir}/gene-info-table.csv',
            gene_length_column='Gene Length',
            sample_group_column='Group'
        )
