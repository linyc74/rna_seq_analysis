from rna_seq_analysis.deseq2 import DESeq2
from .setup import TestCase


class TestDESeq2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        DESeq2(self.settings).main(
            count_table=f'{self.indir}/22_1209_randomize_rna_seq_data_count.csv',
            sample_info_table=f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv',
            gene_id_column='gene_id',
            sample_id_column='sample_id',
            sample_group_column='group'
        )
