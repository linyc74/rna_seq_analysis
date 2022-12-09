import pandas as pd
from rna_seq_analysis.tpm import TPM
from .setup import TestCase


class TestTPM(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = TPM(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count-table.csv', index_col=0),
            gene_info_df=pd.read_csv(f'{self.indir}/gene-info-table.csv', index_col=0),
            gene_length_column='Gene Length'
        )
        expected = pd.read_csv(f'{self.indir}/tpm.csv', index_col=0)
        self.assertDataFrameEqual(expected, actual)
