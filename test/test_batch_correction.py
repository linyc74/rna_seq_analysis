import pandas as pd
from rna_seq_analysis.batch_correction import BatchCorrection
from .setup import TestCase


class TestBatchCorrection(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = BatchCorrection(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
            sample_info_df=pd.read_csv(f'{self.indir}/sample_info_df.csv', index_col=0),
            sample_batch_column='batch',
        )
        expected = pd.read_csv(f'{self.indir}/corrected_count_df.csv', index_col=0)
        self.assertDataFrameEqual(expected, actual)
