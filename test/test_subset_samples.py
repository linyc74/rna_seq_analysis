import pandas as pd

from .setup import TestCase
from rna_seq_analysis.subset_samples import SubsetSamples


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
