import pandas as pd
from os.path import exists
from rna_seq_analysis.pca import PCA
from .setup import TestCase


class TestPCA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        PCA(self.settings).main(
            feature_by_sample_df=pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0),
            sample_info_df=pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0),
            sample_group_column='group',
            colors=[(0.9, 0.2, 0.5, 1.), (0.5, 0.8, 1., 1.)],
            fname='abc'
        )
        for filename in [
            'abc-proportion-explained.csv',
            'abc-sample-coordinate.csv',
            'abc-sample-coordinate.pdf',
            'abc-sample-coordinate.png',
        ]:
            with self.subTest(filename=filename):
                self.assertTrue(exists(f'{self.outdir}/pca/{filename}'))
