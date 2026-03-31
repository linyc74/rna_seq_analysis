import pandas as pd
from os.path import exists
from rna_seq_analysis.heatmap import Heatmap
from .setup import TestCase


class TestHeatmap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Heatmap(self.settings).main(
            feature_by_sample_df=pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0),
            heatmap_read_fraction=0.8,
            fname='abc'
        )
        for filename in [
            'abc.csv',
            'abc.pdf',
            'abc.png',
        ]:
            with self.subTest(filename=filename):
                self.assertTrue(exists(f'{self.outdir}/heatmap/{filename}'))
