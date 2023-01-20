import pandas as pd
from rna_seq_analysis.heatmap import Heatmap
from .setup import TestCase


class TestHeatmap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Heatmap(self.settings).main(
            tpm_df=pd.read_csv(f'{self.indir}/tpm.csv', index_col=0),
            heatmap_read_fraction=0.8
        )
