import pandas as pd
from rna_seq_analysis.gsea import GSEA
from .setup import TestCase


class TestGSEA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        GSEA(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv'),
            gmt=f'{self.indir}/X.gmt'
        )
