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
            tpm_df=pd.read_csv(f'{self.indir}/tpm.csv', index_col=0),
            deseq2_normalized_count_df=pd.read_csv(f'{self.indir}/deseq2_normalized_count.csv', index_col=0),
            heatmap_read_fraction=0.8
        )
        for filename in [
            'heatmap-tpm.csv',
            'heatmap-tpm.pdf',
            'heatmap-tpm.png',
            'heatmap-deseq2.csv',
            'heatmap-deseq2.pdf',
            'heatmap-deseq2.png',
        ]:
            with self.subTest(filename=filename):
                self.assertTrue(exists(f'{self.outdir}/heatmap/{filename}'))

    def test_no_deseq2_input(self):
        Heatmap(self.settings).main(
            tpm_df=pd.read_csv(f'{self.indir}/tpm.csv', index_col=0),
            deseq2_normalized_count_df=None,
            heatmap_read_fraction=0.8
        )
        for filename in [
            'heatmap-tpm.csv',
            'heatmap-tpm.pdf',
            'heatmap-tpm.png',
        ]:
            with self.subTest(filename=filename):
                self.assertTrue(exists(f'{self.outdir}/heatmap/{filename}'))
