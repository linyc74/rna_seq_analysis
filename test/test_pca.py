import pandas as pd
from rna_seq_analysis.pca import PCA
from .setup import TestCase


class TestPCA(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):

        feature_by_sample_df = pd.read_csv(f'{self.indir}/tpm.csv', index_col=0)
        sample_info_df = pd.read_csv(f'{self.indir}/22_1209_randomize_rna_seq_data_sample_info.csv', index_col=0)
        sample_group_column = 'group'

        PCA(self.settings).main(
            feature_by_sample_df=feature_by_sample_df,
            sample_info_df=sample_info_df,
            sample_group_column=sample_group_column,
            output_prefix='pca-tpm'
        )
