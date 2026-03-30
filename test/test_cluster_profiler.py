import pandas as pd
from rna_seq_analysis.cluster_profiler import ClusterProfiler
from .setup import TestCase


class TestClusterProfiler(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()
    
    def test_main(self):
        ClusterProfiler(self.settings).main(
            statistics_df=pd.read_csv(f'{self.indir}/statistics_df.csv', index_col=0),
            organism='human',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_name_column='gene_name',
            gene_q_threshold=0.1,
            pathway_p_threshold=1.0,
            pathway_q_threshold=1.0,
            show_n_pathways=20,
        )
