import pandas as pd
from rna_seq_analysis.deseq2 import DESeq2, volcano_plot
from .setup import TestCase


class TestDESeq2(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        DESeq2(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
            sample_info_df=pd.read_csv(f'{self.indir}/sample_info_df.csv', index_col=0),
            sample_group_column='group',
            control_group_name='normal',
            experimental_group_name='cancer',
            gene_info_df=pd.read_csv(f'{self.indir}/gene_info_df.csv', index_col=0),
            gene_name_column='gene_name',
            gene_description_column='gene_description',
            volcano_plot_label_genes=[
                'FAM238B',
                'RP1L1',
                'BOP1',
                'DRAIC',
                'EIF2B2',
            ],
            colors=[(1.0, 0.3, 0.1, 1.0), (0.2, 0.1, 1.0, 1.0)],
        )

    def test_invalid_group_name(self):
        invalid_group_name = 'X'
        with self.assertRaises(AssertionError):
            DESeq2(self.settings).main(
                count_df=pd.read_csv(f'{self.indir}/count_df.csv', index_col=0),
                sample_info_df=pd.read_csv(f'{self.indir}/sample_info_df.csv', index_col=0),
                sample_group_column='group',
                control_group_name='normal',
                experimental_group_name=invalid_group_name,
                gene_info_df=pd.read_csv(f'{self.indir}/gene_info_df.csv', index_col=0),
                gene_name_column='gene_name',
                gene_description_column=None,
                volcano_plot_label_genes=[
                    'FAM238B',
                    'RP1L1',
                    'BOP1',
                    'DRAIC',
                    'EIF2B2',
                ],
                colors=[(1.0, 0.3, 0.1, 1.0), (0.2, 0.1, 1.0, 1.0)],
            )
        
    def test_volcano_plot(self):
        volcano_plot(
            df=pd.read_csv(f'{self.indir}/statistics_df.csv', index_col=0),
            fold_change_column='log2FoldChange',
            p_value_column='padj',
            gene_name_column='gene_name',
            png=f'{self.outdir}/padj-volcano-plot.png',
            genes_to_label=[
                'FAM238B',
                'RP1L1',
                'BOP1',
                'DRAIC',
                'EIF2B2',
                'IGLV3-9',
                'ZNF702P',
                'IGKV3D-15',
                'SSXP4',
                'MIR486-2',
                'RPL36AP55',
                'PIGM',
                'RPL35AP16',
                'CGB3',
                'XXXXXX',  # invalid gene name
            ],
            up_color=(1.0, 0.3, 0.1, 1.0),
            down_color=(0.2, 0.1, 1.0, 1.0),
        )
