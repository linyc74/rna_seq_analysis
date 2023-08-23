# RNA-Seq Analysis

**Count-based RNA-seq analysis**

## Usage

```bash
git clone https://github.com/linyc74/rna_seq_analysis.git
```

The cloned repository is directly executable as a python script.

```bash
python rna_seq_analysis \
  --count-table PATH/TO/COUNT.CSV \
  --sample-info-table PATH/TO/SAMPLE_INFO.CSV \
  --gene-info-table PATH/TO/GENE_INFO.CSV \
  --gene-sets-gmt PATH/TO/GSEA_GENE_SETS.GMT
```

## Environment

Linux environment dependencies:
- [`GSEA`](https://www.gsea-msigdb.org/gsea/downloads.jsp)

Python:
- `pandas`
- `seaborn`
- `sklearn`

R:
- `DESeq2`
- `sva` (for `ComBat-seq`)
- `goseq`
- `clusterProfiler`
