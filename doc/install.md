# Install

## Conda environment

```bash
conda create -n rna python=3.10  # use a slightly older version
conda activate rna
conda install -c conda-forge r-base=4.4.3

# the following needs to be pre-installed for BiocManager to compile packages
conda install -c conda-forge r-curl=7.0.0  # DESeq2
conda install -c conda-forge r-xml=3.99  # sva
conda install -c conda-forge r-ggplot2=4.0.2  # clusterProfiler
```

Enter the R console to install R packages.

```R
install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("sva")  # sva includes ComBat-seq

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # human
BiocManager::install("org.Mm.eg.db")  # mouse
BiocManager::install("org.Rn.eg.db")  # rat
```

Troubleshoot (workaround) dependency mess for clusterProfiler: `ggtree` needs to be updated directly from the code base under active development.

```R
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
```

After R is installed, install rpy2 and other packages.

```bash
conda install -c conda-forge rpy2=3.6.6
pip install --no-cache-dir numpy==1.26.4 pandas==2.3.3 seaborn==0.13.2 scikit-learn==1.7.2
```

## GSEA

Go to https://www.gsea-msigdb.org/gsea/downloads.jsp (need to provide email) to download GSEA for Linux.

```bash
cd ~/opt
unzip GSEA_Linux_4.3.3.zip
rm GSEA_Linux_4.3.3.zip
```

In the `.bashrc` file add the following line:

```bash
export PATH=$PATH:$HOME/opt/GSEA_Linux_4.3.3
```
