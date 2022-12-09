### Configure `apt` repo

https://cran.r-project.org/bin/linux/ubuntu/

```bash
# update indices
sudo apt update

# install two helper packages we need
sudo apt install --no-install-recommends software-properties-common dirmngr

# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
```

### Install R

Due to some later compilation dependency, need to install `r-base-dev` (not `r-base`).

```bash
sudo apt install --no-install-recommends r-base-dev
```

### Install Packages

Enter the R console with `sudo`, as root access is needed for installation.

```bash
sudo R
```

In the R console:

```R
install.packages("BiocManager")

# upgrade bioconductor
BiocManager::install(version="3.16")

# use bioconductor to install deseq2
BiocManager::install("DESeq2")
```

## Troubleshoot

Previously installed older R packages (e.g. packages installed on R 3.6) can block the installation of DESeq2.
In such cases, completely clean up all R and its associated packages,
and then re-run the whole installation again.

```bash
# this only removes the r-base itself
sudo apt purge r-base

# remove all other r dependencies
sudo apt autoremove

# remove the site-library of R
# this is important, as remaining older packages can block installation
# even after newer version of R is installed 
sudo rm -r /usr/local/lib/R/site-library
```
