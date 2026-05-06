### Install ASCAT with rocker/rstudio
```
sudo apt-get update
sudo apt-get install -y ruby-full zlib1g-dev vim htop libxtst6 texlive-latex-base texlive-latex-recommended texlive-fonts-recommended texlive-latex-extra libuv1-dev liblzma-dev libbz2-dev libcurl4-openssl-dev

install.packages("BiocManager")
install.packages("devtools")
install.packages("remotes")
BiocManager::install(c('GenomicRanges','IRanges'))
devtools::install_github('VanLoo-lab/ascat/ASCAT')
BiocManager::install(c("GenomeInfoDb", "GenomicRanges"))
```
