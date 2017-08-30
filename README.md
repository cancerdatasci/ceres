# ceres
Computational correction of copy-number effect in CRISPR-Cas9 essentiality screens

## Installation instructions

To install CERES, clone the `ceresr` repository and begin an R session. If it is not already installed, make sure to install the package `devtools` by running:

```
install.packages("devtools")
```
You will need several packages available on [Bioconductor](https://bioconductor.org) before installing CERES. To install these, run:

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "Rsamtools", 
            "GenomeInfoDb", "BSgenome", 
            "BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges"))
```

Then, navigate to your local copy of this repository and run:

```
devtools::install("ceresr")
```

Preparing CERES inputs also requires the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) and [`samtools`](http://samtools.sourceforge.net) command line tools. OSX users with `homebrew` installed, these can be installed with the command:

```
brew install bowtie
brew install samtools
```
