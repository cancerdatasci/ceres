# CERES
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

Preparing CERES inputs also requires the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) and [`samtools`](http://samtools.sourceforge.net) command line tools. For OSX users with `homebrew` installed on their machine, these can be installed from the command line:

```
brew tap homebrew/science
brew install bowtie
brew install samtools
```

## Run CERES on example data

Download example data bundle from [depmap.org/ceres](http://depmap.org/ceres) and extract into a directory. (e.g. `./data/download`)

Example data are from screens of 33 cancer cell lines published in [Aguirre et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/27260156) and 14 AML lines s[Wang et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28162770).

Run the example script, ensuring that the data_dir variable points to the directory with the data download.

```
source("./CERES_Vignette.R")
```