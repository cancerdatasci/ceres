# CERES
Computational correction of copy-number effect in CRISPR-Cas9 essentiality screens

## Installation instructions

To install CERES, clone the `ceres` repository and begin an R session. If it is not already installed, make sure to install the package `devtools` by running:

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
devtools::install("ceres")
```

Preparing CERES inputs also requires the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) and [`samtools`](http://samtools.sourceforge.net) command line tools. For OSX users with `homebrew` installed on their machine, these can be installed from the command line:

```
brew tap homebrew/science
brew install bowtie
brew install samtools
```

## Run CERES on example data

Download [example data bundle](https://depmap.org/ceres/data/example_data.zip) from [depmap.org/ceres](http://depmap.org/ceres) and extract into a directory. (e.g. `./data/download`)

Example data are from screens of 33 cancer cell lines published in [Aguirre et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/27260156) and 14 AML lines published in [Wang et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28162770).

Run the example script below, ensuring that the data_dir variable points to the directory with the data download.


```
library(ceres)

data_dir <- "./data/download"

cn_seg_file <- file.path(data_dir, "CCLE_copynumber_2013-12-03.seg.txt")
gene_annot_file <- file.path(data_dir, "CCDS.current.txt")
Sys.setenv(BOWTIE_INDEXES = file.path(data_dir, "bowtie_indexes"))


gecko_dep_file <- file.path(data_dir, "Gecko.gct")
gecko_rep_map <- file.path(data_dir, "Gecko_replicate_map.tsv")

wang_dep_file <- file.path(data_dir, "Wang2017.gct")
wang_rep_map <- file.path(data_dir, "Wang2017_replicate_map.tsv")




gecko_inputs_dir <- file.path("./data/gecko_inputs", Sys.Date())

prepare_ceres_inputs(inputs_dir=gecko_inputs_dir,
                     dep_file=gecko_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=gecko_rep_map,
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

gecko_ceres <-
    wrap_ceres(sg_path=file.path(wang_inputs_dir, "guide_sample_dep.rds"),
               cn_path=file.path(wang_inputs_dir, "locus_sample_cn.rds"),
               guide_locus_path=file.path(wang_inputs_dir, "guide_locus.rds"),
               locus_gene_path=file.path(wang_inputs_dir, "locus_gene.rds"),
               replicate_map_path=file.path(wang_inputs_dir, "replicate_map.rds"),
               run_id="Wang2017",
               params=list(lambda_g=0.68129207))




wang_inputs_dir <- file.path("./data/wang_inputs", Sys.Date())

prepare_ceres_inputs(inputs_dir=wang_inputs_dir,
                     dep_file=wang_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=wang_rep_map,
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

wang_ceres <-
    wrap_ceres(sg_path=file.path(wang_inputs_dir, "guide_sample_dep.rds"),
               cn_path=file.path(wang_inputs_dir, "locus_sample_cn.rds"),
               guide_locus_path=file.path(wang_inputs_dir, "guide_locus.rds"),
               locus_gene_path=file.path(wang_inputs_dir, "locus_gene.rds"),
               replicate_map_path=file.path(wang_inputs_dir, "replicate_map.rds"),
               run_id="Wang2017",
               params=list(lambda_g=0.68129207))
```