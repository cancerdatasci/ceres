# CERES
Computational correction of copy-number effect in CRISPR-Cas9 essentiality screens

## Installation instructions

You will need several packages available on [Bioconductor](https://bioconductor.org) before installing CERES. To install these, run:

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "Rsamtools", 
            "GenomeInfoDb", "BSgenome", 
            "BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges"), type="source")
```

If the `devtools` package is not already installed, install from the R console:

```
install.packages("devtools")
```

To install CERES, either run:

```
devtools::install_github("cancerdatasci/ceres")
```

or clone the `ceres` repository, navigate to the local copy, and run from the R console: 

```
devtools::install("ceres")
```

Note that if C++11 support is not already enabled, you may need to run

```
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
```

prior to running the install command.

Preparing CERES inputs also requires the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) and [`samtools`](http://samtools.sourceforge.net) command line tools. For OSX users with `homebrew` installed on their machine, these can be installed from the command line:

```
brew tap homebrew/science
brew install bowtie
brew install samtools
```

## Run CERES on example data

Download [these zipped files](https://depmap.org/ceres/data/example_data.zip) from [depmap.org/ceres](https://depmap.org/ceres) and extract into a directory. (e.g. `./data/download`). If you haven't already fetched / built them yourself, you may also separately download the necessary bowtie indices [here](https://depmap.org/ceres/data/hg19_bowtie.tar) and place the unzipped files in the `bowtie_indexes` directory of the example folder.

The data in the example files are from screens of 33 cancer cell lines published in [Aguirre et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27260156) and 14 AML lines published in [Wang et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28162770).

Run the example script below, ensuring that the data_dir variable points to the directory with the data download.


```
library(ceres)

### Setup

# Edit this line to point to data directory
data_dir <- "./data/download"

cn_seg_file <- file.path(data_dir, "CCLE_copynumber_2013-12-03.seg.txt")
gene_annot_file <- file.path(data_dir, "CCDS.current.txt")

# Set bowtie index directory. Not needed if $BOWTIE_INDEXES environmental variable is set and includes hg19 index.
Sys.setenv(BOWTIE_INDEXES = file.path(data_dir, "bowtie_indexes"))


gecko_dep_file <- file.path(data_dir, "Gecko.gct")
gecko_rep_map <- file.path(data_dir, "Gecko_replicate_map.tsv")

wang_dep_file <- file.path(data_dir, "Wang2017.gct")
wang_rep_map <- file.path(data_dir, "Wang2017_replicate_map.tsv")



### Run CERES on Gecko data

gecko_inputs_dir <- file.path("./data/gecko_ceres_inputs", Sys.Date())

prepare_ceres_inputs(inputs_dir=gecko_inputs_dir,
                     dep_file=gecko_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=gecko_rep_map,
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

gecko_ceres <-
    wrap_ceres(sg_path=file.path(gecko_inputs_dir, "guide_sample_dep.Rds"),
               cn_path=file.path(gecko_inputs_dir, "locus_sample_cn.Rds"),
               guide_locus_path=file.path(gecko_inputs_dir, "guide_locus.Rds"),
               locus_gene_path=file.path(gecko_inputs_dir, "locus_gene.Rds"),
               replicate_map_path=file.path(gecko_inputs_dir, "replicate_map.Rds"),
               run_id="Gecko",
               params=list(lambda_g=0.68129207))

gecko_ceres_scaled <-
    scale_to_essentials(gecko_ceres$gene_essentiality_results$ge_fit)


### Run CERES on Wang2017 data

wang_inputs_dir <- file.path("./data/wang_ceres_inputs", Sys.Date())

prepare_ceres_inputs(inputs_dir=wang_inputs_dir,
                     dep_file=wang_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=wang_rep_map,
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

wang_ceres <-
    wrap_ceres(sg_path=file.path(wang_inputs_dir, "guide_sample_dep.Rds"),
               cn_path=file.path(wang_inputs_dir, "locus_sample_cn.Rds"),
               guide_locus_path=file.path(wang_inputs_dir, "guide_locus.Rds"),
               locus_gene_path=file.path(wang_inputs_dir, "locus_gene.Rds"),
               replicate_map_path=file.path(wang_inputs_dir, "replicate_map.Rds"),
               run_id="Wang2017",
               params=list(lambda_g=0.68129207))

wang_ceres_scaled <-
    scale_to_essentials(wang_ceres$gene_essentiality_results$ge_fit)
```
