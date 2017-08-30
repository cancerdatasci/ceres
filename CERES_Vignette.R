library(ceresr)

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
