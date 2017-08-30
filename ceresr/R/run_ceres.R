#' @useDynLib ceresr
#' @importFrom Rcpp sourceCpp
#' @import e1071
#' @import readr
#' @import stringr
#' @import magrittr
#' @import GenomicRanges
NULL

#' CERES main routine
#'
#' @param sg_data matrix with sgRNA sequences as row names, samples (i.e. replicates of cell line screens) as column names.
#' @param cn_data matrix with genomic cut sites as row names, cell lines as column names.
#' @param guide_locus data.frame with column names `Guide` and `Locus`. Values must match those corresponding to row names of sg_data and cn_data.
#' @param locus_gene data.frame with column names `Locus` and `Gene`.
#' @param replicate_map data.frame with column names `Replicate` and `CellLine`
#' @param params list of run parameters
#'
#' @return A list of numeric vectors and matrices containing the results of the CERES fit.
#' 
#' @importFrom plyr .
#' @export run_ceres
#'
run_ceres <- function(sg_data, cn_data, guide_locus, locus_gene,
                      replicate_map, params){


  # Set path to log file
  logfile_path <- file.path("log")

  # Checks on data
  sgRNA_data_is_matrix <- "matrix" %in% class(sg_data)
  copy_number_data_is_matrix <- "matrix" %in% class(cn_data)
  guide_locus_map_is_data_frame <- "data.frame" %in% class(guide_locus)
  locus_gene_map_is_data_frame <- "data.frame" %in% class(locus_gene)
  replicate_map_is_data_frame <- "data.frame" %in% class(replicate_map)
  stopifnot(sgRNA_data_is_matrix,
            copy_number_data_is_matrix,
            guide_locus_map_is_data_frame,
            locus_gene_map_is_data_frame,
            replicate_map_is_data_frame)

  #
  cn_data <- cn_data[, unique(replicate_map$CellLine), drop=F]
  sg_data <- sg_data[, unique(replicate_map$Replicate), drop=F]

  # Make sure that sgRNA and CN matrix columns are in the right order and contain same number of cell lines
  replicates_match_columns <- all(replicate_map$Replicate == colnames(sg_data))
  stopifnot(replicates_match_columns)

  # Shrink data, strip data of any columns/rows that are all NA-valued
  non_na_rows <- apply(sg_data, 1, function(x){return(any(!is.na(x)))})
  sg_data <- sg_data[non_na_rows, , drop=F]
  non_na_rows <- apply(cn_data, 1, function(x){return(any(!is.na(x)))})
  cn_data <- cn_data[non_na_rows, , drop=F]

  # Mean-impute missing data in the copy number / sgRNA matrices
  sg_data <- apply(sg_data, 1, function(x){x[is.na(x)] <- mean(x, na.rm=T); return(x)}) %>%
              matrix(., ncol=nrow(sg_data)) %>%
              t() %>%
              magrittr::set_colnames(colnames(sg_data)) %>%
              magrittr::set_rownames(row.names(sg_data))

  cn_data <- apply(cn_data, 1, function(x){x[is.na(x)] <- mean(x, na.rm=T); return(x)}) %>%
              matrix(., ncol=nrow(cn_data)) %>%
              t() %>%
              magrittr::set_colnames(colnames(cn_data)) %>%
              magrittr::set_rownames(row.names(cn_data))

  # Select loci/guides that will be used for fit
  loci_to_use <- guide_locus %>%
                  dplyr::filter(Guide %in% row.names(sg_data)) %>%
                  magrittr::extract2("Locus") %>%
                  unique()

  cn_data <- cn_data[row.names(cn_data) %in% loci_to_use, , drop=F]

  locus_gene %<>% dplyr::filter(Locus %in% row.names(cn_data))
  guide_locus %<>% dplyr::filter(Guide %in% row.names(sg_data))

  # Get numbers of each dimension
  n_genes <- length(unique(locus_gene$Gene))
  n_lines <- ncol(cn_data)
  n_guides <- nrow(sg_data)
  n_loci <- nrow(cn_data)

  # Build gene effect and mapping matrices
  genes <- unique(locus_gene$Gene)

  locus_index <- data.frame(Locus=row.names(cn_data), Locus_Index=(0:(nrow(cn_data)-1)))
  gene_index <- data.frame(Gene=genes, Gene_Index=(0:(length(genes)-1)))
  guide_index <- data.frame(Guide=row.names(sg_data), Guide_Index=(0:(nrow(sg_data)-1)))

  locus_gene %<>%
    dplyr::left_join(locus_index) %>%
    dplyr::left_join(gene_index) %>%
    dplyr::distinct()

  guide_locus %<>%
    dplyr::left_join(locus_index) %>%
    dplyr::left_join(guide_index) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(Locus_Index)) # filter out negative control guides from guide_locus mapping

  guide_locus_triplets <- data.frame(Guide_Index=as.integer(guide_locus$Guide_Index), Locus_Index=as.integer(guide_locus$Locus_Index), Prob=as.integer(guide_locus$Value))
  locus_gene_triplets <- data.frame(Locus_Index=as.integer(locus_gene$Locus_Index), Gene_Index=as.integer(locus_gene$Gene_Index), Map=as.integer(locus_gene$Value))

  guide_locus_triplets <- unique(as.matrix(guide_locus_triplets))
  locus_gene_triplets <- unique(as.matrix(locus_gene_triplets))

  ge_data <- matrix(0, nrow=n_genes, ncol=n_lines) %>%
        as.data.frame() %>%
        dplyr::mutate(Gene_Index=0:(nrow(.)-1)) %>%
        dplyr::left_join(dplyr::select(locus_gene, Gene_Index, Gene) %>% unique()) %>%
        magrittr::set_rownames(., magrittr::extract2(., "Gene")) %>%
        dplyr::select(-Gene, -Gene_Index) %>%
        magrittr::set_colnames(colnames(cn_data)) %>%
        as.matrix()

  ce_data <- as.matrix(rep(0, n_lines)) %>% magrittr::set_names(colnames(cn_data))

  # Read replicate mapping and cell line weighting files
  ColCldf <- data.frame(CellLine=colnames(cn_data), ClIndx=0:(ncol(cn_data)-1))
  ColCl <- replicate_map %>%
              dplyr::left_join(ColCldf) %>%
              magrittr::extract2("ClIndx") %>%
              as.matrix()

  # Cluster cuts to determine CN bins for linear fit
  cat("\nClustering copy number regions...\n")
  quantile_mat <- cn_data %>%
        data.frame(Locus=row.names(cn_data), .) %>%
        reshape2::melt(id.vars="Locus", variable.name="CellLine", value.name="CN") %>%
        dplyr::left_join(guide_locus) %>%
        dplyr::group_by(Guide, CellLine) %>%
        dplyr::summarise(Cuts=sum(CN, na.rm=T)) %>%
        dplyr::ungroup() %>%
        plyr::daply(.(CellLine), function(x){
          cat(paste0("\rClustering cuts for ", unique(x$CellLine), "...                          "))
          return(ceresr::cluster_cn(x$Cuts, n_segments=unique(params[["n_segments"]]) %>% as.integer()))

        }) %>%
        matrix(nrow=n_lines)

  # Final checks
  copy_number_matrix_has_column_per_cell_line <- ncol(cn_data) == n_lines
  gene_effect_matrix_has_column_per_cell_line <- ncol(ge_data) == n_lines
  sgRNA_matrix_has_no_NAs <- !any(is.na(sg_data))
  copy_number_matrix_has_no_NAs <- !any(is.na(cn_data))
  guide_locus_triplets_has_no_NAs <- !any(is.na(guide_locus_triplets))
  locus_gene_triplets_has_no_NAs <- !any(is.na(locus_gene_triplets))
  copy_number_bin_matrix_has_proper_dimension <- (ncol(quantile_mat) == 2*(unique(params[["n_segments"]]) %>% as.integer())) & (nrow(quantile_mat) == n_lines)
  stopifnot(copy_number_matrix_has_column_per_cell_line,
                gene_effect_matrix_has_column_per_cell_line,
                sgRNA_matrix_has_no_NAs,
                copy_number_matrix_has_no_NAs,
                copy_number_bin_matrix_has_proper_dimension)


  # Extract run parameters
  lambda_g <- params[["lambda_g"]] %>% as.numeric()
  lambda_o <- params[["lambda_o"]] %>% as.numeric()
  lambda_s <- params[["lambda_s"]] %>% as.numeric()
  n_segments <- params[["n_segments"]] %>% as.integer()
  validation_set <- params[["validation_set"]] %>% as.integer()
  run_name <- params[["run_name"]] %>% as.character()

  run_mat <- cbind(lambda_g, lambda_o, lambda_s,
                   n_segments, validation_set) %>%
              magrittr::set_rownames(run_name) %>%
              unique()

  # Fit CERES to data
  if(nrow(run_mat) == 1){

    l_g <- run_mat[1, "lambda_g"] %>% as.numeric()
    l_o <- run_mat[1, "lambda_o"] %>% as.numeric()
    l_s <- run_mat[1, "lambda_s"] %>% as.numeric()
    n_s <- run_mat[1, "n_segments"] %>% as.integer()
    v_s <- run_mat[1, "validation_set"] %>% as.integer()
    r_n <- row.names(run_mat)

    res <- fit_ceres(rD = sg_data, rQ = guide_locus_triplets, rM = locus_gene_triplets, rColCl = ColCl,
                             rG = ge_data, rC = cn_data, rTox = ce_data, rQuantileMat = quantile_mat,
                             LAMBDA_G = l_g, LAMBDA_O = l_o, LAMBDA_Smooth = l_s,
                             NSEGMENTS = n_s, MAKE_VALIDATION_SET = v_s,
                             log_file_suffix = r_n, log_file_dir = logfile_path)
    #
    collated_results <- ceresr::collate_ceres_results(res, sg_data, cn_data)
    return(collated_results)

  } else {

    optimization_results <- lapply(row.names(run_mat), function(run){

      l_g <- run_mat[run, "lambda_g"] %>% as.numeric()
      l_o <- run_mat[run, "lambda_o"] %>% as.numeric()
      l_s <- run_mat[run, "lambda_s"] %>% as.numeric()
      n_s <- run_mat[run, "n_segments"] %>% as.integer()
      v_s <- run_mat[run, "validation_set"] %>% as.integer()
      r_n <- run

      res <- fit_ceres(rD = sg_data, rQ = guide_locus_triplets, rM = locus_gene_triplets, rColCl = ColCl,
                               rG = ge_data, rC = cn_data, rTox = ce_data, rQuantileMat = quantile_mat,
                               LAMBDA_G = l_g, LAMBDA_O = l_o, LAMBDA_Smooth = l_s,
                               NSEGMENTS = n_s, MAKE_VALIDATION_SET = v_s,
                               log_file_suffix = r_n, log_file_dir = logfile_path)

      run_df <- data.frame(run_name=run,
                           lambda_g=l_g,
                           lambda_o=l_o,
                           lambda_s=l_s,
                           n_segments=n_s,
                           validation_set=v_s,
                           train_error=res[["train_error"]],
                           val_error=res[["val_error"]])

    }) %>%
      data.table::rbindlist()

    optimal_run <- optimization_results %>%
                    dplyr::filter(val_error == min(val_error))

    l_g <- optimal_run[["lambda_g"]] %>% as.numeric()
    l_o <- optimal_run[["lambda_o"]] %>% as.numeric()
    l_s <- optimal_run[["lambda_s"]] %>% as.numeric()
    n_s <- optimal_run[["n_segments"]] %>% as.integer()
    v_s <- 0
    r_n <- optimal_run[["run_name"]] %>% as.character()

    res <- fit_ceres(rD = sg_data, rQ = guide_locus_triplets, rM = locus_gene_triplets, rColCl = ColCl,
                             rG = ge_data, rC = cn_data, rTox = ce_data, rQuantileMat = quantile_mat,
                             LAMBDA_G = l_g, LAMBDA_O = l_o, LAMBDA_Smooth = l_s,
                             NSEGMENTS = n_s, MAKE_VALIDATION_SET = v_s,
                             log_file_suffix = r_n, log_file_dir = logfile_path)
    #
    collated_results <- ceresr::collate_ceres_results(res, sg_data, cn_data)

    collated_results[["optimization_results"]] <- optimization_results

    return(collated_results)

  }
}
