#' Collate CERES results
#'
#' @param ceres_results List returned by \code{fit_ceres}
#' @param sg_data Matrix with sgRNA sequences as row names, samples (i.e. replicates of cell line screens) as column names
#' @param cn_data Matrix with genomic cut sites as row names, cell lines as column names
#' 
#' @return A list of objects returned by the CERES fit
#' @examples
#' 
#' @export collate_ceres_results
#' 
collate_ceres_results <- function(ceres_results, sg_data, cn_data){
    
  sg_fit <- ceres_results[["sg_fit"]] %>%
              set_rownames(row.names(sg_data)) %>%
              set_colnames(colnames(cn_data))
  
  ge_fit <- ceres_results[["ge_fit"]]
  
  efficacy_fit <- ceres_results[["efficacy_fit"]] %>%
                    as.data.frame() %>%
                    set_colnames(c("guide_index", "locus_index", "efficacy")) %>%
                    dplyr::select(-locus_index) %>%
                    dplyr::distinct() %>%
                    dplyr::left_join(data.frame(sgRNA=row.names(sg_data),
                                                guide_index=0:(nrow(sg_data)-1))) %>%
                    dplyr::select(sgRNA, efficacy) %>%
                    dplyr::distinct()
  
  offset_fit <- ceres_results[["offset_fit"]] %>%
                  data.frame(sgRNA=row.names(sg_data),
                             offset=.)
  
  quantiles <- ceres_results[["quantiles"]] %>%
                set_rownames(colnames(cn_data)) %>%
                set_colnames(paste0("segment", rep(1:(ncol(.)/2), each=2), 
                                              "_", rep(c("start", "end"), (ncol(.)/2)))) %>%
                data.frame(cell_line=row.names(.), .) %>%
                reshape2::melt(id.vars="cell_line", variable.name="quantile", value.name="cn") %>%
                tidyr::separate("quantile", c("segment", "position"), sep="_") %>%
                dplyr::distinct(cell_line, cn, .keep_all=T) %>%
                dplyr::arrange(cell_line)
  
  ce_slopes <- ceres_results[["ce_slopes"]] %>%
                set_rownames(colnames(cn_data)) %>%
                set_colnames(paste0("segment", 1:ncol(.))) %>%
                data.frame(cell_line=row.names(.), .) %>%
                reshape2::melt(id.vars="cell_line", variable.name="segment", value.name="slope")
  
  ce_fit <- ceres_results[["ce_fit"]] %>%
              set_rownames(row.names(sg_data)) %>%
              set_colnames(colnames(cn_data))
  
  intercept_fit <- ceres_results[["intercept_fit"]] %>%
                    data.frame(cell_line=colnames(cn_data), intercept=.)
                
  effective_cuts <- ceres_results[["effective_cuts"]] %>%
                      set_rownames(row.names(sg_data)) %>%
                      set_colnames(colnames(cn_data))
  
  val_set <- ceres_results[["val_set"]] %>%
              data.frame() %>%
              set_colnames(c("sgRNA", "sample_id", "in_set")) %>%
              dplyr::filter(in_set == 1) %>%
              dplyr::select(sgRNA, sample_id) %>%
              dplyr::mutate(sgRNA=row.names(sg_data)[sgRNA+1],
                            sample_id=colnames(sg_data)[sample_id+1])
  
  sgrna_df <- dplyr::left_join(offset_fit, efficacy_fit)
  
  ce_df <- dplyr::left_join(quantiles, ce_slopes) %>%
            dplyr::left_join(intercept_fit)
  
  train_error <- ceres_results[["train_error"]]
  val_error <- ceres_results[["val_error"]]
  
  return(list(cutting_effect_results=list(ce_fit=ce_fit,
                                          effective_cuts=effective_cuts,
                                          ce_df=ce_df),
              gene_essentiality_results=list(ge_fit=ge_fit),
              sgrna_results=list(sgrna_df=sgrna_df,
                                 sg_fit=sg_fit),
              fit_results=list(train_error=train_error,
                               val_error=val_error,
                               val_set=val_set)))
  
}