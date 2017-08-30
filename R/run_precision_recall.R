#' Run precision recall assessment using essential and nonessential genes
#'
#' @param ge_fit Matrix of gene effects by cell line. Expects row names are gene symbols and column names are cell line IDs.
#' 
#' @return List of PR results.
#' 
#' @export run_precision_recall
#' 
run_precision_recall <- function(ge_fit){
  
  
  test_gene_indices <- which(row.names(ge_fit) %in% 
                               union(ceres::hart_essentials[["Gene"]], 
                                     ceres::hart_nonessentials[["Gene"]]))
  
  screen_results <- lapply(colnames(ge_fit), function(x){
                        
                        true_values <- ifelse(row.names(ge_fit[test_gene_indices, x, drop=F]) %in% ceres::hart_essentials[["Gene"]], 1, 0)
                        ceres::evaluate_prediction_auc(-1*ge_fit[test_gene_indices, x, drop=F], true_values)
                        
                      }) %>%
                      set_names(colnames(ge_fit))
  
  screen_plots <- screen_results %>%
                    lapply(function(x){
                      return(x[["performance_plots"]])
                    }) %>%
                    set_names(names(screen_results))
  screen_results_df <- screen_results %>%
                          lapply(function(x){
                            return(x[["performance_df"]])
                          }) %>%
                          set_names(names(screen_results)) %>%
                          data.table::rbindlist(idcol="cell_line")
  
  return(list(screen_plots=screen_plots, screen_results_df=screen_results_df))
  
}