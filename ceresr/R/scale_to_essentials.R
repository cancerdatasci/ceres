#' Normalize CERES gene effects according to essential and nonessential genes
#'
#' @param ge_fit Matrix of gene effects by cell line. Expects row names are gene symbols and column names are cell line IDs.
#' 
#' @return Matrix of scaled gene effects where the median of essential / nonessential genes are -1 and 0, respectively, for all cell lines.
#' 
#' @export scale_to_essentials
#' 
scale_to_essentials <- function(ge_fit){
  
  
  essential_indices <- which(row.names(ge_fit) %in% ceresr::hart_essentials[["Gene"]])
  nonessential_indices <- which(row.names(ge_fit) %in% ceresr::hart_nonessentials[["Gene"]])
  
  scaled_ge_fit <- ge_fit %>%
                    apply(2, function(x){
                      
                      (x - median(x[nonessential_indices], na.rm=T)) %>%
                        divide_by(median(x[nonessential_indices], na.rm=T) - median(x[essential_indices], na.rm=T))
                      
                    })
  
  return(scaled_ge_fit)
  
}