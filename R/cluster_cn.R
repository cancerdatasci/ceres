#' Hierarchically cluster copy number data
#'
#' @description Hierarchically cluster copy number data to create knot points for piecewise linear cutting effect fit
#'
#' @param cn_vec vector of copy number values for each genomic locus in a particular cell line
#' @param min_dpts sets the minimum number of data points within a cluster leaf
#' @param cluster_method the agglomeration method to be used. See \code{\link[stats]{hclust}} for available methods
#' @param n_segments the number of segments in the piecewise linear cutting effect fit
#'
#' @return A vector of knot points of the form c(start_1, end_1, start_2, end_2, start_3, end_3, ..., start_n, end_n)
#' @examples
#'
#' @export cluster_cn
#'
cluster_cn <- function(cn_vec, min_dpts=25, cluster_method="average", n_segments=25){

  res <- cn_vec %>%
    data.frame(cn=.) %>%
    dplyr::arrange(cn) %>%
    dplyr::filter(!is.na(cn)) %>%
    dplyr::mutate(quantile=cut(1:dplyr::n(), breaks=floor(dplyr::n()/min_dpts), labels=F, include.lowest = T)) %>%
    dplyr::group_by(quantile) %>%
    dplyr::summarise(max_cn=max(cn)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster=dist(max_cn) %>%
                            hclust(method=cluster_method) %>%
                            cutree(k=n_segments)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(interval_end=max(max_cn)) %>%
    dplyr::ungroup() %>%
    magrittr::extract2("interval_end") %>%
    rep(., 2) %>%
    sort() %>%
    magrittr::extract(1:(length(.)-1)) %>%
    c(0, .)

  res[length(res)] <- res[length(res)] + 1

  return(res)

}
