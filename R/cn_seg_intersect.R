#' Get a matrix of CN for CRISPR target sites
#' @param guide.gr GRanges object of sgRNA targets
#' @param seg.gr List of GRanges objects of CN segments
#' @param CN.column column to find CN data in seg.gr
#' @param guide.column guide.gr must contain a column with guide sequence
#' @param do_parallel run on multiple cores
#' @return a matrix of CN data, each row is a guide, each column is a sample
#' @importFrom plyr laply
#' @export
#'
intersect_guide_with_copy_number <- function (guide.gr, seg.gr,
                                              CN.column="CN", guide.column = "Guide",
                                              do_parallel = F) {
    stopifnot(is.list(seg.gr))
    single_cell_line <- length(names(seg.gr)) == 1
    laply(seg.gr, function(seg) {
        hits <- GenomicRanges::findOverlaps(guide.gr, seg) %>% as.data.frame %>%
            dplyr::distinct(queryHits, .keep_all = T)
        CN <- rep(NA, length(guide.gr))
        CN[hits$queryHits] <- as.data.frame(seg)[hits$subjectHits,
                                                 CN.column]
        return(CN)
    }, .parallel = do_parallel) %>% 
        {if(single_cell_line) as.matrix(.) else t(.)} %>%
        set_colnames(names(seg.gr)) %>%
        set_rownames(mcols(guide.gr)[, guide.column])
}


#' Get a matrix of CN for genes
#' @param gene.gr GRanges object of genes
#' @param seg.gr List of GRanges objects of CN segments
#' @param CN.column column to find CN data in seg.gr
#' @param gene.column gene.gr must contain a column with guide sequence
#' @param do_parallel run on multiple cores
#' @return a matrix of CN data, each row is a gene, each column is a sample
#' @importFrom plyr laply
#' @export
#'
intersect_gene_with_copy_number <- function (gene.gr, seg.gr,
                                             CN.column="CN", gene.column="Gene",
                                             do_parallel = F) {

    stopifnot(is.list(seg.gr))
    plyr::laply(seg.gr, function(seg) {
        hits <- findOverlaps(gene.gr, seg)
        seg$CN <- mcols(seg)[, CN.column]
        avg.cn <- pintersect(seg[subjectHits(hits)], gene.gr[queryHits(hits)]) %>%
            as.data.frame() %>%
            dplyr::mutate(Query = queryHits(hits),
                          SegCN = CN) %>%
            dplyr::group_by(Query) %>%
            dplyr::summarise(AvgCN = sum(SegCN * width, na.rm = T)/sum(width, na.rm = T))
        CN <- rep(NA, length(gene.gr))
        CN[avg.cn$Query] <- avg.cn$AvgCN
        return(CN)
    }, .parallel = do_parallel) %>% t() %>%
        magrittr::set_colnames(names(seg.gr)) %>%
        magrittr::set_rownames(mcols(gene.gr)[, gene.column])
}



#' Finds overlaps between target segments and data. Calculates result of a function across overlaps.
#' @param seg GRanges object to calculate data over
#' @param gr GRanges object with a data column
#' @param dat.column which mcols name to find data in \code{gr}
#' @param func what function to use (e.g.: \code{median}, \code{mean}, \code{length})
#' @param na.value value to use if segment has no overlaps with gr
#' @param do_parallel
#' @param \dots extra arguments passed to \code{func}
#' @return a vector with one value per segment
#' @importFrom plyr laply
#' @importFrom plyr llply
#' @export
#'
intersect_data_with_segment <- function(seg, gr, dat.column, func=median,
                                        na.value=NA, do_parallel=F, ...) {

    stopifnot(dat.column %in% colnames(mcols(gr)))

    hits <- findOverlaps(gr, seg, type="any")
    hits.list <- split(queryHits(hits), subjectHits(hits))

    dat.vec <- mcols(gr)[, dat.column]

    data.list <- llply(hits.list, function(h) {
        dat.vec[h]
    }, .parallel=do_parallel)

    result <- rep(na.value, length(seg))
    result[as.numeric(names(data.list))] <- laply(data.list, func,
                                                  .parallel=do_parallel, ...)

    return(result)

}
