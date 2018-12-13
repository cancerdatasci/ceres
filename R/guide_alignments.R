# fetch a the corresponding BSGenome object by the genome_id
getBSgenome <- function(genome_id) {
    if(genome_id == "hg19") {
        return (BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    } else if(genome_id == "hg38") {
        return (BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    } else {
        stop(paste0(genome_id, " is not a known genome id"))
    }
}

#' Get a PAM sequence
#' @importFrom BSgenome getSeq
#' @export
#'
fetchPAM <- function(chr, pos, strand, guide_length=20, genome_id="hg19") {
    PAM_start <- ifelse(strand=="+", pos+guide_length, pos-3)
    PAM_end <- ifelse(strand=="+", pos+guide_length+2, pos-1)

    PAM <- getSeq(getBSgenome(genome_id),
                            names=chr,
                            start=PAM_start,
                            end=PAM_end,
                            strand=strand,
                            as.character=TRUE)
    return(PAM)
}




#' Read BAM file report, filter with valid PAM sequence, report # of alignments
#' @param bam.file BAM file path
#' @param max.alns only consider guides with fewer than this many alignments
#' @param pam PAM sequence as regular expression
#' @param include.no.align include guides that do not map to genome
#' @param as.df return data.frame if \code{TRUE}, GRanges object if \code{FALSE}
#' @return data.frame or GRanges object of results
#' @importFrom Rsamtools scanBam
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom BSgenome getSeq
#' @export
#'
guideAlignments <- function(bam.file, max.alns=100, pam="[ACGTN]GG|GG[ACGTN]",
                            genome_id="hg19", chromosomes=paste0("chr", c(as.character(1:22), "X", "Y")),
                            include.no.align=F, as.df=T, guide_length=20) {

    bam <- scanBam(bam.file)
    
    print('in guide alignments')

    ### Bam file to data frame, count number of alignments
    bamdf <- dplyr::as_data_frame(bam[[1]][1:6]) %>%
        dplyr::mutate(Gene = str_extract(qname, "_.+$") %>%
                          str_sub(start=2),
               Guide = str_extract(qname, "^[A-Z]+"),
               Aligned = !is.na(rname)) %>%
        dplyr::filter(Aligned & rname %in% chromosomes) %>%
        dplyr::group_by(qname) %>%
        dplyr::mutate(N.alns = sum(Aligned)) %>%
        dplyr::ungroup()

    print(dim(bamdf))
    bamdf.noAlign <- dplyr::as_data_frame(bam[[1]][1:6]) %>%
        dplyr::mutate(Gene = str_extract(qname, "_.+$") %>%
                         str_sub(start=2),
               Guide = str_extract(qname, "^[A-Z]+"),
               Aligned = !is.na(rname)) %>%
        dplyr::filter(!Aligned) %>%
        dplyr::mutate(N.alns = 0)
    print(dim(bamdf.noAlign))

    ### Limit guides with too many alignments, annotate position from biomart data
    bamdf.filt <- bamdf %>% dplyr::filter(N.alns < max.alns)

    bamdf.pam <- dplyr::mutate(bamdf.filt,
                               Cut.Pos = ifelse(strand=="+", pos+guide_length-4, pos+2),
                               PAM.start = ifelse(strand=="+", pos+guide_length, pos-3),
                               PAM.end = ifelse(strand=="+", pos+guide_length+2, pos-1))

    bamdf.pam$PAM <- getSeq(getBSgenome(genome_id),
                            names=bamdf.pam$rname,
                            start=bamdf.pam$PAM.start,
                            end=bamdf.pam$PAM.end,
                            strand=bamdf.pam$strand,
                            as.character=TRUE)

    bamdf.pam <- bamdf.pam %>%
        dplyr::group_by(qname) %>%
        dplyr::mutate(N.alns = sum(str_detect(PAM, pam))) %>%
        dplyr::ungroup()

    print('bamdf.pam')
    print(bamdf.pam[1:10, ])
    if (any(bamdf.pam$N.alns == 0)) {
        print("At least one guide aligned without corresponding PAM")
        bamdf.noAlign <-
            dplyr::bind_rows(bamdf.noAlign,
                             bamdf.pam %>%
                                 dplyr::filter(N.alns == 0) %>%
                                 dplyr::distinct(qname, Gene, Guide, .keep_all=T) %>%
                                 dplyr::select(qname, Gene, Guide, N.alns) %>%
                                 dplyr::mutate(Aligned = FALSE))
    }
    bamdf.pam <- bamdf.pam %>% dplyr::filter(stringr::str_detect(PAM, pam))


    if (include.no.align) {
        bamdf.final <- dplyr::bind_rows(bamdf.pam, bamdf.noAlign)
    } else {
        bamdf.final <- bamdf.pam
    }
    
    print(dim(bamdf.final))

    if (as.df) {
        return(bamdf.final)
    } else {
        genomeinfo <- Seqinfo(genome=genome_id)[chromosomes]

        return(makeGRangesFromDataFrame(bamdf.final,
                                                       keep.extra.columns=T,
                                                       seqinfo=genomeinfo))
    }
}
