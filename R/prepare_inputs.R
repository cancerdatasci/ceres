#' Map guide to locus
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @export
#'
map_guide_to_locus <- function(guides,
                               genome_id="hg19",
                               bowtie_exe="bowtie",
                               samtools_exe="samtools",
                               temp_dir=tempdir(),
                               write_rds_output_path=NULL, guide_length=20) {

    guides_fa <- file.path(temp_dir, "guides.fa")
    guides_sam <- file.path(temp_dir, "guides.sam")
    guides_bam <- file.path(temp_dir, "guides.bam")

    guides %>% unique %>%
        set_names(.,.) %>%
        Biostrings::DNAStringSet() %>%
        Biostrings::writeXStringSet(guides_fa)

    bowtie_cmd <- paste(bowtie_exe, "-t -p 4 -a -v 0 -f -S", genome_id,
                        guides_fa, guides_sam)
    system(bowtie_cmd)

    samtools_cmd <- paste(samtools_exe, "view -bS -o",
                          guides_bam, guides_sam)
    system(samtools_cmd)

    alns <- guideAlignments(guides_bam, max.alns=100,
                            include.no.align=T, as.df=T, guide_length=guide_length)

    if(!is.null(write_rds_output_path)){
        cat(paste('Writing the mapping of sgRNAs to the genome in', write_rds_output_path, 'csv file'))
        saveRDS(alns, write_rds_output_path)
    }

    return(alns)
}

#' Intersect locus with copy number segment
#' @importFrom GenomeInfoDb Seqinfo
#' @export
#'
intersect_locus_with_cn_seg <-
    function(cn_seg, guide_alns,
             cell_lines=NULL,
             genomeinfo=Seqinfo(genome="hg19"),
             chromosomes=paste0("chr", c(as.character(1:22),"X","Y")),
             do_parallel=F) {

        if (!is.null(cell_lines)) {
            cn_seg <- cn_seg[names(cn_seg) %in% cell_lines]
        }

        cn_seg_gr <- plyr::llply(cn_seg,
                                 makeGRangesFromDataFrame,
                                 seqinfo=genomeinfo, keep.extra.columns=T)

        guide_alns_gr <- guide_alns %>%
            dplyr::filter(!is.na(rname)) %>%
            dplyr::mutate(Chr = rname,
                   Start = Cut.Pos,
                   End = Cut.Pos,
                   AlnID = str_c(Guide, Chr, Start, strand, sep="_")) %>%
            dplyr::distinct(Guide, Chr, Start, End, AlnID) %>%
            dplyr::filter(Chr %in% chromosomes) %>%
            makeGRangesFromDataFrame(seqinfo=genomeinfo, keep.extra.columns=T)

        guide_no_alns <- guide_alns %>%
            dplyr::filter(is.na(rname)) %>%
            dplyr::distinct(Guide, .keep_all=T)

        guide_cn <-
            intersect_guide_with_copy_number(guide_alns_gr, cn_seg_gr,
                                             CN.column="CN",
                                             guide.column="AlnID",
                                             do_parallel=do_parallel) %>%
            rbind(matrix(0, dimnames=list(guide_no_alns$Guide, colnames(.)),
                         nrow=nrow(guide_no_alns), ncol=ncol(.)))
    }

#' @importFrom plyr "."
load_cn_seg_file <-
    function(cn_seg_file,
             chromosomes=paste0("chr", c(as.character(1:22),"X", "Y"))) {
        read_tsv(cn_seg_file,
                 col_types="ccddid") %>%
            set_colnames(c("CellLine", "Chr", "Start", "End",
                           "Num_Probes", "CN")) %>%
            dplyr::mutate(Chr = ifelse(str_detect(Chr, "^chr"), Chr, str_c("chr", Chr)),
                   Start = as.integer(Start),
                   End = as.integer(End)) %>%
            dplyr::filter(Chr %in% chromosomes) %>%
            dplyr::group_by(CellLine) %>%
            dplyr::mutate(CN = if(any(CN < 0)){2 * 2^CN}else{CN}) %>%
            dplyr::ungroup() %>%
            plyr::dlply(.(CellLine))
    }



get_gene_annotations <- function(genes, guide_alns,
                                 chromosomes=paste0("chr",c(as.character(1:22), "X", "Y")),
                                 genomeinfo=Seqinfo(genome="hg19")) {


    gene_annot_grs <- genes %>%
        makeGRangesFromDataFrame(seqinfo=genomeinfo, keep.extra.columns=T)


    guide_aln_grs <- guide_alns %>%
        dplyr::select(Guide, Chr=rname, Start=Cut.Pos, strand) %>%
        dplyr::mutate(End = Start) %>%
        dplyr::filter(Chr %in% chromosomes) %>%
        dplyr::distinct() %>%
        makeGRangesFromDataFrame(seqinfo=genomeinfo, keep.extra.columns=T)

    hits <- findOverlaps(guide_aln_grs, gene_annot_grs, ignore.strand=T) %>%
        as.data.frame

    gene_df <- hits %>%
        dplyr::transmute(Guide = guide_aln_grs$Guide[queryHits],
                         Chr = seqnames(guide_aln_grs)[queryHits] %>% as.character(),
                         Cut.Pos = start(guide_aln_grs)[queryHits] %>% as.integer(),
                         Strand = strand(guide_aln_grs)[queryHits] %>% as.character(),
                         Gene = gene_annot_grs$gene[subjectHits],
                         GeneID = gene_annot_grs$gene_id[subjectHits],
                         CDS_Strand = strand(gene_annot_grs)[subjectHits] %>% as.character(),
                         CDS_Start = start(gene_annot_grs)[subjectHits] %>% as.integer(),
                         CDS_End = end(gene_annot_grs)[subjectHits] %>% as.integer()) %>%
        dplyr::distinct()
}

load_ccds_genes <- function(ccds_file,
                            chromosomes=paste0("chr", c(as.character(1:22),"X","Y"))) {
    ccds <- read_tsv(ccds_file,
                     col_types=cols("#chromosome" = col_character(),
                                    "cds_from" = col_integer(),
                                    "cds_to" = col_integer())) %>%
        dplyr::rename(chromosome=`#chromosome`) %>%
        dplyr::mutate(chromosome = str_c("chr", chromosome)) %>%
        dplyr::filter(ccds_status %in% c("Public", "Reviewed, update pending", "Under review, update"),
                      chromosome %in% chromosomes,
                      !is.na(cds_from), !is.na(cds_to))

    ccds_exon <-  ccds %>%
        dplyr::mutate(cds_interval = str_replace_all(cds_locations, "[\\[\\]]", "") %>%
                   str_split("\\s*,\\s*")) %>%
        tidyr::unnest(cds_interval) %>%
        dplyr::group_by(gene, gene_id, cds_locations) %>%
        dplyr::mutate(exon_code = ifelse(cds_strand=="+", 1:n(), n():1)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(cds_start = str_extract(cds_interval, "^[0-9]+") %>% as.integer,
               cds_end = str_extract(cds_interval, "[0-9]+$") %>% as.integer) %>%
        dplyr::select(gene, gene_id, chromosome, start=cds_start, end=cds_end, strand=cds_strand,
               gene_start=cds_from, gene_end=cds_to, exon_code)
}

#' Generation of the guides from a gct dep_file
#' @param dep_file file path of guide-level dependency data
#' @param write_rds_output_path Optional: Will write guide_dep into a rds file in top of returning the matrix
#'
#' @return Matrix with a description attribute
#' @export
#'
generate_guides <- function(dep_file, write_rds_output_path=NULL){
    guide_dep <- read.gct(dep_file) %>%
        set_rownames(str_extract(rownames(.), "^[ACGT]+")) %>%
        {.[unique(rownames(.)),]} %>%
        remove.rows.all.nas()

    if(!is.null(write_rds_output_path)){
        cat(paste('Writing the generated guides in', write_rds_output_path, 'Rds file'))
        saveRDS(guide_dep, write_rds_output_path)
    }

    return(guide_dep)
}

#' CERES main routine
#'
#' @param inputs_dir directory path to write CERES inputs
#' @param dep_file file path of guide-level dependency data. !!Not necessary if you have generate your guides with generate_guides
#' @param pre_generated_guides Optional: Matrix generated by the generated_guides(dep_file) function. Will fasten this function (since we skip the reading of gct)
#' @param cn_seg file path of segmented copy number data
#' @param gene_annot_file file path of gene annotation
#' @param rep_map file path of replicate map
#' @param guide_alns_file Optional: file path of the guide mapped (use map_guide_to_locus to generate)
#'
#' @return Returns invisibly. Only called for its effects.
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @export
#'
prepare_ceres_inputs <-
    function(inputs_dir,
             dep_file,
             pre_generated_guides_file=NULL,
             cn_seg_file,
             gene_annot_file,
             rep_map_file,
             genome_id="hg19",
             chromosomes=paste0("chr",c(as.character(1:22), "X", "Y")),
             dep_normalize="zmad",
             bowtie_exe="bowtie",
             samtools_exe="samtools",
             do_parallel=F,
             guide_alns_file=NULL) {

        genomeinfo <- Seqinfo(genome=genome_id)[chromosomes]

        dir.create(inputs_dir, recursive=T, showWarnings=F)

        cat("loading dependency data...\n\n")

        guide_dep <- NULL
        if(is.null(pre_generated_guides_file)){
            guide_dep <- generate_guides(dep_file)
        }
        else{

            guide_dep <- readRDS(pre_generated_guides_file)
        }

        guide_length <- nchar(rownames(guide_dep)[2])

        rep_map <- readr::read_tsv(rep_map_file)


        cat("loading copy number data...\n\n")
        cn_seg <- load_cn_seg_file(cn_seg_file, chromosomes=chromosomes)

        guide_alns <- NULL
        if(is.null(guide_alns_file)){
            cat("mapping sgRNAs to the genome...\n\n")
            guide_alns <- map_guide_to_locus(rownames(guide_dep), genome_id=genome_id,
                                         bowtie_exe=bowtie_exe, samtools_exe=samtools_exe,
                                         guide_length=guide_length)
        }
        else{
            cat("reading sgRNAs mapped to the genome...\n\n")
            guide_alns <- readRDS(guide_alns_file)
        }

        cell_lines <- rep_map %>%
            dplyr::filter(Replicate %in% colnames(guide_dep)) %$%
            unique(CellLine)

        cat("getting copy number data per locus...\n\n")
        guide_cn_mat <- intersect_locus_with_cn_seg(cn_seg, guide_alns,
                                                    cell_lines=cell_lines,
                                                    genomeinfo=genomeinfo,
                                                    chromosomes=chromosomes,
                                                    do_parallel=do_parallel)

        locus_cn <- guide_cn_mat[str_detect(rownames(guide_cn_mat), "chr"), , drop=F] %>%
            set_rownames(rownames(.) %>% str_extract("chr.+$")) %>%
            {.[rownames(.),]} %>%
            remove.rows.all.nas()

        non_targeting_cn <- guide_cn_mat[!str_detect(rownames(guide_cn_mat), "chr"), , drop=F]


        guide_locus_df <- guide_alns %>%
            dplyr::transmute(Guide,
                      Locus = str_c(rname, Cut.Pos, strand, sep="_")) %>%
            dplyr::distinct()

        cat("mapping loci to gene coding regions...\n\n")

        ccds <- load_ccds_genes(ccds_file=gene_annot_file,
                                chromosomes=chromosomes)

        gene_df <- get_gene_annotations(ccds, guide_alns,
                                        genomeinfo=genomeinfo,
                                        chromosomes=chromosomes)

        locus_gene_df <- gene_df %>%
            dplyr::transmute(Locus = str_c(Chr, Cut.Pos, Strand, sep="_"),
                             Gene) %>%
            dplyr::distinct()


        cat("normalizing data...\n\n")

        rep_map <- read_tsv(rep_map_file)

        cell_lines_to_use <- intersect(rep_map$CellLine, colnames(locus_cn))

        loci_to_use <- intersect(guide_locus_df$Locus, rownames(locus_cn))

        guides_to_use <-
            intersect(rownames(guide_dep),
                      c(guide_locus_df %>% dplyr::filter(Locus %in% loci_to_use) %$% Guide,
                        rownames(non_targeting_cn)))

        guide_dep <- guide_dep[guides_to_use,
                               rep_map %>%
                                   dplyr::filter(CellLine %in% cell_lines_to_use) %$%
                                   Replicate, drop=F]

        if (dep_normalize=="zmad") {
            guide_dep <- plyr::aaply(guide_dep, 2, function(cl) {
                (cl - median(cl, na.rm=T)) / mad(cl, na.rm=T)
            }) %>% t
        } else if (dep_normalize=="zscore") {
            guide_dep <- plyr::aaply(guide_dep, 2, function(cl) {
                (cl - mean(cl, na.rm=T)) / sd(cl, na.rm=T)
            }) %>% t
        } else if (dep_normalize=="none") {
        } else {
            stop("Error: normalization not recognized")
        }

        guide_locus_df <- guide_locus_df %>%
            dplyr::filter(Guide %in% guides_to_use,
                   Locus %in% loci_to_use) %>%
            dplyr::mutate(Value = 1)

        locus_gene_df <- locus_gene_df %>%
            dplyr::filter(Locus %in% loci_to_use) %>%
            dplyr::mutate(Value = 1)

        cat("writing to disk...\n\n")

        # TODO: Is it useful to copy this to guide_sample_dep if already existing as .Rds via generate_guides?
        saveRDS(guide_dep, file.path(inputs_dir, "guide_sample_dep.Rds"))

        saveRDS(guide_locus_df, file.path(inputs_dir, "guide_locus.Rds"))
        saveRDS(locus_gene_df, file.path(inputs_dir, "locus_gene.Rds"))
        saveRDS(locus_cn, file.path(inputs_dir, "locus_sample_cn.Rds"))
        rep_map %>% dplyr::filter(Replicate %in% colnames(guide_dep)) %>%
            saveRDS(file.path(inputs_dir, "replicate_map.Rds"))

        invisible(NULL)
    }
