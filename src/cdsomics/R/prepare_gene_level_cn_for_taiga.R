#'
#' Prepare gene-level copy number matrix for upload to taiga
#'
#' @param filename: path copy number file
#' @param ccds_filename: path to ccds file
#' @param do_parallel: whether or not to run in parallel
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes input filename is path to a segmented copy number file
#' 
#' @export prepare_gene_level_cn_for_taiga
prepare_gene_level_cn_for_taiga <- function(filename, ccds_filename, do_parallel){
  
  ccds <- load_ccds(ccds_filename)
  
  gene_gr <- ccds %>%
    dplyr::distinct(gene, chromosome, strand, gene_start, gene_end) %>%
    dplyr::group_by(gene, chromosome, strand) %>%
    dplyr::summarise(start = min(gene_start),
                     end = max(gene_end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  seg_gr <- cn_seg %>%
    plyr::llply(GenomicRanges::makeGRangesFromDataFrame, keep.extra.columns=T, .parallel=do_parallel)
  cn_mat <- intersect_gene_with_copy_number(gene_gr, seg_gr,
                                            CN.column="CN",
                                            gene.column="gene",
                                            do_parallel=do_parallel) %>%
                                            {.[!duplicated(rownames(.)),]} %>%
    cdsr::remove_rows_all_na()
  
}

load_ccds <- function(ccds_file){
  
  
  ccds_h37 <- readr::read_tsv(ccds_file,
                       col_types=cols("#chromosome" = col_character(),
                                      "cds_from" = col_integer(),
                                      "cds_to" = col_integer())) %>%
    dplyr::rename(chromosome=`#chromosome`) %>%
    dplyr::mutate(chromosome = str_c("chr", chromosome)) %>%
    dplyr::filter(ccds_status %in% c("Public", "Reviewed, update pending", "Under review, update"),
           chromosome %in% chromosomes,
           !is.na(cds_from), !is.na(cds_to))
  
  ccds_exon_h37 <-  ccds_h37 %>%
    dplyr::mutate(cds_interval = str_replace_all(cds_locations, "[\\[\\]]", "") %>%
             stringr::str_split("\\s*,\\s*")) %>%
    tidyr::unnest(cds_interval) %>%
    dplyr::group_by(gene, gene_id, cds_locations) %>%
    dplyr::mutate(exon_code = ifelse(cds_strand=="+", 1:n(), n():1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cds_start = str_extract(cds_interval, "^[0-9]+") %>% as.integer,
           cds_end = str_extract(cds_interval, "[0-9]+$") %>% as.integer) %>%
    dplyr::select(gene, gene_id, chromosome, start=cds_start, end=cds_end, strand=cds_strand,
           gene_start=cds_from, gene_end=cds_to, exon_code)
  
  return(ccds_exon_h37)
  
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