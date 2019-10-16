#'
#' Prepare RRBS data for upload to Taiga
#' 
#' @param filename1: path to RRBS_1kb file
#' @param filename2: path to RRBS_CpG file
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes that RRBS data arrive in two .txt files (tab-delimited) 
#' The first is 1kb data and the other is CpG cluster data
#' The 1kb dataset has columns TSS_id, gene (which is HGNC symbol), chr, fpos, tpos, strand, avg_coverage, and the CCLE IDs
#' The CpG cluster dataset has columns cluster_id, gene_name (which is HGNC symbol), RefSeq_id, CpG_sites_hg19, avg_coverage and the CCLE IDs                         "avg_coverage"
#' Return four datasets, methylation data with columns of genes (HGNC symbols) and rows of Broad cell line IDs for 1kb data and CpG data
#' and meta data file for the 1kb data and CpG data
#' 
#' @export prepare_RRBS_for_taiga
prepare_RRBS_for_taiga <- function(filename1, filename2){
  
  require(magrittr)
  
  RRBS_1kb <- readr::read_tsv(filename1, guess_max = 10000)
  RRBS_CpG <- readr::read_tsv(filename2, guess_max = 10000)
  
  mapping <- read.csv(url("https://intranet.broadinstitute.org/~datasci/cell_lines/name_mapping.csv"), stringsAsFactors = F)
  mapping <- mapping[-which(duplicated(mapping$ccle_name) ==T),]
  mapping <- mapping %>% magrittr::set_rownames(mapping$ccle_name)
  
  gene_name_1kb <- RRBS_1kb$TSS_id
  gene_name_CpG <- RRBS_CpG$cluster_id
  
  data_1kb <- RRBS_1kb %>% 
    dplyr::select(-TSS_id, -gene, -chr, -fpos, -tpos, -strand, -avg_coverage) %>%
    as.matrix() 
  meta_1kb <- RRBS_1kb %>% 
    dplyr::select(TSS_id, chr, fpos, tpos, strand, avg_coverage) %>%
    as.matrix() 
  
  data_CpG <- RRBS_CpG %>% 
    dplyr::select(-cluster_id, -gene_name, -RefSeq_id, -CpG_sites_hg19, -avg_coverage) %>%
    as.matrix() 
  meta_CpG <- RRBS_CpG %>% 
    dplyr::select(cluster_id, RefSeq_id, CpG_sites_hg19, avg_coverage) %>%
    as.matrix()
  
  broad_id1 <- mapping[colnames(data_1kb), 'broad_id']
  broad_id2 <- mapping[colnames(data_CpG), 'broad_id']
  
  data_1kb <- data_1kb %>% 
    magrittr::set_colnames(broad_id1) %>%
    magrittr::set_rownames(gene_name_1kb) %>%
    t()
  
  meta_1kb <- meta_1kb %>% 
    magrittr::set_rownames(gene_name_1kb)
    
  
  data_CpG <- data_CpG %>% 
    magrittr::set_colnames(broad_id2) %>%
    magrittr::set_rownames(gene_name_CpG) %>%
    t()
  meta_CpG <- meta_CpG %>% 
    magrittr::set_rownames(gene_name_CpG)
  
  return(list(data_1kb, meta_1kb, data_CpG, meta_CpG))
  
}
