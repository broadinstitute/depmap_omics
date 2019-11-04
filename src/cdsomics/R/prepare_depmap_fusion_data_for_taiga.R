#'
#' Prepare DepMap fusion data for upload to Taiga
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @param fusion_data: fusion data matrix
#' 
#' @description Adjust column names and syntax of gene names for fusion data
#' 
#' @export prepare_depmap_fusion_data_for_taiga
prepare_depmap_fusion_data_for_taiga <- function(fusion_data){
  
  left_gene <- fusion_data$LeftGene
  left_genes <- strsplit(left_gene, "\\^")
  left_gene_symbol <- sapply(left_genes, `[`, 1)
  left_gene_ensembl <- gsub("\\..*", "", sapply(left_genes, `[`, 2))
  fusion_data$LeftGene <-  paste0(left_gene_symbol, " (", left_gene_ensembl, ")")

  right_gene <- fusion_data$RightGene
  right_genes <- strsplit(right_gene, "\\^")
  right_gene_symbol <- sapply(right_genes, `[`, 1)
  right_gene_ensembl <- gsub("\\..*", "", sapply(right_genes, `[`, 2))
  fusion_data$RightGene <-  paste0(right_gene_symbol, " (", right_gene_ensembl, ")")
  
  
    return(fusion_data)
  
}

