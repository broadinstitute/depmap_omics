#'
#' Prepare DepMap mutation data for upload to Taiga
#'
#' @param mut_maf: maf file
#' @param gene_symbol_col: name of column containing the gene symbols
#' @param gene_id_col: name of column containing the gene ids
#' @param damaging: initially is set to F, set to T to create binary matrix of damaging mutations 
#' @param other: initially is set to F, set to T to create binary matrix of conserving and non conserving mutations 
#' @param hotspot: initially is set to F, set to T to create binary matrix of hotspot mutations 
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes that mutation calls have been provided for each cell line
#' 
#' @export mutation_maf_to_binary_matrix
mutation_maf_to_binary_matrix <- function(mut_maf,
                                          gene_symbol_col="Hugo_Symbol",
                                          gene_id_col="Entrez_Gene_Id",
                                          damaging=F, other=F, hotspot=F){
  
  mutation_data <- mut_maf %>%
    dplyr::rename_(gene_symbol=gene_symbol_col,
                   gene_id=gene_id_col) %>%
    dplyr::mutate(concat_gene_name=paste0(gene_symbol,
                                          " (", gene_id, ")"))
  
  if(damaging) {
    mutation_data <- mutation_data %>%
                      dplyr::filter(Variant_annotation == 'damaging') %>%
                      dplyr::mutate(value=1) %>%
                      reshape2::acast(DepMap_ID~concat_gene_name, value.var="value",
                                      fun.aggregate = function(x){return(ifelse(any(!is.na(x)), 1, 0))})
    
  } else if(other) {
    mutation_data <- mutation_data %>%
      dplyr::filter(Variant_annotation %in% c('other conserving', 'other non-conserving')) %>%
      dplyr::mutate(value=1) %>%
      reshape2::acast(DepMap_ID~concat_gene_name, value.var="value",
                      fun.aggregate = function(x){return(ifelse(any(!is.na(x)), 1, 0))})
    
  } else if(hotspot) {
    mutation_data <- mutation_data %>%
      dplyr::filter((isTCGAhotspot == TRUE | isCOSMIChotspot == T) & Variant_annotation != 'silent') %>%
      dplyr::mutate(value=1) %>%
      reshape2::acast(DepMap_ID~concat_gene_name, value.var="value",
                      fun.aggregate = function(x){return(ifelse(any(!is.na(x)), 1, 0))})
    
  }
  return(mutation_data)
  
}

