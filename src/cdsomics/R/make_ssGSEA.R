#'
#' Makes single sample GSEA for latest expression data
#'
#' @param data_mat: matrix of gene expression data where columns are HGNC symbol (Entrez ID) and rows are samples
#' @param gmt_path: path to gmt file of gene set object
#' @param gsva_method: method to use in estimation of gene-set enrichment scores per sample, default is ssgsea. 
#' ssGSEA calculates a gene set enrichment score per sample as the normalized difference in empirical cumulative distribution
#' functions of gene epxression ranks inside and outside the gene set
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @descriptison estimate single sample GSEA scores for protein coding TPM dataset
#' 
#' @export make_ssGSEA
make_ssGSEA <- function(data_mat, gmt_path, gsva_method="ssgsea") {
  gsc_obj <- GSEABase::getGmt(gmt_path, 
                              collectionType = GSEABase::BroadCollection(),
                              geneIdType = GSEABase::EntrezIdentifier())
  
  entrez_ids <- stringr::str_match(colnames(data_mat), '\\(([0-9]+)\\)')[,2]
  stopifnot(sum(duplicated(entrez_ids)) == 0)
  colnames(data_mat) <- entrez_ids
  
  ssGSEA <- GSVA::gsva(t(data_mat), gsc_obj, method = gsva_method)
  
  # returns gene-set by sample matrix of ssGSEA enrichment scores
  return(t(ssGSEA))
}



