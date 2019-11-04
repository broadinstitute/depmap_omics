#' Replace cell lines with embargoed Broad WES CN data with chordoma, Sanger or SNP CN data
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @param gene_cn_broad: gene level copy number data containing data from internal and blacklisted cell lines
#' @param gene_cn_chordoma: gene level copy number data for the chordoma lines
#' @param gene_cn_sanger: gene level copy number data from Sanger
#' @param gene_cn_SNP: gene level copy number data using SNP data
#' @param replace_CLs: cell lines that are not allowed for this dataset (have embargoed WES)
#' 
#' @description 
#' 
#' @export subset_gene_level_cn
subset_gene_level_cn <- function(gene_cn_broad, gene_cn_chordoma, gene_cn_sanger, gene_cn_SNP, replace_CLs) {
  for(CL in replace_CLs) {
    if(CL %in% rownames(gene_cn_chordoma)) {
      gene_cn_broad[CL,] <- gene_cn_chordoma[CL,colnames(gene_cn_broad)] %>% log2()
    } else if(CL %in% rownames(gene_cn_sanger)) {
      gene_cn_broad[CL,] <- gene_cn_sanger[CL,colnames(gene_cn_broad)] %>% log2()
    } else if(CL %in% rownames(gene_cn_SNP)) {
      gene_cn_broad[CL,] <- gene_cn_SNP[CL,colnames(gene_cn_broad)] %>% log2()
    } else {
      gene_cn_broad <- gene_cn_broad[-CL,]
    }
  }
  
  return(gene_cn_broad)
  
}

