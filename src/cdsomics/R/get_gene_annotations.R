#'
#' Pair Ensembl gene ids with the correct HGNC gene symbol
#'
#' @param filename: path to gencode gene annotations file, possible files are gencode.v19.annotation.gtf or gencode.v29.annotation.gff3
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Parses the GENCODE gene annotations file to get the Ensembl ID and HGNC symbols for genes
#' 
#' @export ensembl_hgnc_gene_pairing
ensembl_hgnc_gene_pairing <- function(filename) {
  
  hgnc_gene_dataset <- taigr::load.from.taiga(data.name='hgnc-87ab', data.version=3, data.file='hgnc_complete_set-2018q3')
  
  gencode_gene_data <- rtracklayer::readGFF(filename, tags=c("gene_id", "gene_name", "transcript_id"))
  
  gencode_gene_data$gene_id <- gsub("\\..*", "", gencode_gene_data$gene_id)
  remove_id <- which(duplicated(gencode_gene_data$gene_id)==T)
  gencode_gene_data <- gencode_gene_data[-remove_id,]
  gencode_gene_data <- gencode_gene_data[,c("gene_id", "gene_name")]
  colnames(gencode_gene_data) <- c("ensembl_gene_id", "hgnc_id")
  
  gencode_gene_data <- plyr::join(gencode_gene_data, hgnc_gene_dataset[,c('symbol', 'ensembl_gene_id', 'locus_group', 'entrez_id')], by='ensembl_gene_id', type='left')
  gencode_gene_data$old_symbol <- gencode_gene_data$hgnc_id
  
  gencode_gene_data[which(!is.na(gencode_gene_data$symbol)), 'hgnc_id'] <- gencode_gene_data$symbol[which(!is.na(gencode_gene_data$symbol))]
  
  duplicated_genes <- which(duplicated(gencode_gene_data$ensembl_gene_id)==T)
  if(length(duplicated_genes) > 0) {
    print(paste(length(duplicated_genes), "duplicated ensembl ID(s)!"))
    gencode_gene_data <- gencode_gene_data[-duplicated_genes,]
    
  }
  
  return(gencode_gene_data)
}

#'
#' Test that GENCODE Ensembl ID-HGNC symbol pairings match HGNC dataswt
#'
#' @param gene_map: matrix containing the Ensembl ID and HGNC symbol pairings from GENCODE
#' @param version: dataset from HGNC with gene annotations
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description confirms that gene ID pairings are correct
#' 
#' @export test_gencode_ids
test_gencode_ids <- function(gene_map, hgnc_gene_dataset) {
  gene_map <- filter(gene_map, ensembl_gene_id %in% hgnc_gene_dataset$ensembl_gene_id)
  print(paste("comparing", nrow(gene_map), "genes"))
  
  gene_comparison <- merge(gene_map, hgnc_gene_dataset, by = 'ensembl_gene_id')
  
  print(length(which(gene_comparison$hgnc_id.x == gene_comparison$symbol)))
  
  print(gene_comparison[which(gene_comparison$hgnc_id.x != gene_comparison$symbol), c("hgnc_id.x", "symbol", "ensembl_gene_id")])
  
}



#'
#' Identify the gene names for protein coding genes
#'
#' @param gene_map: matrix containing the Ensembl ID and HGNC symbol pairings, the Ensembl ID column must be named ensembl_gene_id
#' @param hgnc_gene_dataset: dataset from HGNC with gene annotations
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description gets the Ensembl ID, HGNC symbol, and Entrez ID for protein coding genes
#' 
#' @export protein_coding_genes
protein_coding_genes <- function(gene_map, hgnc_gene_dataset) {
  gene_map <- dplyr::filter(gene_map, ensembl_gene_id %in% hgnc_gene_dataset$ensembl_gene_id)
  gene_dataset <- plyr::join(gene_map, hgnc_gene_dataset[,c(2,4,19,20)], by='ensembl_gene_id')
  protein_genes <- dplyr::filter(gene_dataset, locus_group=="protein-coding gene")
  
  return(protein_genes)
}



