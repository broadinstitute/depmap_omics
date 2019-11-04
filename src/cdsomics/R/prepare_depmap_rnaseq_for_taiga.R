#'
#' Prepare DepMap TPM RNAseq data for upload to Taiga
#'
#' @param rnaseq_data: rnaseq file, assumes rows are genes, gene_id and transcript_id columns, and rest of columns are samples
#' @param log_transform: whether or not to log2 transform the data
#' @param just_protein_coding: whether or not to just subset to protein coding genes and use entrez ids
#' @param gencode_file: path to the gencode file with Ensembl ID and HGNC symbol pairings
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Takes in a matrix of TPM values then log2 + 1 transforms it and returns gene expression matrix of
#' cell lines by genes
#' 
#' @export prepare_depmap_TPM_for_taiga
prepare_depmap_TPM_for_taiga <- function(rnaseq_data, 
                                            log_transform=F,
                                            just_protein_coding=F,
                                            gencode_file = "~/data_files/gene_annotations/gencode.v19.annotation.gtf"){
  
  require(magrittr)
  
  gene_map <- ensembl_hgnc_gene_pairing(gencode_file)
  row_metadata <- rnaseq_data %>% 
    dplyr::select(gene_id, `transcript_id(s)`) %>%
    dplyr::mutate(ensembl_gene_id=gsub("\\..*", "", gene_id)) %>% as.data.frame(.)
  
  rownames(row_metadata) <- row_metadata$ensembl_gene_id
  rownames(gene_map) <- gene_map$ensembl_gene_id
  
  intersect_genes <- row_metadata$ensembl_gene_id %in% gene_map$ensembl_gene_id
  row_metadata <- row_metadata[intersect_genes,]
  
  rnaseq_data <- rnaseq_data[intersect_genes,]
  gene_map <- gene_map[row_metadata$ensembl_gene_id,]
  
  if((nrow(gene_map)!=nrow(row_metadata)) & length(which(gene_map$ensembl_gene_id == row_metadata$ensembl_gene_id)) != nrow(row_metadata)) {
    print("ERROR! Incorrectly pairing gene names")
  }
  row_metadata$HGNC_symbol <- gene_map$hgnc_id
  row_metadata$concat_gene_name <- paste0(row_metadata$HGNC_symbol, " (", row_metadata$ensembl_gene_id, ")")
   
  data_matrix <- rnaseq_data %>% 
    dplyr::select(-gene_id, -`transcript_id(s)`) %>%
    as.matrix() %>%
    magrittr::set_rownames(row_metadata[["concat_gene_name"]]) %>%
    t() 
  
  if(log_transform) {
    data_matrix <- log2(data_matrix + 1)
  }
  
  # broad_ids <- gsub(".*\\((.*)\\).*", "\\1", rownames(data_matrix))
  # 
  # data_matrix <- data_matrix %>%
  #   magrittr::set_rownames(broad_ids)  
  
  if(just_protein_coding) {
      hgnc_gene_dataset <- taigr::load.from.taiga(data.name='hgnc-87ab', data.version=3, data.file='hgnc_complete_set-2018q3')
      pc_genes <- protein_coding_genes(row_metadata, hgnc_gene_dataset)
      pc_genes <- pc_genes %>%
        dplyr::mutate(concat_entrez_ID=paste0(symbol,
                                              " (",
                                              gsub("[.].*", "", entrez_id),
                                              ")"))
      if(length(unique(pc_genes$symbol)) < nrow(pc_genes)) {
        print("ERROR! HGNC symbols not unique")
      }
      
      print(paste("using", nrow(pc_genes), "protein coding genes"))
    data_matrix <- data_matrix[,pc_genes$concat_gene_name] %>% 
      magrittr::set_colnames(pc_genes$concat_entrez_ID)
  }
  
  
  return(data_matrix)
  
}

#'
#' Prepare DepMap TPM transcript RNAseq data for upload to Taiga
#'
#' @param rnaseq_data: rnaseq file, assumes rows are genes, gene_id and transcript_id columns, and rest of columns are samples
#' @param gencode_file: path to the gencode file with Ensembl ID and HGNC symbol pairings
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Takes in a matrix of transcript data and returns transcript data with gene names
#' 
#' @export prepare_depmap_transcripts_for_taiga
prepare_depmap_transcripts_for_taiga <- function(rnaseq_data, 
                                         gencode_file = "~/data_files/gene_annotations/gencode.v19.annotation.gtf"){
  
  require(magrittr)
  
  gene_map <- rtracklayer::readGFF(gencode_file, tags=c("gene_id", "gene_name", "transcript_id"))
  gene_map$transcript_id <- gsub("\\..*", "", gene_map$transcript_id)
  duplicate_or_missing_genes <- which(duplicated(gene_map$transcript_id)==T | is.na(gene_map$transcript_id)==T)
  if(length(duplicate_or_missing_genes) >0 ){
    gene_map <- gene_map[-duplicate_or_missing_genes,]
    
  }

  row_metadata <- rnaseq_data %>% 
    dplyr::select(gene_id, `transcript_id(s)`) %>%
    dplyr::mutate(transcript_id=gsub("\\..*", "", `transcript_id(s)`)) %>% as.data.frame(.)
  
  rownames(row_metadata) <- row_metadata$transcript_id
  rownames(gene_map) <- gene_map$transcript_id
  
  intersect_genes <- row_metadata$transcript_id %in% gene_map$transcript_id
  row_metadata <- row_metadata[intersect_genes,]
  
  rnaseq_data <- rnaseq_data[intersect_genes,]
  gene_map <- gene_map[row_metadata$transcript_id,]
  
  if((nrow(gene_map)!=nrow(row_metadata)) & length(which(gene_map$transcript_id == row_metadata$transcript_id)) != nrow(row_metadata)) {
    print("ERROR! Incorrectly pairing gene names")
  }
  
  row_metadata$HGNC_symbol <- gene_map$gene_name
  row_metadata$concat_gene_name <- paste0(row_metadata$HGNC_symbol, " (", row_metadata$transcript_id, ")")
  
  data_matrix <- rnaseq_data %>% 
    dplyr::select(-gene_id, -`transcript_id(s)`) %>%
    as.matrix() %>%
    magrittr::set_rownames(row_metadata[["concat_gene_name"]]) %>%
    t() 
  
  return(data_matrix)
  
}



#'
#' Prepare DepMap RNAseq data for upload to Taiga
#'
#' @param rnaseq_data: expression matrix
#' @param log_transform: whether or not to log2 transform the data
#' @param noisy_floor: if using a noisy floor, what value should be used
#' @param just_protein_coding: whether or not to just subset to protein coding genes and use entrez ids
#' @param gencode_file: path to the gencode file with Ensembl ID and HGNC symbol pairings
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes that RNAseq data arrive in .gct file (tab-delimited) 
#' with `Name`, `Description`, and CCLE cell line names as columns. The `Name` field contains
#' the Ensembl ID (with version number appended) for each gene, and the `Description` field 
#' contains the HUGO gene symbol.
#' 
#' @export prepare_depmap_rnaseq_for_taiga
prepare_depmap_rnaseq_for_taiga <- function(rnaseq_data, 
                                            log_transform=F,
                                            noisy_floor=NULL,
                                            just_protein_coding=F,
                                            gencode_file = "~/data_files/gene_annotations/gencode.v19.annotation.gtf"){
  
  require(magrittr)
  set.seed(4)
  gene_map <- ensembl_hgnc_gene_pairing(gencode_file)

  row_metadata <- rnaseq_data %>% 
    dplyr::select(Name, Description) %>%
    dplyr::mutate(concat_gene_name=paste0(Description, 
                                          " (", 
                                          gsub("[.].*", "", Name), 
                                          ")")) %>%
    dplyr::mutate(ensembl_gene_id=gsub("[.].*", "", Name))
  
  
  data_matrix <- rnaseq_data %>% 
    dplyr::select(-Name, -Description) %>%
    as.matrix() %>%
    magrittr::set_rownames(row_metadata[["concat_gene_name"]]) %>%
    t() 
  
  
  if(log_transform) {
    data_matrix <- log2(data_matrix + 1)
  }
  # broad_ids <- gsub(".*\\((.*)\\).*", "\\1", rownames(data_matrix))
  # 
  # data_matrix <- data_matrix %>%
  #   magrittr::set_rownames(broad_ids)  
  
  # if(!is.null(noisy_floor)){
  #   
  #   data_matrix %<>%
  #     magrittr::inset(is.infinite(.), noisy_floor)
  #   
  #   data_matrix[data_matrix <= noisy_floor] <- noisy_floor
  #   data_matrix[data_matrix == noisy_floor & !is.na(data_matrix)] <- data_matrix[data_matrix == noisy_floor & !is.na(data_matrix)] %>% magrittr::add(rnorm(n=length(.), sd=0.1))
  #   
  # }
  
  
    if(just_protein_coding) {
      hgnc_gene_dataset <- taigr::load.from.taiga(data.name='hgnc-87ab', data.version=3, data.file='hgnc_complete_set-2018q3')
      pc_genes <- protein_coding_genes(row_metadata, hgnc_gene_dataset)
      pc_genes <- pc_genes %>%
        dplyr::mutate(concat_entrez_ID=paste0(symbol,
                                              " (",
                                              gsub("[.].*", "", entrez_id),
                                              ")"))
      if(length(unique(pc_genes$symbol)) < nrow(pc_genes)) {
        print("ERROR! HGNC symbols not unique")
      }
      data_matrix <- data_matrix[,pc_genes$concat_gene_name] %>% 
        magrittr::set_colnames(row_metadata_subset$concat_entrez_ID)
    }
  
  return(data_matrix)
  
}
