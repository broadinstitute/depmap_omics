#'
#' Prepare miRNA data for upload to Taiga
#' 
#' @param filename: path to miRNA file
#' @param log_transform: boolean, whether or not log2(.+50) transform the data
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes that miRNA data arrive in .gct file (tab-delimited) 
#' with `Name`, `Description`, and CCLE cell line names as columns. The `Name` field contains
#' the Accession number for each gene, and the `Description` field 
#' contains the miRNA ID.
#' 
#' @export prepare_miRNA_for_taiga
prepare_miRNA_for_taiga <- function(filename, 
                                            log_transform=T){
  
  require(magrittr)
  
  miRNA_data <- readr::read_tsv(filename, skip=2, guess_max = 10000)
  
  mapping <- read.csv(url("https://intranet.broadinstitute.org/~datasci/cell_lines/name_mapping.csv"), stringsAsFactors = F)
  mapping <- mapping[-which(duplicated(mapping$ccle_name) ==T),]
  mapping <- mapping %>% magrittr::set_rownames(mapping$ccle_name)
  
  miRNA_IDs <- miRNA_data$Description
  
  data_matrix <- miRNA_data %>% 
    dplyr::select(-Name, -Description) %>%
    as.matrix() 
  
  broad_id <- mapping[colnames(data_matrix), 'broad_id']
    
  data_matrix <- data_matrix %>% 
    magrittr::set_colnames(broad_id) %>%
    magrittr::set_rownames(miRNA_IDs) %>%
    t() %>%
    {if(!log_transform) . else log2(. + 50)}
  
  
  return(data_matrix)
  
}
