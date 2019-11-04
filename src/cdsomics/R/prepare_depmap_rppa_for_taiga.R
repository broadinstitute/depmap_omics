#'
#' Prepare DepMap RPPA data for upload to Taiga
#'
#' @param rppa_filename: path to rppa file
#' @param antibody_info_filename: path to antibody info file
#' @param antibody_id_col: name of column containing the antibody names
#' @param gene_id_col: name of column containing the gene names
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Assumes that RPPA data arrive as a .csv file with CCLE cell line
#' names on the rows and antibody identifiers on the columns. The antibody info file
#' is assumed to be a .csv file with a mapping from antibody identifier to HUGO symbol.
#' 
#' @export prepare_depmap_rppa_for_taiga
prepare_depmap_rppa_for_taiga <- function(rppa_filename, antibody_info_filename,
                                          antibody_id_col="Antibody_Name",
                                          gene_id_col="Target_Genes"){
  
  rppa_data <- readr::read_csv(rppa_filename) %>%
                magrittr::set_rownames(magrittr::extract2(., "X1")) %>%
                dplyr::select(-X1) %>%
                as.matrix()
  
  antibody_info <- readr::read_csv(antibody_info_filename)
  
  return(list(rppa_data=rppa_data, antibody_info=antibody_info))
  
}