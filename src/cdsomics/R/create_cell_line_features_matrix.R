#'
#' Prepare DepMap features data for upload to Taiga
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Updates the cell line feature matrix as new data is added
#' 
#' @export create_cell_line_features_matrix
create_cell_line_features_matrix <- function(){
  masterfile <- taigr::load.from.taiga(data.name='master-cell-line-export-0306', data.version=55)
  
}

