#'
#' Map one set of row and column identifiers to another
#'
#' @importFrom magrittr "%>%"
#' 
#' @param data_matrix: input data frame
#' @param row_map: data frame containing mapping for row identifiers
#' @param col_map: data frame containing mapping for column identifiers
#' @param row_from: parameter in row_map containing the original row identifiers
#' @param row_to: parameter in row_map contaning the new row identifiers
#' @param col_from: parameter in col_map containing the original column identifiers
#' @param col_to: parameter in col_map containing the new column identifiers
#' 
#' @description Maps matrix row and column names to new identifiers given in `row_map` and `col_map`
#' 
#' @export map_dim_names
map_dim_names <- function(data_matrix, 
                          row_map, col_map,
                          row_from, row_to,
                          col_from, col_to){
  
  
  new_row_names <- map_ids(input_ids = row.names(data_matrix),
                                     map_df = row_map,
                                     from = row_from,
                                     to = row_to)
  
  new_col_names <- map_ids(input_ids = colnames(data_matrix),
                                     map_df = col_map,
                                     from = col_from,
                                     to = col_to)
  
  data_matrix %>% 
    magrittr::set_rownames(new_row_names) %>%
    magrittr::set_colnames(new_col_names)
  
  
}
