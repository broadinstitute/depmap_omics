#'
#' Map one set of identifiers to another
#'
#' @importFrom magrittr "%>%"
#' 
#' @param input_ids: the input identifiers
#' @param map_df: data frame containing the identifier mapping
#' @param from: parameter in map_df that you are mapping from
#' @param to: parameter in map_df that you are mapping to
#' 
#' @description Maps identifiers given in `input_ids` to new identifiers contained in `map_df`
#' 
#' @export map_ids
map_ids <- function(input_ids, map_df, from, to){
  
  
  map_df %<>% dplyr::filter(!is.na(to))
  
  if(!all(c(from, to) %in% colnames(map_df))){
    stop("map_df does not contain both `from` and `to` fields...")
  }
  
  if(!all(input_ids %in% map_df[[from]])){
      stop("Not all input_ids mapped in map_df")
  }
  
  if(any(duplicated(map_df[[from]]))){
    stop("User supplied duplicate mappings (IDs in `from` map to multiple IDs in `to`)")
  }
  
  data.frame(input_ids) %>%
    magrittr::set_colnames(from) %>%
    dplyr::left_join(map_df, by = from) %>%
    magrittr::extract2(to)
  
}
