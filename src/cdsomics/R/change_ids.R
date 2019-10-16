#'
#' Map one set of identifiers to another
#'
#' @param input_df: input data frame
#' @param map_df: input data frame of mapping ids
#' @param from: parameter in map_df containing the original ids
#' @param to: parameter in map_df containing the new ids
#' @param dim: whether you are changing rows (1) or columns (2)
#' 
#' @importFrom magrittr "%>%"
#' 
#' @description Maps identifiers given in `input_ids` to new identifiers contained in `map_df`
#' 
#' @export change_ids
#' 
change_ids <- function(input_df, map_df, from, to, dim){
  
  
  map_df %<>% dplyr::filter(!is.na(to))
  
  input_ids <- ""
  if(dim==1) {
    input_ids <- rownames(input_df)
  } else {
    input_ids <- colnames(input_df)
  }
  
  if(!all(c(from, to) %in% colnames(map_df))){
    stop("map_df does not contain both `from` and `to` fields...")
  }
  
  if(!all(input_ids %in% map_df[[from]])){
    user_question <- readline("Not all input_ids contained in map_df, would you like to subset the data to just entries where a mapping exists (y/n)? ")
    if(regexpr(user_question, 'y', ignore.case = TRUE) == 1) {
      missing_in_maps <- which(!input_ids %in% map_df[[from]])
      input_ids <- input_ids[-missing_in_maps]
      if(dim==1) {
        input_df <- input_df[-missing_in_maps,]
      } else {
        input_df <- input_df[,-missing_in_maps]
      }
    } else {
      stop("Not all input_ids contained in map_df")
    }
  }
  
  if(any(duplicated(map_df[[from]]))){
    stop("User supplied duplicate mappings (IDs in `from` map to multiple IDs in `to`)")
  }
  
  
  rownames(map_df) <- map_df[,from]
  if(length(which(is.na(map_df[input_ids,to])))>0){
    user_question <- readline("Not all new mapping ids contained in map_df, would you like to subset the data to just entries where a mapping exists (y/n)? ")
    if(regexpr(user_question, 'y', ignore.case = TRUE) == 1) {
      missing_out_maps <- which(is.na(map_df[input_ids,to]))
      input_ids <- input_ids[-missing_out_maps]
      if(dim==1) {
        input_df <- input_df[-missing_out_maps,]
      } else {
        input_df <- input_df[,-missing_out_maps]
      }
    } else {
      stop("Not all input_ids contained in map_df")
    }
  }
  
  
  rownames(map_df) <- map_df[,from]
  new_ids <- map_df[input_ids, to]
  if(dim == 1) {
    input_df <- input_df %>% 
      magrittr::set_rownames(new_ids)
  } else {
    input_df <- input_df %>% 
      magrittr::set_colnames(new_ids)
  }
  
  return(input_df)
}
