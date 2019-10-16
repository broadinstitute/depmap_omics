#'
#' Add a new taiga nickname
#' 
#' @param nickname: nickname for the taiga file
#' @param data_name: taiga data name
#' @param data_file: taiga data file
#' @param data_version: taiga data version
#' @param filename: name of the taiga nickname file if it exists
#' 
#' @importFrom magrittr "%>%"
#' 
#' @description Adds a new nickname to the taiga_nicknames file. Creates a new taiga_nicknames file if none found.
#' 
#' @export add_taiga_nickname
add_taiga_nickname <- function(nickname, data_name, data_file, data_version=NA, filename=NULL){
  
  if(!is.null(filename)){
    to_write <- filename
    if(file.exists(filename)){
      append <- T
    } else{
      append <- F
    }
  } else if(file.exists("~/.taiga/taiga_nicknames")){
    to_write <- "~/.taiga/taiga_nicknames"
    append <- T
  } else if(file.exists("./taiga_nicknames")){
    to_write <- "./taiga_nicknames"
    append <- T
  } else {
    message("Nicknames file not found or does not exist; new file will be created at `~/.taiga/taiga_nicknames`")
    to_write <- "~/.taiga/taiga_nicknames"
    append <- F
  }
  
  readr::write_tsv(data.frame(taiga_nickname = nickname,
                              data_name = data_name,
                              data_file = data_file,
                              data_version = data_version,
                              time_stamp = Sys.time()),
                   path = to_write,
                   append = append)
  
}