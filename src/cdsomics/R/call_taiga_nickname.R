#'
#' Load taiga data using a nickname
#' 
#' @param nickname: nickname for taiga file you are referencing
#' @param filename: path to file of taiga nicknames if not just taiga_nicknames
#' 
#' @importFrom magrittr "%>%"
#' 
#' @description Loads data from taiga given most recent nickname pre-specified by the user.
#' 
#' @export call_taiga_nickname
call_taiga_nickname <- function(nickname, filename=NULL){
  
  require(magrittr)
  require(taigr)
  
  # Check for tg-nicknames in ~/.taiga
  if(!is.null(filename)){
    tg_nicknames <- readr::read_tsv(filename)
  } else if(file.exists("~/.taiga/taiga_nicknames")){
    tg_nicknames <- readr::read_tsv("~/.taiga/taiga_nicknames")
  } else if(file.exists("./taiga_nicknames")){
    tg_nicknames <- readr::read_tsv("./taiga_nicknames")
  } else {
    stop("Please specify location of nickname file...")
  }
  
  has_nickname_columns <- all(c("taiga_nickname", "data_name", "data_file", "data_version", "time_stamp") %in% colnames(tg_nicknames))
  if(!has_nickname_columns){
    stop("Nickname file must have colnames `taiga_nickname`, `data_name`, `data_file`, `data_version` and `time_stamp`") 
  }
  
  tg_nicknames %<>% 
    dplyr::filter(taiga_nickname == nickname) %>%
    plyr::arrange(dplyr::desc(time_stamp)) %>%
    magrittr::extract(1, 1:ncol(.))
  
  if(is.na(tg_nicknames[["data_version"]])){
    tg_nicknames[["data_version"]] <- NULL
  }
  
  taigr::load.from.taiga(data.name=tg_nicknames[["data_name"]],
                         data.file=tg_nicknames[["data_file"]],
                         data.version=tg_nicknames[["data_version"]])
  
}