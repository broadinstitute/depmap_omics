#'
#' Test the new data to make sure it mostly matches the previous dataset
#'
#' @param new_data: newest version of the RNAseq dataset
#' @param previous_data: previous release of RNAseq data
#' @param data_type: either RNAseq data or CN for copy number data
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Look at the correlation between cell lines in common between the new data and the previous release of data
#' 
#' @export test_new_data
test_new_data <- function(new_data, previous_data, data_type='RNAseq') {
  common_CLs <- intersect(rownames(new_data), rownames(previous_data))
  common_genes <- intersect(colnames(new_data), colnames(previous_data))
  
  print(paste(length(common_genes), "common genes out of", ncol(new_data), "genes in the new data and", ncol(previous_data), 
              "genes in the old data"))
  print(paste(length(common_CLs), "common cell lines out of", nrow(new_data), "cell lines in the new data and", nrow(previous_data), 
            "cell lines in the old data"))
  
  print(setdiff(rownames(previous_data), rownames(new_data)))
  
  new_data <- new_data[common_CLs, common_genes]
  previous_data <- previous_data[common_CLs, common_genes]
 
  print("correlation of overlap: ")
  
  data_cor <- unlist(sapply(seq(1,length(common_CLs)), function(i) cor(new_data[i,], previous_data[i,], use='pairwise'))) 
  
  if(data_type=='CN') {
    gene_level_diff <- unlist(sapply(seq(1,length(common_CLs)), function(i) length(which(abs(new_data[i,] - previous_data[i,]) > .5))))
    print(summary(gene_level_diff))
    #print(plot(data_cor, gene_level_diff))
  }
  print(summary(data_cor))
}

#'
#' Test the new data to make sure it contains only allowed cell lines
#'
#' @param df: matrix of data, assumes rownames are cell lines unless cell_line_column is not null
#' @param restricted_cell_lines:embargoed and blacklisted cell lines
#' @param data_type: either RNAseq or WES
#' @param dataset: either internal, DMC, or public
#' @param cell_line_column: defaults to assuming cell lines are rownames, otherwise input the name of the column containing the cell lines
#' @param is_CN: is the dataset copy number data, T or F (only matters for public CN)
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Check that the data does not include embargoed or blacklisted cell lines
#' 
#' @export test_released_cell_lines
test_released_cell_lines <- function(df, restricted_cell_lines, data_type, dataset, cell_line_column='rows', is_CN=F) {
  not_allowed_cell_lines <- restricted_cell_lines$blacklist$Blacklist
  if(data_type == 'RNAseq') {
    if(dataset == 'DMC') {
      not_allowed_cell_lines <- c(not_allowed_cell_lines, restricted_cell_lines$RNAseq_embargo$IBM)  
    } else if(dataset == 'public') {
      not_allowed_cell_lines <- c(not_allowed_cell_lines, restricted_cell_lines$RNAseq_embargo$IBM, restricted_cell_lines$RNAseq_embargo$DMC)  
      }
  } else {
    if(dataset == 'public' & is_CN == F) {
      not_allowed_cell_lines <- c(not_allowed_cell_lines, restricted_cell_lines$WES_embargo$DMC)
    } else if (dataset == 'public' & is_CN == T) {
      possible_error <- intersect(restricted_cell_lines$WES_embargo$DMC, rownames(df))
      if(length(possible_error)>0) {
        print(paste('public CN data contains', length(possible_error), "cell line(s) with restricted broad WES data, possible error"))
        print(possible_error)
      }
    }
  }
  
  if(cell_line_column=='rows') {
    error_CLs <- intersect(not_allowed_cell_lines, rownames(df))
    if(length(error_CLs) > 0) {
      print(paste('data includes', length(error_CLs), 'cell line(s) that should not be included in this dataset'))
      print(error_CLs)
    }
  } else {
    cur_CLs <- df %>%
      dplyr::select(cell_line_column)
    error_CLs <- intersect(not_allowed_cell_lines, cur_CLs$DepMap_ID)
    if(length(error_CLs) > 0) {
      print(paste('data includes', length(error_CLs), 'cell line(s) that should not be included in this dataset'))
      print(error_CLs)
    }
  }
  
}

