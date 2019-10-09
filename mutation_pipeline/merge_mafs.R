# Usage: Rscript merge_mafs.R <mutect1_filename> <mutect2_filename> <output prefix>

install.packages("data.table", type = "source",
                 repos = "https://Rdatatable.github.io/data.table")
library(data.table)

args<-commandArgs(TRUE)
mutect1_file=as.character(args[[1]]) 
mutect2_file=as.character(args[[2]])
prefix=as.character(args[[3]])

mafs <- c(mutect1_file,mutect2_file)

merge_mafs = function(mafs, MAFobj = FALSE, ...){
  
  maf = lapply(mafs, data.table::fread, stringsAsFactors = FALSE, fill = TRUE,
               showProgress = TRUE, header = TRUE, skip = "Hugo_Symbol")
  names(maf) = gsub(pattern = "\\.maf$", replacement = "", x = basename(path = mafs), ignore.case = TRUE)
  maf = data.table::rbindlist(l = maf, fill = TRUE, idcol = "sample_id", use.names = TRUE)
  
  if(MAFobj){
    maf = read.maf(maf = maf, ...)
  }
  
  maf
}

merged_maf <- merge_mafs(mafs)
merged_maf$sample_id <- NULL
write.table(merged_maf, file=paste(prefix,".merged.maf",sep = ""), sep = "\t", quote = FALSE,row.names = FALSE)
