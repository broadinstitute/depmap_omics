args <- commandArgs(TRUE)
inMAFfn=as.character(args[[1]]) 

maf <- read.delim(args[[1]], sep="\t", header=TRUE, stringsAsFactors = FALSE, comment.char = "#", quote="")
maf$i_ExAC_AF = maf$exac_af
write.table(maf, file= 'CLEANED.maf', sep="\t", row.names = FALSE, quote = FALSE)