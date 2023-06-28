args<-commandArgs(TRUE)

countfile <- args[1];
gmtfile <- args[2];
method <- args[3]

library(GSEABase)
library(GSVA)
counts <- read.csv(countfile, row.names=1)
mat <- data.matrix(counts, rownames.force = T)
colnames(mat) <- colnames(counts)
gsc_obj <- GSEABase::getGmt(gmtfile,
collectionType = GSEABase::BroadCollection(),
geneIdType = GSEABase::EntrezIdentifier())
gsea <- GSVA::gsva(mat, gsc_obj, method = method)
write.table(gsea, file = "/tmp/res_genepy_ssGSEA.tsv", sep = '\t', quote = F)
