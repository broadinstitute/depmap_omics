countfile = args[1];
gmtfile = args[2];

library(GSEABase)
library(GSVA)
counts <- read.csv(countfile, row.names=1)
mat <- data.matrix(counts, rownames.force = T)
colnames(mat) <- colnames(counts)
gsc_obj <- GSEABase::getGmt(gmtfile,
collectionType = GSEABase::BroadCollection(),
geneIdType = GSEABase::EntrezIdentifier())
gsea <- GSVA::gsva(mat, gsc_obj, method = 'ssgsea')
write.table(gsea, file = "temp/res_ssGSEA.csv", sep = ',', quote = F)
