# Input: RPKM file from rnaseqc_counts, cell line sample_id
# Output: File containing columns for Name, Description, RPKM, and TPM (where TPM was converted from RPKM)

# Extract command line arguments: 1st is the RPKM file and 2nd is the cell line sample_id
args<-commandArgs(TRUE)
reads_file=as.character(args[[1]])
sample_id = as.character(args[[2]])

# Read in the RPKM file 
"rpkm_data" <- data.frame(read.delim(reads_file, header=TRUE))
colnames(rpkm_data) = c("Name", "Description", "RPKM")

# Conversion to TPM
rpkm_sum <- sum(rpkm_data$RPKM)
rpkm_data$TPM <- (rpkm_data$RPKM / rpkm_sum) * 1000000

# Write output file 
write.table(rpkm_data,file=paste(sample_id,"gene_rpkm_tpm.gct",sep="."), quote = FALSE, sep = "\t", row.names = FALSE)
