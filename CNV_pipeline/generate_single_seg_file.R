library(dplyr)

## reads selected fields from maf files and aggregates them 
args <- commandArgs(TRUE)
print(args)

output_file_name = args[1]
input_file_names = args[2]

# Load files and ensure that we are only reading in files that exist
fns = as.character(read.csv(input_file_names, header=FALSE)[,])
fe = sapply(fns, file.exists)
samples = fns[which(fe)]
nsamples = length(samples)

# Iterate across samples, read tsv, add the DepMap_ID column and keep the rest of them
for(i in 1:nsamples){
  f <- samples[i] 
  sample <- f %>% gsub('\\..*', '', .) %>% gsub('.*/', '', .) # Remove file paths info
  
  segs <- read_tsv(f, comment = '@', col_names = TRUE, cols(
    CONTIG = col_character(),
    START = col_double(),
    END = col_double(),
    NUM_POINTS_COPY_RATIO = col_double(),
    MEAN_LOG2_COPY_RATIO = col_double(),
    CALL = col_character()
  )) %>%
    mutate(Sample=sample) %>% dplyr::select(Sample, everything())
  
  # Write to the output file
  write.table(segs,  file=output_file_name,  sep='\t', row.names = FALSE, quote = FALSE, append = (i>1), col.names = !(i>1))
}
