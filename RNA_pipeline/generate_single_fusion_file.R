library(methods)
library(readr)
library(stringr)
library(plyr)
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
  
  # Extract the arxspan id from the file name
  # sample <- stringr::str_extract(string = f, pattern = 'ACH\\-[0-9]+')
  
  # In this workspace, samples are not indexed by ARXSPAND ID but by CCLE_name. Will need to update script in the future
  sample <- gsub('\\.fusions.annotated$', '', gsub('.*/', '', f))
  segs <- read_tsv(f, col_names = TRUE, col_types = cols(.default = "c")) %>%
    mutate(DepMap_ID=sample) %>% 
    dplyr::select(DepMap_ID, everything())
  
  # Write to the output file
  write.table(segs,  file=output_file_name,  sep='\t', row.names = FALSE, quote = FALSE, append = (i>1), col.names = !(i>1))
}
