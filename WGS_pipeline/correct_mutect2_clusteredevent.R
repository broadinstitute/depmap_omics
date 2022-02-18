library(vcfR)
library(fuzzyjoin)
library(tools)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

filepath <- args[1]
vcf <- vcfR::read.vcfR(filepath)
df <- vcf@fix %>% as_tibble() %>% transmute(chromosome = CHROM, start = POS,ref = REF,alt = ALT,
                                                 filter = FILTER) %>% type_convert()
df <- df %>% mutate(germline = str_detect(filter,"germline|panel_of_normals"),pos = start,start = pos - 50,end = pos + 50)
neighbor_df <- df %>% genome_join(df %>% transmute(chromosome,start,end,neighbor_germline = germline),
                              by = c("chromosome","start","end"))
neighbor_df <- neighbor_df %>% rename(chromosome = chromosome.x)
count_df <- neighbor_df %>% group_by(chromosome,pos,ref,alt,filter,germline) %>% 
  summarise(n_somatic_neighbor = sum(!neighbor_germline)) %>% ungroup()
count_df <- count_df %>% mutate(filter = ifelse(filter %in% c("clustered_events","clustered_events;haplotype") & 
                                      n_somatic_neighbor <= 2,"PASS",filter))
vcf@fix[,"FILTER"] <- count_df$filter

write.vcf(vcf,str_c(file_path_sans_ext(filepath),".cluster-corrected.vcf.gz"))