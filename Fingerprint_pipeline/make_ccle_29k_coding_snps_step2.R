library(tidyverse)
# These steps are required for lifting over the map to hg38

# Combine chromosome map files
# cat 1.filtered.map 2.filtered.map 3.filtered.map 4.filtered.map 5.filtered.map 6.filtered.map 7.filtered.map 8.filtered.map 9.filtered.map 10.filtered.map 11.filtered.map 12.filtered.map 13.filtered.map 14.filtered.map 15.filtered.map 16.filtered.map 17.filtered.map 18.filtered.map 19.filtered.map 20.filtered.map 21.filtered.map 22.filtered.map X.filtered.map > ccle_29k_coding_snps_hg19

# Convert to map to vcf by running these lines
map <- read_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg19",
                col_names = c("chrom","pos","name","major","minor","maf","anchor"),col_types = "cdcccdc")

map %>% mutate(anchor = ifelse(is.na(anchor),"none",anchor)) %>% 
  transmute(chrom,pos,name,major,minor,qual = "100",filter = "PASS",info = str_c("AF=",maf,";","ANCHOR=",anchor)) %>% 
  write_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg19.vcf",col_names = F)

# Add hg19 header to ccle_29k_coding_snps_hg19.vcf

# Run java -Xmx6g -jar ~/Library/gatk/gatk-package-4.2.0.0-local.jar LiftoverVcf --INPUT ccle_29k_coding_snps_hg19.vcf --OUTPUT ccle_29k_coding_snps_hg38.vcf  --CHAIN ~/Data/VCFs/Liftover/b37ToHg38.over.chain -REFERENCE_SEQUENCE ~/Data/VCFs/Liftover/hg38.fa --REJECT reject.vcf --RECOVER_SWAPPED_REF_ALT

# Filter out snps that don't liftover
rejected_vcf <- read_tsv("~/Data/Fingerprinting/output_75_20/reject.vcf",
                                 col_names = c("chrom","pos","name","major","minor","qual","filter","info"))

map %<>% filter(!name %in% rejected_vcf$name & !anchor %in% rejected_vcf$name)

map %>% mutate(anchor = iflese(is.na(anchor),"",anchor)) %>% 
  write_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg19",col_names = F)

# Convert updated map to vcf
map %>% mutate(anchor = ifelse(is.na(anchor),"none",anchor)) %>% 
  transmute(chrom,pos,name,major,minor,qual = "100",filter = "PASS",info = str_c("AF=",maf,";","ANCHOR=",anchor)) %>% 
  write_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg19.vcf",col_names = F)

# Add hg19 header to ccle_29k_coding_snps_hg19.vcf

# Run java -Xmx6g -jar ~/Library/gatk/gatk-package-4.2.0.0-local.jar LiftoverVcf --INPUT ccle_29k_coding_snps_hg19.vcf --OUTPUT ccle_29k_coding_snps_hg38.vcf  --CHAIN ~/Data/VCFs/Liftover/b37ToHg38.over.chain -REFERENCE_SEQUENCE ~/Data/VCFs/Liftover/hg38.fa --REJECT reject.vcf --RECOVER_SWAPPED_REF_ALT

# Convert hg38 vcf back to map format
vcf <- read_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg38.vcf",
                col_names = c("chrom","pos","name","major","minor","qual","filter","info"))

map_38 <- vcf %>% transmute(chrom,pos,name,major,minor,maf = str_split_fixed(info,"=|;",5)[,2],
                            anchor = str_split_fixed(info,"=|;",5)[,4])

map_38 %>% mutate(anchor = ifelse(anchor == "none","",anchor)) %>% write_tsv("~/Data/Fingerprinting/output_75_20/ccle_29k_coding_snps_hg38",col_names = F)

                  