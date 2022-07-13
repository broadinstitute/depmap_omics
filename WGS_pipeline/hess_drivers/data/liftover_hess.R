library(tidyverse)

hess <- read_csv("~/Data/Mutect2/reference/mmc3.csv")
hess %<>% arrange(chr,pos) %>% mutate(id = row_number()) 
hess %>% transmute(CHROM = chr, POS = pos,ID = id,REF = str_sub(base_change,4,4),ALT = str_sub(base_change,7,7),
                   QUAL = ".", FILTER = "PASS",INFO = ".") %>% write_tsv("~/Desktop/hg19.vcf")

system(str_c("java -jar ~/Library/picard.jar LiftoverVcf -I ~/Desktop/hg19.vcf -O ~/Desktop/hg38.vcf",
             " -CHAIN ~/Data/VCFs/Liftover/b37ToHg38.over.chain -REJECT $TMPDIR/rejected_variants.vcf", 
             " -R ~/Data/VCFs/Liftover/hg38.fa --MAX_RECORDS_IN_RAM 10000"), intern=TRUE)

complement <- tibble(c = c("A","T","C","G"),a = c("T","A","G","C"))
hess %>% transmute(CHROM = chr, POS = pos,ID = id,ref = str_sub(base_change,4,4),alt = str_sub(base_change,7,7),
                   QUAL = ".", FILTER = "PASS",INFO = ".") %>% 
  left_join(complement %>% transmute(ref = a,REF = c)) %>% 
  left_join(complement %>% transmute(alt = a,ALT = c)) %>% 
  dplyr::select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO) %>% 
  write_tsv("~/Desktop/hg19_complement.vcf")

system(str_c("java -jar ~/Library/picard.jar LiftoverVcf -I ~/Desktop/hg19_complement.vcf -O ~/Desktop/hg38_complement.vcf",
             " -CHAIN ~/Data/VCFs/Liftover/b37ToHg38.over.chain -REJECT $TMPDIR/rejected_variants.vcf", 
             " -R ~/Data/VCFs/Liftover/hg38.fa --MAX_RECORDS_IN_RAM 10000"), intern=TRUE)

hg38 <- read_tsv("~/Desktop/hg38.vcf",comment = "##") %>% rename(CHROM = `#CHROM`)
hg38_complement <- read_tsv("~/Desktop/hg38_complement.vcf",comment = "##") %>% rename(CHROM = `#CHROM`)
hg38 %<>% bind_rows(hg38_complement)

test <- hess %>% left_join(hg38 %>% dplyr::select(CHROM,POS,REF,ALT,ID),by = c("id" = "ID"))
