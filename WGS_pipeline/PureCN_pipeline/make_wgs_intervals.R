library(magrittr)
library(tidyverse)

chr_lens <- read_tsv("~/Data/Copy_Number/example_wgs/chromosome_lengths.txt",col_names = c("chromosome","length"))
targets <- list()
for (chr in chr_lens$chromosome) {
  targets %<>% append(list(tibble(chromosome = chr,start = seq(1, 500000000, by=10000)) %>% mutate(end = start + 199)))
}
targets %<>% bind_rows()
targets %<>% left_join(chr_lens,by = "chromosome")
targets %<>% filter(end < length)

targets %>% write_tsv("~/Data/Copy_Number/example_wgs/wgs_hg38_intervals.bed",col_names  = F)

                   
