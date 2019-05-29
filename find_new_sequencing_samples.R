library(googlesheets)

# Pull the MF omics tab
gap <- gs_title("Master_DepMap_PV")
omics <- gap %>% gs_read('Omics')

## WES samples in DepMap (CCLE and Sanger included)

latest_wes_samples_all <- read_tsv('~/Downloads/sample_metadata_broad-firecloud-ccle_DepMap_WES_CN_hg38.tsv', col_types = cols(.default = 'c')) %>% 
  dplyr::select(sample=`entity:sample_id`, GCP_path=WES_bam)

# We only want to ID new Broad samples, so remove Sanger samples
broad_only_samples <- latest_wes_samples_all %>%
  filter(!grepl('^sanger_', sample)) %>%
  mutate(DepMap_ID=stringr::str_extract(pattern='ACH\\-[0-9]+', string=sample))

# Find samples that are not in the DepMap data
omics %>%
  dplyr::rename(DepMap_ID=`DepMap ID`) %>%
  filter(`WES Status` == 'Complete', `WES Owner` != 'Sanger') %>%
  left_join(., broad_only_samples, by='DepMap_ID') %>%
  dplyr::select(`DepMap_ID`, CCLE_name, `WES Status`, `WES Owner`, `WES Data Links`, GCP_path, GCP_ID=sample) %>%
  filter(is.na(GCP_ID))
