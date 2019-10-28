# test new omics data to see that it looks reasonable/mostly matches the previous release

library(celllinemapr)
library(magrittr)
library(taigr)
library(cdsomics)
library(readr)
data_release <- '19Q3'
check_local <- T

RNAseq_embargo <- read_csv(paste0("~/data_files/", data_release, "/RNAseq_Embargo.csv"))
WES_embargo <- read_csv(paste0("~/data_files/", data_release, "/WES_Embargo.csv"))
blacklist <- read_csv(paste0("~/data_files/", data_release, "/blacklist.csv"))
restricted_cell_lines <- list(RNAseq_embargo = RNAseq_embargo, WES_embargo = WES_embargo, blacklist = blacklist)


# RNAseq ------------------------------------------
read_counts_previous <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=16, data.file='CCLE_depMap_19Q2_RNAseq_reads')
TPM_previous <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=16, data.file='CCLE_depMap_19Q2_TPM')
TPM_protein_coding_previous <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=16, data.file='CCLE_depMap_19Q2_TPM_ProteinCoding')
transcripts_previous <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=16, data.file='CCLE_depMap_19Q2_TPM_transcripts')

if(check_local) {
  TPM <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  TPM_protein_coding <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"),row.names=1, check.names = F, stringsAsFactors=FALSE))
  TPM_transcripts <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM_transcripts.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  read_counts <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_RNAseq_reads.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  read_counts <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=13, data.file='CCLE_depMap_19Q1_RNAseq_reads')
  TPM <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=13, data.file='CCLE_depMap_19Q1_TPM')
  TPM_protein_coding <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=13, data.file='CCLE_depMap_19Q1_TPM_ProteinCoding')
  TPM_transcripts <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=13, data.file='CCLE_depMap_19Q1_TPM_transcripts')
}
test_new_data(TPM_protein_coding, TPM_protein_coding_previous) 
test_new_data(TPM, TPM_previous) 
test_new_data(TPM_transcripts, transcripts_previous) 
test_new_data(read_counts, read_counts_previous) 

test_released_cell_lines(TPM_protein_coding, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal')
test_released_cell_lines(TPM, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal')
test_released_cell_lines(TPM_transcripts, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal')
test_released_cell_lines(read_counts, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal')

if(check_local) {
  TPM_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  TPM_protein_coding_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"),row.names=1, check.names = F, stringsAsFactors=FALSE))
  #TPM_transcripts_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM_transcripts.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  read_counts_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_RNAseq_reads.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  read_counts_DMC <- load.from.taiga(data.name='depmap-rnaseq-expression-data-80ef', data.version=1, data.file='CCLE_depMap_19Q1_RNAseq_reads')
  TPM_DMC <- load.from.taiga(data.name='depmap-rnaseq-expression-data-80ef', data.version=1, data.file='CCLE_depMap_19Q1_TPM')
  TPM_protein_coding_DMC <- load.from.taiga(data.name='depmap-rnaseq-expression-data-80ef', data.version=1, data.file='CCLE_depMap_19Q1_TPM_ProteinCoding')
  TPM_transcripts_DMC <- load.from.taiga(data.name='depmap-rnaseq-expression-data-80ef', data.version=1, data.file='CCLE_depMap_19Q1_TPM_transcripts')
}
test_new_data(TPM_DMC, TPM_previous) 
test_new_data(TPM_protein_coding_DMC, TPM_protein_coding_previous) 
test_new_data(TPM_transcripts_DMC, transcripts_previous) 
test_new_data(read_counts_DMC, read_counts_previous) 

test_released_cell_lines(TPM_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC')
test_released_cell_lines(TPM_protein_coding_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC')
test_released_cell_lines(TPM_transcripts_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC')
test_released_cell_lines(read_counts_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC')

TPM_previous_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=13, data.file='CCLE_depMap_19Q2_TPM')
TPM_protein_coding_previous_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=13, data.file='CCLE_depMap_19Q2_TPM_ProteinCoding')
transcripts_previous_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=13, data.file='CCLE_depMap_19Q2_TPM_transcripts')
read_counts_previous_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=13, data.file='CCLE_depMap_19Q2_RNAseq_reads')

if(check_local) {
  TPM_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  TPM_protein_coding_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"),row.names=1, check.names = F, stringsAsFactors=FALSE))
  #TPM_transcripts_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM_transcripts.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  read_counts_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_RNAseq_reads.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  read_counts_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=10, data.file='CCLE_depMap_19Q1_RNAseq_reads')
  TPM_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=10, data.file='CCLE_depMap_19Q1_TPM')
  TPM_protein_coding_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=10, data.file='CCLE_depMap_19Q1_TPM_ProteinCoding')
  TPM_transcripts_public <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=10, data.file='CCLE_depMap_19Q1_TPM_transcripts')
}
test_new_data(TPM_public, TPM_previous_public) 
test_new_data(TPM_protein_coding_public, TPM_protein_coding_previous_public) 
test_new_data(TPM_transcripts_public, transcripts_previous_public) 
test_new_data(read_counts_public, read_counts_previous_public) 

test_released_cell_lines(TPM_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public')
test_released_cell_lines(TPM_protein_coding_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public')
test_released_cell_lines(TPM_transcripts_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public')
test_released_cell_lines(read_counts_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public')


# CN ---------------------------------------------------------------

CN_previous <- load.from.taiga(data.name='depmap-wes-cn-data-81a7', data.version=11, data.file='internal_19Q2_gene_cn')

WES_CN <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/internal_", data_release, "_gene_cn.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE)) 
CN_segs <- read_csv(paste0("~/data_files/", data_release, "/internal/internal_", data_release, "_cn_segs.csv"))

#WES_CN <- load.from.taiga(data.name='depmap-wes-cn-data-81a7', data.version=8, data.file='internal_19Q1_gene_cn')

test_new_data(WES_CN, CN_previous, data_type='CN') 
test_released_cell_lines(WES_CN, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'internal')
length(intersect(blacklist$Blacklist, CN_segs$DepMap_ID))

WES_CN_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/DMC_", data_release, "_gene_cn.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
CN_segs_DMC <- read_csv(paste0("~/data_files/", data_release, "/DMC/DMC_", data_release, "_cn_segs.csv"))

#WES_CN_DMC <- load.from.taiga(data.name='depmap-cn-data-9b9d', data.version=1, data.file='DMC_19Q1_gene_cn')

test_new_data(WES_CN_DMC, CN_previous, data_type='CN') 
test_released_cell_lines(WES_CN_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'DMC')
length(intersect(blacklist$Blacklist, CN_segs_DMC$DepMap_ID))

CN_previous_public <- load.from.taiga(data.name='depmap-wes-cn-data-97cc', data.version=16, data.file='public_19Q2_cn_segs')
CN_previous_internal <- load.from.taiga(data.name='depmap-wes-cn-data-81a7', data.version=9, data.file='internal_19Q1_gene_cn')

WES_CN_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/public_", data_release, "_gene_cn.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
CN_segs_public <- read_csv(paste0("~/data_files/", data_release, "/public/public_", data_release, "_cn_segs.csv"))
#WES_CN_public <- load.from.taiga(data.name='depmap-wes-cn-data-97cc', data.version=13, data.file='public_19Q1_gene_cn')

test_new_data(WES_CN_public, CN_previous_public, data_type='CN') 
test_new_data(WES_CN_public, CN_previous_internal, data_type='CN') 
test_released_cell_lines(WES_CN_public, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'public', is_CN=T)
length(intersect(blacklist$Blacklist, CN_segs_public$DepMap_ID))
length(intersect(WES_embargo$DMC, CN_segs_public$DepMap_ID))


# mutation ------------------------------------------------------------
mutations_previous <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='depmap_19Q2_mutation_calls')
damaging.mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='damaging_mutation')
hotspot.mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='hotspot_mutation')
other.mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='other_mutation')


if(check_local) {
  mutations <- as.data.frame(read.csv(paste0("~/data_files/", data_release, "/internal/depmap_", data_release, "_mutation_calls.csv")))
  damaging_mutation <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/damaging_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  other_mutation <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/other_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  hotspot_mutation <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/internal/hotspot_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  damaging_mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='damaging_mutation')
  mutations <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='depmap_19Q2_mutation_calls')
  hotspot_mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='hotspot_mutation')
  other_mutation <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=13, data.file='other_mutation')
  
}

length(unique(mutations$DepMap_ID))
length(unique(mutations_previous$DepMap_ID))

test_new_data(damaging_mutation, damaging.mutation) 
test_new_data(other_mutation, other.mutation) 
test_new_data(hotspot_mutation, hotspot.mutation) 

test_released_cell_lines(mutations, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'internal', cell_line_column = 'DepMap_ID')
test_released_cell_lines(damaging_mutation, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'internal')
test_released_cell_lines(other_mutation, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'internal')
test_released_cell_lines(hotspot_mutation, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'internal')

if(check_local) {
  mutations_DMC <- as.data.frame(read.csv(paste0("~/data_files/", data_release, "/DMC/depmap_", data_release, "_mutation_calls.csv")))
  damaging_mutation_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/damaging_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  other_mutation_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/other_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  hotspot_mutation_DMC <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/DMC/hotspot_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  damaging_mutation_DMC <- load.from.taiga(data.name='depmap-mutation-calls-dfce', data.version=6, data.file='damaging_mutation')
  mutations_DMC <- load.from.taiga(data.name='depmap-mutation-calls-dfce', data.version=6, data.file='depmap_19Q2_mutation_calls')
  hotspot_mutation_DMC <- load.from.taiga(data.name='depmap-mutation-calls-dfce', data.version=6, data.file='hotspot_mutation')
  other_mutation_DMC <- load.from.taiga(data.name='depmap-mutation-calls-dfce', data.version=6, data.file='other_mutation')
}

length(unique(mutations_DMC$DepMap_ID))

test_new_data(damaging_mutation_DMC, damaging.mutation) 
test_new_data(other_mutation_DMC, other.mutation) 
test_new_data(hotspot_mutation_DMC, hotspot.mutation) 

test_released_cell_lines(mutations_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'DMC', cell_line_column = 'DepMap_ID')
test_released_cell_lines(damaging_mutation_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'DMC')
test_released_cell_lines(other_mutation_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'DMC')
test_released_cell_lines(hotspot_mutation_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'DMC')

mutations_previous_public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='depmap_19Q2_mutation_calls')
damaging.mutation.public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='damaging_mutation')
hotspot.mutation.public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='hotspot_mutation')
other.mutation.public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='other_mutation')

if(check_local) {
  mutations_public <- read.csv(paste0("~/data_files/", data_release, "/public/depmap_", data_release, "_mutation_calls.csv"))
  damaging_mutation_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/damaging_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  other_mutation_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/other_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
  hotspot_mutation_public <- as.matrix(read.csv(paste0("~/data_files/", data_release, "/public/hotspot_mutation.csv"), row.names=1, check.names = F, stringsAsFactors=FALSE))
} else {
  damaging_mutation_public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='damaging_mutation')
  mutations_public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='depmap_19Q2_mutation_calls')
  hotspot_mutation_public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='hotspot_mutation')
  other_mutation_public <- load.from.taiga(data.name='depmap-mutation-calls-9a1a', data.version=10, data.file='other_mutation')
}

length(unique(mutations_public$DepMap_ID))
length(unique(mutations_previous_public$DepMap_ID))

test_new_data(damaging_mutation_public, damaging.mutation.public) 
test_new_data(other_mutation_public, other.mutation.public) 
test_new_data(hotspot_mutation_public, hotspot.mutation.public) 


test_released_cell_lines(mutations_public, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'public', cell_line_column = 'DepMap_ID')
test_released_cell_lines(damaging_mutation_public, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'public')
test_released_cell_lines(other_mutation_public, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'public')
test_released_cell_lines(hotspot_mutation_public, restricted_cell_lines = restricted_cell_lines, data_type = 'WES', dataset = 'public')

# fusions ----------------------------------------------------------------
filtered_fusions_previous_public <- load.from.taiga(data.name='gene-fusions-6212', data.version=6, data.file='filtered_fusions_19Q2')
unfiltered_fusions_previous_public <- load.from.taiga(data.name='gene-fusions-6212', data.version=6, data.file='unfiltered_fusions_19Q2')

fusions <- load.from.taiga(data.name='gene-fusions-8b7a', data.version=6, data.file='filtered_fusions_19Q2')
unfusions <- load.from.taiga(data.name='gene-fusions-8b7a', data.version=6, data.file='unfiltered_fusions_19Q2')
# fusions <- read_csv(paste0("~/data_files/", data_release, "/internal/filtered_fusions_19Q1.csv"))
# unfusions <- read_csv(paste0("~/data_files/", data_release, "/internal/unfiltered_fusions_19Q1.csv")) 

length(unique(fusions$DepMap_ID))
test_released_cell_lines(fusions, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal', cell_line_column = 'DepMap_ID')
test_released_cell_lines(unfusions, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'internal', cell_line_column = 'DepMap_ID')

#fusions_DMC <- read_csv(paste0("~/data_files/", data_release, "/DMC/filtered_fusions_19Q1.csv")) 
filtered_fusions_DMC <- load.from.taiga(data.name='gene-fusions-375f', data.version=1, data.file='filtered_fusions_19Q1')
unfiltered_unfusions_DMC <- load.from.taiga(data.name='gene-fusions-375f', data.version=1, data.file='unfiltered_fusions_19Q1')

length(unique(filtered_fusions_DMC$DepMap_ID))
test_released_cell_lines(filtered_fusions_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC', cell_line_column = 'DepMap_ID')
test_released_cell_lines(unfiltered_fusions_DMC, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'DMC', cell_line_column = 'DepMap_ID')

fusions_public <- read_csv(paste0("~/data_files/", data_release, "/public/filtered_fusions_", data_release, ".csv")) 
unfiltered_fusions_public <- read_csv(paste0("~/data_files/", data_release, "/public/unfiltered_fusions_", data_release, ".csv")) 
unfiltered_fusions_public <- load.from.taiga(data.name='gene-fusions-6212', data.version=2, data.file='filtered_fusions_19Q1')
unfusions_public <- load.from.taiga(data.name='gene-fusions-6212', data.version=2, data.file='unfiltered_fusions_19Q1')

length(unique(unfiltered_fusions_public$DepMap_ID))
test_released_cell_lines(filtered_fusions_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public', cell_line_column = 'DepMap_ID')
test_released_cell_lines(unfiltered_fusions_public, restricted_cell_lines = restricted_cell_lines, data_type = 'RNAseq', dataset = 'public', cell_line_column = 'DepMap_ID')

length(intersect(unfiltered_fusions_previous_public$DepMap_ID, unfiltered_fusions_public$DepMap_ID))
