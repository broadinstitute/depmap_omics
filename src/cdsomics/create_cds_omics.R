
library(celllinemapr)
library(magrittr)
library(taigr)
library(cdsomics)
library(readr)
data_release <- '19Q3'

## public
# cell lines that were internal 6 months ago - in this case the 19Q1 release
RNAseq_internal_past_tpm <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=15, data.file='CCLE_depMap_19Q1_TPM')
RNAseq_internal_past_tpm_protein_coding <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=15, data.file='CCLE_depMap_19Q1_TPM_ProteinCoding')
RNAseq_internal_past_reads <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=15, data.file='CCLE_depMap_19Q1_RNAseq_reads')
RNAseq_internal_past_transcripts <- load.from.taiga(data.name='depmap-rnaseq-expression-data-363a', data.version=15, data.file='CCLE_depMap_19Q1_TPM_transcripts')

CN_gene_internal_past <- load.from.taiga(data.name='depmap-wes-cn-data-81a7', data.version=9, data.file='internal_19Q1_gene_cn')
CN_seg_internal_past <- load.from.taiga(data.name='segmented-cn-wes-prioritzed-7fe1', data.version=20, data.file='wes_priority_cn_seg_profiles')

mutation_internal_past <- load.from.taiga(data.name='depmap-mutation-calls-9be3', data.version=12, data.file='depmap_19Q1_mutation_calls')


RNAseq_embargo <- read_csv(paste0("~/data_files/", data_release, "/RNAseq_Embargo.csv"))
WES_embargo <- read_csv(paste0("~/data_files/", data_release, "/WES_Embargo.csv"))
blacklist <- read_csv(paste0("~/data_files/", data_release, "/blacklist.csv"))

RNAseq_public_cell_lines <- rownames(RNAseq_internal_past_transcripts)
RNAseq_public_cell_lines <- setdiff(RNAseq_public_cell_lines, RNAseq_embargo$IBM)
RNAseq_public_cell_lines <- setdiff(RNAseq_public_cell_lines, RNAseq_embargo$DMC)
RNAseq_public_cell_lines <- setdiff(RNAseq_public_cell_lines, blacklist$Blacklist)

WES_public_cell_lines <- rownames(CN_gene_internal_past)
WES_public_cell_lines <- setdiff(WES_public_cell_lines, blacklist$Blacklist)
WES_public_replace <-  WES_embargo$DMC

mutation_public_cell_lines <- unique(mutation_internal_past$DepMap_ID)
mutation_public_cell_lines <- setdiff(mutation_public_cell_lines, WES_embargo$DMC)
mutation_public_cell_lines <- setdiff(mutation_public_cell_lines, blacklist$Blacklist)

# Expression - internal and public is hg38

## internal

### TPM

TPM_mat <- load.from.taiga(data.name='depmap-expression-87f8', data.version=10, data.file='expression.19Q3.genes')
print(colnames(TPM_mat)[1:10])
TPM_mat <- TPM_mat[,2:ncol(TPM_mat)]
tpm_gene_ids <- gsub("\\..*", "", TPM_mat$gene_id)
if(nrow(TPM_mat) != length(unique(tpm_gene_ids))) {
  print("Duplicate ensembl ids")
  print(nrow(TPM_mat) - length(unique(tpm_gene_ids)))
  TPM_mat <- TPM_mat[-which(duplicated(tpm_gene_ids)==T),]
}

TPM <- prepare_depmap_TPM_for_taiga(TPM_mat, log_transform=T, just_protein_coding=F, 
                                    gencode_file = "~/data_files/gene_annotations/gencode.v29.annotation.gff3") 

remove_CL_from_internal <- intersect(rownames(TPM), blacklist$Blacklist)

TPM_protein_coding <- prepare_depmap_TPM_for_taiga(TPM_mat, log_transform=T, just_protein_coding=T,
                                                   gencode_file = "~/data_files/gene_annotations/gencode.v29.annotation.gff3") 

transcripts_mat <- load.from.taiga(data.name='depmap-expression-87f8', data.version=10, data.file='expression.19Q3.transcripts')
print(colnames(transcripts_mat)[1:10])
transcripts_mat <- transcripts_mat[,2:ncol(transcripts_mat)]
colnames(transcripts_mat)[1] <- "transcript_id(s)"
tpm_transcripts <- gsub("\\..*", "", transcripts_mat$`transcript_id(s)`)
if(nrow(transcripts_mat) != length(unique(tpm_transcripts))) {
  print("Duplicate transcript ids")
  print(nrow(transcripts_mat) - length(unique(tpm_transcripts)))
  transcripts_mat <- transcripts_mat[-which(duplicated(tpm_transcripts)==T),]
}
TPM_transcripts <- prepare_depmap_transcripts_for_taiga(transcripts_mat, gencode_file = "~/data_files/gene_annotations/gencode.v29.annotation.gff3") 


### counts
counts_mat <- load.from.taiga(data.name='depmap-expression-87f8', data.version=10, data.file='expression.19Q3.counts')
print(colnames(counts_mat)[1:10])
counts_mat <- counts_mat[,2:ncol(counts_mat)]
counts_gene_ids <- gsub("\\..*", "", counts_mat$gene_id)
if(nrow(counts_mat) != length(unique(counts_gene_ids))) {
  print("Duplicate ensembl ids")
  print(length(which(duplicated(counts_gene_ids)==T)))
  counts_mat <- counts_mat[-which(duplicated(counts_gene_ids)==T),]
}
read_counts <- prepare_depmap_TPM_for_taiga(counts_mat, gencode_file = "~/data_files/gene_annotations/gencode.v29.annotation.gff3") 

if(length(remove_CL_from_internal) > 0) {
  TPM <- TPM[-which(rownames(TPM) %in% remove_CL_from_internal),]
  TPM_protein_coding <- TPM_protein_coding[-which(rownames(TPM_protein_coding) %in% remove_CL_from_internal),]
  TPM_transcripts <- TPM_transcripts[-which(rownames(TPM_transcripts) %in% remove_CL_from_internal),]
  read_counts <- read_counts[-which(rownames(read_counts) %in% remove_CL_from_internal),]
  
}
# public
length(setdiff(RNAseq_public_cell_lines, rownames(TPM)))

TPM_public <- TPM[RNAseq_public_cell_lines,]
TPM_transcripts_public <- TPM_transcripts[RNAseq_public_cell_lines,]
TPM_protein_coding_public <- TPM_protein_coding[RNAseq_public_cell_lines,]
read_counts_public <- read_counts[RNAseq_public_cell_lines,]


# DMC
DMC_cell_lines <- rownames(TPM_protein_coding)
DMC_cell_lines <- setdiff(DMC_cell_lines, RNAseq_embargo$IBM)
DMC_cell_lines <- setdiff(DMC_cell_lines, blacklist$Blacklist)

TPM_DMC <- TPM[DMC_cell_lines,]
TPM_transcripts_DMC <- TPM_transcripts[DMC_cell_lines,]
TPM_protein_coding_DMC <- TPM_protein_coding[DMC_cell_lines,]
read_counts_DMC <- read_counts[DMC_cell_lines,]



# Copy Number - internal is hg38, public is hg38 liftover

## internal
WES_CN_seg <- load.from.taiga(data.name='segmented-cn-wes-prioritzed-7fe1', data.version=33, data.file='wes.19Q3.segmented')
WES_CN_gene <- load.from.taiga(data.name='segmented-cn-wes-prioritzed-7fe1', data.version=33, data.file='wes.19Q3.gene')
WES_CN_gene <- log2(WES_CN_gene + 1)

CN_internal_remove_lines <- intersect(rownames(WES_CN_gene), blacklist$Blacklist)
if(length(CN_internal_remove_lines) > 0) {
  WES_CN_seg <- filter(WES_CN_seg, !DepMap_ID %in%  CN_internal_remove_lines)
  WES_CN_gene <- WES_CN_gene[-which(rownames(WES_CN_gene) %in% CN_internal_remove_lines),]
  
}

## public
# need to pull from here b/c we have changed transformation
CN_gene_internal_past <- load.from.taiga(data.name='segmented-cn-wes-prioritzed-7fe1', data.version=20, data.file='wes_priority_cn_gene_matrix')
CN_gene_internal_past <- log2(CN_gene_internal_past + 1)
WES_CN_gene_public <- CN_gene_internal_past[WES_public_cell_lines,]
WES_CN_seg_public <- filter(CN_seg_internal_past, DepMap_ID %in% WES_public_cell_lines)

WES_public_replace <- intersect(WES_public_replace, rownames(WES_CN_gene_public))
CN_CL_ind <- which(WES_CN_seg_public$DepMap_ID %in% WES_public_replace & WES_CN_seg_public$Source =='Broad WES')
CN_public_replace <- unique(WES_CN_seg_public$DepMap_ID[CN_CL_ind])
paste("replacing CN data for public dataset for", length(CN_public_replace), "cell lines")


if(length(CN_public_replace) > 0) {
  # CN priority (1) Chordoma, (2) Sanger, (3) SNP
  Chordoma.HG38.genes <- load.from.taiga(data.name='snp-sanger-chordoma-d62a', data.version=3, data.file='Chordoma.HG38.genes')
  Sanger.HG38.genes <- load.from.taiga(data.name='snp-sanger-chordoma-d62a', data.version=3, data.file='Sanger.HG38.genes')
  SNP.HG38.genes <- load.from.taiga(data.name='snp-sanger-chordoma-d62a', data.version=3, data.file='SNP.HG38.genes')
  
  rownames(Chordoma.HG38.genes) <- paste0(Chordoma.HG38.genes$SYMBOL, " (", Chordoma.HG38.genes$EGID, ")")
  rownames(Sanger.HG38.genes) <- paste0(Chordoma.HG38.genes$SYMBOL, " (", Chordoma.HG38.genes$EGID, ")")
  rownames(SNP.HG38.genes) <- paste0(Chordoma.HG38.genes$SYMBOL, " (", Chordoma.HG38.genes$EGID, ")")
  
  Chordoma.HG38.genes <- Chordoma.HG38.genes %>% 
    dplyr::select(-EGID, -SYMBOL, -CHR) %>%
    as.matrix() %>%
    t() 
  
  Sanger.HG38.genes <- Sanger.HG38.genes %>% 
    dplyr::select(-EGID, -SYMBOL, -CHR) %>%
    as.matrix() %>%
    t() 
  
  SNP.HG38.genes <- SNP.HG38.genes %>% 
    dplyr::select(-EGID, -SYMBOL, -CHR) %>%
    as.matrix() %>%
    t() 
  
  Chordoma.HG38.genes <- Chordoma.HG38.genes[which(rownames(Chordoma.HG38.genes) %in% celllinemapr::arxspan.to.ccle(rownames(WES_CN))),]
  Sanger.HG38.genes <- Sanger.HG38.genes[which(rownames(Sanger.HG38.genes) %in% celllinemapr::arxspan.to.ccle(rownames(WES_CN))),]
  SNP.HG38.genes <- SNP.HG38.genes[which(rownames(SNP.HG38.genes) %in% celllinemapr::arxspan.to.ccle(rownames(WES_CN))),]
  
  rownames(Chordoma.HG38.genes) <- celllinemapr::ccle.to.arxspan(rownames(Chordoma.HG38.genes))
  rownames(Sanger.HG38.genes) <- celllinemapr::ccle.to.arxspan(rownames(Sanger.HG38.genes))
  SNP.HG38.genes <- SNP.HG38.genes[-which(rownames(SNP.HG38.genes) == "UMRC6_KIDNEY"),]
  rownames(SNP.HG38.genes) <- celllinemapr::ccle.to.arxspan(rownames(SNP.HG38.genes))
  
  length(intersect(colnames(WES_CN), colnames(Chordoma.HG38.genes)))
  length(intersect(colnames(WES_CN), colnames(Sanger.HG38.genes)))
  length(intersect(colnames(WES_CN), colnames(SNP.HG38.genes)))
  
  WES_CN_public <- subset_gene_level_cn(WES_CN_public, Chordoma.HG38.genes, Sanger.HG38.genes, SNP.HG38.genes, CN_public_replace)
}

#DMC

DMC_WES_cell_lines <- rownames(WES_CN_gene)
DMC_WES_cell_lines <- setdiff(DMC_WES_cell_lines, blacklist$Blacklist)

WES_CN_gene_DMC <- WES_CN_gene[DMC_WES_cell_lines,]

WES_CN_seg <- WES_CN_seg[,2:ncol(WES_CN_seg)]
WES_CN_seg_DMC <- filter(WES_CN_seg, DepMap_ID %in% DMC_WES_cell_lines)


# Mutation - hg19 still

## internal
mutation_maf <- load.from.taiga(data.name='depmap-mutations-maf-35fe', data.version=9, data.file='mutations.19Q3')
mutations <- maf_add_variant_annotations(mutation_maf)
mutations <- dplyr::filter(mutations, !(DepMap_ID %in% blacklist$Blacklist))
damaging_mutation <- mutation_maf_to_binary_matrix(mutations, damaging =  TRUE)
other_mutation <- mutation_maf_to_binary_matrix(mutations, other = TRUE)
hotspot_mutation <- mutation_maf_to_binary_matrix(mutations, hotspot = TRUE)

## public
mutations_public <- filter(mutations, DepMap_ID %in% mutation_public_cell_lines)
damaging_mutation_public <- mutation_maf_to_binary_matrix(mutations_public, damaging =  TRUE)
other_mutation_public <- mutation_maf_to_binary_matrix(mutations_public, other = TRUE)
hotspot_mutation_public <- mutation_maf_to_binary_matrix(mutations_public, hotspot = TRUE)

## DMC
DMC_mutation_cell_lines <- unique(mutations$DepMap_ID)
DMC_mutation_cell_lines <- setdiff(DMC_mutation_cell_lines, blacklist$Blacklist)

mutations_DMC <- filter(mutations, DepMap_ID %in% DMC_mutation_cell_lines)
damaging_mutation_DMC <- mutation_maf_to_binary_matrix(mutations_DMC, damaging =  TRUE)
other_mutation_DMC <- mutation_maf_to_binary_matrix(mutations_DMC, other = TRUE)
hotspot_mutation_DMC <- mutation_maf_to_binary_matrix(mutations_DMC, hotspot = TRUE)


# Fusion data - internal and public are hg38

#internal
fusions_filtered <- load.from.taiga(data.name='depmap-fusions-7990', data.version=6, data.file='fusions.19Q3.filtered')
fusions_unfiltered <- load.from.taiga(data.name='depmap-fusions-7990', data.version=6, data.file='fusions.19Q3.unfiltered')

filtered_fusions <- prepare_depmap_fusion_data_for_taiga(fusions_filtered)
unfiltered_fusions <- prepare_depmap_fusion_data_for_taiga(fusions_unfiltered)

filtered_fusions <- filtered_fusions[,2:ncol(filtered_fusions)]
unfiltered_fusions <- unfiltered_fusions[,2:ncol(unfiltered_fusions)]

filtered_fusions <- dplyr::filter(filtered_fusions, !(DepMap_ID %in% blacklist$Blacklist))
unfiltered_fusions <- dplyr::filter(unfiltered_fusions, !(DepMap_ID %in% blacklist$Blacklist))

#public
filtered_fusions_public <- dplyr::filter(filtered_fusions, DepMap_ID %in% RNAseq_public_cell_lines)
unfiltered_fusions_public <- dplyr::filter(unfiltered_fusions, DepMap_ID %in% RNAseq_public_cell_lines)

#DMC
fusion_CLs <- union(unique(filtered_fusions$DepMap_ID), unique(unfiltered_fusions$DepMap_ID))
DMC_fusion_cell_lines <- setdiff(fusion_CLs, RNAseq_embargo$IBM)
DMC_fusion_cell_lines <- setdiff(DMC_fusion_cell_lines, blacklist$Blacklist)

filtered_fusions_DMC <- dplyr::filter(filtered_fusions, DepMap_ID %in% DMC_fusion_cell_lines)
unfiltered_fusions_DMC <- dplyr::filter(unfiltered_fusions, DepMap_ID %in% DMC_fusion_cell_lines)


# ssGSEA
#ssGSEA <- load.from.taiga(data.name='ssgsea-c39a', data.version=1, data.file='depMap_19Q1_ssGSEA')

ssGSEA <- cdsomics::make_ssGSEA(TPM_protein_coding_public, "~/data_files/19Q1/msigdb.v6.2.entrez.gmt")

DMC_cell_lines <- setdiff(DMC_cell_lines, RNAseq_embargo$IBM)
ssGSEA_DMC <- ssGSEA[DMC_cell_lines,]

public_gsea_lines <- intersect(rownames(ssGSEA), RNAseq_public_cell_lines)

ssGSEA_public <- ssGSEA[public_gsea_lines,]

# Write files
write.csv(as.data.frame(TPM), paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM.csv"))
write.csv(as.data.frame(TPM_protein_coding), paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"))
write.csv(as.data.frame(TPM_transcripts), paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_TPM_transcripts.csv"))
write.csv(as.data.frame(read_counts), paste0("~/data_files/", data_release, "/internal/CCLE_depMap_", data_release, "_RNAseq_reads.csv"))

write.csv(as.data.frame(TPM_public), paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM.csv"))
write.csv(as.data.frame(TPM_protein_coding_public), paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"))
write.csv(as.data.frame(TPM_transcripts_public), paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_TPM_transcripts.csv"))
write.csv(as.data.frame(read_counts_public), paste0("~/data_files/", data_release, "/public/CCLE_depMap_", data_release, "_RNAseq_reads.csv"))

write.csv(as.data.frame(TPM_DMC), paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM.csv"))
write.csv(as.data.frame(TPM_protein_coding_DMC), paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM_ProteinCoding.csv"))
write.csv(as.data.frame(TPM_transcripts_DMC), paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_TPM_transcripts.csv"))
write.csv(as.data.frame(read_counts_DMC), paste0("~/data_files/", data_release, "/DMC/CCLE_depMap_", data_release, "_RNAseq_reads.csv"))

write.csv(as.data.frame(damaging_mutation), paste0("~/data_files/", data_release, "/internal/damaging_mutation.csv"))
write.csv(as.data.frame(other_mutation), paste0("~/data_files/", data_release, "/internal/other_mutation.csv"))
write.csv(as.data.frame(hotspot_mutation), paste0("~/data_files/", data_release, "/internal/hotspot_mutation.csv"))
write_csv(as.data.frame(mutations), paste0("~/data_files/", data_release, "/internal/depmap_", data_release, "_mutation_calls.csv"))

write.csv(as.data.frame(damaging_mutation_public), paste0("~/data_files/", data_release, "/public/damaging_mutation.csv"))
write.csv(as.data.frame(other_mutation_public), paste0("~/data_files/", data_release, "/public/other_mutation.csv"))
write.csv(as.data.frame(hotspot_mutation_public), paste0("~/data_files/", data_release, "/public/hotspot_mutation.csv"))
write_csv(as.data.frame(mutations_public), paste0("~/data_files/", data_release, "/public/depmap_", data_release, "_mutation_calls.csv"))

write.csv(as.data.frame(damaging_mutation_DMC), paste0("~/data_files/", data_release, "/DMC/damaging_mutation.csv"))
write.csv(as.data.frame(other_mutation_DMC), paste0("~/data_files/", data_release, "/DMC/other_mutation.csv"))
write.csv(as.data.frame(hotspot_mutation_DMC), paste0("~/data_files/", data_release, "/DMC/hotspot_mutation.csv"))
write_csv(as.data.frame(mutations_DMC), paste0("~/data_files/", data_release, "/DMC/depmap_", data_release, "_mutation_calls.csv"))

write.csv(as.data.frame(filtered_fusions), paste0("~/data_files/", data_release, "/internal/filtered_fusions_", data_release, ".csv"))
write.csv(as.data.frame(unfiltered_fusions), paste0("~/data_files/", data_release, "/internal/unfiltered_fusions_", data_release, ".csv"))

write.csv(as.data.frame(filtered_fusions_public), paste0("~/data_files/", data_release, "/public/filtered_fusions_", data_release, ".csv"))
write.csv(as.data.frame(unfiltered_fusions_public), paste0("~/data_files/", data_release, "/public/unfiltered_fusions_", data_release, ".csv"))

write.csv(as.data.frame(filtered_fusions_DMC), paste0("~/data_files/", data_release, "/DMC/filtered_fusions_", data_release, ".csv"))
write.csv(as.data.frame(unfiltered_fusions_DMC), paste0("~/data_files/", data_release, "/DMC/unfiltered_fusions_", data_release, ".csv"))


write.csv(as.data.frame(WES_CN_gene), paste0("~/data_files/", data_release, "/internal/internal_", data_release, "_gene_cn.csv"))
write.csv(as.data.frame(WES_CN_gene_public), paste0("~/data_files/", data_release, "/public/public_", data_release, "_gene_cn.csv"))
write.csv(as.data.frame(WES_CN_gene_DMC), paste0("~/data_files/", data_release, "/DMC/DMC_", data_release, "_gene_cn.csv"))

write_csv(as.data.frame(WES_CN_seg), paste0("~/data_files/", data_release, "/internal/internal_", data_release, "_cn_segs.csv"))
write_csv(as.data.frame(WES_CN_seg_public), paste0("~/data_files/", data_release, "/public/public_", data_release, "_cn_segs.csv"))
write_csv(as.data.frame(WES_CN_seg_DMC), paste0("~/data_files/", data_release, "/DMC/DMC_", data_release, "_cn_segs.csv"))

write.csv(as.data.frame(ssGSEA), paste0("~/data_files/", data_release, "/internal/depMap_", data_release, "_ssGSEA.csv"))
write.csv(as.data.frame(ssGSEA_DMC), paste0("~/data_files/", data_release, "/DMC/depMap_", data_release, "_ssGSEA.csv"))
write.csv(as.data.frame(ssGSEA_public), paste0("~/data_files/", data_release, "/public/depMap_", data_release, "_ssGSEA.csv"))


write.csv(as.data.frame(RRBS[1]), paste0("~/data_files/", data_release, "/CCLE_RRBS_TSS_1kb_20180614.csv"))
write.csv(as.data.frame(RRBS[2]), paste0("~/data_files/", data_release, "/CCLE_RRBS_TSS_1kb_info_20180614.csv"))
write.csv(as.data.frame(RRBS[3]), paste0("~/data_files/", data_release, "/CCLE_RRBS_TSS_CpG_clusters_20180614.csv"))
write.csv(as.data.frame(RRBS[4]), paste0("~/data_files/", data_release, "/CCLE_RRBS_TSS_CpG_clusters_info_20180614.csv"))
write.csv(as.data.frame(miRNA), paste0("~/data_files/", data_release, "/CCLE_miRNA_20180525.csv"))


# stats
#mutation_stats <- c(length(unique(mutations$Depmap_ID)), )

#miRNA <- prepare_miRNA_for_taiga(paste0("~/data_files/", data_release, "/CCLE_miRNA_20180525_raw.gct"))
RRBS <- prepare_RRBS_for_taiga("~/data_files/18Q3/internal/CCLE_RRBS_TSS_1kb_20180614_raw.txt", paste0("~/data_files/", data_release, "/internal/CCLE_RRBS_TSS_CpG_clusters_20180614_raw.txt"))


# metabolomics
metabolomics <- load.from.taiga(data.name='metabolomics-cd0c', data.version=2)
rownames(metabolomics) <- ccle.to.arxspan(rownames(metabolomics))

write.csv(as.data.frame(data), paste0("~/data_files/", data_release, "/metabolomics.csv"))

gene_mapping <- as.data.frame(read_csv("~/HGNC_to_Entrez.csv"))
#transcription activity
TF_activity_CLs <- as.data.frame(read_csv("~/TF_activity_in_cell_lines.csv"))
TF_activity_tumors <- as.data.frame(read_csv("~/TF_activity_in_tumors.csv"))
Sanger_CL_mapping <- load.from.taiga(data.name='sanger-cell-line-mapping-adba', data.version=1)

CL_mapping <- as.data.frame(colnames(TF_activity_CLs)[2:ncol(TF_activity_CLs)])
CCLE_map <- character()
Broad_IDs <- character()
for(i in 1:nrow(CL_mapping)) {
  if(CL_mapping[i,1] %in% Sanger_CL_mapping$`GDSC1000 cosmic id`) {
    CCLE_name <- Sanger_CL_mapping$`CCLE name`[grep(CL_mapping[i,1], Sanger_CL_mapping$`GDSC1000 cosmic id`)]
    if(CCLE_name == "N/A") {
      CCLE_name <- NA
    }
    CCLE_map <- c(CCLE_map, CCLE_name)
    mapping_ind <- grep(CCLE_name, mapping$ccle_name)
    if(length(mapping_ind) < 1) {
      mapping_ind <- grep(CCLE_name, mapping$canonical_ccle_name)
    }
    if(length(mapping_ind) > 1) {
      mapping_ind <- mapping_ind[1]
    }
    broad_id <- mapping$broad_id[mapping_ind]
    if(length(broad_id) < 1) {
      broad_id <- NA
    }
    Broad_IDs <- c(Broad_IDs, broad_id)
    
  } else {
    CCLE_map <- c(CCLE_map, NA)
    Broad_IDs <- c(Broad_IDs, NA)
  }
} 

CL_mapping$CCLE_name <- CCLE_map
CL_mapping$Broad_ID <- Broad_IDs
CL_mapping$CCLE_ID[924] <- "RKN_SOFT_TISSUE"
colnames(CL_mapping)[1:2] <- c("Sanger_ID", "CCLE_ID")

entrez_gene_ids <- gene_mapping[TF_activity_CLs$Gene,]$Entrez_ID
fusion_gene <- TF_activity_CLs$Gene[which(is.na(entrez_gene_ids)==T)]
fusion_entrez_id <- paste0(gene_mapping$Entrez_ID[which(gene_mapping$HGNC_symbol == "EWSR1")], "-", gene_mapping$Entrez_ID[which(gene_mapping$HGNC_symbol == "FLI1")])
entrez_gene_ids[which(is.na(entrez_gene_ids)==T)] <- fusion_entrez_id
TF_gene_names <- paste(TF_activity_CLs$Gene, paste0("(", entrez_gene_ids, ")"))
rownames(TF_activity_CLs) <- TF_gene_names

entrez_gene_ids2 <- gene_mapping[TF_activity_tumors$Gene,]$Entrez_ID
entrez_gene_ids2[which(is.na(entrez_gene_ids2)==T)] <- fusion_entrez_id
TF_gene_names2 <- paste(TF_activity_tumors$Gene, paste0("(", entrez_gene_ids2, ")"))
rownames(TF_activity_tumors) <- TF_gene_names2

TF_activity_CLs <- TF_activity_CLs[,2:ncol(TF_activity_CLs)]
TF_activity_tumors <- TF_activity_tumors[,2:ncol(TF_activity_tumors)]
write.csv(as.data.frame(t(TF_activity_CLs)), paste0("~/data_files/", data_release, "/cell_line_TF_activity.csv"))
write.csv(as.data.frame(t(TF_activity_tumors)), paste0("~/data_files/", data_release, "/tumor_TF_activity.csv"))
write_csv(as.data.frame(CL_mapping), paste0("~/data_files/", data_release, "/Sanger_CCLE_cell_line_mapping.csv"), col_names = T)

