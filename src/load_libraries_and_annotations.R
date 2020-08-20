### SET UP FOR MOST RMARKDOWN ###

install_all <- function(){
  install.packages(c('plyr','dplyr','tidyverse','magrittr','reshape2','useful','ggplot2','ggthemes',
    'ggrepel','gridExtra','ggridges','GGally','plotly','VennDiagram','RColorBrewer','extrafont','cowplot',
    'network','data.table','DT','readr','readxl','clues','mclust','pheatmap','Rtsne','NMF','hash',
     "stringr", "irr", "zoo", "devtools", "scales", "rlang", "rmarkdown","lsa" ))
  font_import()
  loadfonts()
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(c("GSEABase","limma","org.Hs.eg.db","GenomicRanges","DESeq2",'rtracklayer'))
  print("if can't use broad intranet, install from source with [R CMD INSTALL .] for 'taigr',\
    'cdsr','svacd', cell_line_mapping/celllinemapr")
}



# DATA FRAME MANIPULATION
library(plyr)
library(dplyr)
library(tidyverse)
library(magrittr)
library(reshape2)
library(useful) # For 'corner' function
# Other
#library(sva) # For ComBat
library(limma)
library(stringr)
library(irr)
library(zoo)
library(scales)
library(rlang)
library(rmarkdown)
library(knitr)
library(lsa)
# Plotting
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(gridExtra)
library(ggridges)
library(GGally)
library(plotly)
library(VennDiagram)
library(RColorBrewer)
library(extrafont)
library(cowplot)
# library(networkD3)
loadfonts()

# For tables
library(data.table)
library(DT)

# Reading in data
library(readr)
library(readxl)
library(GenomicRanges)

# GSEA
library(GSEABase)

# CDS
library(taigr)
library(cdsr)
library(celllinemapr)

# Clustering
library(clues)
library(mclust)
library(pheatmap)
library(Rtsne)
#library(NMF)

zmad <- function(x) {
  med_x <- median(x, na.rm = T)
  mad_x <- mad(x, constant = 0.6744897501960817, na.rm = T)
  
  return((x - med_x)/(mad_x))
}


# FIt-SNE implementation (functon to use is fftRtsne. Use fast_tsne_path indicated below)
# This was downloaded from here https://github.com/KlugerLab/FIt-SNE. 
# source('~/Documents/Analysis/RScripts/FIt-SNE/fast_tsne.R')
# fast_tsne_path <- '~/Documents/Analysis/RScripts/FIt-SNE/bin/fast_tsne'

# Define repeated kmeans function here
# ds: sample x feature matrix
# pca: whether or not to run PCA first
# num_pcs: number of PCs to keep
# iterations: number of iterations to run (number of kmeans)
# n_clusters: number of clusters to separate into to generate the similarity score
# plus_minus: if you want to vary the n_clusters to be within a random range (e.g, in range 190-210, then set n_clusers to 200 and plus_minus to 10)
repeated_k_means <- function(ds, pca=TRUE, num_pcs=200, iterations=1000, n_clusters=2, plus_minus=0) {
  # Run PCA
  for_clustering <- ds
  if (pca) {
    # Run PCA
    num_pcs <- min(num_pcs, ncol(ds))
    pca_res <- prcomp(ds)$x
    
    for_clustering <- pca_res[,paste0('PC', seq(1, num_pcs))]
  }
  
  # Run the clustering
  keep_count <- NULL
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  
  for (r in seq(1, iterations)) {
    increment_clust <- sample(seq((n_clusters-plus_minus), n_clusters+plus_minus), 1)
    clusters <- suppressWarnings(kmeans(as.matrix(for_clustering), increment_clust)$cluster)
    
    if (is.null(keep_count)) {
      keep_count <- data.frame(run_1=clusters)
    } else {
      keep_count[names(clusters),paste0('run_', r)] <- clusters[names(clusters)]
    }
    setTxtProgressBar(pb, r)
  }
  print('Completed clustering iterations. Aggregating results...')
  
  # Now calculate the similarity score
  repeated_kmeans_similarity <- apply(keep_count, 1, FUN = function(x) apply(keep_count, 1, FUN = function(y) length(which(x==y))))
  repeated_kmeans_similarity <- repeated_kmeans_similarity/iterations

  return(repeated_kmeans_similarity)
}



# Load MSIGDB gene sets and others
#gsc_data <- readr::read_rds(taigr::download.raw.from.taiga(data.name='msigdb-gene-set-collections-8453', data.version=2))

# Pull annotations (and later can format in one place if needed)
# This formats them based on the release master file
pull_mf <- function(sample.info.df) {
  ped_types <- c("Ewing", 'Medulloblastoma', 'Neuroblastoma', 'Osteosarcoma', 'Rhabdoid', 'Rhabdomyosarcoma', 'Synovial Sarcoma', 'B-cell ALL', 'Retinoblastoma') 

  sample.info.df %<>%
    set_colnames(gsub('disease_sutype', 'disease_subtype', colnames(.))) %>%
    filter(!is.na(DepMap_ID) & DepMap_ID != '') %>%
    # For some reason, this cell lines has no disease information...
    mutate(disease=ifelse(DepMap_ID=='ACH-000008' & (trimws(disease)=='' | is.na(disease)), 'Skin', disease)) %>%
    mutate(Type=case_when(
      disease_subtype %in% c('Ewing_sarcoma') ~ 'Ewing',
      disease_subtype %in% c('synovial_sarcoma') ~ 'Synovial Sarcoma',
      disease_subtype == 'B-cell_ALL' ~ 'B-cell ALL',
      toupper(disease_subtype) %in% toupper(ped_types) ~ stringr::str_to_title(gsub('_', ' ', disease_subtype)),
      disease_subtype %in% c('chordoma') ~ 'Chordoma',
      disease_subtype %in% c('uterus_endometrium') ~ 'Endometrium',
      TRUE ~ stringr::str_to_title(gsub('_', ' ', disease))
    )) %>%
    # Custom relabelling of cell lines
    mutate(Type=ifelse(DepMap_ID %in% c("ACH-000772", "ACH-000051", "ACH-000689", "ACH-001750"), 'Other', Type)) %>%
    mutate(Type=ifelse(DepMap_ID %in% c('ACH-001033'), "Rhabdoid", Type)) %>% # This might not be entirely accurate
    mutate(Type=ifelse(Type=='Medullo', 'Medulloblastoma', Type)) %>%
    mutate(T2=ifelse(
      DepMap_ID %in% c("ACH-001741", "ACH-001743", "ACH-000100", "ACH-001189", "ACH-001050", "ACH-001740", "ACH-000833", "ACH-001765", "ACH-001745", "ACH-001184"), 'alveolar', ifelse(
        DepMap_ID %in% c("ACH-001790", "ACH-001751", "ACH-001096", "ACH-000169","ACH-001196"), 'embryonal', disease_subtype
      ))) %>%
    mutate(PvA=ifelse(
      Type %in% ped_types,
      'Pediatric', 
      'Other'
    ))
  
  return(sample.info.df)
}

# This loads the mapping of genes to the genome
# the default is to use hg38
load_gene_mapping <- function(genome_version='hg38') {
  library(org.Hs.eg.db) # This is using hg38 
  
  ## Filter coordinates to exclude multiple transcripts that are too spread out
  filter.coordinates <- function(loc, coordinate.end, max.spread=2000000) {
    ## Only keep coordinates on propper chromosomes
    loc <- as.numeric(loc[names(loc) %in% c(as.character(1:22), "X", "Y")])
    if (length(loc)==0) {
      return(NA)
    } else if (coordinate.end=="start") {
      ## start coordinate
      if (diff(range(loc)) > max.spread) {
        # print(loc)
      }
      return(ifelse(diff(range(loc))>max.spread, NA, min(abs(loc))))
    } else {
      ## end coordinate
      return(ifelse(diff(range(loc))>max.spread, NA, max(abs(loc))))
    }
  }
  
  ## Only keep propper chromosomes
  filter.chromosomes <- function(chrom) {
    chrom <- intersect(chrom, c(as.character(1:22), "X", "Y"))
    ifelse(length(chrom)>0, chrom, NA)
  }
  
  allENTREZG <- data.frame(array(NA, dim=c(length(mappedkeys(org.Hs.egSYMBOL)),5)))
  colnames(allENTREZG) <- c("EGID", "SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")
  
  allENTREZG$EGID <- mappedkeys(org.Hs.egSYMBOL)
  allENTREZG$SYMBOL <- unlist(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egSYMBOL))
  
  allENTREZG$CHR <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHR), filter.chromosomes))
  
  allENTREZG$CHRLOC <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHRLOC), filter.coordinates, coordinate.end="start" ))
  allENTREZG$CHRLOCEND <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHRLOCEND), filter.coordinates, coordinate.end="end"))
  
  allENTREZG <- subset(allENTREZG, !is.na(CHR) & !is.na(CHRLOC) & !is.na(CHRLOCEND))
  
  if (genome_version=='hg19') {
    allENTREZG <- load.from.taiga(data.name='depmap-wes-cn-data--08f3', data.version=12, data.file='WES_CN_gene.depMap_18q4b') %>%
      dplyr::select(EGID, SYMBOL, CHR, CHRLOC, CHRLOCEND)
  }
  
  # Remove chromosome Y genes, as currently we do not call CN on the Y chromosome
  # allENTREZG %<>% filter(gsub('chr', '', CHR) != 'Y')
  return(allENTREZG)
}

# Map genes to segments
generate_gene_level_matrix_from_segments <- function(gene_mapping, segments) {
  # However we decide to get the gene annotations... we now do the below to map genes to segments to get gene level calls
  # This is very fast (~2 minutes)
  allENTREZG_as_granges <- gene_mapping %>% 
    dplyr::select(start=CHRLOC, end=CHRLOCEND, seqnames=CHR, everything()) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # Now overlap segments and genes to get gene level data
  # unique_samples <- segments_gaps_filled$DepMap_ID %>% sample(size = 100)
  segments_gaps_filled_as_granges <- GenomicRanges::makeGRangesFromDataFrame(segments, keep.extra.columns = T)
  segments_gaps_filled_as_granges_list <- split(segments_gaps_filled_as_granges, segments_gaps_filled_as_granges$DepMap_ID)
  
  # Shamelessly stolen from CERES
  single_cell_line <- length(names(segments_gaps_filled_as_granges_list)) == 1
  gene_level_data <- laply(segments_gaps_filled_as_granges_list, function(seg) {
    hits <- GenomicRanges::findOverlaps(allENTREZG_as_granges, seg) %>% as.data.frame %>%
      dplyr::distinct(queryHits, .keep_all = T)
    CN <- rep(NA, length(allENTREZG_as_granges))
    CN[hits$queryHits] <- as.data.frame(seg)[hits$subjectHits,'Segment_Mean']
    return(CN)
  }, .parallel = FALSE) %>% 
    {if(single_cell_line) as.matrix(.) else t(.)} %>%
    set_colnames(names(segments_gaps_filled_as_granges_list)) %>%
    as.data.frame()
  
  # Add the annotations back
  gene_level_mat_columns_to_add <- c('SYMBOL', 'EGID', 'seqnames', 'start', 'end')
  gene_level_data[,gene_level_mat_columns_to_add] <- allENTREZG_as_granges %>% 
    GenomicRanges::as.data.frame() %>%
    dplyr::select(gene_level_mat_columns_to_add)
  gene_level_data %<>% dplyr::rename(Chromosome=seqnames, CHRLOC=start, CHRLOCEND=end) 
  
  return(gene_level_data)
}

## Publication theme ggplot
theme_Publication <- function(base_size=14, base_family="Arial") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.justification = 'top',
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            # legend.key.size= unit(0.2, "cm"),
            # legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# Color palette (keep adding as necessary)
cols_to_use_for_groups <- c(
  "#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33",
  # My added colors
  "orange", "cyan", "magenta","green","blue","yellow"
)

scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication",manual_pal(values = cols_to_use_for_groups), ...)
}

scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication",manual_pal(values = cols_to_use_for_groups), ...)
}

scale_color_Publication <- function(...){
  discrete_scale("colour","Publication",manual_pal(values = cols_to_use_for_groups), ...)
}





