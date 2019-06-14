## Guillaume Kugener
## for BroadInsitute
## in 2017

# Jérémie Kalfon edited 2019
####
#
# HELPER FUNC  ######################################
#
#########


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

## Only keep proper chromosomes
filter.chromosomes <- function(chrom) {
  chrom <- intersect(chrom, c(as.character(1:22), "X", "Y"))
  ifelse(length(chrom)>0, chrom, NA)
}


library(org.Hs.eg.db) # This is using hg38 
library(plyr)
library(dplyr)
library(readr)
library(stringr)

fusionFusions <- function(input_file_names,output_file_name){
  ## reads selected fields from maf files and aggregates them 
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
}

AddSamplesTo <- function(dataset, listOfSamples){
  # listOfSamples:
  #   $genes_count data 
  #   $trancripts data
  # dataset:
  #   $genes_count list of file path containing sample depmap id
  #   $trancripts list of file path containing sample depmap id
  #   $genes_tpm list of file path containing sample depmap id


  #### TPM (genes)
  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  if (length(listOfSamples$genes_count) > 0) {
    for (f in listOfSamples$genes_count) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(t)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      dataset$genes_tpm %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }

  ##### counts
  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  if (length(listOfSamples$genes_count) > 0) {
    for (f in listOfSamples$genes_count) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(dataset$genes_count)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, expected_count) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      dataset$genes_count %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }

  #### transcripts
  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  if (length(listOfSamples$transcripts) > 0) {
    for (f in listOfSamples$transcripts) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, transcript_id, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id', s_id))
      
      dataset_transcripts$transcripts %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id'))
    }
  }
  return(dataset)
}

########
#
# Copy Number #########################
#
###########


process_segments <- function(segment_file) {


  new_copy_number <- readr::read_tsv(segment_file, col_types = cols(
    Sample = col_character(),
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  ))
  
  new_copy_number %<>% dplyr::select(Sample, Chromosome=CONTIG, Start=START, End=END, 
  Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    
    dplyr::mutate(Source=case_when(
      grepl('^(dm|ccle2|ibm)_', Sample) ~ 'Broad WES',
      grepl('^chordoma_', Sample) ~ 'Chordoma WES',
      TRUE ~ "Other WES"))

  # Untransform data if it is log2 transformed (which it will be if from the FireCloud pipeline)
  if (min(new_copy_number$Segment_Mean) < 0) {
    new_copy_number %<>% mutate(Segment_Mean=2^Segment_Mean)
  }
  return(new_copy_number)
}

filter_for_CCLE <- function(new_copy_number){

  # We shouldn't have anything labelled other, so check this
  if (nrow(new_copy_number %>% filter(Source=='Other WES')) > 0) {
    print('ERROR. THERE IS A SAMPLE NOT FROM BROAD OR THE CHORDOMA FOUNDATION')
  }

  # Make sure there are no duplicates
  counts_unique_by_source <- new_copy_number %>%
    mutate(DepMap_ID=stringr::str_extract(string=Sample, pattern='ACH\\-[0-9]+')) %>%
    distinct(Source, DepMap_ID, Sample) %>%
    group_by(Source, DepMap_ID) %>%
    dplyr::summarise(count=n()) %>%
    filter(count > 1)

  if (nrow(counts_unique_by_source) != 0) {
    print('ERROR. THERE IS A DUPLICATE SAMPLE IN THE SET FOR A SINGLE SOURCE')
  }
  # If it passes above, then can remove prefix from DepMap IDs
  new_copy_number %<>% dplyr::mutate(DepMap_ID=stringr::str_extract(pattern='ACH\\-[0-9]+', string = Sample))
  return(new_copy_number)
}


######################
# The function below fills in the gaps in the segmented level data
#
# In this release, the following changes were made:
#
# * Number of new cell lines: `r length(new_cell_lines)`
# * Number of cell lines moving from Sanger WES/Broad SNP -> Broad WES CN: `r length(replaced_cell_lines)`
# * Total cell lines with CN this release: `r length(unique(combined_new_prioritized_dataset$DepMap_ID))`
#
# ## Interpolate segments
#
# In this section we perform to operations on the segment level data
#
# 1. Fill in the gaps: there may be gaps between segments in the copy number data, 
#  leading to the possibility genes mapping to these gaps and being NAed further downstream.
# 2. Extend the ends: there are genes that map outside the targeting regions (in WES).
#  To address these cases, we can extend the ends of the segments so that these genes are not NAed.
# 
# @args:
#   - segments: 
#       data.frame with DepMap_ID (what samples are separated on), 
#       Chromosome, Start, End, Segment_Mean, Num_Probes
#
# @returns: a segments data.frame of the same size with gaps in ranges filled
######################

interpolate_gaps_in_segmented <- function(segments) {
  segments_as_granges <- GenomicRanges::makeGRangesFromDataFrame(segments, keep.extra.columns = T)
  segments_as_granges_list <- split(segments_as_granges, segments_as_granges$DepMap_ID)
  
  # Determine which ones have too many gaps
  segments_gaps_filled <- NULL
  count <- 0
  disjoin_different_cls <- c()
  for (cl in names(segments_as_granges_list)) {
    if (count %% 100 == 0) {
      print(count/length(names(segments_as_granges_list))* 100)
    }
    count <- count + 1
    
    if (length(segments_as_granges_list[[cl]] %>% disjoin()) != length(segments_as_granges_list[[cl]])) {
      disjoin_different_cls %<>% c(., cl)
      source_using_current <- unique(segments_as_granges_list[[cl]]$Source)
      
      if (source_using_current != 'Broad SNP') {
        print(paste0('Disjoin different for: ', cl, ' source: ', source_using_current))
      }
      
      ttt <- segments_as_granges_list[[cl]] %>%
        as.data.frame() %>%
        dplyr::select(-Source, -width, -DepMap_ID)
      
      ddd <- segments_as_granges_list[[cl]] %>% 
        disjoin() %>%
        as.data.frame()
      
      get_missing_in_disjoined <- ddd %>%
        left_join(., ttt, by=c('seqnames', 'start', 'end', 'strand')) %>%
        dplyr::rename(SM_start_end=Segment_Mean, Num_Probes_start_end=Num_Probes) %>%
        left_join(., ttt %>% dplyr::select(-end), by=c('seqnames', 'start', 'strand')) %>%
        dplyr::rename(SM_start=Segment_Mean, Num_Probes_start=Num_Probes) %>%
        left_join(., ttt %>% dplyr::select(-start), by=c('seqnames','end', 'strand')) %>%
        dplyr::rename(SM_end=Segment_Mean, Num_Probes_end=Num_Probes) %>%
        filter(!(is.na(SM_start_end) & !is.na(SM_start) & !is.na(SM_end))) %>%
        mutate(SM_final=case_when(!is.na(SM_start_end) ~ SM_start_end,
          !is.na(SM_start) ~ SM_start,
          !is.na(SM_end) ~ SM_end,
          TRUE ~ 1)) %>%
        mutate(Num_Probes_final=case_when(!is.na(Num_Probes_start_end) ~ Num_Probes_start_end,
          !is.na(Num_Probes_start) ~ Num_Probes_start,
          !is.na(Num_Probes_end) ~ Num_Probes_end,
          TRUE ~ as.integer(1))) %>%
        mutate(DepMap_ID=cl, Source=source_using_current) %>%
        dplyr::select(seqnames, start, end, Num_Probes=Num_Probes_final, 
          Segment_Mean=SM_final, DepMap_ID, Source) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
      
      segments_as_granges_list[[cl]] <- get_missing_in_disjoined
    }
    # We want to extend the ranges to fill in any gaps. We use GenomicRanges gaps function to identify
    # the gaps in the segmented data that are present per cell line. If there is a gap between
    # two segments, we extend the ends to meet halfway between the gap. This is what the code snippet below
    # is performing
    segments_gaps_filled %<>% rbind(segments_as_granges_list[[cl]] %>%
      GenomicRanges::as.data.frame() %>%
      dplyr::select(-strand, -width) %>%
      left_join(., 
        gaps(segments_as_granges_list[[cl]]) %>% 
          GenomicRanges::as.data.frame() %>%
          dplyr::select(-strand, -width) %>%
          dplyr::mutate(right_start=start, right_end=end) %>%
          dplyr::mutate(end=start-1) %>%
          dplyr::select(-start),
        by=c('seqnames', 'end')) %>%
      left_join(., 
          gaps(segments_as_granges_list[[cl]]) %>% 
          GenomicRanges::as.data.frame() %>%
          dplyr::select(-strand, -width) %>%
          dplyr::mutate(left_start=start, left_end=end) %>%
          dplyr::mutate(start=end+1) %>%
          dplyr::select(-end),
        by=c('seqnames', 'start')) %>%
      # Edit the range values if needed
      dplyr::mutate(
        start=case_when(
          !is.na(left_end) ~ (left_start + floor((left_end-left_start)/2)+1),
          TRUE ~ start
        ),
        end=case_when(
          !is.na(right_end) ~ (right_start + floor((right_end-right_start)/2)),
          TRUE ~ end
      )) %>% dplyr::select(DepMap_ID, seqnames, start, end, Num_Probes, Segment_Mean, Source))
  } 
  return(list(segs=segments_gaps_filled, disjoin_different_cls=disjoin_different_cls))
}

# We use this method to extend the ends of the segments
extend_ends_of_segments <- function(segments, hg38_cyto_band_reference='data/hg38_cytoband.gz') {
  # Starts to 1
  segments %<>%
    group_by(DepMap_ID, seqnames) %>%
    dplyr::mutate(m=min(start)) %>%
    dplyr::mutate(start=ifelse(start==m, 1, start)) %>%
    ungroup() %>%
    dplyr::select(-m)
  
  # Ends to max of chromosome value
  chromosome_maxes <- read_tsv(hg38_cyto_band_reference, col_names = F) %>%
    group_by(X1) %>% dplyr::summarise(max=max(X3)) %$% setNames(max, X1)
  
  segments %<>%
    group_by(DepMap_ID, seqnames) %>%
    dplyr::mutate(m=max(end)) %>%
    dplyr::mutate(end=ifelse(end==m, chromosome_maxes[paste0(seqnames)], end)) %>%
    ungroup() %>%
    dplyr::select(-m)
  return(segments)
}


# This chunk is kept in here for reference, but does not need to be run for the release
# This was used to combine the copy number calls from Sanger and the Broad from WES
# For the Broad samples, we want the calls that use the ICE/AGILENT PON for chr1-22,X
# For Sanger samples we want the calls that use the SANGER specific AGILENT PON for chr1-22 and then the AGILENT PONT for X
# This is what this chunk accomplished and then saved to a tsv for upload to taiga

# Format for upload to taiga
combine_broad_sanger <- function(Sanger_filename= 'Downloads/Sanger.called.seg',
  Broad_filename = '~/Downloads/DM19Q2_COMPLETE.called.seg', outfile="'~/Downloads/all_wes_data_19q2_v2.tsv'"){
  agilent_ice_pon_based_calls <- readr::read_tsv(filename, col_types = cols(
    Sample = col_character(), 
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  )) %>% dplyr::select(Sample, Chromosome=CONTIG, Start=START, End=END, 
  Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    dplyr::mutate(Source=ifelse(grepl('sanger', Sample), 'Sanger WES', 
      ifelse(grepl('chordoma', Sample), 'Chordoma WES', 'Broad WES'))) %>%
    # We only keep the X calls from here 
    filter((Source %in% c('Chordoma WES', 'Broad WES')) | (Source == 'Sanger WES' & gsub('chr', '', Chromosome) == 'X'))

  sanger_pon_based_calls <- readr::read_tsv(, col_types = cols(
    Sample = col_character(), 
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  )) %>% dplyr::select(Sample, Chromosome=CONTIG, Start=START, End=END, 
  Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    dplyr::mutate(Source='Sanger WES') %>%
    # We remove the sex chromosomes from here as this was a mixed panel
    filter(gsub('chr', '', Chromosome) %in% seq(1,22))

  # Combine two datasets above, validate that only one cell line per source, 
  # change the DepMap_ID column to have only the DepMap_ID
  combined_copynumber_calls <- rbind(sanger_pon_based_calls, agilent_ice_pon_based_calls)

  counts_unique_by_source <- combined_copynumber_calls %>%
    mutate(DepMap_ID=stringr::str_extract(string=Sample, pattern='ACH\\-[0-9]+')) %>%
    distinct(Source, DepMap_ID, Sample) %>%
    group_by(Source, DepMap_ID) %>%
    dplyr::summarise(count=n()) %>%
    filter(count > 1)

  if (nrow(counts_unique_by_source) != 0) {
    print('ERROR. THERE IS A DUPLICATE SAMPLE IN THE SET FOR A SINGLE SOURCE')
  }

  combined_copynumber_calls %<>% mutate(DepMap_ID=stringr::str_extract(pattern='ACH\\-[0-9]+', string=Sample)) %>%
    dplyr::select(-Sample) %>%
    mutate(Chromosome=gsub('chr', '', Chromosome)) %>%
    arrange(DepMap_ID, Chromosome)

  write.table(combined_copynumber_calls, file = outfile, sep = '\t', quote = F, row.names = F) 
  # What we upload to taiga
}


# ## Reprioritize processed data

# Here, we combine the previous release segments with the additional 
# segments we have processed this quarter. We replace Sanger WES and Broad SNP prioritized cell 
# lines with the Broad WES version if it is now available in the newly downloaded data.

reprioritize_data <- function(new_copy_number, wes.priority.cn.seg.profiles){
  
  new_copy_number %<>% magrittr::set_colnames(
    c('DepMap_ID','Chromosome','Start','End','Num_Probes','Segment_Mean','Source'))
  print(new_copy_number)
  if (min(wes.priority.cn.seg.profiles$Segment_Mean) < 0) {
    wes.priority.cn.seg.profiles %<>% mutate(Segment_Mean=2^Segment_Mean)
  }

  broad_wes_cell_lines_in_new <- new_copy_number %>% filter(Source=='Broad WES') %$% unique(DepMap_ID)
  replaced_cell_lines <- wes.priority.cn.seg.profiles %>%
    distinct(DepMap_ID, Source) %>%
    dplyr::filter(Source != 'Broad WES') %$%
    intersect(DepMap_ID, broad_wes_cell_lines_in_new)

  new_cell_lines <- broad_wes_cell_lines_in_new %>%
    setdiff(., wes.priority.cn.seg.profiles$DepMap_ID)

  # Use from previous release
  sanger_lines_using <- wes.priority.cn.seg.profiles %>% filter(Source=='Sanger WES') %$% DepMap_ID %>% unique()

  # So we only keep the Sanger lines that were previously used
  new_copy_number %<>% dplyr::filter(Source=='Broad WES' | (Source=='Sanger WES' & (DepMap_ID %in% sanger_lines_using)) | (Source == 'Chordoma WES' & !(DepMap_ID %in% broad_wes_cell_lines_in_new)))
  new_copy_number %<>% filter(!(Source=='Sanger WES' & DepMap_ID %in% broad_wes_cell_lines_in_new))

  # Use this snippet if simply adding new Broad WES samples to the previous release dataset
  # Add the Broad WES data to the new dataset
  combined_new_prioritized_dataset <- rbind(
    wes.priority.cn.seg.profiles %>% dplyr::filter(!(DepMap_ID %in% new_copy_number$DepMap_ID)),
    new_copy_number
  )

  # Check that all cell lines only represented once (this should be empty)
  combined_new_prioritized_dataset %>%
    distinct(Source, DepMap_ID) %>%
    group_by(DepMap_ID) %>%
    dplyr::mutate(count=n()) %>%
    # filter(Source=='Broad SNP')
    filter(count > 1)

  # Check the combined matrix contains all same samples
  if (length(setdiff(unique(wes.priority.cn.seg.profiles$DepMap_ID), unique(combined_new_prioritized_dataset$DepMap_ID))) > 0) {
    stop('THERE IS A CELL LINE MISSING IN THE NEW RELEASE SET THAT WAS PREVIOUSLY INCLUDED')
  }

  # Quick check that we added all the new cell lines
  if ((length(unique(combined_new_prioritized_dataset$DepMap_ID)) - length(unique(wes.priority.cn.seg.profiles$DepMap_ID))) != length(new_cell_lines)) {
    stop('FAILED COMBINING NEW AND PREVIOUS RELEASE DATA')
  }
  return(combined_new_prioritized_dataset)
}

######################
# The function below generates the gene level matrix from a gene mapping
#
# @args:
#   - gene_mapping: 
#       data.frame with EGID, SYMBOL, CHR, CHRLOC, CHRLOCEND
#   - segments: 
#       data.frame with DepMap_ID (what samples are separated on), 
#       Chromosome, Start, End, Segment_Mean, Num_Probes
#
# @returns: a data.frame with EGID, SYMBOL, CHR, CHRLOC, CHRLOCEND followed 
#           by one column for each sample containing the Segment_Mean of the 
#           gene
######################
generate_gene_level_matrix_from_segments <- function(gene_mapping, segments) {
  # However we decide to get the gene annotations... we now do the below to map genes 
  # to segments to get gene level calls
  # This is very fast (~2 minutes)
  allENTREZG_as_granges <- gene_mapping %>% 
    dplyr::select(start=CHRLOC, end=CHRLOCEND, seqnames=CHR, everything()) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # Now overlap segments and genes to get gene level data
  # unique_samples <- segments_gaps_filled$DepMap_ID %>% sample(size = 100)
  segments_gaps_filled_as_granges <- GenomicRanges::makeGRangesFromDataFrame(segments, keep.extra.columns = T)
  segments_gaps_filled_as_granges_list <- split(segments_gaps_filled_as_granges, 
    segments_gaps_filled_as_granges$DepMap_ID)
  
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
  gene_level_data_hg38 <- gene_level_data %>%
  dplyr::mutate(rn=paste0(SYMBOL, ' (', EGID, ')')) %>%
  dplyr::select(c('rn', unique(segments$DepMap_ID))) %>%
  column_to_rownames(var='rn') %>% t()
  return(gene_level_data_hg38)
}

generate_entrez_genes <- function(genome_version='hg38'){
  # This is how Achilles maps genes to chromosomes. However, for CCLE, we use a different reference.
  # I am leaving this snippet here just for reference but we do not use it
  # Pull the gene coordindates from ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human
  cds_gene_paths <- list(
    hg19='ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs105/CCDS.current.txt',
    hg38='ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.20180614.txt'
  )
  gene_mapping_ccds <- read_tsv(url(cds_gene_paths[[genome_version]]), col_types = cols(.default = 'c')) 
  # This is how Achilles gets the genes. There are far fewer genes than using org.Hs.eg.db
  allENTREZG <- data.frame(array(NA, dim=c(length(mappedkeys(org.Hs.egSYMBOL)),5)))
  colnames(allENTREZG) <- c("EGID", "SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")

  allENTREZG$EGID <- mappedkeys(org.Hs.egSYMBOL)
  allENTREZG$SYMBOL <- unlist(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egSYMBOL))

  allENTREZG$CHR <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHR), filter.chromosomes))

  allENTREZG$CHRLOC <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHRLOC), 
    filter.coordinates, coordinate.end="start" ))
  allENTREZG$CHRLOCEND <- unlist(lapply(mget(mappedkeys(org.Hs.egSYMBOL), org.Hs.egCHRLOCEND), 
    filter.coordinates, coordinate.end="end"))

  allENTREZG <- subset(allENTREZG, !is.na(CHR) & !is.na(CHRLOC) & !is.na(CHRLOCEND))

  if (genome_version=='hg19') {
    allENTREZG <- load.from.taiga(data.name='depmap-wes-cn-data--08f3', data.version=12,
      data.file='WES_CN_gene.depMap_18q4b') %>%
      dplyr::select(EGID, SYMBOL, CHR, CHRLOC, CHRLOCEND)
  }

  # Remove chromosome Y genes, as currently we do not call CN on the Y chromosome
  allENTREZG %<>% filter(gsub('chr', '', CHR) != 'Y')
  return(allENTREZG)
}

filter_segs <- function(blacklist_url='', segments_gaps_filled){
  # Remove the internally embargoed lines
  black_listed_lines <- read_xlsx(, sheet = 'Blacklist') %>%
    magrittr::set_colnames(c('Blacklist', 'Type')) %$% Blacklist
  # Save the segmented level data (untransformed)
  segments_gaps_filled %<>%
    filter(!(DepMap_ID %in% black_listed_lines)) %>%
    dplyr::select(DepMap_ID, Chromosome=seqnames, Start=start, End=end, Num_Probes, Segment_Mean, Source)
  return(segments_gaps_filled)
}



#######
#
# Expression # ####################
#
########

read_counts <- function(filepath, counts_samples_to_add_one_off=c()){
  # Process counts
  counts_genes <- read_tsv(
    file = filepath
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  if (length(counts_samples_to_add_one_off) > 0) {
    for (f in counts_samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(counts_genes)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, expected_count) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      counts_genes %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }
  return(counts_genes)
}

read_tpm <- function(filepath, samples_to_add_one_off=c()){
  # TPM (genes)
  tpm_genes <- read_tsv(
    file = filepath
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)

  if (length(samples_to_add_one_off) > 0) {
    for (f in samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(counts_genes)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      tpm_genes %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }
  return(tpm_genes)
}

read_transcripts <- function(filepath, transcripts_samples_to_add_one_off=c()){
  # Transcripts (genes)
  transcripts <- read_tsv(
    file = filepath,
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  if (length(transcripts_samples_to_add_one_off) > 0) {
    for (f in transcripts_samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, transcript_id, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id', s_id))
      
      transcripts %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id'))
    }
  }
  return(transcripts)
}

compare_releases <- function(tpm_genes, previous_release_tpm){
  # Compare the previous release to the current release (using correlations)
  tentative_new_release_tpm <- tpm_genes %>% 
    mutate(gene=stringr::str_extract(string=gene_id, pattern='ENSG[0-9]+')) %>%
    filter(!grepl('PAR_Y', gene_id)) %>% # Not sure what these are
    filter(!grepl('ERCC', gene_id)) %>% # These are the ERCC spike ins
    dplyr::select(c('gene', colnames(.)[grepl('ACH\\-[0-9]+', colnames(.))])) %>%
    column_to_rownames(var='gene') %>%
    t() %>%
    # If we have two samples with the same arx span id (we shouldn't if the sample set 
    # is designed correctly) then it will complain below
    set_rownames(stringr::str_extract(string=row.names(.), pattern = 'ACH\\-[0-9]+'))

  # Now do correlations (of log2+1 TPM)
  overlap_genes <- intersect(colnames(previous_release_tpm), colnames(tentative_new_release_tpm))
  overlap_cell_lines <- intersect(row.names(previous_release_tpm), row.names(tentative_new_release_tpm))

  # Check to see if any cell lines from the previous release are not present in this dataset 
  # (this should not be the case unless there is a known processing error, so this list should be empty)
  row.names(previous_release_tpm) %>% setdiff(row.names(tentative_new_release_tpm))

  # Correlations of samples (could also just look at the set of most variable )
  # Intersect the top 2000 most variable in both
  tpm_19q1_most_variable <- apply(previous_release_tpm, 2, sd) %>% .[order(-.)] %>% names() %>% .[1:2000]
  tpm_19q2_most_variable <- apply(log2(tentative_new_release_tpm+1), 2, sd) %>% .[order(-.)] %>% names() %>% .[1:2000]

  # 95% overlap
  most_variable_for_correlations <- intersect(tpm_19q1_most_variable, tpm_19q2_most_variable)
  length(most_variable_for_correlations)/2000

  correlation_rnaseq_data_releases <- cor(
    t(previous_release_tpm[overlap_cell_lines, most_variable_for_correlations]),
    t(log2(tentative_new_release_tpm[overlap_cell_lines, most_variable_for_correlations]+1)),
  )
  return(correlation_rnaseq_data_releases)
}

# TODO: process the exons in future releaes...

rename_function <- function(columns) {
  columns_new <- ifelse(columns %in% c('Name', 'Description', 'gene_id', 'transcript_id',
    "transcript_id(s)"), columns, ifelse(grepl('ACH\\-[0-9]+$', columns), 
    stringr::str_extract(string=columns, pattern='ACH\\-[0-9]+'), ccle.to.arxspan(columns, ignore.problems = T)
  ))
  
  return(columns_new)
}


# Comparison function
gene_level_comparison_matrix_gen <- function(new_mat, old_mat, cell_lines, genes_using) {
  # Get overlapping genes
  overlapping_genes <- intersect(colnames(new_mat), colnames(old_mat)) %>% intersect(genes_using)
  overlapping_cell_lines <- intersect(row.names(new_mat), row.names(old_mat)) %>% intersect(cell_lines)
  
  # Calculate differences (cell line level, gene level)
  differences <- log2(new_mat[overlapping_cell_lines, overlapping_genes]) - log2(
    old_mat[overlapping_cell_lines, overlapping_genes])
  
  # cell overview
  cl_diffs <- apply(differences, 1, FUN = function(x) length(x[!is.na(x) & abs(x) > 0.5])/length(x[!is.na(x)]))
  gene_diffs <- apply(differences, 2, FUN = function(x) length(x[!is.na(x) & abs(x) > 0.5]))
  
  return(list(gene=gene_diffs, cl=cl_diffs))
}


####
#
# FUSIONS ###############################
#
####

read_fusions <- function(filepath){
  tentative_new_release_unfiltered_fusions <- read_tsv(filepath, col_types = cols(.default = 'c'))
  ccle_totals <- tentative_new_release_unfiltered_fusions %>%
    group_by(LeftBreakpoint, RightBreakpoint) %>%
    dplyr::summarise(count=n()) %>%
    arrange(-count)
  # Add number of times observed in CCLE
  tentative_new_release_unfiltered_fusions %<>% left_join(., ccle_totals %>% dplyr::rename(CCLE_count=count),
    by=c('LeftBreakpoint', 'RightBreakpoint'))
  return(tentative_new_release_unfiltered_fusions)
}


filter_fusions <- function(fusions){
  # Filters we use
  tentative_new_release_filtered_fusions <- fusions %>%
    # (1) Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes, 
    filter(!grepl('^HLA\\-', `#FusionName`)) %>%
    filter(!grepl('chrM', LeftBreakpoint), !grepl('chrM', RightBreakpoint)) %>%
    # (2) Remove red herring fusions
    filter(!grepl('GTEx_recurrent', annots, ignore.case = T)) %>%
    filter(!grepl('DGD_PARALOGS', annots)) %>%
    filter(!grepl('HGNC_GENEFAM', annots)) %>%
    filter(!grepl('Greger_Normal', annots)) %>%
    filter(!grepl('Babiceanu_Normal', annots)) %>%
    filter(!grepl('ConjoinG', annots)) %>%
    filter(!grepl('NEIGHBORS', annots)) %>%
    # (3) Remove recurrent in this dataset (>= 25 samples)
    filter(CCLE_count < 25) %>%
    # (4) Removed fusion with (SpliceType=" INCL_NON_REF_SPLICE" and LargeAnchorSupport="No" and minFAF<0.02), or 
    filter(!(SpliceType=="INCL_NON_REF_SPLICE" & LargeAnchorSupport=="NO_LDAS" & FFPM < 0.1)) %>%
    filter(FFPM > 0.05) 
    # STAR-Fusion suggests using 0.1, but after looking at the 
    # translocation data, this looks like it might be too aggressive
  return(tentative_new_release_filtered_fusions)
}


###########
#
# MUTATIONS #########################################################
#
###########
#


# We will replace these functions with cleaned versions later
varstr = function(maf){
  return(paste(gsub('-Tumor', '', gsub('fh_', '', maf$Tumor_Sample_Barcode)), 
    maf$Chromosome, maf$Start_position, maf$End_position, maf$Reference_Allele, maf$Tumor_Seq_Allele2, sep='_'))
}

varlocstr = function(maf){return(paste(maf$Chromosome,maf$Start_position, sep='_'))}

varcellocstr = function(maf){
  return(paste(gsub('-Tumor', '', gsub('fh_', '', maf$Tumor_Sample_Barcode)), maf$Chromosome, 
    maf$Start_position, sep='_'))
}

varonlystr = function(maf){
  return(paste( maf$Chromosome, maf$Start_position, maf$End_position, maf$Reference_Allele,
    maf$Tumor_Seq_Allele2, sep='_'))
}

# This should not make a new maf
fixSampleName = function(Tumor_Sample_Barcode) { 
  Tumor_Sample_Barcode %<>% gsub('(dm|ccle2|fh|SANGER|chordoma)_', '', .) %>%
    gsub('\\-Tumor', '', .) %>%
    CleanCellLineName(.) %>%
    ccle.to.arxspan(check.unique.mapping = F)
  return(Tumor_Sample_Barcode)
}

removeRec = function(maf, thrCCLErat = 0.05, tcgathr=5, 
  rescueListFN='/Volumes/xchip_cle/Resources/variantFilter/snp_indels_rescue_list.txt') {
  
  vlstr = varlocstr(maf)
  vclstr = varcellocstr(maf)
  
  frq = table(vlstr[!duplicated(vclstr)])
#  hsi = which(maf$TCGAhsCnt<=20)
#  thr = max(frq[match(vlstr[hsi],names(frq))]) +1; ## at this frequency, all the tcga_hs mutations have freq >20 
#  hsi = which((maf$TCGAhsCnt<=20)&(maf$TCGAhsCnt>=5))
#  thr = max(frq[match(vlstr[hsi],names(frq))]) +1; ## at this frequency, all the tcga_hs mutations have freq >20 
#  blacklist = names(which(frq>thr))
#  blacklist = setdiff(blacklist, vlstr[which(maf$TCGAhsCnt>0)])
  
  N = length(unique(maf$Tumor_Sample_Barcode))
  thr = thrCCLErat * N; 
  blacklist = names(which(frq>thr))
  blacklist = setdiff(blacklist, vlstr[which(maf$TCGAhsCnt>=tcgathr)])
  if(rescueListFN!=''){
    mafRescue = read.delim(rescueListFN, sep='\t', header=TRUE, stringsAsFactors = FALSE, comment.char = '#', colClasses='character'); 
    vlstrRescue = varlocstr(mafRescue)
    blacklist = setdiff(blacklist,vlstrRescue)
  } 
  
  pass = !(vlstr %in% blacklist)
  #sort(table(maf$Hugo_Symbol[!pass]))
  #subset(maf, !pass & Hugo_Symbol=='AKT1')
  maf$pass = pass; # we will just mark the duplicates and not remove them 
  
  #  res= list(res=subset(maf,pass),pass=pass, blacklist=blacklist, thr=thr)
  res= list(res=maf,pass=pass, blacklist=blacklist, thr=thr)
}

removeLowQ = function(maf){
  maf = subset(maf, i_qual+i_read_depth > 50)  ## filter GATK indels based on quality
}

intersectmafs = function(maf1, maf2, retainTCGAhs=TRUE){
  cls = intersect(unique(maf1$Tumor_Sample_Barcode), unique(maf2$Tumor_Sample_Barcode))
  maf1 = subset(maf1, Tumor_Sample_Barcode %in% cls)
  maf2 = subset(maf2, Tumor_Sample_Barcode %in% cls)
  
  v1 = varstr(maf1)
  v2 = varstr(maf2)
  
  m1 = match(v1, v2)  
  m2 = match(v2, v1)  
  
  maf2 = maf2[,match(colnames(maf1), colnames(maf2))] 
  
  if(retainTCGAhs){
    maf3 = rbind(maf1[maf1$isTCGAhotspot | !is.na(m1), ], maf2[ maf2$isTCGAhotspot & is.na(m2),])
  }else{
    maf3 = maf1[!is.na(m1), ]
  }
  maf3=maf3
}

makeBlacklistFile = function(wessnp,wessnp_blcklist,current_blacklist_FN, new_blcklist_FN){
  blacklist = read.delim(current_blacklist_FN, sep='\t', stringsAsFactors = FALSE)
  a=subset(wessnp, !is.na(match(paste(wessnp$Chromosome, wessnp$Start_position, sep='_'), wessnp_blcklist)))
  a=subset(a, !duplicated(a$Genome_Change))
  a= a[,c('Hugo_Symbol','Chromosome','Start_position','End_position','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele2','Protein_Change','Tumor_Sample_Barcode')]
  colnames(a)=colnames(blacklist)
  write.table(rbind(blacklist,a), 
              sep='\t', quote=FALSE, row.names=FALSE, file = new_blcklist_FN)
}

addAC = function(snp, snpc, maf, mafnam){
  mafcls=unique(maf$Tumor_Sample_Barcode)
  mafc=varstr(maf); 
    maf = subset(maf, mafc%in%snpc)
    mafc=varstr(maf); 
  res = rep('',length(snpc))
  res[match(mafc, snpc)]= paste(maf$t_alt_count, maf$t_ref_count, sep=':')
  res[!(snp$Tumor_Sample_Barcode %in% mafcls)]=NA
  snp[[mafnam]]=res
  snp
}

addAC_cleaned <- function(given, to_add, given_name){
  res <- given %>%
    left_join(.,
      to_add %>% 
        dplyr::select(c(distinct_columns_to_keep, 't_alt_count', 't_ref_count')) %>%
        mutate(!!given_name := paste0(t_alt_count, ':', t_ref_count)) %>%
        dplyr::select(-t_alt_count, -t_ref_count),
      by = distinct_columns_to_keep
    )
  return(res)
}

polish = function(maf){
  mafc = varonlystr(maf)
  umafc = unique(mafc[duplicated(mafc)])
  mtch = match(mafc,umafc)
  maf[is.na(maf)]='NA'
  colnames(maf)[match("Tumor_Seq_Allele2", colnames(maf))]='Tumor_Seq_Allele1'
  maf
}

