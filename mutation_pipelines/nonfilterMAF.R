# mghandi, chrislo, 7/19/16
#Applies several filters on the MAF file, mainly to remove germline and artifacts

filterMaf = function(inMAFfn,
                     outMAFfn,
                     outMAFannotatedFN="",
                     minAF=0, 
                     maxExAC_AF = 0, 
                     maxExAC_AF_Cosmic = 1, 
                     maxPon_loglike = 0,
                     onlyCoding = FALSE,
                     minCoverage = 3,
                     minAltReads = 2, 
                     TCGAhotspotMutFN = "/xchip/cle/Resources/variantFilter/all_hotspots_for_Mahmoud_083116.tsv",
                     TCGAhotspotMinCnt = 3,
                     COSMIChotspotMinCnt = 5,
                     #blackListFN = "/xchip/cle/Resources/variantFilter/pancan_mutation_blacklist.v14.hg19.txt", 
                     intersectWithMafFN = ""
                 
                     ){
  
#  inMAFfn: input MAF 
#  outMAFfn: output MAF (filtered)
#  outMAFannotatedFN: output maf file, only annotated but not filtered 
  
#  minAF: min allelic Fraction 
#  maxExAC_AF: max ExAC allelic fraction   
#  maxPon_loglike : max panel of normal (PoN) filter logLikelihood score
#  onlyCoding: filte out variants in non coding regions (introns, IGR, ..) 
#  minCoverage: min number of total reads covering the variant
#  minAltReads: min number of reads supporting alt allele 
#  TCGAhotspotMutFN: list of hotspot mutations to rescue (also see TCGAhotspotMinCnt)  
#  TCGAhotspotMinCnt: minimum number of counts needed to be rescued
#  blackListFN: list of mutations to be removed 
#  intersectWithMafFN: a second maf file name to be used to intersect with the first MAF. This was added to be used when merging outputs of two different filters. 
  
  
  convertToNumeric = function(ExAC_AF){  
  ## converts string to numeric, if there are multiple values separated by |, it returns max value 
    pickMax = function(x){
      max(as.numeric(unlist(strsplit(x, '|', fixed = TRUE))))
    }
    #ExAC_AF = maf$i_ExAC_AF
    if(typeof(ExAC_AF)=='character' ){
      ii = grep('|', ExAC_AF, fixed = TRUE)
      if(length(ii)>0){
        ExAC_AF[ii] = sapply(ExAC_AF[ii], pickMax)
      }
      ExAC_AF = as.numeric(ExAC_AF)
    }
    if(typeof(ExAC_AF)!='double' ){
      ExAC_AF = as.numeric(ExAC_AF)
    }
    return(ExAC_AF)
  }
  
  
  if (!(toupper(intersectWithMafFN) %in% c("","EMPTY","NA","NULL"))){
    maf = read.delim(intersectWithMafFN, sep='\t', header=TRUE, stringsAsFactors = FALSE, comment.char = '#', quote=''); 
    OtherMafcp = paste(maf$Chromosome,maf$Start_position,maf$End_position,maf$Reference_Allele, maf$Tumor_Seq_Allele2, sep='_')
  }
  #inMAFfn = '~/Downloads/NCIH1975_LUNG-Tumor.snp_wes_gfannot.maf';
  #inMAFfn = '~/CCLE/tmp.txt'
  #inMAFfn = '~/CCLE/tmp5.txt'
  #intersectWithMafFN= '~/CCLE/tmp2.txt'
  maf = read.delim(inMAFfn, sep='\t', header=TRUE, stringsAsFactors = FALSE, comment.char = '#', quote=''); 

  filtOther = rep(NA, nrow(maf))
  if (!(toupper(intersectWithMafFN) %in% c("","EMPTY","NA","NULL"))){
    mafcp = paste(maf$Chromosome,maf$Start_position,maf$End_position,maf$Reference_Allele, maf$Tumor_Seq_Allele2, sep='_')
    filtOther = is.na(match(mafcp,OtherMafcp)); 
  }
  
  # filter: TRUE means failed (remove variant)
  
  # Cosmic cnts (overlapping muts) 
  
  COSMICcnt = unlist(sapply(as.character(maf$COSMIC_overlapping_mutations),
    function(x){sum(as.numeric(gsub(')','', sapply(unlist(strsplit(x, '|', fixed=TRUE)),
      function(x){ unlist(strsplit(x,'(', fixed=TRUE))[2] }))))})); 
  if(is.null(maf$COSMIC_overlapping_mutations)){
    COSMICcnt=  maf$COSMIC_n_overlapping_mutations
    if(is.null(COSMICcnt)){
      COSMICcnt=  maf$i_COSMIC_n_overlapping_mutations
      
    }
  }
  COSMICcnt= convertToNumeric(COSMICcnt)
  isCOSMIChotspot = (COSMICcnt>=COSMIChotspotMinCnt); 
 
  # hotspot mutations, note: only TCGA for naw, uses exact position matching (may not cover some indels properly) 
  if (toupper(TCGAhotspotMutFN) %in% c("","EMPTY","NA","NULL")){
    isTCGAhotspot = rep(NA, nrow(maf)); 
    TCGAhsCnt = isTCGAhotspot; 
  }else{
    hs = read.delim(TCGAhotspotMutFN, sep='\t', stringsAsFactors = FALSE)
    hscp = paste(hs$chr, hs$pos, sep='_')
    mafcp = paste(maf$Chromosome, maf$Start_position, sep='_')
    TCGAhsCnt = hs$count[match(mafcp,hscp)]; #number of times observed in TCGA (>2)
    isTCGAhotspot = !is.na(TCGAhsCnt); 
    isTCGAhotspot[which(TCGAhsCnt<TCGAhotspotMinCnt)]=FALSE
  }
    
  if(length(maf$i_ExAC_AF)>0){
    maf$i_ExAC_AF=convertToNumeric(maf$i_ExAC_AF)
    filtExAC = (maf$i_ExAC_AF >maxExAC_AF)
    filtExAC[which(((COSMICcnt>0)|(TCGAhsCnt>0)) & (maf$i_ExAC_AF <= maxExAC_AF_Cosmic))]=FALSE; 
  }else{
    maf$ExAC_AF=convertToNumeric(maf$ExAC_AF)
    filtExAC = (maf$ExAC_AF >maxExAC_AF)
    filtExAC[which(((COSMICcnt>0)|(TCGAhsCnt>0)) & (maf$ExAC_AF <= maxExAC_AF_Cosmic))]=FALSE; 
  }

  af=as.numeric(maf$t_alt_count)/(as.numeric(maf$t_ref_count)+as.numeric(maf$t_alt_count))
  filtAF = (af<minAF) | is.na(af); 

  if (maxPon_loglike<0){
    filtPon_loglike = (maf$pon_loglike>maxPon_loglike)
  }else{
    filtPon_loglike = rep(NA, nrow(maf)) 
  }
  
  filtNonCoding = !(maf$Variant_Classification %in% c('De_novo_Start_OutOfFrame','Nonstop_Mutation',
    'Nonsense_Mutation','Missense_Mutation','Silent','Splice_Site','Start_Codon_SNP',
    'Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Start_Codon_Del',
    'Start_Codon_Ins','Stop_Codon_Del','Stop_Codon_Ins'))

  isDeleterious = (maf$Variant_Classification %in% c('De_novo_Start_OutOfFrame','Nonsense_Mutation',
    'Nonstop_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins','Start_Codon_Del',
    'Start_Codon_Ins','Stop_Codon_Del','Stop_Codon_Ins'))
  
  filtCoverage = ((as.numeric(maf$t_ref_count)+as.numeric(maf$t_alt_count)) < minCoverage)| is.na(af); 

  # blacklist 
  filtBlacklist = rep(NA, nrow(maf)); 
  
  
  pass = !(filtExAC | filtAF | filtPon_loglike| filtCoverage| filtOther);  ## assumes that R understands that NA | TRUE = TRUE |NA = TRUE !
  pass[is.na(pass)]=TRUE
  pass[which(isTCGAhotspot)] = TRUE
  pass[which(filtBlacklist)] = FALSE
  pass[which(onlyCoding&filtNonCoding)] = FALSE
  
  maf$isDeleterious = isDeleterious
  maf$isTCGAhotspot = isTCGAhotspot
  maf$TCGAhsCnt = TCGAhsCnt
  maf$isCOSMIChotspot = isCOSMIChotspot
  maf$COSMIChsCnt = COSMICcnt
  
  maf$PASS = pass
  
  write.table(subset(maf,PASS), file= outMAFfn, sep='\t', row.names = FALSE, quote = FALSE)

  if (!(toupper(outMAFannotatedFN) %in% c("","EMPTY","NA","NULL"))){
    maf$filtExAC = filtExAC
    maf$filtAF = filtAF
    maf$filtPon_loglike = filtPon_loglike
    maf$filtNonCoding = filtNonCoding
    maf$filtCoverage = filtCoverage
    maf$filtBlacklist = filtBlacklist
    maf$filtOther = filtOther; 
    write.table(maf, file= outMAFannotatedFN, sep='\t', row.names = FALSE, quote = FALSE)
  }    
}   
  
  
args = c('/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/FilterMAF/broadinstitute.org/cancer.genome.analysis/11274/18/filterMAF.R',
         '/fh/subscription-ccle/maf_pon_filter/fh_786O_KIDNEY/25580760/fh_786O_KIDNEY.pon_annotated.maf',
         'fh_786O_KIDNEY.indel_wgs_gfpass.maf',
         'fh_786O_KIDNEY.indel_wgs_gfannot.maf',
         '0.1',
         '0',
         '1E-4',
         '-2.5',
         'TRUE',
         '4',
         '2',
         '/xchip/cle/Resources/variantFilter/all_hotspots_for_Mahmoud_083116_ExAC_filt.tsv',
         '3',
         '10',
         '/xchip/cle/Resources/variantFilter/pancan_mutation_blacklist.v14.hg19.txt',
         '/fh/subscription-ccle/maf_pon_filter/fh_786O_KIDNEY/25580586/fh_786O_KIDNEY.pon_annotated.pass.maf')[-1]

args<-commandArgs(TRUE)
inMAFfn=as.character(args[[1]]) 
outMAFfn =as.character(args[[2]]) 
outMAFannotatedFN = as.character(args[[3]]) 
minAF = as.numeric(args[[4]]) 
maxExAC_AF = as.numeric(args[[5]]) 
maxExAC_AF_Cosmic = as.numeric(args[[6]]) 
maxPon_loglike = as.numeric(args[[7]]) 
onlyCoding = toupper(as.character(args[[8]])) %in% c('T','TRUE')
minCoverage = as.numeric(args[[9]]) 
minAltReads =as.numeric(args[[10]]) 
TCGAhotspotMutFN = as.character(args[[11]]) 
TCGAhotspotMinCnt= as.numeric(args[[12]]) 
COSMIChotspotMinCnt =as.numeric(args[[13]]) 
blackListFN =as.character(args[[14]]) 
intersectWithMafFN = as.character(args[[15]]) 

filterMaf(inMAFfn=inMAFfn ,
          outMAFfn = outMAFfn,
          outMAFannotatedFN = outMAFannotatedFN ,
          minAF = minAF,
          maxExAC_AF = maxExAC_AF ,
          maxExAC_AF_Cosmic = maxExAC_AF_Cosmic,
          maxPon_loglike = maxPon_loglike,
          onlyCoding = onlyCoding, 
          minCoverage = minCoverage,
          minAltReads =minAltReads,
          TCGAhotspotMutFN = TCGAhotspotMutFN,
          TCGAhotspotMinCnt=TCGAhotspotMinCnt,
          COSMIChotspotMinCnt = COSMIChotspotMinCnt,
          blackListFN = blackListFN,
          intersectWithMafFN = intersectWithMafFN)

