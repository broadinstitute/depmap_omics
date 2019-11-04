#'
#' Add mutation annotations
#'
#' @param mut_maf: matrix MAF matrix
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang "!!"
#' 
#' @description Add columns to maf file describing each mutation as either a silent, damaging
#' non-conserving, or other-conserving mutations
#' 
#' @export maf_add_variant_annotations
maf_add_variant_annotations <- function(mut_maf){
  
  mut_maf$Variant_annotation <- NA
  variant_map <- c(Silent = 'silent',
                   Splice_Site = 'damaging',
                   Missense_Mutation = 'other non-conserving',
                   Nonsense_Mutation = 'damaging',
                   De_novo_Start_OutOfFrame = 'damaging',
                   Nonstop_Mutation = 'other non-conserving',
                   Frame_Shift_Del = 'damaging',
                   Frame_Shift_Ins = 'damaging',
                   In_Frame_Del = 'other non-conserving',
                   In_Frame_Ins = 'other non-conserving',
                   Stop_Codon_Del = 'other non-conserving',
                   Stop_Codon_Ins = 'other non-conserving',
                   Start_Codon_SNP = 'damaging',
                   Start_Codon_Del = 'damaging',
                   Start_Codon_Ins = 'damaging',
                   `5'Flank` = 'other conserving',
                   Intron = 'other conserving',
                   IGR = 'other conserving',
                   `3'UTR` = 'other conserving',
                   `5'UTR` = 'other conserving')
  
  mut_maf %<>%
    dplyr::mutate(Variant_annotation = plyr::revalue(Variant_Classification, variant_map)) 
  
  mut_maf$DepMap_ID <- mut_maf$Tumor_Sample_Barcode
  
  return(mut_maf)
  
}

