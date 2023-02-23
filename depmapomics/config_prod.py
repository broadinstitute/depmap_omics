from depmapomics.config_global import *

########################## GENERIC PARAMS

isCCLE = True
doCleanup = True

## google storage
BAM_GCS_BUCKET = "gs://cclebams"

RNA_GCS_PATH_HG38 = BAM_GCS_BUCKET + "/rnasq_hg38/"
RNA_GCS_PATH = BAM_GCS_BUCKET + "/rna/"
WGS_GCS_PATH = BAM_GCS_BUCKET + "/wgs/"
WGS_GCS_PATH_HG38 = BAM_GCS_BUCKET + "/wgs_hg38/"
WES_GCS_PATH = BAM_GCS_BUCKET + "/wes/"

## TAIGA specific

TAIGA_ETERNAL = "depmap-a0ab"
TAIGA_ETERNAL_UPLOAD = TAIGA_ETERNAL

DEPMAP_TAIGA = "arxspan-cell-line-export-f808"

TAIGA_MUTATION = "mutations-latest-ed72"
TAIGA_CN = "cn-latest-d8d4"
TAIGA_CN_ACHILLES = "cn-achilles-version-06ca"
TAIGA_EXPRESSION = "expression-d035"
TAIGA_FUSION = "fusions-95c9"
TAIGA_LEGACY_CN = "copy-number-5f61"

VIRTUAL_FOLDER = "8d9c4c0691154a1f86b1b6e67c3fb683"

VIRTUAL = {
    "internal": "",
    "dmc": "",
    "public": "",
}

## our working workspace (reference)
RNAWORKSPACE = "broad-firecloud-ccle/DepMap_hg38_RNAseq"


WGSWORKSPACE = "broad-firecloud-ccle/DepMap_WGS_CN"
WESCNWORKSPACE = "broad-firecloud-ccle/DepMap_WES_CN_hg38"
WESMUTWORKSPACE = "broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline"


FPWORKSPACE = "broad-firecloud-ccle/CCLE_SNP_QC"

TAIGA_FP = "ccle-bam-fingerprints-6f30"
TAIGA_FP_FILENAME = "fingerprint_lod_matrix"


# upload mapping, taiga latest to file name dicts
LATEST2FN_NUMMAT_MODEL = {
    TAIGA_CN: VIRTUAL_FILENAMES_NUMMAT_CN_MODEL,
    TAIGA_EXPRESSION: VIRTUAL_FILENAMES_NUMMAT_EXP_MODEL,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_NUMMAT_MUT_MODEL,
}

LATEST2FN_TABLE_MODEL = {
    TAIGA_CN: VIRTUAL_FILENAMES_TABLE_CN_MODEL,
    TAIGA_FUSION: VIRTUAL_FILENAMES_TABLE_FUSION_MODEL,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_TABLE_MUT_MODEL,
}

LATEST2FN_NUMMAT_PR = {
    TAIGA_CN: VIRTUAL_FILENAMES_NUMMAT_CN_PR,
    TAIGA_EXPRESSION: VIRTUAL_FILENAMES_NUMMAT_EXP_PR,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_NUMMAT_MUT_PR,
}

LATEST2FN_TABLE_PR = {
    TAIGA_CN: VIRTUAL_FILENAMES_TABLE_CN_PR,
    TAIGA_FUSION: VIRTUAL_FILENAMES_TABLE_FUSION_PR,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_TABLE_MUT_PR,
}
