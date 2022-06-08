import numpy as np
from depmapomics.config_global import *

########################## GENERIC PARAMS

SHEETNAME = "ccle sample tracker sandbox"

isCCLE = True
doCleanup = True

## google storage
BAM_GCS_BUCKET = "gs://cclebams-sandbox"

RNA_GCS_PATH_HG38 = BAM_GCS_BUCKET + "/rnasq_hg38/"
RNA_GCS_PATH = BAM_GCS_BUCKET + "/rna/"
WGS_GCS_PATH = BAM_GCS_BUCKET + "/wgs/"
WGS_GCS_PATH_HG38 = BAM_GCS_BUCKET + "/wgs_hg38/"
WES_GCS_PATH = BAM_GCS_BUCKET + "/wes/"


REFSHEET_URL = "https://docs.google.com/spreadsheets/d/161mmmHF5nc3nFhgpp_7zub5Erod8xDbXnNZ-u-7BWKI"

DEPMAP_PV = "https://docs.google.com/spreadsheets/d/1uqCOos-T9EMQU7y2ZUw4Nm84opU5fIT1y7jet1vnScE"

POTENTIAL_LIST = "https://docs.google.com/spreadsheets/d/1YuKEgZ1pFKRYzydvncQt9Y_BKToPlHP-oDB-0CAv3gE"

SAMPLES_FOUND_NAME = "depmap ALL samples found"

SAMPLES_NOT_FOUND_NAME = "depmap samples not found"

SAMPLES_NOT_FOUND_URL = "https://docs.google.com/spreadsheets/d/1yC3brpov3JELvzNoQe3eh0W196tfXzvpa0jUezMAxIg"

SAMPLES_MISSING_ARXSPAN_NAME = "depmap samples missing arxspan"

SAMPLES_MISSING_ARXSPAN_URL = "https://docs.google.com/spreadsheets/d/1htfgpXrMvXDlqbcZltpq6vOE_Vo2YZ3-3mdaXV-Irzk"

## TAIGA specific

TAIGA_ETERNAL = "depmap-a0ab"
TAIGA_ETERNAL_UPLOAD = "eternal-74b2"

DEPMAP_TAIGA = "arxspan-cell-line-export-f808"

TAIGA_MUTATION = "mutations-latest-f263"
TAIGA_CN = "cn-latest-8bea"
TAIGA_CN_ACHILLES = "cn-achilles-version-43ea"
TAIGA_EXPRESSION = "expression-869e"
TAIGA_FUSION = "fusions-64c4"
TAIGA_LEGACY_CN = "copy-number-5f61"

VIRTUAL_FOLDER = "aee7ec053d434091a670bc64a9d7a3c1"

VIRTUAL = {
    "internal": "",
    "ibm": "",
    "dmc": "",
    "public": "",
}

## our working workspace (reference)
RNAWORKSPACE = "broad-firecloud-ccle/DEV_DepMap_hg38_RNAseq"


WGSWORKSPACE = "broad-firecloud-ccle/DEV_DepMap_WGS_CN"
WESCNWORKSPACE = "broad-firecloud-ccle/DepMap_WES_CN_hg38-sandbox"
WESMUTWORKSPACE = "broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline-sandbox"


FPWORKSPACE = "broad-firecloud-ccle/CCLE_SNP_QC-copy"

TAIGA_FP = "ccle-bam-fingerprints-4f4a"
TAIGA_FP_FILENAME = "fingerprint_lod_matrix"


########### Gumbo configs #############
GUMBO_SHEET = "https://docs.google.com/spreadsheets/d/10Lg0xkT5OHLYgJ9VKpkh8VR64TXfxPVJXRVAckU8uBg"
GUMBO_SHEETNAME = "Backfilled profile IDs"


# upload mapping, taiga latest to file name dicts
LATEST2FN_NUMMAT = {
    TAIGA_CN: VIRTUAL_FILENAMES_NUMMAT_CN,
    TAIGA_EXPRESSION: VIRTUAL_FILENAMES_NUMMAT_EXP,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_NUMMAT_MUT,
}

LATEST2FN_TABLE = {
    TAIGA_CN: VIRTUAL_FILENAMES_TABLE_CN,
    TAIGA_FUSION: VIRTUAL_FILENAMES_TABLE_FUSION,
    TAIGA_MUTATION: VIRTUAL_FILENAMES_TABLE_MUT,
}
