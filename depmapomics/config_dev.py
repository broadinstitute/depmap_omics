SHEETNAME = 'ccle sample tracker sandbox'

BAM_GCS_BUCKET = 'gs://cclebams-sandbox'

RNA_GCS_PATH_HG38 = 'gs://cclebams-sandbox/rnasq_hg38/'

RNA_GCS_PATH = 'gs://cclebams-sandbox/rna/'

WGS_GCS_PATH = 'gs://cclebams-sandbox/wgs/'

WGS_GCS_PATH_HG38 = 'gs://cclebams-sandbox/wgs_hg38/'

WES_GCS_PATH = 'gs://cclebams-sandbox/wes/'

REFSHEET_URL = 'https://docs.google.com/spreadsheets/d/161mmmHF5nc3nFhgpp_7zub5Erod8xDbXnNZ-u-7BWKI'

SAMPLES_FOUND_NAME = 'depmap ALL samples found'

SAMPLES_NOT_FOUND_NAME = 'depmap samples not found'

SAMPLES_NOT_FOUND_URL = 'https://docs.google.com/spreadsheets/d/1yC3brpov3JELvzNoQe3eh0W196tfXzvpa0jUezMAxIg'

SAMPLES_MISSING_ARXSPAN_NAME = 'depmap samples missing arxspan'

SAMPLES_MISSING_ARXSPAN_URL = 'https://docs.google.com/spreadsheets/d/1htfgpXrMvXDlqbcZltpq6vOE_Vo2YZ3-3mdaXV-Irzk'

TAIGA_ETERNAL_UPLOAD = 'eternal-74b2'

TAIGA_MUTATION = 'mutations-latest-f263'

TAIGA_CN = 'cn-latest-8bea'

TAIGA_CN_ACHILLES = 'cn-achilles-version-43ea'

TAIGA_EXPRESSION = 'expression-869e'

TAIGA_FUSION = 'fusions-64c4'

VIRTUAL_FOLDER = 'aee7ec053d434091a670bc64a9d7a3c1'

RNAWORKSPACE = 'broad-firecloud-ccle/DEV_DepMap_hg38_RNAseq'

WGSWORKSPACE = 'broad-firecloud-ccle/DEV_DepMap_WGS_CN'

WESCNWORKSPACE = 'broad-firecloud-ccle/DepMap_WES_CN_hg38-sandbox'

WESMUTWORKSPACE = 'broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline-sandbox'

FPWORKSPACE = 'broad-firecloud-ccle/CCLE_SNP_QC-copy'

TAIGA_FP = 'ccle-bam-fingerprints-4f4a'

GUMBO_SHEET = 'https://docs.google.com/spreadsheets/d/10Lg0xkT5OHLYgJ9VKpkh8VR64TXfxPVJXRVAckU8uBg'

GUMBO_SHEETNAME = 'Backfilled profile IDs'

LATEST2FN_NUMMAT_MODEL = {'cn-latest-8bea': {'merged_gene_cn_profile': 'OmicsCNGene', 'merged_absolute_gene_cn_profile': 'OmicsAbsoluteCNGene', 'merged_loh_profile': 'OmicsLoH', 'globalGenomicFeatures_profile': 'OmicsSignatures'}, 'expression-869e': {'proteinCoding_genes_tpm_logp1_profile': 'OmicsExpressionProteinCodingGenesTPMLogp1', 'gene_set_enrichment_profile': 'OmicsExpressionGeneSetEnrichment'}, 'mutations-latest-f263': {'somaticMutations_genotypedMatrix_hotspot_profile': 'OmicsSomaticMutationsMatrixHotspot', 'somaticMutations_genotypedMatrix_damaging_profile': 'OmicsSomaticMutationsMatrixDamaging', 'somaticMutations_genotypedMatrix_driver_profile': 'OmicsSomaticMutationsMatrixDriver'}}

LATEST2FN_TABLE_MODEL = {'cn-latest-8bea': {}, 'fusions-64c4': {'fusions_filtered_profile': 'OmicsFusionFiltered'}, 'mutations-latest-f263': {'somaticMutations_profile': 'OmicsSomaticMutations', 'structuralVariants_profile': 'OmicsStructuralVariants'}}

LATEST2FN_NUMMAT_PR = {'cn-latest-8bea': {'globalGenomicFeatures_profile': 'OmicsSignaturesProfile'}, 'expression-869e': {'genes_expectedCount_profile': 'OmicsExpressionGenesExpectedCountProfile', 'transcripts_expectedCount_profile': 'OmicsExpressionTranscriptsExpectedCountProfile', 'gene_set_enrichment_profile': 'OmicsExpressionGeneSetEnrichmentProfile'}, 'mutations-latest-f263': {}}

LATEST2FN_TABLE_PR = {'cn-latest-8bea': {'merged_segments_profile': 'OmicsCNSegmentsProfile', 'merged_absolute_segments_profile': 'OmicsAbsoluteCNSegmentsProfile'}, 'fusions-64c4': {'fusions_unfiltered_profile': 'OmicsFusionUnfilteredProfile'}, 'mutations-latest-f263': {'somaticMutations_profile': 'OmicsSomaticMutationsProfile', 'structuralVariants_profile': 'OmicsStructuralVariantsProfile'}}

