BAM_GCS_BUCKET = "gs://cclebams"

RNA_GCS_PATH_HG38 = "gs://cclebams/rnasq_hg38/"

RNA_GCS_PATH = "gs://cclebams/rna/"

WGS_GCS_PATH = "gs://cclebams/wgs/"

WGS_GCS_PATH_HG38 = "gs://cclebams/wgs_hg38/"

WES_GCS_PATH = "gs://cclebams/wes/"

TAIGA_ETERNAL_UPLOAD = "depmap-a0ab"

TAIGA_MUTATION = "mutations-latest-ed72"

TAIGA_CN = "cn-latest-d8d4"

TAIGA_CN_ACHILLES = "cn-achilles-version-06ca"

TAIGA_EXPRESSION = "expression-d035"

TAIGA_FUSION = "fusions-95c9"

VIRTUAL_FOLDER = "8d9c4c0691154a1f86b1b6e67c3fb683"

RNAWORKSPACE = "broad-firecloud-ccle/DepMap_hg38_RNAseq"

WGSWORKSPACE = "broad-firecloud-ccle/DepMap_WGS_CN"

WESCNWORKSPACE = "broad-firecloud-ccle/DepMap_WES_CN_hg38"

FPWORKSPACE = "broad-firecloud-ccle/CCLE_SNP_QC"

TAIGA_FP = "ccle-bam-fingerprints-6f30"

LATEST2FN_NUMMAT_MODEL = {
    "cn-latest-d8d4": {
        "merged_gene_cn_profile": "OmicsCNGene",
        "merged_absolute_gene_cn_profile": "OmicsAbsoluteCNGene",
        "merged_loh_profile": "OmicsLoH",
        "globalGenomicFeatures_profile": "OmicsSignatures",
    },
    "expression-d035": {
        "proteinCoding_genes_tpm_logp1_profile": "OmicsExpressionProteinCodingGenesTPMLogp1",
        "gene_set_enrichment_profile": "OmicsExpressionGeneSetEnrichment",
    },
    "mutations-latest-ed72": {
        "somaticMutations_genotypedMatrix_hotspot_profile": "OmicsSomaticMutationsMatrixHotspot",
        "somaticMutations_genotypedMatrix_damaging_profile": "OmicsSomaticMutationsMatrixDamaging",
        "somaticMutations_genotypedMatrix_driver_profile": "OmicsSomaticMutationsMatrixDriver",
    },
}

LATEST2FN_TABLE_MODEL = {
    "cn-latest-d8d4": {},
    "fusions-95c9": {"fusions_filtered_profile": "OmicsFusionFiltered"},
    "mutations-latest-ed72": {
        "somaticMutations_profile": "OmicsSomaticMutations",
        "structuralVariants_profile": "OmicsStructuralVariants",
    },
}

LATEST2FN_NUMMAT_PR = {
    "cn-latest-d8d4": {"globalGenomicFeatures_profile": "OmicsSignaturesProfile"},
    "expression-d035": {
        "genes_expectedCount_profile": "OmicsExpressionGenesExpectedCountProfile",
        "transcripts_expectedCount_profile": "OmicsExpressionTranscriptsExpectedCountProfile",
        "transcripts_tpm_logp1_profile": "OmicsExpressionTranscriptsTPMLogp1Profile",
        "gene_set_enrichment_profile": "OmicsExpressionGeneSetEnrichmentProfile",
    },
    "mutations-latest-ed72": {},
}

LATEST2FN_TABLE_PR = {
    "cn-latest-d8d4": {
        "merged_segments_profile": "OmicsCNSegmentsProfile",
        "merged_absolute_segments_profile": "OmicsAbsoluteCNSegmentsProfile",
    },
    "fusions-95c9": {"fusions_unfiltered_profile": "OmicsFusionUnfilteredProfile"},
    "mutations-latest-ed72": {
        "somaticMutations_profile": "OmicsSomaticMutationsProfile",
        "structuralVariants_profile": "OmicsStructuralVariantsProfile",
    },
}
