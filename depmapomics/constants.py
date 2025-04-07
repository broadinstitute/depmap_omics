SAMPLESETNAME = "25Q2"

RELEASE = "25q2"

### !!! for the sake of testing only, need to replace once finalized!!!
OMICS_ID_MAPPING_TABLE_TAIGAID = 'test-files-for-simone-postproc-ac2c'

OMICS_ID_MAPPING_TABLE_NAME = "merged_table.20250327_1239selcols"

WORKING_DIR = "output/"

GCS_PAYER_PROJECT = "broad-firecloud-ccle"

SAMPLEID = "DepMap_ID"

SOURCE_PRIORITY = ["BROAD", "DEPMAP", "IBM", "CCLE2", "SANGER", "CHORDOMA", "PRISM"]

ENSEMBL_SERVER_V = "http://may2021.archive.ensembl.org/biomart"

HG38BAMCOL = ["bam_filepath", "bai_filepath"]

LEGACY_BAM_COLNAMES = ["hg19_bam_filepath", "hg19_bai_filepath"]

HG38_CRAM_COLNAMES = ["hg38_cram_filepath", "hg38_crai_filepath"]

WGSSETENTITY = "sample_set"

WESSETENTITY = "sample_set"

FPALLBATCHPAIRSETS = "all"

PREV_VIRTUAL = {
    "public": "public-22q1-305b",
    "ibm": "ibm-22q1-cce1",
    "dmc": "dmc-22q1-d00a",
    "internal": "internal-22q1-1778",
}

DATASETS = ["internal", "dmc", "public"]

PROCQC = [
    "allelic_counts_tumor",
    "delta_MAD_tumor",
    "denoised_MAD_tumor",
    "scaled_delta_MAD_tumor",
    "denoised_copy_ratios_lim_4_plot_tumor",
    "denoised_copy_ratios_plot_tumor",
    "modeled_segments_plot_tumor",
    "gatk_cnv_all_plots",
    "lego_plotter_pngs",
    "copy_number_qc_report",
    "ffpe_OBF_figures",
    "mut_legos_html",
    "oxoG_OBF_figures",
    "tumor_bam_base_distribution_by_cycle_metrics",
    "tumor_bam_converted_oxog_metrics",
]

BAMQC = [
    "duplication_metrics",
    "bqsr_report",
    "tumor_bam_alignment_summary_metrics",
    "tumor_bam_bait_bias_summary_metrics",
    "tumor_bam_gc_bias_summary_metrics",
    "tumor_bam_hybrid_selection_metrics",
    "tumor_bam_insert_size_histogram",
    "tumor_bam_insert_size_metrics",
    "tumor_bam_pre_adapter_summary_metrics",
    "tumor_bam_quality_by_cycle_metrics",
    "tumor_bam_quality_distribution_metrics",
    "tumor_bam_quality_yield_metrics",
]

STR_BED = "data/str_hg38.bed"
STR_COLNAME = "str"
STR_METHOD = "hipstr"

STR_LOCI_14 = [
    "CSF1PO",
    "D13S317",
    "D16S539",
    "D18S51",
    "D21S11",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "FGA",
    "PentaD",
    "PentaE",
    "TH01",
    "TPOX",
]
STR_LOCI_13 = [
    "CSF1PO",
    "D13S317",
    "D16S539",
    "D18S51",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "FGA",
    "PentaD",
    "PentaE",
    "TH01",
    "TPOX",
]

GUIDESBED = {
    "avana": "data/avana_guides.bed",
    "humagne": "data/humagne_guides.bed",
    "ky": "data/ky_score_guides.bed",
    "tkov3": "data/tkov_25Q2.bed",
    "brunello": "data/brunello_25Q2.bed"
}

BINARY_MUT_COLNAME_DICT = {
    "avana": "avana_binary_mut",
    "humagne": "humagne_binary_mut",
    "ky": "ky_binary_mut",
    "brunello": "brunello_binary_mut",
    "tkov3": "tkov3_binary_mut"
}

GENECHANGETHR = 0.025

SEGMENTSTHR = 1500

MAXYCHROM = 150

COLRENAMING = {
    "CONTIG": "Chromosome",
    "START": "Start",
    "END": "End",
    "end": "End",
    "seqnames": "Chromosome",
    "start": "Start",
    "Sample": "DepMap_ID",
    "NUM_POINTS_COPY_RATIO": "NumProbes",
    "MEAN_LOG2_COPY_RATIO": "SegmentMean",
    "CALL": "Status",
}

PURECN_MIN_GOF = 70

PURECN_MAX_PLOIDY = 5

PURECN_COLRENAMING = {
    "start": "Start",
    "end": "End",
    "chr": "Chromosome",
    "Sampleid": "DepMap_ID",
    "type": "LoHStatus",
    "C": "SegmentAbsoluteCN",
    "M": "MinorAlleleAbsoluteCN",
}

PURECN_LOH_COLNAME = "PureCN_loh_merged"

PURECN_FAILED_COLNAME = "PureCN_failed"

PURECN_SAMPLESET = "PureCN"

PURECN_LOHVALUES = [
    "LOH",
    "COPY-NEUTRAL LOH",
    "WHOLE ARM COPY-NEUTRAL LOH",
    "WHOLE ARM LOH",
]

SIGTABLE_TERRACOLS = {
    "PureCN_ploidy",
    "PureCN_wgd",
    "PureCN_loh_fraction",
    "PureCN_cin_allele_specific_ploidy_robust",
}

MISC_SIG_TERRACOLS = {"msisensor2_score"}

SIGTABLE_BINARYCOLS = ["PureCN_wgd"]

SIGTABLE_RENAMING = {
    "PureCN_loh_fraction": "LoHFraction",
    "PureCN_wgd": "WGD",
    "PureCN_cin_allele_specific_ploidy_robust": "CIN",
    "PureCN_ploidy": "Ploidy",
    "msisensor2_score": "MSIScore",
}

HG38_CENTROMERE = "data/Modeled_regions_for_GRCh38.tsv"

CNA_ARMS = [
    "1p",
    "1q",
    "2p",
    "2q",
    "3p",
    "3q",
    "4p",
    "4q",
    "5p",
    "5q",
    "6p",
    "6q",
    "7p",
    "7q",
    "8p",
    "8q",
    "9p",
    "9q",
    "10p",
    "10q",
    "11p",
    "11q",
    "12p",
    "12q",
    "13q",
    "14q",
    "15q",
    "16p",
    "16q",
    "17p",
    "17q",
    "18p",
    "18q",
    "19p",
    "19q",
    "20p",
    "20q",
    "21q",
    "22q",
]


GENEMASKTHRESH = 0.5

SEGDUP_BED = "data/segDup_majorAllele_withAltContigs_98pcFracMatch_merged.bed"

RM_BED = "data/repeatMasker_max10_noAlt_merged.bed"

BEDTOOLSPATH = "~/bin/bedtools2/bin/"

MINFREQTOCALL = 0.25

SV_COLNAME = "expanded_filtered_sv_bedpe"

SV_FILENAME = "all_sv.csv"

SV_MAT_FILENAME = "all_sv_matrix.csv"

SV_INTERNAL_AF_CUTOFF = 0.5

SV_HEADER = [
    "CHROM_A",
    "START_A",
    "END_A",
    "ID",
    "STRAND_A",
    "TYPE",
    "FILTER",
    "REF_A",
    "ALT_A",
    "SVLEN_A",
    "MATEID_A",
    "SVINSLEN_A",
    "BND_DEPTH_A",
    "MATE_BND_DEPTH_A",
    "SYMBOL_A",
    "GENEID_A",
    "vep_SV_overlap_name_A",
    "vep_SV_overlap_AF_A",
    "CHROM_B",
    "START_B",
    "END_B",
    "STRAND_B",
    "REF_B",
    "ALT_B",
    "SVLEN_B",
    "MATEID_B",
    "SVINSLEN_B",
    "BND_DEPTH_B",
    "MATE_BND_DEPTH_B",
    "SYMBOL_B",
    "GENEID_B",
    "vep_SV_overlap_name_B",
    "vep_SV_overlap_AF_B",
    "DEL_SYMBOLS",
    "DUP_SYMBOLS",
    "PR",
    "SR",
    "Rescue",
    SAMPLEID,
]

MAF_COL = "depmap_maf_25q2"

HGNC_MAPPING = "data/new_and_old_hgnc_symbols.txt"

DNA_CHANGE_COL = "DNAChange"

CHROM_COL = "Chrom"

POS_COL = "Pos"

IMMORTALIZED_THR = 0.05

HUGO_COL = "HugoSymbol"

LIKELY_LOF_COL = "LikelyLoF"

HOTSPOT_COL = "Hotspot"

CCLE_DELETERIOUS_COL = "CCLEDeleterious"

HESS_COL = "hess_driver"

ONCOKB_HOTSPOT_COL = "oncokb_hotspot"

COSMIC_TIER_COL = "cosmic_tier"

PROTEIN_DICT = {
    "Cys": "C",
    "Asp": "D",
    "Ser": "S",
    "Gln": "Q",
    "Lys": "K",
    "Ile": "I",
    "Pro": "P",
    "Thr": "T",
    "Phe": "F",
    "Asn": "N",
    "Gly": "G",
    "His": "H",
    "Leu": "L",
    "Arg": "R",
    "Trp": "W",
    "Ala": "A",
    "Val": "V",
    "Glu": "E",
    "Tyr": "Y",
    "Met": "M",
}

# see https://github.com/broadinstitute/cds-omics-wip/blob/main/projects/pipeline_improvement_mutation/update_oncokb_gene_lists.ipynb 
# for details on how the oncokb gene lists are generated
ONCOKB_ONCOGENE_LIST = "depmapomics/tasks/oncokb_dm/data/oncogene_oncokb_20250205.txt"

# generated by mapping hugo symbols in depmapomics/tasks/oncokb_dm/data/oncogene_oncokb_20250205.txt to ENSG ids
# using: hgnc_complete_set = tc.get(name='hgnc-gene-table-e250', version=3, file='hgnc_complete_set')
ONCOKB_ONCOGENE_ENSG_LIST = "depmapomics/tasks/oncokb_dm/data/oncogene_and_ts_oncokb_ensg_20250205.txt"

ONCOKB_TM_LIST = "depmapomics/tasks/oncokb_dm/data/tumor_suppressor_oncokb_20250205.txt"

MUTCOL_DEPMAP = {
    "chrom": "Chrom",
    "pos": "Pos",
    "ref": "Ref",
    "alt": "Alt",
    "af": "AF",
    "dp": "DP",
    "ref_count": "RefCount",
    "alt_count": "AltCount",
    "gt": "GT",
    "ps": "PS",
    "variant_type": "VariantType",
    "variant_info": "VariantInfo",
    "dna_change": "DNAChange",
    "protein_change": "ProteinChange",
    "hugo_symbol": "HugoSymbol",
    "exon": "Exon",
    "intron": "Intron",
    "ensembl_gene_id": "EnsemblGeneID",
    "ensembl_feature_id": "EnsemblFeatureID",
    "hgnc_name": "HgncName",
    "hgnc_family": "HgncFamily",
    "uniprot_id": "UniprotID",
    "dbsnp_rs_id": "DbsnpRsID",
    "gc_content": "GcContent",
    "lof_gene_name": "LofGeneName",
    "lof_gene_id": "LofGeneId",
    "lof_number_of_transcripts_in_gene": "LofNumberOfTranscriptsInGene",
    "lof_percent_of_transcripts_affected": "LofPercentOfTranscriptsAffected",
    "nmd": "NMD",
    # 'clnsig': 'ClnSig',
    "molecular_consequence": "MolecularConsequence",
    # 'af_exac': 'AfExac',
    # 'af_esp': 'AfEsp',
    # 'af_tgp': 'AfTgp',
    "vep_impact": "VepImpact",
    "vep_biotype": "VepBiotype",
    "vep_hgnc_id": "VepHgncID",
    "vep_existing_variation": "VepExistingVariation",
    "vep_mane_select": "VepManeSelect",
    "vep_ensp": "VepENSP",
    "vep_swissprot": "VepSwissprot",
    "sift": "Sift",
    "polyphen": "Polyphen",
    "gnomade_af": "GnomadeAF",
    "gnomadg_af": "GnomadgAF",
    "vep_clin_sig": "VepClinSig",
    "vep_somatic": "VepSomatic",
    "vep_pli_gene_value": "VepPliGeneValue",
    "vep_loftool": "VepLofTool",
    "oncogene_high_impact": "OncogeneHighImpact",
    "tumor_suppressor_high_impact": "TumorSuppressorHighImpact",
    # 'achilles_top_genes': 'AchillesTopGenes',
    # 'structural_relation': 'StructuralRelation',
    # 'associated_with': 'AssociatedWith',
    # 'lof': 'Lof',
    # 'driver': 'Driver',
    # 'likely_driver': 'LikelyDriver',
    "transcript_likely_lof": "TranscriptLikelyLof",
    "brca1_func_score": "Brca1FuncScore",
    "civic_id": "CivicID",
    "civic_description": "CivicDescription",
    "civic_score": "CivicScore",
    # 'popaf': 'Popaf',
    # 'likely_gof': 'LikelyGof',
    "likely_lof": "LikelyLoF",
    "hess_driver": "HessDriver",
    "hess_signture": "HessSignature",
    "revel_score": "RevelScore",
    "pharmgkb_id": "PharmgkbId",
    "dida_id": "DidaID",
    "dida_name": "DidaName",
    "gwas_disease": "GwasDisease",
    "gwas_pmid": "GwasPmID",
    "gtex_gene": "GtexGene",
    # 'cosmic_tier': 'CosmicTier',
    # 'oncokb_effect': 'OncokbEffect',
    # 'oncokb_hotspot': 'OncokbHotspot',
    # 'oncokb_oncogenic': 'OncokbOncogenic',
    "provean_prediction": "ProveanPrediction",
    # 'segdup': 'Segdup',
    # 'rm': 'Rm',
    "am_class": "AMClass",
    "am_pathogenicity": "AMPathogenicity",
    "rescue": "Rescue",
    "rescue_reason": "RescueReason",
    "DepMap_ID": "DepMap_ID",
}

# columns generated in postprocessing, not in terra
MUTCOL_ADDITIONAL = {"hotspot": "Hotspot"}

MUTCOL_STANDARDMAF = {
    "Chromosome": "Chromosome",
    "Start_Position": "Start_Position",
    "End_Position": "End_Position",
    "Variant_Type": "Variant_Type",
    "Reference_Allele": "Reference_Allele",
    "Tumor_Seq_Allele1": "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2": "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode": "Tumor_Sample_Barcode",
    "Variant_Classification": "Variant_Classification",
    "Protein_Change": "Protein_Change",
    "Hugo_Symbol": "Hugo_Symbol",
    "NCBI_Build": "NCBI_Build",
    "Strand": "Strand",
    "InFrame": "InFrame",
}

FUSION_COLNAME = [
    "FusionName",
    "JunctionReadCount",
    "SpanningFragCount",
    "SpliceType",
    "LeftGene",
    "LeftBreakpoint",
    "RightGene",
    "RightBreakpoint",
    "LargeAnchorSupport",
    "FFPM",
    "LeftBreakDinuc",
    "LeftBreakEntropy",
    "RightBreakDinuc",
    "RightBreakEntropy",
    "annots",
]

FUSION_RED_HERRING = [
    "GTEx_recurrent",
    "DGD_PARALOGS",
    "HGNC_GENEFAM",
    "Greger_Normal",
    "Babiceanu_Normal",
    "ConjoinG",
    "NEIGHBORS",
]

FUSION_MAXFREQ = 0.1

FUSION_MINFFPM = 0.05

FUSION_MAXFFPM = 0.1

SAMPLESETNAME_STRANDED = "all_stranded"

STARBAMCOLTERRA = ["internal_bam_filepath", "internal_bai_filepath"]

RSEM_TRANSCRIPTS = ["rsem_transcripts_expected_count", "rsem_transcripts_tpm"]

RSEMFILENAME_GENE = ["genes_tpm", "genes_expected_count", "genes_effective_length"]

RSEMFILENAME_GENE_STRANDED = [
    "genes_tpm_stranded",
    "genes_expected_count_stranded",
    "genes_effective_length_stranded",
]

PROTEINEFILENAMES = ["proteincoding_genes_tpm", "proteincoding_genes_expected_count"]

RSEMFILENAME_TRANSCRIPTS = ["transcripts_tpm", "transcripts_expected_count"]

RSEMFILENAME_TRANSCRIPTS_STRANDED = [
    "transcripts_tpm_stranded",
    "transcripts_expected_count_stranded",
]

RNASEQC2_GENE_COUNT_COL = "rnaseqc2_gene_counts"

RNASEQC2_GENE_COUNT_COL_STRANDED = "rnaseqc2_gene_counts_stranded"

RSEMFILENAME = [
    "genes_tpm",
    "genes_expected_count",
    "transcripts_tpm",
    "transcripts_expected_count",
]

SSGSEAFILEPATH = "data/genesets/msigdb.v7.2.symbols.gmt"

RNAMINSIMI = 0.95

RNASEQC_THRESHOLDS_LOWQUAL = {
    "minmapping": 0.85,
    "minendmapping": 0.75,
    "minefficiency": 0.75,
    "maxendmismatch": 0.02,
    "maxmismatch": 0.02,
    "minhighqual": 0.8,
    "minexon": 0.7,
    "maxambiguous": 0.05,
    "maxsplits": 0.1,
    "maxalt": 0.2,
    "maxchim": 0.05,
    "minreads": 20000000,
    "minlength": 80,
    "maxgenes": 35000,
    "mingenes": 12000,
}

RNASEQC_THRESHOLDS_FAILED = {
    "minmapping": 0.7,
    "minendmapping": 0.66,
    "minefficiency": 0.6,
    "maxendmismatch": 0.02,
    "maxmismatch": 0.02,
    "minhighqual": 0.7,
    "minexon": 0.66,
    "maxambiguous": 0.1,
    "maxsplits": 0.1,
    "maxalt": 0.5,
    "maxchim": 0.2,
    "minreads": 20000000,
    "minlength": 80,
    "maxgenes": 35000,
    "mingenes": 10000,
}

README_folder = "../depmap-release-readmes/"

README_currentfolder = "../depmap-release-readmes/release-22q4/"

README_changes = "\n\n\n"

SKIP_UPLOADING_README = True

FUSIONreadme = "\n# Fusions\n\nPORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.\n\n/!\\ This is the most up to date version of the CCLE CN data.\n\n## Annotations\n\nDescription: Gene fusions derived from RNAseq data.\n\nRows: cell lines, IDs contained in the column DepMap_ID\n\nUnfiltered data contains all output fusions, while the filtered data uses the filters suggested by the star fusion docs. These filters are:\n- FFPM > 0.1 -  a cutoff of 0.1 means&nbsp;at least 1 fusion-supporting RNAseq fragment per 10M total reads\n- Remove known false positives, such as GTEx recurrent fusions and certain paralogs\n- Genes that are next to each other\n- Fusions with mitochondrial breakpoints\n- Removing fusion involving mitochondrial chromosomes or HLA genes\n- Removed common false positive fusions (red herring annotations as described in the STAR-Fusion docs)\n- Recurrent fusions observed in CCLE across cell lines (in more than 10% of our samples)\n- Removed fusions where SpliceType='INCL_NON_REF_SPLICE' and LargeAnchorSupport='NO_LDAS' and FFPM < 0.1\n- FFPM < 0.05\n"

RNAseqreadme = "\n# RNAseq\n\nPORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.\n\n/!\\ This is the most up to date version of the CCLE RNA data.\n\n## Annotations:\n\ntranscriptions (Transcripts rpkm):\n\ngenes (gene rpkm):\n__Rows__:\n__Columns__:\nCounts (gene counts):\n__Rows__:\n__Columns__:\nGene level CN data:\n__Rows__:\n__Columns__:\n DepMap cell line IDs\n gene names in the format HGNC\\_symbol (Entrez\\_ID)\nDepMap\\_ID, Chromosome, Start, End, Num\\_Probes, Segment\\_Mean\n"

CNreadme = "\n# Copy Number\n\nPORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.\n\n/!\\ This is the most up to date version of the CCLE CN data.\n\n# Notations:\n\nall: everything\n\nallWES: all data comes from the WExomeS samples we posses\n\nallWGS: all data comes from the WGenomeS samples we posses\n\nwithreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.\n\nmerged: everything from both WGS and WES\n\nlatest: only the latest sequencing versions of the samples were kept\n\n\nGene level CN data:\n\n__Rows__: cell line IDs\n\n__Columns__: gene names in the format HGNC\\_symbol (Entrez\\_ID)\n\nSegment level data:\n\n__Columns__: DepMap\\_ID, Chromosome, Start, End, Segment\\_Mean, Num\\_Probes, Calls"

Achillesreadme = "\n# Copy Number\n\nCombined segment and gene-level CN calls from Broad WES, Sanger WES, and Broad SNP. Relative CN, log2(x+1) transformed.\n\nPORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal. Must use subsetted dataset instead. These data will not make it on the portal starting 19Q1. With the DMC portal, there is new cell line release prioritization as to which lines can be included, so a new taiga dataset will be created containing CN for the portal.\n\nThese data are generated for Achilles to pull from to run CERES.\n\n\nGene level CN data:\n\n__Rows__: DepMap cell line IDs\n\n__Columns__: gene names in the format HGNC\\_symbol (Entrez\\_ID)\n\nSegment level data:\n\n__Columns__: DepMap\\_ID, Chromosome, Start, End, Num\\_Probes, Segment\\_Mean"

Mutationsreadme = "\n# Mutations\n\nPORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.\n\n/!\\ This is the most up to date version of the CCLE Mutatios data.\nThe data is most likely of a better quality that what is on other folder. It is however in beta version as not all changes have either been confirmed or accepted by the DepMap Ops and the DepMap Portal Team.\n\n# Notations:\n\nall: every cell lines we have\n\nWES: all data comes from the WExomeS samples we posses\n\nWGS: all data comes from the WGenomeS samples we posses\n\nwithreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.\n\nmerged: everything from both WGS and WES\n\nlatest: only the latest sequencing versions of the samples were kept\n\ngenes (gene rpkm):\n__Rows__:\n__Columns__:\nCounts (gene counts):\n__Rows__:\n__Columns__:\nGene level CN data:\n__Rows__:\n__Columns__:\n DepMap cell line IDs\n gene names in the format HGNC\\_symbol (Entrez\\_ID)\nDepMap\\_ID, Chromosome, Start, End, Num\\_Probes, Segment\\_Mean\n "

GUMBO_CLIENT_USERNAME = "szhang"

DATE_COL_DICT = {
    "internal": "InternalReleaseDate",
    "dmc": "ConsortiumReleaseDate",
    "public": "PublicReleaseDate",
}

MODEL_TABLE_NAME = "model"

MODEL_TABLE_INDEX = "ModelID"

MC_TABLE_NAME = "model_condition"

MC_TABLE_INDEX = "ModelConditionID"

PR_TABLE_NAME = "omics_profile"

PR_TABLE_INDEX = "ProfileID"

SEQ_TABLE_NAME = "omics_sequencing"

SEQ_TABLE_INDEX = "SequencingID"

SAMPLE_TABLE_NAME = "sample"

ACH_CHOICE_TABLE_COLS = ["ModelConditionID", "ProfileID", "ProfileType"]

ACH_CHOICE_TABLE_NAME = "achilles_choice_table"

DEFAULT_TABLE_COLS = ["ModelID", "ProfileID", "ProfileType"]

DEFAULT_TABLE_NAME = "default_model_table"

PROFILE_TABLE_COLS = ["ModelCondition", "ModelID", "Datatype"]

RELEASE_PR_TABLE_NAME = "OmicsProfiles"

VIRTUAL_FILENAMES_NUMMAT_EXP_MODEL = {
    "proteinCoding_genes_tpm_logp1_profile": "OmicsExpressionProteinCodingGenesTPMLogp1",
    "gene_set_enrichment_profile": "OmicsExpressionGeneSetEnrichment",
}

VIRTUAL_FILENAMES_NUMMAT_CN_MODEL = {
    "merged_gene_cn_profile": "OmicsCNGene",
    "merged_absolute_gene_cn_profile": "OmicsAbsoluteCNGene",
    "merged_loh_profile": "OmicsLoH",
    "globalGenomicFeatures_profile": "OmicsSignatures",
}

VIRTUAL_FILENAMES_NUMMAT_MUT_MODEL = {
    "somaticMutations_genotypedMatrix_hotspot_profile": "OmicsSomaticMutationsMatrixHotspot",
    "somaticMutations_genotypedMatrix_damaging_profile": "OmicsSomaticMutationsMatrixDamaging",
    "somaticMutations_genotypedMatrix_driver_profile": "OmicsSomaticMutationsMatrixDriver",
}

VIRTUAL_FILENAMES_GERMLINE_MODEL = {
    "binary_germline_mutation": "OmicsGuideMutationsBinaryAvana"
}

VIRTUAL_FILENAMES_TABLE_FUSION_MODEL = {
    "fusions_filtered_profile": "OmicsFusionFiltered"
}

VIRTUAL_FILENAMES_TABLE_CN_MODEL = {}

VIRTUAL_FILENAMES_TABLE_MUT_MODEL = {
    "somaticMutations_profile": "OmicsSomaticMutations",
    "structuralVariants_profile": "OmicsStructuralVariants",
}

VIRTUAL_FILENAMES_GUIDEMUT = {
    "binary_mutation_avana": "OmicsGuideMutationsBinaryAvana",
    "binary_mutation_ky": "OmicsGuideMutationsBinaryKY",
    "binary_mutation_humagne": "OmicsGuideMutationsBinaryHumagne",
    "binary_mutation_brunello": "OmicsGuideMutationsBinaryBrunello",
    "binary_mutation_tkov3": "OmicsGuideMutationsBinaryTKOv3"
}

VIRTUAL_FILENAMES_NUMMAT_EXP_PR = {
    "genes_expectedCount_profile": "OmicsExpressionGenesExpectedCountProfile",
    "transcripts_expectedCount_profile": "OmicsExpressionTranscriptsExpectedCountProfile",
    "gene_set_enrichment_profile": "OmicsExpressionGeneSetEnrichmentProfile",
    "rnaseqc_count_mat_profile": "OmicsExpressionRNASeQCGeneCountProfile",
}

VIRTUAL_FILENAMES_NUMMAT_CN_PR = {
    "globalGenomicFeatures_profile": "OmicsSignaturesProfile"
}

VIRTUAL_FILENAMES_NUMMAT_MUT_PR = {}

VIRTUAL_FILENAMES_GERMLINE_PR = {}

VIRTUAL_FILENAMES_TABLE_FUSION_PR = {
    "fusions_unfiltered_profile": "OmicsFusionUnfilteredProfile"
}

VIRTUAL_FILENAMES_TABLE_CN_PR = {
    "merged_segments_profile": "OmicsCNSegmentsProfile",
    "merged_absolute_segments_profile": "OmicsAbsoluteCNSegmentsProfile",
}

VIRTUAL_FILENAMES_TABLE_MUT_PR = {
    "somaticMutations_profile": "OmicsSomaticMutationsProfile",
    "structuralVariants_profile": "OmicsStructuralVariantsProfile",
}

DEPMAP_PV = "https://docs.google.com/spreadsheets/d/1uqCOos-T9EMQU7y2ZUw4Nm84opU5fIT1y7jet1vnScE"

TAIGA_FP_FILENAME = "fingerprint_lod_matrix"

TAIGA_LEGACY_CN = "copy-number-5f61"

POTENTIAL_LIST = "https://docs.google.com/spreadsheets/d/1YuKEgZ1pFKRYzydvncQt9Y_BKToPlHP-oDB-0CAv3gE"

DEPMAP_TAIGA = "arxspan-cell-line-export-f808"

TAIGA_ETERNAL = "depmap-a0ab"

HGNC_MAPPING_TABLE_TAIGAID = "hgnc-gene-table-e250"

HGNC_MAPPING_TABLE_NAME = "hgnc_complete_set"

HGNC_MAPPING_TABLE_VERSION = 3


# OMICS_ID_MAPPING_TABLE_VERSION = 1

# generated by copynumbers.maskGenes()
MASKED_GENE_LIST = "data/genes_to_mask_minus_rescue_ensg.txt"
