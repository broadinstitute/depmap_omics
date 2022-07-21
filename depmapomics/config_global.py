import numpy as np


LINES_TO_RELEASE = []
SAMPLESETNAME = "CCLF"
RELEASE = SAMPLESETNAME.lower()

CACHE_PATH = "~/.depmapomics/"
TMP_PATH = "/tmp/"

### GCP credentials
SHEETCREDS = "../../.credentials.json"
MY_ID = "~/.client_secret.json"
MYSTORAGE_ID = "~/.storage.json"
GCS_PAYER_PROJECT = "broad-firecloud-ccle"

SAMPLEID = "DepMap_ID"

# from most prioritized to least prioritized:
SOURCE_PRIORITY = ["BROAD", "DEPMAP", "IBM", "CCLE2", "SANGER", "CHORDOMA"]

ENSEMBL_SERVER_V = "http://nov2020.archive.ensembl.org/biomart"


############## TERRA

HG38BAMCOL = ["internal_bam_filepath", "internal_bai_filepath"]
LEGACY_BAM_COLNAMES = ["legacy_bam_filepath", "legacy_bai_filepath"]

############## LOADING

TORAISE = (["ACH-001195"],)

TO_UPDATE = {
    "primary_disease": ["Primary Disease"],
    "sex": ["CCLF Gender"],
    "collection_site": ["Sample Collection Site"],
    "subtype": ["lineage_subtype"],
    "subsubtype": ["lineage_sub_subtype"],
    "lineage": ["lineage"],
    # 'source': ['Flagship'],
    "parent_cell_line": ["Parental ID"],
    "comments": ["Comments"],
    # "mediatype": ["Culture Medium", "Culture Type"],
    "stripped_cell_line_name": ["Stripped Cell Line Name"],
    "cellosaurus_id": ["RRID"],
    "age": ["CCLF Age"],
}

MAXAGE = "2021-01-01"

EXTRACT_TO_CHANGE = {"from_arxspan_id": "participant"}

REPLACE = {"T": "Tumor", "N": "Normal", "m": "Unknown", "L": "Unknown"}

values = (["legacy_bam_filepath", "legacy_bai_filepath"],)

filetypes = ["bam", "bai"]

MATCH = ["ACH-", "CDS-"]

# default values in the GP workspaces and our sample tracker (to change if you use another workspace/
# sample tracker)
EXTRACT_DEFAULTS = {
    "name": "sample_alias",
    "bai": "crai_or_bai_path",
    "bam": "cram_or_bam_path",
    "ref_bam": "hg19_bam_filepath",
    "ref_type": "Datatype",
    "ref_bai": "hg19_bai_filepath",
    "version": "version",
    "primary_disease": "primary_disease",
    "ref_arxspan_id": "arxspan_id",
    "ref_name": "stripped_cell_line_name",
    "source": "source",
    "size": "size",
    "legacy_size": "legacy_size",
    "from_arxspan_id": "individual_alias",
    "ref_id": "sample_id",
    "PDO_id_terra": "PDO",
    "PDO_id_gumbo": "PDO-ID",
    "update_time": "update_time",
    "from_patient_id": "individual_alias",
    "patient_id": "participant_id",
    "ref_date": "date_sequenced",
    "hs_hs_library_size": "hs_hs_library_size",
    "hs_het_snp_sensitivity": "hs_het_snp_sensitivity",
    "hs_mean_bait_coverage": "hs_mean_bait_coverage",
    "hs_mean_target_coverage": "hs_mean_target_coverage",
    "hs_on_target_bases": "hs_on_target_bases",
    "total_reads": "total_reads",
    "release_date": "sequencing_date",
    "hash": "crc32c_hash",
    "legacy_hash": "legacy_crc32c_hash",
    "mean_depth": "mean_depth",
    "root_sample_id": "root_sample_id",
    "sm_id": "SM-ID",
    "profile_id": "ProfileID",
    "expected_type": "ExpectedType",
}

# minimum bam file size in bytes for each sequencing type
MINSIZES = {
    "rna": 2_000_000_000,
    "wes": 3_000_000_000,
    "wgs": 50_000_000_000,
}

# known cell lines that are from the same patient (not called?)
# samepatient = [
#     ["ACH-000635", "ACH-000717", "ACH-000864", "ACH-001042", "ACH-001547"],
#     ["ACH-002291", "ACH-001672"],
#     ["ACH-001706", "ACH-001707"],
# ]

# known duplicate arxspan-ids
DUP_ARXSPANS = {"ACH-001620": "ACH-001605", "ACH-001621": "ACH-001606"}

# rename ccle_name TODO: ask becky what to do
# (not called?)
# rename = {"PEDS117": "CCLFPEDS0009T"}

## old GP storage buckets
# rnaworkspace2 = "broad-firecloud-ccle/CCLE_DepMap_RNAseq"
# rnaworkspace4 = "broad-genomics-delivery/Cancer_Cell_Line_Factory_CCLF_RNAseq"
# rnaworkspace5 = "nci-mimoun-bi-org/CCLF_RNA_2_0"
# rnaworkspace3 = "broad-genomics-delivery/CCLE_DepMap_RNAseq"
# rnaworkspace1 = "broad-genomics-delivery/Getz_IBM_CellLines_RNASeqData"

## curent GP buckets
rnaworkspace6 = "terra-broad-cancer-prod/CCLE_DepMap_RNAseq"
rnaworkspace7 = "terra-broad-cancer-prod/Getz_IBM_CellLines_RNASeqData"

## and their correesponding sample source
# rnasource1 = "ibm"
# rnasource2 = "ccle"
# rnasource3 = "ccle"
# rnasource4 = "cclf"
# rnasource5 = "cclf"

rnasource6 = "DEPMAP"
rnasource7 = "IBM"

RNAWORKSPACES = {
    "DEPMAP": "terra-broad-cancer-prod/CCLE_DepMap_RNAseq",
    "IBM": "terra-broad-cancer-prod/Getz_IBM_CellLines_RNASeqData",
}

## curent WGS GP buckets
wgsworkspace1 = "terra-broad-cancer-prod/DepMap_WGS"
wgsworkspace2 = "terra-broad-cancer-prod/Getz_IBM_CellLines_WGS"

## and their corresponding sample source
wgssource1 = "DEPMAP"
wgssource2 = "IBM"

WGSWORKSPACES = {
    "DEPMAP": "terra-broad-cancer-prod/DepMap_WGS",
    "IBM": "terra-broad-cancer-prod/Getz_IBM_CellLines_WGS",
}

WGSSETENTITY = "sample_set"
WESSETENTITY = "pair_set"

############## Fingerprinting

FPALLBATCHPAIRSETS = "all"

# Local directory to save intermediate files to
WORKING_DIR = "output/"

PREV_VIRTUAL = {}

DATASETS = ["internal", "ibm", "dmc", "public"]
RUN_NOTEBOOKS = ["WGS_CCLE.ipynb", "RNA_CCLE.ipynb"]
# 20Q3
# PREV_VIRTUAL={}
# PREV_VIRTUAL['public'] = 'public-20q3-3d35'
# PREV_VIRTUAL['dmc'] = 'dmc-20q3-deprecated-never-released--5f55'
# PREV_VIRTUAL['internal'] = 'internal-20q3-00d0'

# 20Q4
# PREV_VIRTUAL={}
# PREV_VIRTUAL['public'] = 'public-20q4-a4b3'
# PREV_VIRTUAL['dmc'] = 'dmc-20q4-fcf4'
# PREV_VIRTUAL['ibm'] = 'ibm-20q4-269f'
# PREV_VIRTUAL['internal'] = 'internal-20q4-2540'

# 21Q1
# PREV_VIRTUAL['public'] = 'public-21q1-4b39'
# PREV_VIRTUAL['ibm'] = 'ibm-21q1-abd9'
# PREV_VIRTUAL['dmc'] = 'dmc-21q1-0e11'
# PREV_VIRTUAL['internal'] = 'internal-21q1-4fc4'

# 21Q2
# PREV_VIRTUAL['public'] = 'public-21q2-110d'
# PREV_VIRTUAL['ibm'] = 'ibm-21q2-9ed1'
# PREV_VIRTUAL['dmc'] = 'dmc-21q2-27e1'
# PREV_VIRTUAL['internal'] = 'internal-21q2-9d16'

# 21Q3
# PREV_VIRTUAL["public"] = "public-21q3-bf1e"
# PREV_VIRTUAL["ibm"] = "ibm-21q3-179f"
# PREV_VIRTUAL["dmc"] = "dmc-21q3-482c"
# PREV_VIRTUAL["internal"] = "internal-21q3-fe4c"

# 21Q4
# PREV_VIRTUAL["public"] = "public-21q4-a0d6"
# PREV_VIRTUAL["ibm"] = "ibm-21q4-4e18"
# PREV_VIRTUAL["dmc"] = "dmc-21q4-5725"
# PREV_VIRTUAL["internal"] = "internal-21q4-ac0a"

# 21Q4v2
# PREV_VIRTUAL["public"] = "public-21q4v2-103d"
# PREV_VIRTUAL["ibm"] = "ibm-21q4v2-2d92"
# PREV_VIRTUAL["dmc"] = "dmc-21q4v2-0e7c"
# PREV_VIRTUAL["internal"] = "internal-21q4v2-403b"

# 22Q1
PREV_VIRTUAL["public"] = "public-22q1-305b"
PREV_VIRTUAL["ibm"] = "ibm-22q1-cce1"
PREV_VIRTUAL["dmc"] = "dmc-22q1-d00a"
PREV_VIRTUAL["internal"] = "internal-22q1-1778"


############## DNAseq

WGS_methods = [
    "gatk/PreProcessingForVariantDiscovery_GATK4/8",
    "GP-TAG/Manta_SomaticSV/9",
    "gkugener/ArrayOfFilesToTxt/1",
    "vdauwera/BamToUnmappedRGBams/4",
    "gatk/CNV_Somatic_Pair_Workflow/9",
    "gkugener/Aggregate_CN_seg_files/2",
]

CNWES_methods = [
    "gatk/PreProcessingForVariantDiscovery_GATK4/8",
    "GP-TAG/Manta_SomaticSV/9",
    "gkugener/ArrayOfFilesToTxt/1",
    "vdauwera/BamToUnmappedRGBams/4",
    "gatk/CNV_Somatic_Pair_Workflow/9",
    "gkugener/Aggregate_CN_seg_files/2",
]

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

KNOWN_DROP = ["CDS-R22IHj", "CDS-xMnTwN", "CDS-2FC7DW", "CDS-ToOF9G"]

# germline matrix
WESVCFDIR = "/tmp/vcfs_wes/"
WGSVCFDIR = "/tmp/vcfs_wgs/"
GUIDESBED = "../data/avana_guides.bed"

# rescue certain lines that are blacklisted in the tracker but we want them for MUTATION ONLY
RESCUE_FOR_MUTATION_WES = {
    "CDS-Rl87Z1": "ACH-001956",
    "CDS-mys9Dm": "ACH-001955",
    "CDS-TzQAjG": "ACH-001957",
    "CDS-TuKZau": "ACH-001709",
    "CDS-4lWqEA": "ACH-000859",
    "CDS-Fyjj8I": "ACH-000116",
}

RESCUE_FOR_MUTATION_WGS = {"CDS-AqZLna": "ACH-002512"}

############## CN

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
    "Sample": SAMPLEID,
    "NUM_POINTS_COPY_RATIO": "Num_Probes",
    "MEAN_LOG2_COPY_RATIO": "Segment_Mean",
    "CALL": "Status",
}

PURECN_COLRENAMING = {
    "start": "Start",
    "end": "End",
    "chr": "Chromosome",
    "Sampleid": SAMPLEID,
    "seg.mean": "Segment_Mean",
    "type": "LOH_status",
    "C": "Absolute_CN",
}

PURECN_LOH_COLNAME = "PureCN_loh_merged"
PURECN_FAILED_COLNAME = "PureCN_failed"

# if the loh function outputs one of these, record as 1 in the loh bool matrix
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
    "PureCN_curated",
    "PureCN_curated_solution",
    "PureCN_cin_allele_specific_ploidy_robust",
}
MISC_SIG_TERRACOLS = {"msisensor2_score"}
PURECN_FAILED_COLNAME = "PureCN_failed"
SIGTABLE_BINARYCOLS = [
    "PureCN_wgd",
    "PureCN_curated",
]

SIGTABLE_RENAMING = {
    "PureCN_loh_fraction": "loh_fraction",
    "PureCN_wgd": "wgd",
    "PureCN_cin_allele_specific_ploidy_robust": "cin",
    "PureCN_ploidy": "ploidy",
}

SOURCE_RENAME = {
    "CCLF": "Broad",
    "CHORDOMA": "Chordoma",
    "SANGER": "Sanger",
    "IBM": "Broad",
    "CCLE2": "Broad",
    np.nan: "Broad",
    "DEPMAP": "Broad",
    "IBM WES": "Broad WES",
    "Broad CCLF": "Broad WES",
}

wrongwes_arxspan = {
    "ACH-001189",
    "ACH-002303",
    "ACH-002315",
    "ACH-002341",
    "ACH-001011",
    "ACH-001108",
    "ACH-001187",
    "ACH-002875",
    "ACH-002874",
    "ACH-001955",  # chordoma lines
    "ACH-001956",
    "ACH-001957",
}

CYTOBANDLOC = "data/hg38_cytoband.gz"

toreprocess = ["CDS-C2RlCj", "CDS-8GqFo5"]

############## Mutations

MINFREQTOCALL = 0.25

SV_COLNAME = "filtered_annotated_sv"
SV_FILENAME = "all_sv.bedpe"

MUTATION_GROUPS = {
    "other conserving": ["5'Flank", "Intron", "IGR", "3'UTR", "5'UTR"],
    "other non-conserving": [
        "In_Frame_Del",
        "In_Frame_Ins",
        "Stop_Codon_Del",
        "Stop_Codon_Ins",
        "Missense_Mutation",
        "Nonstop_Mutation",
    ],
    "silent": ["Silent"],
    "damaging": [
        "De_novo_Start_OutOfFrame",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Splice_Site",
        "Start_Codon_Del",
        "Start_Codon_Ins",
        "Start_Codon_SNP",
        "Nonsense_Mutation",
    ],
}

MUTCOL_DEPMAP = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "NCBI_Build",
    "Chromosome",
    "Start_position",
    "End_position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Allele",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Genome_Change",
    "Annotation_Transcript",
    SAMPLEID,
    "cDNA_Change",
    "Codon_Change",
    "Protein_Change",
    "isDeleterious",
    "isTCGAhotspot",
    "TCGAhsCnt",
    "isCOSMIChotspot",
    "COSMIChsCnt",
    "ExAC_AF",
    "Variant_annotation",
    "CGA_WES_AC",
]


############## FUSION

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

############## EXPRESSION

STARBAMCOLTERRA = ["internal_bam_filepath", "internal_bai_filepath"]

RSEM_TRANSCRIPTS = ["rsem_transcripts_expected_count", "rsem_transcripts_tpm"]

RSEMFILENAME_GENE = ["genes_tpm", "genes_expected_count"]
PROTEINEFILENAMES = ["proteincoding_genes_tpm", "proteincoding_genes_expected_count"]

RSEMFILENAME_TRANSCRIPTS = ["transcripts_tpm", "transcripts_expected_count"]

RSEMFILENAME = RSEMFILENAME_GENE + RSEMFILENAME_TRANSCRIPTS

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

PREVIOUS_QC_FAIL = [
    "CDS-12DTEw",
    "CDS-9hv1zM",
    "CDS-A6GSeQ",
    "CDS-aWlMRt",
    "CDS-B1ywOH",
    "CDS-BixxtG",
    "CDS-DRM3l2",
    "CDS-jOlYT4",
    "CDS-KMhiT9",
    "CDS-M6mnMA",
    "CDS-pYwECX",
    "CDS-v6E624",
    "CDS-vxTqNJ",
    "CDS-YxtmkI",
    "CDS-fk564T",
    "CDS-kU30H5",
    "CDS-G0F5f5",
    "CDS-ABH0uZ",
]

###################### READMEs output

README_folder = "../depmap-release-readmes/"

README_currentfolder = README_folder + "release-" + RELEASE + "/"

README_changes = """


"""

SKIP_UPLOADING_README = True

FUSIONreadme = """
# Fusions

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE CN data.

## Annotations

Description: Gene fusions derived from RNAseq data.

Rows: cell lines, IDs contained in the column DepMap_ID

Unfiltered data contains all output fusions, while the filtered data uses the filters suggested by the star fusion docs. These filters are:
- FFPM > 0.1 -  a cutoff of 0.1 means&nbsp;at least 1 fusion-supporting RNAseq fragment per 10M total reads
- Remove known false positives, such as GTEx recurrent fusions and certain paralogs
- Genes that are next to each other
- Fusions with mitochondrial breakpoints
- Removing fusion involving mitochondrial chromosomes or HLA genes
- Removed common false positive fusions (red herring annotations as described in the STAR-Fusion docs)
- Recurrent fusions observed in CCLE across cell lines (in more than 10% of our samples)
- Removed fusions where SpliceType='INCL_NON_REF_SPLICE' and LargeAnchorSupport='NO_LDAS' and FFPM < 0.1
- FFPM < 0.05
"""


RNAseqreadme = """
# RNAseq

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE RNA data.

## Annotations:

transcriptions (Transcripts rpkm):

genes (gene rpkm):
__Rows__:
__Columns__:
Counts (gene counts):
__Rows__:
__Columns__:
Gene level CN data:
__Rows__:
__Columns__:
 DepMap cell line IDs
 gene names in the format HGNC\_symbol (Entrez\_ID)
DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean
"""


CNreadme = """
# Copy Number

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE CN data.

# Notations:

all: everything

allWES: all data comes from the WExomeS samples we posses

allWGS: all data comes from the WGenomeS samples we posses

withreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.

merged: everything from both WGS and WES

latest: only the latest sequencing versions of the samples were kept


Gene level CN data:

__Rows__: cell line IDs

__Columns__: gene names in the format HGNC\_symbol (Entrez\_ID)

Segment level data:

__Columns__: DepMap\_ID, Chromosome, Start, End, Segment\_Mean, Num\_Probes, Calls"""


Achillesreadme = """
# Copy Number

Combined segment and gene-level CN calls from Broad WES, Sanger WES, and Broad SNP. Relative CN, log2(x+1) transformed.

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal. Must use subsetted dataset instead. These data will not make it on the portal starting 19Q1. With the DMC portal, there is new cell line release prioritization as to which lines can be included, so a new taiga dataset will be created containing CN for the portal.

These data are generated for Achilles to pull from to run CERES.


Gene level CN data:

__Rows__: DepMap cell line IDs

__Columns__: gene names in the format HGNC\_symbol (Entrez\_ID)

Segment level data:

__Columns__: DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean"""


Mutationsreadme = """
# Mutations

PORTAL TEAM SHOULD NOT USE THIS: There are lines here that should not make it even to internal.

/!\ This is the most up to date version of the CCLE Mutatios data.
The data is most likely of a better quality that what is on other folder. It is however in beta version as not all changes have either been confirmed or accepted by the DepMap Ops and the DepMap Portal Team.

# Notations:

all: every cell lines we have

WES: all data comes from the WExomeS samples we posses

WGS: all data comes from the WGenomeS samples we posses

withreplicates: if we have two different sequencing from a sample, we kept both, see the depmap sample tracker for annotations [https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE](https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE). this dataset is more geared toward QC or in-depth analysis of a particular cell line.

merged: everything from both WGS and WES

latest: only the latest sequencing versions of the samples were kept

genes (gene rpkm):
__Rows__:
__Columns__:
Counts (gene counts):
__Rows__:
__Columns__:
Gene level CN data:
__Rows__:
__Columns__:
 DepMap cell line IDs
 gene names in the format HGNC\_symbol (Entrez\_ID)
DepMap\_ID, Chromosome, Start, End, Num\_Probes, Segment\_Mean
 """

########### Gumbo configs #############

# names of columns that contain release dates
DATE_COL_DICT = {
    "internal": "InternalReleaseDate",
    "ibm": "IBMReleaseDate",
    "dmc": "ConsortiumReleaseDate",
    "public": "PublicReleaseDate",
}

MODEL_TABLE_NAME = "Model"
MC_TABLE_NAME = "ModelCondition"
PR_TABLE_NAME = "OmicsProfile"
SEQ_TABLE_NAME = "OmicsSequencing"
SAMPLE_TABLE_NAME = "Sample"
ACH_CHOICE_TABLE_COLS = ["ModelConditionID", "ProfileID", "ProfileType"]
ACH_CHOICE_TABLE_NAME = "achilles_choice_table"
DEFAULT_TABLE_COLS = ["ModelID", "ProfileID", "ProfileType"]
DEFAULT_TABLE_NAME = "default_model_table"

# dictionaries mapping names of pr-level matrices on latest to names on virtual
# NumericMatrixCSV matrices:
VIRTUAL_FILENAMES_NUMMAT_EXP = {
    "genes_expectedCount_profile": "Omics_expression_genes_expectedCount",
    "genes_tpm_logp1_profile": "Omics_expression_genes_tpm_logp1",
    "transcripts_expectedCount_profile": "Omics_expression_transcripts_expectedCount",
    "transcripts_tpm_logp1_profile": "Omics_expression_transcripts_tpm_logp1",
    "proteinCoding_genes_expectedCount_profile": "Omics_expression_proteinCoding_genes_expectedCount",
    "proteinCoding_genes_tpm_logp1_profile": "Omics_expression_proteinCoding_genes_tpm_logp1",
}

VIRTUAL_FILENAMES_NUMMAT_EXP_INTERNAL = {
    "gene_set_enrichment_profile": "Omics_gene_set_enrichment",
}

VIRTUAL_FILENAMES_NUMMAT_CN = {
    "merged_gene_cn_profile": "Omics_merged_gene_cn",
}

VIRTUAL_FILENAMES_NUMMAT_MUT = {
    "somaticMutations_boolMatrix_hotspot_profile": "Omics_somaticMutations_boolMatrix_hotspot",
    "somaticMutations_boolMatrix_otherNonconserving_profile": "Omics_somaticMutations_boolMatrix_otherNonconserving",
    "somaticMutations_boolMatrix_damaging_profile": "Omics_somaticMutations_boolMatrix_damaging",
    "somaticMutations_boolMatrix_otherConserving_profile": "Omics_somaticMutations_boolMatrix_otherConserving",
}

VIRTUAL_FILENAMES_GERMLINE = {
    "binary_germline_mutation": "Omics_binary_germline_mutation",
}

# TableCSV matrices:
VIRTUAL_FILENAMES_TABLE_FUSION = {
    "fusions_filtered_profile": "Omics_fusions_filtered",
    "fusions_unfiltered_profile": "Omics_fusions_unfiltered",
}

VIRTUAL_FILENAMES_TABLE_CN = {
    "merged_segments_profile": "Omics_merged_segments",
}

VIRTUAL_FILENAMES_TABLE_MUT = {
    "somaticMutations_profile": "Omics_somaticMutations",
}
