import numpy as np

########################## GENERIC PARAMS

CACHE_PATH = '~/.depmapomics/'
TMP_PATH = '/tmp/'
ENSEMBL_SERVER_V = "http://nov2020.archive.ensembl.org/biomart"

SHEETCREDS = '../.credentials.json'
MY_ID = '~/.client_secret.json',
MYSTORAGE_ID = "~/.storage.json"

SHEETNAME='ccle sample tracker'

TAIGA_ETERNAL = 'depmap-a0ab'

REFSHEET_URL = "https://docs.google.com/spreadsheets/d/1Pgb5fIClGnErEqzxpU7qqX6ULpGTDjvzWwDN8XUJKIY"

SAMPLEID="DepMap_ID"

############## TERRA

HG38BAMCOL = ['internal_bam_filepath',
              "internal_bai_filepath"]


############## LOADING



############## CN

COLRENAMING = {'CONTIG': 'Chromosome',
               'START': 'Start',
               'END': 'End',
               'end': 'End',
               'seqnames': 'Chromosome',
               'start': 'Start',
               'Sample': SAMPLEID,
               'NUM_POINTS_COPY_RATIO': 'Num_Probes',
               'MEAN_LOG2_COPY_RATIO': 'Segment_Mean',
               'CALL': 'Status'}

SOURCE_RENAME = {'CCLF': 'Broad WES', 'CHORDOMA': 'Chordoma WES',
                'SANGER': 'Sanger WES', 'IBM': 'Broad WES',
                np.nan: 'Broad WES', 'DEPMAP': 'Broad WES',
                'IBM WES': "Broad WES", 'Broad CCLF': "Broad WES"}

PROCQC = ["allelic_counts_tumor", "delta_MAD_tumor", "denoised_MAD_tumor",
          "scaled_delta_MAD_tumor", "denoised_copy_ratios_lim_4_plot_tumor",
          "denoised_copy_ratios_plot_tumor", "modeled_segments_plot_tumor",
          "gatk_cnv_all_plots", "lego_plotter_pngs", "copy_number_qc_report",
          "ffpe_OBF_figures", "mut_legos_html", "oxoG_OBF_figures",
          "tumor_bam_base_distribution_by_cycle_metrics",
          "tumor_bam_converted_oxog_metrics"]

BAMQC = ["duplication_metrics", "bqsr_report",
         "tumor_bam_alignment_summary_metrics",
         "tumor_bam_bait_bias_summary_metrics",
         "tumor_bam_gc_bias_summary_metrics",
         "tumor_bam_hybrid_selection_metrics",
         "tumor_bam_insert_size_histogram",
         "tumor_bam_insert_size_metrics",
         "tumor_bam_pre_adapter_summary_metrics",
         "tumor_bam_quality_by_cycle_metrics",
         "tumor_bam_quality_distribution_metrics",
         "tumor_bam_quality_yield_metrics"]

############## Mutations

MUTATION_GROUPS = {
    "other conserving": ["5'Flank", "Intron", "IGR", "3'UTR", "5'UTR"],
    "other non-conserving": ["In_Frame_Del", "In_Frame_Ins", "Stop_Codon_Del",
                             "Stop_Codon_Ins", "Missense_Mutation", "Nonstop_Mutation"],
    'silent': ['Silent'],
    "damaging": ['De_novo_Start_OutOfFrame', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                 'Splice_Site', 'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP', 'Nonsense_Mutation']
}

MUTCOL_DEPMAP = ['Hugo_Symbol', 'Entrez_Gene_Id', 'NCBI_Build', 'Chromosome',
                 'Start_position', 'End_position', 'Strand', 'Variant_Classification',
                 'Variant_Type', 'Reference_Allele', 'Tumor_Allele', 'dbSNP_RS',
                 'dbSNP_Val_Status', 'Genome_Change', 'Annotation_Transcript',
                 SAMPLEID, 'cDNA_Change', 'Codon_Change', 'Protein_Change', 'isDeleterious',
                 'isTCGAhotspot', 'TCGAhsCnt', 'isCOSMIChotspot', 'COSMIChsCnt',
                 'ExAC_AF', "Variant_annotation", 'CGA_WES_AC', 'HC_AC',
                 'RD_AC', 'RNAseq_AC', 'SangerWES_AC', 'WGS_AC']


############## FUSION

FUSION_COLNAME = ['FusionName', 'JunctionReadCount',
                  'SpanningFragCount', 'SpliceType', 'LeftGene', 'LeftBreakpoint',
                  'RightGene', 'RightBreakpoint', 'LargeAnchorSupport', 'FFPM',
                  'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc',
                  'RightBreakEntropy', 'annots']

FUSION_RED_HERRING = ['GTEx_recurrent', 'DGD_PARALOGS', 'HGNC_GENEFAM',
                      'Greger_Normal', 'Babiceanu_Normal', 'ConjoinG', 'NEIGHBORS']

############## EXPRESSION

RNAGSPATH38="gs://cclebams/rnasq_hg38/"

STARBAMCOLTERRA = [ "star_bam_file", 'star_bam_index']

RSEM_TRANSCRIPTS = ['rsem_transcripts_expected_count',
                    'rsem_transcripts_tpm']

RSEMFILENAME_GENE=["genes_tpm",
"genes_expected_count"]

RSEMFILENAME_TRANSCRIPTS=["transcripts_tpm", "transcripts_expected_count"]

RSEMFILENAME = RSEMFILENAME_GENE+RSEMFILENAME_TRANSCRIPTS

SSGSEAFILEPATH = "data/genesets/msigdb.v7.2.symbols.gmt"

PATHTOGENEPY = "../"

RNASEQC_THRESHOLDS_LOWQUAL = {'minmapping': 0.85, 'minendmapping': 0.75, 'minefficiency': 0.75,
                              'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.8,
                              'minexon': 0.7, "maxambiguous": 0.05, "maxsplits": 0.1,
                              "maxalt": 0.2, "maxchim": 0.05, "minreads": 20000000,
                              "minlength": 80, "maxgenes": 35000, "mingenes": 12000}


RNASEQC_THRESHOLDS_FAILED = {'minmapping': 0.7, 'minendmapping': 0.66, 'minefficiency': 0.6,
                             'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.7,
                             'minexon': 0.66, "maxambiguous": 0.1, "maxsplits": 0.1,
                             "maxalt": 0.5, "maxchim": 0.2, "minreads": 20000000,
                             "minlength": 80, "maxgenes": 35000, "mingenes": 10000}


###################### README

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
