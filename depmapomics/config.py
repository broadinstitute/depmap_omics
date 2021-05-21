CACHE_PATH = '~/.depmapomics/'
TMP_PATH = '/tmp/'
ENSEMBL_SERVER_V = "http://nov2020.archive.ensembl.org/biomart"

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
