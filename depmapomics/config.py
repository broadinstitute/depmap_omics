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
