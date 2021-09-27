# there are some issues in the older versions of omics data on virtual that this flag deals with
# some columns got renamed in the mutation file and some data was uploaded as tsv instead of csv
from gsheets.api import Sheets


LEGACY_PATCH_FLAGS = {'rename_column': False, 'tsv2csv': False}

# release ids on taiga
TENTATIVE_VIRTUAL = {'name': 'tentative-virtual-d84e', 'version': None}
PORTAL = 'internal'

VIRTUAL_RELEASES = {'21q2': {'public': {'name': 'public-21q2-110d', 'version': 13},
                             'ibm': {'name': 'ibm-21q2-9ed1', 'version': 15},
                             'dmc': {'name': 'dmc-21q2-27e1', 'version': 14},
                             'internal': {'name': 'internal-21q2-9d16', 'version': 17}},
                    '21q3': {'internal': {'name': 'internal-21q3-fe4c', 'version': 12},
                             'ibm': {'name': 'ibm-21q3-179f', 'version': 8},
                             'dmc': {'name': 'dmc-21q3-482c', 'version': 7},
                             'public': {'name': 'public-21q3-bf1e', 'version': 7}}
                    } # release ids on taiga

PREV_RELEASE = VIRTUAL_RELEASES['21q3'][PORTAL]
NEW_RELEASE = TENTATIVE_VIRTUAL


LINES_TO_DROP_COMMON = {'ACH-000010', 'ACH-001078', 'ACH-001146', 'ACH-001173',
                        'ACH-001741', 'ACH-001790', 'ACH-002022', 'ACH-002184',
                        'ACH-002260'}
LINES_TO_DROP_DNA = LINES_TO_DROP_COMMON | {'ACH-002475'}
LINES_TO_DROP_RNA = LINES_TO_DROP_COMMON | {'ACH-000658', 'ACH-001212', 'ACH-001316'}
LINES_TO_DROP = {'DNA': LINES_TO_DROP_DNA, 'RNA': LINES_TO_DROP_RNA}

# these are the columns that if merged with an older release (assuming that old data was not altered),
# should uniquely identify each row of the file to find equal values in each column
FUSIONS_MERGE_COLS = ['DepMap_ID', 'LeftGene', 'RightGene', 'LeftBreakpoint', 'RightBreakpoint']
SEGMENT_CN_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start', 'End']
MUTATIONS_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start_position', 'End_position', 'Tumor_Seq_Allele1']

# if there are new files that were added in the previous release, add them here
FILES_RELEASED_BEFORE = ['CCLE_expression', 'CCLE_expression_proteincoding_genes_expected_count',
                         'CCLE_RNAseq_transcripts', 'CCLE_expression_transcripts_expected_count',
                         'CCLE_expression_full', 'CCLE_RNAseq_reads', 'CCLE_fusions', 'CCLE_fusions_unfiltered',
                         'CCLE_gene_cn', 'CCLE_segment_cn', 'CCLE_mutations']

# correlation thresholds above which we consider two releases as 'similar'
CORRELATION_THRESHOLDS = {'CCLE_gene_cn': 0.99, 'all_expressions': 0.99999}

SKIP_ARXSPAN_COMPARISON = True # set to False if you want to test whether some arxspans were added/removed

PLOTS_OUTPUT_FILENAME_PREFIX = '/tmp/plots_' # location/prefix for saving output plots

LINES_TO_RELEASE_SHEET = 'https://docs.google.com/spreadsheets/d/1-Iz_TrLy2DT2hFZKr6m-GJsVeQbKMA06Uf2RmGRlUdA/edit?usp=sharing'
sheets_obj = Sheets.from_files('~/.client_secret.json', '~/.storage.json')
sheets = sheets_obj.get(LINES_TO_RELEASE_SHEET).sheets
LINES_TO_RELEASE = sheets[0].to_frame(header=0, index_col=None)
LINES_TO_RELEASE.columns = LINES_TO_RELEASE.columns.str.lower()

# all the file attributes
FILE_ATTRIBUTES = [
    {'file': 'CCLE_expression', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_proteincoding_genes_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_transcripts', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_transcripts_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_full', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_reads', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_fusions', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS, 'expected_changed_cols':['CCLE_count']},
    {'file': 'CCLE_fusions_unfiltered', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS, 'expected_changed_cols':['CCLE_count']},
    # {'file': 'CCLE_ssGSEA', 'ismatrix': True, 'hasNA': False,'omicssource':'RNA', 'gene_id': None},
    {'file': 'CCLE_gene_cn', 'ismatrix': True, 'hasNA': True, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_segment_cn', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': SEGMENT_CN_MERGE_COLS, 'expected_changed_cols':[]},
    {'file': 'CCLE_mutations', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': MUTATIONS_MERGE_COLS, 'expected_changed_cols':[]},
    {'file': 'CCLE_mutations_bool_damaging', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_hotspot', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_otherconserving', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_nonconserving', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'}
]

# comment/uncomment to use all/subset of files for testing
# FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x['file'] in ['CCLE_expression', 'CCLE_expression_full'])]
# FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x['file'] in ['CCLE_mutations'])]
FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x['omicssource'] in ['RNA']) and x['ismatrix']]
# FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x['file'] in ['CCLE_fusions', 'CCLE_fusions_unfiltered'])]

# the following information is used to create a tentative virtual
MUTATIONS_TAIGA_ID = 'mutations-latest-ed72'
FUSIONS_TAIGA_ID = 'fusions-95c9'
EXPRESSION_TAIGA_ID = 'expression-d035'
CN_TAIGA_ID = 'cn-achilles-version-06ca'
# CN_TAIGA_ID = 'cn-latest-d8d4'

TAIGA_IDS_LATEST = {
    MUTATIONS_TAIGA_ID:[
        ('CCLE_mutations', 'merged_somatic_mutations_withlegacy'),
        ('CCLE_mutations_bool_damaging', 'merged_somatic_mutations_boolmatrix_fordepmap_damaging'),
        ('CCLE_mutations_bool_nonconserving', 'merged_somatic_mutations_boolmatrix_fordepmap_othernoncons'),
        ('CCLE_mutations_bool_otherconserving', 'merged_somatic_mutations_boolmatrix_fordepmap_othercons'),
        ('CCLE_mutations_bool_hotspot', 'merged_somatic_mutations_boolmatrix_fordepmap_hotspot')
    ],
    FUSIONS_TAIGA_ID:[
        ('CCLE_fusions_unfiltered', 'fusions_latest'),
        ('CCLE_fusions', 'filteredfusions_latest')
    ],
    EXPRESSION_TAIGA_ID:[
        ('CCLE_expression_full', 'genes_tpm_logp1'),
        ('CCLE_RNAseq_transcripts', 'transcripts_tpm_logp1'),
        ('CCLE_RNAseq_reads', 'genes_expected_count'),
        ('CCLE_expression', 'proteincoding_genes_tpm_logp1'),
        ('CCLE_expression_proteincoding_genes_expected_count', 'proteincoding_genes_expected_count'),
        ('CCLE_expression_transcripts_expected_count', 'transcripts_expected_count'),
        # ('CCLE_ssGSEA', 'gene_sets_all')
    ],
    CN_TAIGA_ID:[
        ('CCLE_gene_cn', 'achilles_gene_cn'),
        ('CCLE_segment_cn', 'achilles_segment')
        # ('CCLE_gene_cn', 'merged_genecn_all'),
        # ('CCLE_segment_cn', 'merged_segments_all')
    ]
}
