TEMP_VIRTUAL_TAIGA_ID = 'temp-virtual-2ec2'
TEMP_VIRTUAL_TAIGA_VERSION = '2'

MUTATIONS_TAIGA_ID = 'mutations-latest-ed72'
FUSIONS_TAIGA_ID = 'fusions-95c9'
EXPRESSION_TAIGA_ID = 'expression-d035'
CN_TAIGA_ID = 'cn-achilles-version-06ca'

TAIGA_IDS_LATEST = {
    MUTATIONS_TAIGA_ID:[
        ('CCLE_mutations', 'wes_somatic_mutations_all_21Q1'),
        ('CCLE_mutations_bool_damaging', 'all_somatic_mutations_boolmatrix_fordepmap_damaging'),
        ('CCLE_mutations_bool_othernoncons', 'all_somatic_mutations_boolmatrix_fordepmap_othernoncons'),
        ('CCLE_mutations_bool_othercons', 'all_somatic_mutations_boolmatrix_fordepmap_othercons'),
        ('CCLE_mutations_bool_hotspot', 'all_somatic_mutations_boolmatrix_fordepmap_hotspot')
    ],
    FUSIONS_TAIGA_ID:[
        ('CCLE_fusions_unfiltered', 'fusions_21Q2'),
        ('CCLE_fusions', 'filtered_fusions_21Q2')
    ],
    EXPRESSION_TAIGA_ID:[
        ('CCLE_expression_full', 'expression_21Q2_genes_tpm'),
        ('CCLE_RNAseq_transcripts', 'expression_21Q2_transcripts_tpm'),
        ('CCLE_RNAseq_reads', 'expression_21Q2_genes_expected_count'),
        ('CCLE_expression', 'expression_21Q2_proteincoding_genes_tpm'),
        ('CCLE_expression_proteincoding_genes_expected_count', 'expression_21Q2_proteincoding_genes_expected_count'),
        ('CCLE_expression_transcripts_expected_count', 'expression_21Q2_transcripts_expected_count')
    ],
    CN_TAIGA_ID:[
        ('CCLE_gene_cn', 'all_21Q2_gene_cn'),
        ('CCLE_segment_cn', 'all_21Q2_segment')
    ]
}

FILE_ATTRIBUTES = [
    {'file': 'CCLE_expression', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_expression_proteincoding_genes_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_mutations_bool_damaging', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_mutations_bool_hotspot', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_mutations_bool_othercons', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_mutations_bool_othernoncons', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez'},
    {'file': 'CCLE_RNAseq_transcripts', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst'},
    {'file': 'CCLE_expression_transcripts_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst'},
    {'file': 'CCLE_expression_full', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg'},
    {'file': 'CCLE_RNAseq_reads', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg'},
    {'file': 'CCLE_gene_cn', 'ismatrix': True, 'hasNA': True, 'gene_id': 'entrez'},
    {'file': 'CCLE_fusions', 'ismatrix': False},
    {'file': 'CCLE_fusions_unfiltered', 'ismatrix': False},
    {'file': 'CCLE_segment_cn', 'ismatrix': False},
    {'file': 'CCLE_mutations', 'ismatrix': False}
]
