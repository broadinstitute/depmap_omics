VIRTUAL_RELEASE = {'name': 'tentative-virtual-d84e', 'version': 8}
REFERENCE_RELEASE = {'name': 'internal-21q1-4fc4', 'version': 39}

FUSIONS_MERGE_COLS = ['DepMap_ID', 'LeftBreakpoint', 'RightBreakpoint']
SEGMENT_CN_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start', 'End']
MUTATIONS_MERGE_COLS = ['DepMap_ID', 'Chromosome', 'Start_position', 'End_position', 'Tumor_Seq_Allele1']
FILE_ATTRIBUTES = [
    {'file': 'CCLE_expression', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_proteincoding_genes_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_transcripts', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_transcripts_expected_count', 'ismatrix': True, 'hasNA': False, 'gene_id': 'enst', 'omicssource':'RNA'},
    {'file': 'CCLE_expression_full', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_RNAseq_reads', 'ismatrix': True, 'hasNA': False, 'gene_id': 'ensg', 'omicssource':'RNA'},
    {'file': 'CCLE_fusions', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS},
    {'file': 'CCLE_fusions_unfiltered', 'ismatrix': False, 'omicssource':'RNA', 'merge_cols': FUSIONS_MERGE_COLS},
    {'file': 'CCLE_ssGSEA', 'ismatrix': True, 'hasNA': False,'omicssource':'RNA', 'gene_id': None},
    {'file': 'CCLE_gene_cn', 'ismatrix': True, 'hasNA': True, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_segment_cn', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': SEGMENT_CN_MERGE_COLS},
    {'file': 'CCLE_mutations', 'ismatrix': False, 'omicssource':'DNA', 'merge_cols': MUTATIONS_MERGE_COLS},
    {'file': 'CCLE_mutations_bool_damaging', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_hotspot', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_othercons', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'},
    {'file': 'CCLE_mutations_bool_othernoncons', 'ismatrix': True, 'hasNA': False, 'gene_id': 'entrez', 'omicssource':'DNA'}
]


# the following information is used to create a tentative virtual
MUTATIONS_TAIGA_ID = 'mutations-latest-ed72'
FUSIONS_TAIGA_ID = 'fusions-95c9'
EXPRESSION_TAIGA_ID = 'expression-d035'
CN_TAIGA_ID = 'cn-achilles-version-06ca'

TAIGA_IDS_LATEST = {
    MUTATIONS_TAIGA_ID:[
        ('CCLE_mutations', 'all_somatic_mutations_all_21Q2_depmapversion'),
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
        ('CCLE_expression_full', 'expression_21Q2_genes_tpm_logp1'),
        ('CCLE_RNAseq_transcripts', 'expression_21Q2_transcripts_tpm_logp1'),
        ('CCLE_RNAseq_reads', 'expression_21Q2_genes_expected_count'),
        ('CCLE_expression', 'expression_21Q2_proteincoding_genes_tpm_logp1'),
        ('CCLE_expression_proteincoding_genes_expected_count', 'expression_21Q2_proteincoding_genes_expected_count'),
        ('CCLE_expression_transcripts_expected_count', 'expression_21Q2_transcripts_expected_count'),
        ('CCLE_ssGSEA', 'gene_sets_21Q2_all')
    ],
    CN_TAIGA_ID:[
        ('CCLE_gene_cn', 'all_21Q2_gene_cn'),
        ('CCLE_segment_cn', 'all_21Q2_segment')
    ]
}
