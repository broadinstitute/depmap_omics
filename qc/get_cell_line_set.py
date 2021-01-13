from taigapy import TaigaClient
tc = TaigaClient()

def tcget(name=None, version=None, file=None):
    # get arxspan ids from the Taiga file 
    # this function avoids overloading the memory because the entire variable is local to it
    tc_df = tc.get(name=name, version=version, file=file)
    
    if 'DepMap_ID' in tc_df.columns:
        arxspan_ids = set(tc_df['DepMap_ID'])
    else:
        arxspan_ids = set(tc_df.index)
    assert all([x.startswith('ACH-') for x in all_21Q1_gene_cn])    
    return arxspan_ids

def get_all_arxspans(name='internal-20q3-00d0', version=9):
    # get the set of all the arxspans from a specific Taiga upload
    CCLE_expression = tcget(name=name, version=version, file='CCLE_expression')
    CCLE_expression_full = tcget(name=name, version=version, file='CCLE_expression_full')
    CCLE_fusions = tcget(name=name, version=version, file='CCLE_fusions')
    CCLE_fusions_unfiltered = tcget(name=name, version=version, file='CCLE_fusions_unfiltered')
    CCLE_gene_cn = tcget(name=name, version=version, file='CCLE_gene_cn')
    CCLE_mutations = tcget(name=name, version=version, file='CCLE_mutations')
    CCLE_RNAseq_reads = tcget(name=name, version=version, file='CCLE_RNAseq_reads')
    CCLE_RNAseq_transcripts = tcget(name=name, version=version, file='CCLE_RNAseq_transcripts')
    CCLE_segment_cn = tcget(name=name, version=version, file='CCLE_segment_cn')
    
    arxspan_ids = (CCLE_expression | CCLE_expression_full | CCLE_fusions | 
                   CCLE_fusions_unfiltered | CCLE_gene_cn | CCLE_mutations |
                   CCLE_RNAseq_reads | CCLE_RNAseq_transcripts | CCLE_segment_cn)
    return arxspan_ids

