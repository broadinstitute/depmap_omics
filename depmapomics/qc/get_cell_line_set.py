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
    assert all([x.startswith('ACH-') for x in arxspan_ids])    
    return arxspan_ids

def get_all_arxspans(name=None, version=None, verbose=False,
        files=['CCLE_expression', 'CCLE_expression_full', 'CCLE_fusions',
               'CCLE_fusions_unfiltered', 'CCLE_gene_cn', 'CCLE_mutations',
               'CCLE_RNAseq_reads', 'CCLE_RNAseq_transcripts', 'CCLE_segment_cn']):
    # get the set of all the arxspans from a specific Taiga upload
    arxspan_ids = set()
    for file in files:
        if verbose:
            print('Fetching {}'.format(file))
        arxspan_ids = arxspan_ids | tcget(name=name, version=version, file=file)
        
    return arxspan_ids
