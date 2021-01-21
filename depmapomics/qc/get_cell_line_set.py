from taigapy import TaigaClient
tc = TaigaClient()

DEFAULT_FILENAMES=['CCLE_expression', 'CCLE_expression_full', 'CCLE_fusions',
                   'CCLE_fusions_unfiltered', 'CCLE_gene_cn', 'CCLE_mutations',
                   'CCLE_RNAseq_reads', 'CCLE_RNAseq_transcripts', 'CCLE_segment_cn']

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

# def get_all_arxspans_pooled(name=None, version=None, verbose=False, files=DEFAULT_FILENAMES):
#     # get the set of all the arxspans from a specific Taiga upload
#     arxspan_ids = set()
#     for file in files:
#         if verbose:
#             print('\tfetching {}'.format(file))
#         arxspan_ids = arxspan_ids | tcget(name=name, version=version, file=file)
        
#     return arxspan_ids

def applyfunc_to_json(json_dict, func,  verbose=False):
    arxspans = {}
    if type(json_dict)==dict:
        output_dict = {}
        for k, v in json_dict.items():
            if verbose:
                print(k)
            output_dict[k] = applyfunc_to_json(v, func, verbose=verbose) 
        return output_dict
    else:
        return func(json_dict)

# def rename_dict(json_dict, verbose=False):
#     arxspans = {}
#     if type(json_dict)==dict:
#         output_dict = {}
#         for k, v in json_dict.items():
#             if verbose:
#                 print(k)
#             if k == 'CCLE_unfiltered_fusions':
#                 output_dict['CCLE_fusions_unfiltered'] = rename_dict(v, verbose=verbose) 
#             else:
#                 output_dict[k] = rename_dict(v, verbose=verbose) 
#         return output_dict
#     else:
#         return json_dict    

def get_release_diffs(arxspan_dict, lines_to_release, quarters = ['20q3', '21q1']):
    release_diff = {}
    release_diff_reverse = {}
    for portal, portal_dict in arxspan_dict[quarters[0]].items():
        release_diff[portal] = {}
        release_diff_reverse[portal] = {}
        for file in portal_dict.keys():
            expected_lines = set(lines_to_release[portal].dropna()) | arxspan_dict[quarters[0]][portal][file]
            
            release_diff[portal][file] = arxspan_dict[quarters[1]][portal][file] - expected_lines
            release_diff_reverse[portal][file] = expected_lines - arxspan_dict[quarters[1]][portal][file]
    return release_diff, release_diff_reverse  

def pool_arxspans_per_portal(arxspan_dict):
    arxspan_dict_per_portal = {}
    for quarter, qurter_dict in arxspan_dict.items():
        portal_arxspans_dict = {}
        for portal, portal_dict in qurter_dict.items():
            portal_arxspans = set()
            for file in portal_dict.keys():
                portal_arxspans = portal_arxspans | arxspan_dict[quarter][portal][file]
            portal_arxspans_dict[portal] = portal_arxspans
        arxspan_dict_per_portal[quarter] = portal_arxspans_dict
    return arxspan_dict_per_portal


def get_release_diff_pooled(arxspan_dict, lines_to_release, quarters):
    arxspans = pool_arxspans_per_portal(arxspan_dict)
    arxspan_diff = {}
    arxspan_revdiff = {}
    for portal in arxspans[quarters[0]].keys():
        arxspans_expected = (arxspans[quarters[0]][portal] | set(lines_to_release[portal].dropna()))
        arxspan_diff[portal] = arxspans[quarters[1]][portal] - arxspans_expected
        arxspan_revdiff[portal] = arxspans_expected - arxspans[quarters[1]][portal]
    return arxspan_diff, arxspan_revdiff
    
def get_all_arxspans(taiga_dict_expanded, verbose=False):
    arxspan_dict = applyfunc_to_json(taiga_dict_expanded, lambda x: tcget(name=x[0], version=x[1], file=x[2]), verbose=True)
    return arxspan_dict

def propagate_taiga_dict_with_filenames(taiga_dict, files=DEFAULT_FILENAMES):
    return applyfunc_to_json(taiga_dict, lambda x:{file:x+[file] for file in files})

def pretty_print(release_diffs):
    def _convert_to_text(x):
        if x == set():
            return 'None'
        else:
            return ', '.join(x)
    text = ''
    release_diffs_text = applyfunc_to_json(release_diffs, lambda x: _convert_to_text(x))
    for portal, portal_dict in release_diffs_text.items():
        text += portal + ':\n'
        if type(portal_dict) == dict:
            for file in portal_dict.keys():
                text+= '\t{}: {}\n'.format(file, release_diffs_text[portal][file])
            text+= '\n'
        else:
            text+= '\t{}\n'.format(release_diffs_text[portal])
    return(text)


def pretty_print_diff(arxspan_dict, lines_to_release, quarters = ['20q3', '21q1'], savefile=False):
    release_diffs, release_diffs_reverse = get_release_diffs(arxspan_dict, lines_to_release, quarters = quarters)
    release_diffs_pooled, release_diffs_reverse_pooled = get_release_diff_pooled(arxspan_dict, lines_to_release, quarters)
    text= 'lines added ({} compared to {})\n'.format(*quarters[::-1])
    text+= pretty_print(release_diffs)
    text += '\n\tlines pooled across portals\n'
    text += pretty_print(release_diffs_pooled)
    
    text += '\n'+'_'*20 + '\n'

    text += 'lines dropped ({} compared to {})\n'.format(*quarters[::-1])
    text += pretty_print(release_diffs_reverse)
    text += '\n\tlines pooled across portals\n'
    text += pretty_print(release_diffs_reverse_pooled)
    if savefile:
        with open('diff_{}_vs_{}.txt'.format(*quarters[::-1]), 'w') as f:
            f.write(text)
    
    return text

# def get_arxspans_across_portals_pooled(taiga_dict, files=DEFAULT_FILENAMES, verbose=False):
#     arxpans = {}
#     for key, val in taiga_dict.items():
#         print('fetching {}/{}'.format(*val))
#         arxpans[key] = get_all_arxspans_pooled(name=val[0], version=val[1], verbose=verbose, files=files)
#     return arxpans

# def get_release_diff_pooled(arxspans, lines_to_release, quarters, portals = ['public', 'internal', 'dmc', 'ibm']):
#     arxspan_diff = {}
#     arxspan_revdiff = {}
#     for portal in portals:
#         arxspans_expected = (arxspans[quarters[0]][portal] | set(lines_to_release[portal].dropna()))
#         arxspan_diff[portal] = arxspans[quarters[1]][portal] - arxspans_expected
#         arxspan_revdiff[portal] = arxspans_expected - arxspans[quarters[1]][portal]
        
#     def printx(x):
#         def replace_empty(v):
#             if v == set():
#                 return {'None'}
#             else:
#                 return v
        
#         if type(x)==dict:
#             text='\n'.join(['\t{}: {}'.format(k, ', '.join(replace_empty(v))) for k,v in x.items()])
#         elif type(x)==set:
#             text = ', '.join(replace_empty(x))
#         print(text)
#     print('extra arxspans per portal...')
#     printx(arxspan_diff)
#     print('\nmissing arxspans per portal...')
#     printx(arxspan_revdiff)
#     print('\nextra arxspans across portals...')
#     printx(set().union(*arxspan_diff.values()))
#     print('\nmissing arxspans across portals...')
#     printx(set().union(*arxspan_revdiff.values()))