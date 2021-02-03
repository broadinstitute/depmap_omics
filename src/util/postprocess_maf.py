import pandas as pd
from util.util import gsutil_cp

def download_maf_from_workspace(refwm, sample_set_ids = ['all_ice', 'all_agilent'],
                       output_maf='/tmp/mutation_filtered_terra_merged.txt'):
    sample_sets = refwm.get_sample_sets()
    dfs = []
    for sample_set_id in sample_sets.index.intersection(sample_set_ids):
        gsutil_cp(sample_sets.loc[sample_set_id, 'filtered_CGA_MAF_aggregated'], 
                  '/tmp/tmp.txt', payer_project_id='broad-firecloud-ccle', verbose=False);
        df = pd.read_csv('/tmp/tmp.txt', sep='\t', low_memory=False)
        dfs.append(df)
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(output_maf, index=False, sep='\t')
    return dfs_concat