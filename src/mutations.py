import pandas as pd
from genepy.google import gcp

def download_maf_from_workspace(refwm, sample_set_ids = ['all_ice', 'all_agilent'],
                       output_maf='/tmp/mutation_filtered_terra_merged.txt'):
    sample_sets = refwm.get_sample_sets()
    dfs = []
    for sample_set_id in sample_sets.index.intersection(sample_set_ids):
        gcp.cpFiles(sample_sets.loc[sample_set_id, 'filtered_CGA_MAF_aggregated'], 
                  '/tmp/tmp.txt', payer_project_id='broad-firecloud-ccle', verbose=False);
        df = pd.read_csv('/tmp/tmp.txt', sep='\t', low_memory=False)
        dfs.append(df)
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(output_maf, index=False, sep='\t')
    return dfs_concat


def removeDuplicates(a, loc, prepended=['dm', 'ibm', 'ccle']):
  """
  This function is used to subset a df to only the columns with the most up to date names

  We consider a naming convention preprended_NAME_version and transform it into NAME with latest NAMES

  Args:
  ----
    a: the dataframe where loc contain the names
    loc: the location of names
    prepended: the set of possible prepended values

  Returns:
  -------
    a: the subsetted dataframe
  """
  values = []
  if len(prepended) > 0:
    for i in a[loc]:
      i = i.split('_')
      if i[0] in prepended:
        values.append('_'.join(i[1:]))
      else:
        values.append('_'.join(i))
      if i[-1] == '2':
        print(i)
    a[loc] = values
  a = a.sort_values(by=[loc])
  todrop = []
  for i in range(len(a[loc]) - 1):
    e = a[loc][i + 1].split('_')
    if len(e[-1]) == 1:
      if int(e[-1]) > 1 and e[0] == a[loc][i].split('_')[0]:
        todrop.append(a[loc][i])
        print(a[loc][i])
        print(e)
  a = a.set_index(loc).drop(todrop).reset_index()
  return a


def annotateLikelyImmortalized(maf, sample_col="DepMap_ID", 
genome_change_col="Genome_Change", TCGAlocs=['TCGAhsCnt', 'COSMIChsCnt'], 
max_recurrence=0.05, min_tcga_true_cancer=5):
  maf['is_likely_immortalization'] = False
  leng = len(set(maf[sample_col]))
  tocheck = []
  for k, v in Counter(maf[genome_change_col].tolist()).items():
    if v > max_recurrence*leng:
      tocheck.append(k)
  for val in list(set(tocheck)-set([np.nan])):
    if np.nan_to_num(maf[maf[genome_change_col] == val][TCGAlocs], 0).max() < min_tcga_true_cancer:
      maf.loc[maf[maf[genome_change_col]
                              == val].index, 'is_likely_immortalization'] = True
  return maf


def addAnnotation(maf, NCBI_Build='37', Strand="+"):
  """
  adds NCBI_Build and Strand annotaation on the whole maf file
  """
  maf['NCBI_Build'] = NCBI_Build
  maf['Strand'] = Strand
  return maf
