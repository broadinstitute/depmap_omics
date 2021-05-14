import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
import numpy as np
from genepy.mutations import filterAllelicFraction, filterCoverage
from collections import Counter


def download_maf_from_workspace(refwm, sample_set_ids=['all_ice', 'all_agilent'],
                                output_maf='/tmp/mutation_filtered_terra_merged.txt'):
    sample_sets = refwm.get_sample_sets()
    dfs = []
    for sample_set_id in sample_sets.index.intersection(sample_set_ids):
        cpFiles(sample_sets.loc[sample_set_id, 'filtered_CGA_MAF_aggregated'],
                    '/tmp/tmp.txt', payer_project_id='broad-firecloud-ccle', verbose=False)
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


def add_variant_annotation_column(maf):
    mutation_groups={
        "other conserving": ["5'Flank", "Intron", "IGR", "3'UTR", "5'UTR"],
        "other non-conserving":["In_Frame_Del", "In_Frame_Ins", "Stop_Codon_Del",
            "Stop_Codon_Ins", "Missense_Mutation", "Nonstop_Mutation"],
        'silent': ['Silent'],
        "damaging":['De_novo_Start_OutOfFrame','Frame_Shift_Del','Frame_Shift_Ins',
            'Splice_Site', 'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP','Nonsense_Mutation']
    }

    rename = {}
    for k,v in mutation_groups.items():
        for e in v:
            rename[e] = k
    maf['Variant_annotation'] = [rename[i] for i in maf['Variant_Classification'].tolist()]
    return maf

def postprocess_mutations_filtered_wes(refworkspace, sample_set_name = 'all',
                                       output_file='/tmp/wes_somatic_mutations.csv'):
    refwm = dm.WorkspaceManager(refworkspace).disable_hound()
    filtered = refwm.get_sample_sets().loc[sample_set_name]['filtered_CGA_MAF_aggregated']
    print('copying aggregated filtered mutation file')
    cpFiles([filtered], "/tmp/mutation_filtered_terra_merged.txt")
    print('reading the mutation file')
    mutations = pd.read_csv('/tmp/mutation_filtered_terra_merged.txt', sep='\t', low_memory=False)
    mutations = mutations.rename(columns={"i_ExAC_AF":"ExAC_AF",
                                          "Tumor_Sample_Barcode":'DepMap_ID',
                                          "Tumor_Seq_Allele2":"Tumor_Allele"}).\
    drop(columns=['Center','Tumor_Seq_Allele1'])
    # mutations = annotate_likely_immortalized(mutations, TCGAlocs = ['TCGAhsCnt', 'COSMIChsCnt'], max_recurrence=0.05, min_tcga_true_cancer=5)
    print('writing CGA_WES_AC column')
    mutations['CGA_WES_AC'] = [str(i[0]) + ':' + str(i[1]) for i in np.nan_to_num(mutations[['t_alt_count','t_ref_count']].values,0).astype(int)]
    # apply version:
    # mutations['CGA_WES_AC'] = mutations[['t_alt_count', 't_ref_count']].fillna(0).astype(int).apply(lambda x: '{:d}:{:d}'.format(*x), raw=True, axis=1)
    print('filtering coverage')
    mutations = filterCoverage(mutations, loc=['CGA_WES_AC'], sep=':', cov=2)
    print('filtering allelic fractions')
    mutations = filterAllelicFraction(mutations, loc=['CGA_WES_AC'], sep=':', frac=0.1)
    print('adding NCBI_Build and strand annotations')
    mutations = addAnnotation(mutations, NCBI_Build='37', Strand="+")
    print('adding the Variant_annotation column')
    mutations = add_variant_annotation_column(mutations)
    print('saving results to output file')
    mutations.to_csv('/tmp/wes_somatic_mutations.csv', index=False)
    return mutations
