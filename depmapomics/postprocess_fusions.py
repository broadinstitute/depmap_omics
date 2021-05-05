import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
import re


def addToMainFusion(input_filenames, main_filename, DepMap_ID="DepMap_ID"):
  """
  Given a tsv fusion files from RSEM algorithm, merge it to a tsv set of fusion data

  Args:
  ----
    input_filenames: a set of filepath to input the files should be c|tsv from Terra fusion pipeline
    main_filename: a filepath to input the files should be c|tsv from Terra aggregation pipeline
  """
  maindata = pd.read_csv(main_filename, sep='\t')
  if '.' in maindata[DepMap_ID][0]:
    maindata[DepMap_ID] = [i[0]
                           for i in maindata[DepMap_ID].str.split('.').tolist()]
  samples = set(maindata[DepMap_ID].tolist())
  with open(main_filename, 'a') as f:
    for input_filename in input_filenames:
      df = pd.read_csv(input_filename, sep='\t')
      input_filename = input_filename.split('/')[-1].split('.')[0]
      if input_filename in samples:
        print(input_filename + " is Already in main fusions")
      df[DepMap_ID] = pd.Series(
        [input_filename] * len(df.index.tolist()), index=df.index)
      cols = df.columns.tolist()
      cols = cols[-1:] + cols[: -1]
      df = df[cols]
      df.to_csv(f, header=False, sep='\t', index=False)


def filterFusions(fusions, maxfreq=0.1, minffpm=0.05, CCLE_count="CCLE_count", red_herring=['GTEx_recurrent', 'DGD_PARALOGS', 'HGNC_GENEFAM', 'Greger_Normal', 'Babiceanu_Normal', 'ConjoinG', 'NEIGHBORS']):
  """
  Given a fusion file from star fusion, filters it (will also filter Mitochrondria and HLA genes)

  Args:
  -----
    fusions: dataframe the fusion as a dataframe should contain: LeftBreakpoint, RightBreakpoint, FusionName, annots, SpliceType, LargeAnchorSupport, FFPM columns
    maxfreq: int the max allowed frequency of that fusion across our samples
    DepMap_ID: str: colname for the sample ids
    CCLE_count: str: colname where are stored counts of that fusion name across our samples
    minffpm: int minimum ffpm freq to filter on
    red_herring: list[str] of flags to filter on
  """
  fusions = fusions.copy()
  # remove recurrent
  fusions = fusions[fusions[CCLE_count] <
                    len(set(fusions["DepMap_ID"]))*maxfreq]
  # (1) Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes,
  fusions = fusions[~(fusions['LeftBreakpoint'].str.contains(
      'chrM') & fusions['RightBreakpoint'].str.contains('chrM'))]
  fusions = fusions[~fusions['FusionName'].str.contains('^HLA\\-')]
  # (2) Remove red herring fusions
  fusions = fusions[~fusions['annots'].str.contains(
      '|'.join(red_herring), case=False)]
  # (4) Removed fusion with (SpliceType=" INCL_NON_REF_SPLICE" and
  # LargeAnchorSupport="No" and minFAF<0.02), or
  fusions = fusions[~((fusions['SpliceType'] == "INCL_NON_REF_SPLICE") & (
      fusions['LargeAnchorSupport'] == "NO_LDAS") & (fusions['FFPM'] < 0.1))]
  # STAR-Fusion suggests using 0.1, but after looking at the
  # translocation data, this looks like it might be too aggressive
  fusions = fusions[fusions['FFPM'] > minffpm]
  return fusions


def renameFusionGene(a):
  """
  Given a fusion name from star-fusion, renames it
  """
  return [str(i.split('^')).replace(', ', ' (').replace("'", "")[1:-1]+')' for i in a]


def standardize_leftgene_rightgene(fusions):
    """
    converts [GENE_NAME]^[ENSG] --> [GENE_NAME] ([ENSG])
    Example: "SMAD4^ENSG00000141646.14" --> "SMAD4 (ENSG00000141646.14)"
    """
    fusions[['LeftGene', 'RightGene']] = fusions[['LeftGene', 'RightGene']]\
        .applymap(lambda x: '{} ({})'.format(*x.split(r'^')))
    return fusions


def postprocess_fusions(refworkspace, sample_id_name='DepMap_ID', sample_set_name = 'all',
                        output_fusion_file='/tmp/fusions.csv',
                        output_fusion_filtered_file='/tmp/filtered_fusions.csv'):
    refwm = dm.WorkspaceManager(refworkspace).disable_hound()
    aggregated = refwm.get_sample_sets().loc[sample_set_name]['fusions_star']
    cpFiles([aggregated], "/tmp/expression.fusion.tsv")
    fusions = pd.read_csv("/tmp/expression.fusion.tsv",
                          names=[sample_id_name, 'FusionName', 'JunctionReadCount',
                                 'SpanningFragCount', 'SpliceType', 'LeftGene', 'LeftBreakpoint',
                                 'RightGene', 'RightBreakpoint', 'LargeAnchorSupport', 'FFPM',
                                 'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc',
                                 'RightBreakEntropy', 'annots'], skiprows=1, sep='\t')

    CCLE_count = fusions[['LeftBreakpoint', 'RightBreakpoint']]\
        .value_counts()\
        .to_frame(name='CCLE_count')
    fusions = pd.merge(fusions, CCLE_count, on=['LeftBreakpoint', 'RightBreakpoint'])

    fusions_filtered = filterFusions(fusions, maxfreq=0.1, sample_id_name=sample_id_name)

    fusions = standardize_leftgene_rightgene(fusions)
    fusions_filtered = standardize_leftgene_rightgene(fusions_filtered)

    fusions.to_csv(output_fusion_file, index=False)
    fusions_filtered.to_csv(output_fusion_filtered_file, index=False)

    return fusions, fusions_filtered
