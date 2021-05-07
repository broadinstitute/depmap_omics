import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
import re
import os.path
import seaborn as sns
from depmapomics.config import FUSIONreadme

from depmapomics import tracker

from taigapy import TaigaClient
tc = TaigaClient()


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


def filterFusions(fusions, maxfreq=0.1, sampleCol, minffpm=0.05, countCol="CCLE_count",
                  red_herring=['GTEx_recurrent', 'DGD_PARALOGS', 'HGNC_GENEFAM', 'Greger_Normal', 'Babiceanu_Normal', 'ConjoinG', 'NEIGHBORS']):
  """
  Given a fusion file from star fusion, filters it (will also filter Mitochrondria and HLA genes)

  Args:
  -----
    fusions: dataframe the fusion as a dataframe should contain: LeftBreakpoint, RightBreakpoint, FusionName, annots, SpliceType, LargeAnchorSupport, FFPM columns
    maxfreq: int the max allowed frequency of that fusion across our samples
    DepMap_ID: str: colname for the sample ids
    countCol: str: colname where are stored counts of that fusion name across our samples
    minffpm: int minimum ffpm freq to filter on
    red_herring: list[str] of flags to filter on
  """
  fusions = fusions.copy()
  # remove recurrent
  fusions = fusions[fusions[countCol] <
                    len(set(fusions[sampleCol]))*maxfreq]
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
  Given a list of fusion names from star-fusion, renames and returns them
  """
  return [str(i.split('^')).replace(', ', ' (').replace("'", "")[1:-1]+')' for i in a]


def standardizeGeneNames(fusions):
    """
    TODO: todocument
    converts [GENE_NAME]^[ENSG] --> [GENE_NAME] ([ENSG])
    Example: "SMAD4^ENSG00000141646.14" --> "SMAD4 (ENSG00000141646.14)"
    """
    fusions[['LeftGene', 'RightGene']] = fusions[['LeftGene', 'RightGene']]\
        .applymap(lambda x: '{} ({})'.format(*x.split(r'^')))
    return fusions


def postProcessFusions(refworkspace, samplesetname, sampleCol='DepMap_ID', samplesetToLoad = 'all',
                        colnames=['FusionName', 'JunctionReadCount',
                        'SpanningFragCount', 'SpliceType', 'LeftGene', 'LeftBreakpoint',
                        'RightGene', 'RightBreakpoint', 'LargeAnchorSupport', 'FFPM',
                        'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc',
                        'RightBreakEntropy', 'annots'], doplot=True,
                        countCol="CCLE_count", folder="temp", rnFunc=None,
                        **kwargs ):
    """
    TODO: todocument
    """
    refwm = dm.WorkspaceManager(refworkspace)
    print("loading fusions")
    aggregated = refwm.get_sample_sets().loc[samplesetToLoad]['fusions_star']
    fusions = pd.read_csv(aggregated,
                          names=[sampleCol]+colnames, skiprows=1, sep='\t')

    print("postprocessing fusions")
    fusions.RightGene = renameFusionGene(fusions.RightGene)
    fusions.LeftGene = renameFusionGene(fusions.LeftGene)
    fusions = standardizeGeneNames(fusions)
    
    count = fusions[['LeftBreakpoint', 'RightBreakpoint']]\
        .value_counts()\
        .to_frame(name=countCol)
    fusions = pd.merge(fusions, count, on=[
                       'LeftBreakpoint', 'RightBreakpoint'])

    fusions_filtered = filterFusions(
        fusions, sampleCol=sampleCol, countCol=countCol, **kwargs)
    if doplot:
        sns.kdeplot(fusions[countCol])
    
    print("saving")
    fusions.to_csv(os.path.join(folder,'/fusions_unRenamed_'+samplesetname+'.csv'), index=False)
    if rnFunc:
        print('renaming')
        renaming = rnFunc(set(fusions[sampleCol]))
        fusions = fusions[fusions[sampleCol].isin(renaming.keys())].replace(
            {sampleCol: renaming}).reset_index(drop=True)
        fusions.to_csv(os.path.join(folder, 'fusions_' +
                                    samplesetname+'.csv'), index=False)
    
    fusions_filtered.to_csv(os.path.join(
        folder, 'filtered_fusions_'+samplesetname+'.csv'), index=False)

    return fusions, fusions_filtered


def CCLEPostProcessFusions(refworkspace, samplesetname, fusionSamplecol="depmap_id", 
                           taiga_dataset="fusions-95c9", dataset_description=FUSIONreadme,
                           **kwargs):
    """
    TODO: todocument
    """
    refwm = dm.WorkspaceManager(refworkspace)
    fusions, _ = postProcessFusions(refworkspace=refworkspace, samplesetname=samplesetname,
        rnFunc=lambda x: tracker.removeOlderVersions(
            names=x, refsamples=refwm.get_samples(),
            arxspan_id="arxspan_id", version="version",
            ),
        **kwargs)
    
    print('')
    prev = tc.get(name='depmap-a0ab', file='CCLE_fusions_unfiltered')
    print('new')
    print(set(fusions[fusionSamplecol]) - set(prev[fusionSamplecol]))

    print('removed')
    print(set(prev[fusionSamplecol]) - set(fusions[fusionSamplecol]))

    print("changes in fusion names")
    pf = prev.copy()
    pf["id"] = pf[fusionSamplecol]+"_"+pf["FusionName"]
    f = fusions.copy()
    f["id"] = f[fusionSamplecol]+"_"+f["FusionName"]
    print(len(set(pf[~pf.id.isin(f.id.tolist())][fusionSamplecol])))

    print("changes in junction readd counts")
    f["sid"] = f[fusionSamplecol]+"_"+f["FusionName"] + \
        "_" + f["JunctionReadCount"].astype(str)
    pf["sid"] = pf[fusionSamplecol]+"_"+pf["FusionName"] + \
        "_" + pf["JunctionReadCount"].astype(str)
    print(len(set(pf[~pf.sid.isin(f.sid.tolist())][fusionSamplecol])))
    tc.update_dataset(dataset_permaname=taiga_dataset,
                      changes_description="new "+samplesetname+" release!",
                      upload_files=[
                          {
                              "path": 'temp/fusions_'+samplesetname+'.csv',
                              "format": "TableCSV",
                              "encoding": "utf-8"
                          },
                          {
                              "path": 'temp/filtered_fusions_'+samplesetname+'.csv',
                              "format": "TableCSV",
                              "encoding": "utf-8"
                          },
                          {
                              "path": "temp/fusions_withReplicates_"+samplesetname+".csv",
                              "format": "TableCSV",
                              "encoding": "utf-8"
                          },
                      ],
                      dataset_description=dataset_description)
