import os.path

import dalmatian as dm
import pandas as pd
import numpy as np
from scipy.stats import zscore

from genepy.utils import helper as h
from genepy import rna, terra

from depmapomics import terra as myterra
from depmapomics.qc import rna as myQC
from depmapomics import tracker as track

from depmapomics.config import *
from gsheets import Sheets


def addSamplesRSEMToMain(input_filenames, main_filename):
  """
  given a tsv RNA files from RSEM algorithm, merge it to a tsv set of RNA data

  Args:
  ----
    input_filenames: a list of dict like file path in Terra gs://, outputs from the rsem pipeline
    main_filename: a dict like file paths in Terra gs://, outputs from rsem aggregate
  """
  genes_count = pd.read_csv('temp/' + main_filename['rsem_genes_expected_count'].split('/')[-1],
                            sep='\t', compression='gzip')
  transcripts_tpm = pd.read_csv('temp/' + main_filename['rsem_transcripts_tpm'].split('/')[-1],
                                sep='\t', compression='gzip')
  genes_tpm = pd.read_csv('temp/' + main_filename['rsem_genes_tpm'].split('/')[-1],
                          sep='\t', compression='gzip')

  for input_filename in input_filenames:
    name = input_filename['rsem_genes'].split(
        '/')[-1].split('.')[0].split('_')[-1]
    rsem_genes = pd.read_csv(
        'temp/' + input_filename['rsem_genes'].split('/')[-1], sep='\t')
    rsem_transcripts = pd.read_csv(
        'temp/' + input_filename['rsem_isoforms'].split('/')[-1], sep='\t')
    genes_count[name] = pd.Series(
        rsem_genes['expected_count'], index=rsem_genes.index)
    transcripts_tpm[name] = pd.Series(
        rsem_transcripts['TPM'], index=rsem_transcripts.index)
    genes_tpm[name] = pd.Series(rsem_genes['TPM'], index=rsem_genes.index)

  genes_count.to_csv('temp/' + main_filename['rsem_genes_expected_count'].split('/')[-1], sep='\t',
                     index=False, index_label=False, compression='gzip')
  transcripts_tpm.to_csv('temp/' + main_filename['rsem_transcripts_tpm'].split('/')[-1], sep='\t',
                         index=False, index_label=False, compression='gzip')
  genes_tpm.to_csv('temp/' + main_filename['rsem_genes_tpm'].split('/')[-1], sep='\t',
                   index=False, index_label=False, compression='gzip')


def solveQC(tracker, failed, save="", newname="arxspan_id"):
  """create a renaming dict to rename the columns of the QC file

  based on which samples have failed QC and which are blacklisted int he sample tracker

  Args:
    tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
    failed (list): list of samples that failed QC
    save (str): save the renaming dict to a file
    newname (str): name of the column in the tracker to rename to
  Returns:
    dict: a dict to rename the samples to
  """
  newfail = []
  rename = {}
  # finding other replicates to solve failed ones
  for val in failed:
    if val not in tracker:
      continue
    a = tracker.loc[val][newname]
    res = tracker[(tracker.datatype == 'rna')
                  & (tracker[newname] == a)]
    if len(res) > 1:
      for k in res.index:
        if k not in failed:
          rename[val] = k
    else:
      newfail.append(val)
  print("samples that failed:")
  print(newfail)
  if save:
    h.listToFile(newfail, save+"_rnafailed.txt")
  return rename


def updateTracker(selected, failed, lowqual, tracker, samplesetname, refworkspace=RNAWORKSPACE,
                  sheetname=SHEETNAME, sheetcreds=SHEETCREDS,
                  onlycol=STARBAMCOLTERRA, newgs=RNA_GCS_PATH_HG38,
                  dry_run=False, qcname="star_logs", match=".Log.final.out", todrop=[]):
  """updates the sample tracker with the new samples and the QC metrics

  Args:
    tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
    selected (list[str]): which samples were selected in the release of the analysis
    samplesetname (str): the name of the sample set or of the current analysis
    samplesinset (list[str]): list of samples in the analysis.
    lowqual (list[str]): list of samples that failed QC
    newgs (str, optional): google storage path where to move the files. Defaults to ''.
    sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
    sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
    qcname (str, optional): Terra column containing QC files. Defaults to "star_logs".
    refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
    onlycol (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to STARBAMCOLTERRA.
    todrop (list, optional): list of samples to be dropped. Defaults to []
  """
  refwm = dm.WorkspaceManager(refworkspace)
  samplesinset = [i['entityName'] for i in refwm.get_entities(
      'sample_set').loc[samplesetname].samples]
  starlogs = myterra.getQC(workspace=refworkspace, only=samplesinset,
                           qcname=qcname, match=match)
  for k, v in starlogs.items():
    if k == 'nan':
      continue
    a = tracker.loc[k, 'processing_qc']
    a = '' if a is np.nan else a
    tracker.loc[k, 'processing_qc'] = str(v) + ',' + a
    if tracker.loc[k, 'bam_qc'] != v[0]:
      tracker.loc[k, 'bam_qc'] = v[0]
  tracker.loc[tracker[tracker.datatype.isin(['rna'])].index, samplesetname]=0
  track.update(tracker, selected, samplesetname, failed, lowqual, newgs,
                     sheetcreds, sheetname, refworkspace, onlycol,  dry_run, todrop=todrop)


def loadFromRSEMaggregate(refworkspace, todrop=[], filenames=RSEMFILENAME,
                          sampleset="all", renamingFunc=None):
  """Load the rsem aggregated files from Terra

  Args:
    refworkspace (str): the workspace where to load the files from
    todrop (list[str], optional): list of samples to drop. Defaults to [].
    filenames (list[str], optional): the filenames to load. Defaults to RSEMFILENAME.
    sampleset (str, optional): the sample set to load. Defaults to 'all'.
    renamingFunc (function, optional): the function to rename the samples 
      (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
  
  Returns:
    dict(str: pd.df): the loaded dataframes
    dict: the renaming dict used to rename the dfs columns

  """
  files = {}
  renaming = {}
  refwm = dm.WorkspaceManager(refworkspace)
  res = refwm.get_sample_sets().loc[sampleset]
  for val in filenames:
    file = pd.read_csv(res[val], compression='gzip', header=0,
                       sep='\t', quotechar='"', error_bad_lines=False)
    if renamingFunc is not None:
      # removing failed version
      renaming = renamingFunc(file.columns[2:], todrop)
    else:
      renaming.update({i: i for i in file.columns[2:] if i not in todrop})
    renaming.update({'transcript_id(s)': 'transcript_id'})
    # we remove the failed samples where we did not found anything else to replace them with
    files[val] = file[file.columns[:2].tolist()+[i for i in file.columns[2:]
                                                 if i in set(renaming.keys())]].rename(columns=renaming)
  return files, renaming


def subsetGenes(files, gene_rename, filenames=RSEM_TRANSCRIPTS,
                drop=[], index="transcript_id"):
  """
  Subset the rsem transcripts file to keep only the genes of interest

  Args:
    files (dict(str: pd.dfs)): the rsem transcripts dfs to subset samples x genes
    gene_rename (dict): the gene renaming dict (here we expect a dict of ensembl transcript ids: gene names)
    filenames (list[str], optional): the dict dfs to look at. Defaults to RSEM_TRANSCRIPTS.
    drop (list[str], optional): the genes to drop. Defaults to [].
    index (str, optional): the index to use. Defaults to 'transcript_id'.

  Returns:
    dict(str: pd.df): the subsetted dfs
  """
  print('subsetting '+index+' columns')
  rename_transcript = {}
  missing = []
  for val in filenames:
    if len(rename_transcript) == 0 and index == "transcript_id":
      for _, v in files[val].iterrows():
        if v['gene_id'].split('.')[0] in gene_rename:
          rename_transcript[v['transcript_id'].split('.')[0]] = gene_rename[
              v['gene_id'].split('.')[0]].split(' (')[0] + ' (' + v.transcript_id.split('.')[0] + ')'
        else:
          missing.append(v.gene_id.split('.')[0])
      print('missing: '+str(len(missing))+' genes')
    file = files[val].drop(columns=drop).set_index(index)
    file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
    r = [i.split('.')[0] for i in file.index]
    dup = h.dups(r)
    if len(dup) > 0:
      print(dup)
      raise ValueError('duplicate '+index)
    file.index = r
    file = file.rename(index=rename_transcript if len(
        rename_transcript) != 0 else gene_rename).T
    files[val] = file
  return files


def extractProtCod(files, mybiomart, protcod_rename,
                   filenames=RSEMFILENAME_GENE, dropNonMatching=False,
                   rep=('genes', 'proteincoding_genes')):
  """extracts protein coding genes from a merged RSEM gene dataframe and a biomart dataframe

  Args:
    files (dict(str: pd.dfs)): the rsem transcripts dfs to subset samples x genes
    mybiomart (pd.df): the biomart dataframe should contain the following columns: 
      'ensembl_gene_id', 'entrezgene_id', 'gene_biotype'
    protcod_rename (dict(str, str)): the protein coding gene renaming dict 
      (here we expect a dict of ensembl transcript ids: gene names)
    filenames (list[str], optional): the dict dfs to look at. Defaults to RSEMFILENAME_GENE.
    rep (tuple, optional): how to rename the protein gene subseted df copies in the dict. Defaults to ('genes', 'proteincoding_genes').

  Raises:
    ValueError: if the biomart dataframe does not contain the required columns

  Returns:
    dict(str: pd.df): the subsetted dfs
  """
  for val in filenames:
    name = val.replace(rep[0], rep[1])
    files[name] = files[val].drop(columns='transcript_id').set_index('gene_id')
    files[name] = files[name][(files[name].sum(1) != 0) & (files[name].var(1) != 0)]
    r = [i.split('.')[0] for i in files[name].index]
    dup = h.dups(r)
    if len(dup) > 0:
        print(dup)
        raise ValueError('duplicate genes')
    files[name].index = r
    files[name] = files[name][files[name].index.isin(set(mybiomart.ensembl_gene_id))].rename(index=protcod_rename)
    # removing genes that did not match.. pretty unfortunate
    if dropNonMatching:
      files[name] = files[name].loc[[i for i in files[name].index if ' (' in i]]
    # check: do we have any duplicates?
    # if we do, managing duplicates
    if len(set(h.dups(files[name].index.tolist()))) > 0:
      print("we have duplicate gene names!!")
      for dup in h.dups(files[name].index):
        a = files[name].loc[dup].sum()
        files[name].drop(index=dup)
        files[name].loc[dup] = a
    files[name] = files[name].T

  return files


async def ssGSEA(tpm_genes, geneset_file=SSGSEAFILEPATH, recompute=True):
  """the way we run ssGSEA on the CCLE dataset

  Args:
    tpm_genes (pd.df): the tpm genes dataframe
    geneset_file (str, optional): the path to the geneset file. Defaults to SSGSEAFILEPATH.

  Returns:
    pd.df: the ssGSEA results
  """
  tpm_genes = tpm_genes.copy()
  tpm_genes.columns = [i.split(' (')[0] for i in tpm_genes.columns]

  # summing the different exons/duplicates
  for i in h.dups(tpm_genes.columns):
    val = tpm_genes[i].sum(1)
    tpm_genes = tpm_genes.drop(columns=i)
    tpm_genes[i] = val

  # total size of data
  print("total size of data")
  print(len(set([val for val in tpm_genes.columns if '.' not in val])))
  tpm_genes = pd.DataFrame(data=zscore(np.log2(
      tpm_genes+1), nan_policy="omit"), columns=tpm_genes.columns, index=tpm_genes.index)

  # MAYBE NOT NEEDED
  #### merging splicing variants into the same gene
  #counts_genes_merged, _, _= h.mergeSplicingVariants(counts_genes.T, defined='.')

  enrichments = (await rna.gsva(tpm_genes.T,
                                geneset_file=geneset_file, method='ssgsea', recompute=recompute)).T
  enrichments.index = [i.replace('.', '-') for i in enrichments.index]
  return enrichments


def saveFiles(files, folder=TMP_PATH, rep=('rsem', 'expression')):
  """
  saves the files in the dict to the folder

  Args:
    files (dict(str: pd.df)): the dfs to save
    folder (str, optional): the folder to save the files. Defaults to TMP_PATH.
    rep (tuple, optional): how to rename (parts of) the files. Defaults to ('rsem', 'expression').
  """
  print('storing files in {}'.format(folder))
  for k, val in files.items():
    val.to_csv(os.path.join(folder, k.replace(
        rep[0], rep[1])+'.csv'))
    if 'tpm' in k:
      val.apply(lambda x: np.log2(x + 1)).to_csv(os.path.join(folder,
                                                              k.replace(rep[0], rep[1]) +
                                                              '_logp1.csv'))


async def postProcess(refworkspace, samplesetname,
                      save_output="", doCleanup=False,
                      colstoclean=[], ensemblserver=ENSEMBL_SERVER_V,
                      todrop=[], samplesetToLoad="all", priority=[],
                      geneLevelCols=RSEMFILENAME_GENE,
                      trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
                      ssGSEAcol="genes_tpm", renamingFunc=None, useCache=False,
                      dropNonMatching=False, recompute_ssgsea=True,
                      compute_enrichment=True,
                      ):
  """postprocess a set of aggregated Expression table from RSEM in the CCLE way

  (usually using the aggregate_RSEM terra worklow)

  Args:
    refworkspace (str): terra workspace where the ref data is stored
    sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
    save_output (str, optional): whether to save our data. Defaults to "".
    doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
    colstoclean (list, optional): the columns to clean in the terra workspace. Defaults to [].
    ensemblserver (str, optional): ensembl server biomart version . Defaults to ENSEMBL_SERVER_V.
    todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
    priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
    useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
    samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
    geneLevelCols (list, optional): the columns that contain the gene level 
      expression data in the workspace. Defaults to RSEMFILENAME_GENE.
    trancriptLevelCols (list, optional): the columns that contain the transcript 
      level expression data in the workspacce. Defaults to RSEMFILENAME_TRANSCRIPTS.
    ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
    renamingFunc (function, optional): the function to use to rename the sample columns
      (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
    compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
    dropNonMatching (bool, optional): whether to drop the non matching genes 
      between entrez and ensembl. Defaults to False.
    recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
  """
  if not samplesetToLoad:
    samplesetToLoad = samplesetname
  refwm = dm.WorkspaceManager(refworkspace)
  if save_output:
    terra.saveWorkspace(refworkspace, save_output+'terra/')
  print("load QC and generate QC report")
  samplesinset = [i['entityName'] for i in refwm.get_entities(
      'sample_set').loc[samplesetname].samples]

  _, lowqual, failed = myQC.plot_rnaseqc_results(refworkspace, samplesinset,
                                                 save=bool(save_output),
                                                 output_path=save_output+"rna_qcs/")

  failed = failed.index.tolist()
  print('those samples completely failed qc: ', failed)
  print("rescuing from whitelist list: ", set(failed) & set(priority))
  failed = list(set(failed) - set(priority))
  failed.extend(todrop)

  if doCleanup:
    print("cleaninp up data")
    res = refwm.get_samples()
    for val in colstoclean:
      if val in res.columns.tolist():
        refwm.disable_hound().delete_entity_attributes(
            'sample', res[val], delete_files=True)
      else:
        print(val+' not in the workspace\'s data')

  print("generating gene names")

  mybiomart = h.generateGeneNames(
      ensemble_server=ensemblserver, useCache=useCache)
  # creating renaming index, keeping top name first
  gene_rename = {}
  for _, i in mybiomart.iterrows():
    if i.ensembl_gene_id not in gene_rename:
      gene_rename.update({i.ensembl_gene_id: i.hgnc_symbol+' ('+i.ensembl_gene_id+')'})
  protcod_rename = {}
  for _, i in mybiomart[(~mybiomart.entrezgene_id.isna()) &
                            (mybiomart.gene_biotype == 'protein_coding')].iterrows():
    if i.ensembl_gene_id not in protcod_rename:
      protcod_rename.update({i.ensembl_gene_id: i.hgnc_symbol+' ('+str(int(i.entrezgene_id))+')'})

  print("loading files")
  files, renaming = loadFromRSEMaggregate(refworkspace, todrop=failed, filenames=trancriptLevelCols+geneLevelCols,
                                          sampleset=samplesetToLoad, renamingFunc=renamingFunc)
  if save_output:
    h.dictToFile(renaming, save_output+"rna_sample_renaming.json")
    lowqual.to_csv(save_output+"rna_lowqual_samples.csv")
    h.listToFile(failed, save_output+"rna_failed_samples.txt")
  print("renaming files")
  # gene level
  if len(geneLevelCols) > 0:
    #import pdb; pdb.set_trace()
    files = extractProtCod(files, mybiomart[mybiomart.gene_biotype == 'protein_coding'],
                           protcod_rename, dropNonMatching=dropNonMatching,
                           filenames=geneLevelCols)
    # assert {v.columns[-1] for k,v in files.items()} == {'ACH-000052'}
    files = subsetGenes(files, gene_rename, filenames=geneLevelCols,
                        index="gene_id", drop="transcript_id")
    # assert {v.columns[-1] for k,v in files.items()} == {'ACH-000052'}
  if len(trancriptLevelCols) > 0:
    files = subsetGenes(
        files, gene_rename, filenames=trancriptLevelCols, drop="gene_id", index="transcript_id")

  if compute_enrichment:
    print("doing ssGSEA")
    enrichments = await ssGSEA(files[ssGSEAcol], recompute=recompute_ssgsea)
    print("saving files")
    enrichments.to_csv(save_output+'gene_sets_all.csv')
  saveFiles(files, save_output)
  print("done")

  return files, failed, samplesinset, renaming, lowqual, enrichments if compute_enrichment else None


async def _CCLEPostProcessing(refworkspace=RNAWORKSPACE, samplesetname=SAMPLESETNAME, refsheet_url=REFSHEET_URL,
                              colstoclean=['fastq1', 'fastq2',
                                            'recalibrated_bam', 'recalibrated_bam_index'],
                              ensemblserver=ENSEMBL_SERVER_V, doCleanup=True,
                              my_id=MY_ID, mystorage_id=MYSTORAGE_ID, samplesetToLoad="all",
                              tocompare={"genes_expected_count": "CCLE_RNAseq_reads",
                                          "genes_tpm": "CCLE_expression_full",
                                          "proteincoding_genes_tpm": "CCLE_expression"},
                              sheetname=SHEETNAME, sheetcreds=SHEETCREDS, todrop=KNOWN_DROP,
                              prevcounts='ccle',
                              taiga_dataset=TAIGA_EXPRESSION, minsimi=0.95, dropNonMatching=True,
                              dataset_description=RNAseqreadme, **kwargs):
  """the full CCLE Expression post processing pipeline (used only by CCLE)

  @see postprocessing() to reproduce our analysis and for parameters

  Args:
    refworkspace (str): terra workspace where the ref data is stored
    sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
    save_output (str, optional): whether to save our data. Defaults to "".
    doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
    colstoclean (list, optional): the columns to clean in the terra workspace. Defaults to [].
    ensemblserver (str, optional): ensembl server biomart version . Defaults to ENSEMBL_SERVER_V.
    todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
    priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
    useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
    samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
    geneLevelCols (list, optional): the columns that contain the gene level 
      expression data in the workspace. Defaults to RSEMFILENAME_GENE.
    trancriptLevelCols (list, optional): the columns that contain the transcript 
      level expression data in the workspacce. Defaults to RSEMFILENAME_TRANSCRIPTS.
    ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
    renamingFunc (function, optional): the function to use to rename the sample columns
      (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
    compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
    dropNonMatching (bool, optional): whether to drop the non matching genes 
      between entrez and ensembl. Defaults to False.
    recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
    prevcounts (str, optional): the previous counts to use to QC the data for the release. Defaults to 'ccle'.
    taiga_dataset (str, optional): the taiga dataset path to use for uploading results. Defaults to TAIGA_EXPRESSION.
    minsimi (float, optional): the minimum similarity to use for comparison to previous dataset. Defaults to 0.95.
    dataset_description (str, optional): the taiga dataset description to use. Defaults to RNAseqreadme.
    sheetname (str, optional): the sheet name to use for updating the sample tracker 
      (should be an actual google spreadsheet). Defaults to SHEETNAME.
    sheetcreds (str, optional): path to the google sheet credentials file to use. Defaults to SHEETCREDS.
    refsheet_url (str, optional): the url of the google sheet containing the data. Defaults to REFSHEET_URL.
    tocompare (dict, optional): the columns to compare. Defaults to {"genes_expected_count": "CCLE_RNAseq_reads", "genes_tpm": "CCLE_expression_full", "proteincoding_genes_tpm": "CCLE_expression"}.
    my_id (str, optional): path to the id containing file for google sheet. Defaults to MY_ID.
    mystorage_id (str, optional): path to the id containing file for google storage. Defaults to MYSTORAGE_ID.
  """
  from taigapy import TaigaClient
  tc = TaigaClient()

  if prevcounts == "ccle":
    prevcounts = tc.get(name=TAIGA_ETERNAL,
           file='CCLE_RNAseq_reads')

  sheets = Sheets.from_files(my_id, mystorage_id)
  ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

  todrop += ccle_refsamples[ccle_refsamples.blacklist == 1].index.tolist()
  priority = ccle_refsamples[ccle_refsamples.prioritized == 1].index.tolist()

  def rn(r, todrop):
    renaming = track.removeOlderVersions(
        names=r, refsamples=ccle_refsamples[ccle_refsamples.datatype == "rna"],
        priority="prioritized")
    # if we have a replaceable failed version in our dataset
    rename = solveQC(ccle_refsamples, todrop)
    for k, _ in renaming.items():
      if k in rename:
        renaming[rename[k]] = renaming.pop(k)
    return renaming

  folder = os.path.join("temp", samplesetname, "")
  h.createFoldersFor(folder)
  files, failed, _, renaming, lowqual, _ = await postProcess(refworkspace, samplesetname,
                                                    save_output=folder, doCleanup=doCleanup, priority=priority,
                                                    colstoclean=colstoclean, ensemblserver=ensemblserver,
                                                    todrop=todrop, samplesetToLoad=samplesetToLoad,
                                                    geneLevelCols=RSEMFILENAME_GENE,
                                                    trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS, #compute_enrichment=False,
                                                    ssGSEAcol="genes_tpm", renamingFunc=rn,
                                                    dropNonMatching=dropNonMatching, **kwargs)

  print("doing validation")
  nonoverlap = set(prevcounts.columns) ^ set(
      files.get('genes_expected_count').columns)
  print("number of non overlaping genes:")
  print(len(nonoverlap))
  # have we lost any samples compared to last release?
  lost = set(prevcounts.index) - set(files.get('genes_expected_count').index)
  print("of which, lost genes:")
  print(lost)
  # do we have samples that are missanotated compared to previous releases (replicate level)
  #notindataset, missannotated, unmatched = findMissAnnotatedReplicates(replevel, prevcounts, renaming)
  #for k,v in unmatched.items():
  #    if ccle_refsamples.loc[k].arxspan_id!=v:
  #        print(k,v)
  # do we have samples that are missanotated compared to previous releases (sample level)
  unmatched = rna.getDifferencesFromCorrelations(
      files.get('genes_expected_count'), prevcounts, minsimi=minsimi)
  print("differences in correlations against the previous release")
  print(unmatched)
  # Is it because of  duplicate version?
  print('do we see it as a duplicate in the tracker?')
  rnasamples = ccle_refsamples[ccle_refsamples.datatype == 'rna']
  for i, _ in unmatched:
    print(len(rnasamples[rnasamples.arxspan_id == i]))

  #CCLE_expression, CCLE_expression_full, ,
  print("comparing to previous release")
  #h.compareDfs(files["rsem_transcripts_tpm"], tc.get(name=TAIGA_ETERNAL, file='CCLE_RNAseq_transcripts'))
  #h.compareDfs(files["rsem_transcripts_expected_count"], tc.get(name=TAIGA_ETERNAL, file='CCLE_expression_transcripts_expected_count'))
  # h.compareDfs(enrichments, tc.get(name=TAIGA_ETERNAL, file='CCLE_fusions_unfiltered'))
  for key, val in tocompare.items():
    _, omissmatchCols, _, omissmatchInds, newNAs, new0s = h.compareDfs(
        files[key], tc.get(name=TAIGA_ETERNAL, file=val))
    print(key)
    assert len(omissmatchCols) == 0
    assert len(omissmatchInds) == 0
    assert newNAs == 0
    print("New 0s: ", new0s)

  print("updating the tracker")

  updateTracker(set(renaming.keys()) - set(['transcript_id(s)']), failed,
                lowqual[lowqual.sum(1) > 3].index.tolist(),
                ccle_refsamples, samplesetname, refworkspace,
                sheetname=sheetname, sheetcreds=sheetcreds, todrop=todrop)

  print("uploading to taiga")
  tc.update_dataset(changes_description="new "+samplesetname+" release!",
                    dataset_permaname=taiga_dataset,
                    upload_files=[
                      {
                        "path": folder+"proteincoding_genes_tpm_logp1.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"transcripts_tpm_logp1.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"genes_tpm_logp1.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"genes_tpm.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"transcripts_tpm.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"proteincoding_genes_tpm.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"transcripts_expected_count.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"proteincoding_genes_expected_count.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      },
                      {
                        "path": folder+"genes_expected_count.csv",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8"
                      }
                      # {
                      #     "path": folder+'gene_sets_all.csv',
                      #     "format": "NumericMatrixCSV",
                      #     "encoding": "utf-8"
                      # },
                  ],
                  upload_async=False,
                  dataset_description=dataset_description)
  print("done")