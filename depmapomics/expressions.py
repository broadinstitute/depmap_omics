import os.path
import asyncio

import dalmatian as dm
import pandas as pd
import numpy as np
from scipy.stats import zscore 

from genepy.google import gcp
from genepy.utils import helper as h
from genepy.google.google_sheet import dfToSheet
from genepy import rna, terra

from depmapomics import terra as myterra
from depmapomics.qc import rna as myQC
from depmapomics import utils
from depmapomics import tracker as track

from depmapomics.config import *

from gsheets import Sheets
from taigapy import TaigaClient
tc = TaigaClient()


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
  

def solveQC(tracker, failed, save=""):
  """
  # TODO todocument
  """
  newfail = []
  rename = {}
  # finding other replicates to solve failed ones
  for val in failed:
    a = tracker.loc[val].arxspan_id
    res = tracker[(tracker.datatype == 'rna')
                        & (tracker.arxspan_id == a)]
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
    

def updateTracker(refworkspace, selected, failed, lowqual, tracker, samplesinset, samplesetname, 
                  sheetname=SHEETNAME, sheetcreds=SHEETCREDS, 
                  onlycol=STARBAMCOLTERRA, newgs=RNAGSPATH38,
                  dry_run=False, keeppath=False, qcname="star_logs", match=".Log.final.out"):
  """
  # TODO: to document
  """
  
  starlogs = myterra.getQC(workspace=refworkspace, only=samplesinset,
                            qcname=qcname, match=match)
  for k, v in starlogs.items():
    if k == 'nan':
      continue
    if tracker.loc[k, 'bam_qc'] != v[0]:
      tracker.loc[k, 'bam_qc'] = v[0]
  
  ## copy star bam file to our cclebams/rnasq_hg38/ bucket
  renamed, _ = terra.changeGSlocation(workspacefrom=refworkspace, newgs=newgs, 
                                      onlysamples=samplesinset, onlycol=onlycol, 
                                      entity="sample", keeppath=keeppath, dry_run=dry_run)

  tracker.loc[samplesinset, ['legacy_size', 'legacy_crc32c_hash']
                      ] = tracker.loc[samplesinset][['size', 'crc32c_hash']].values
  tracker.loc[samplesinset, HG38BAMCOL] = renamed[onlycol[:2]].values
  tracker.loc[samplesinset, 'size'] = [gcp.extractSize(
    i)[1] for i in gcp.lsFiles(renamed[onlycol[0]].tolist(), '-l')]
  tracker.loc[samplesinset, 'crc32c_hash'] = [gcp.extractHash(
    i) for i in gcp.lsFiles(renamed[onlycol[0]].tolist(), '-L')]
  tracker.loc[samplesinset, 'md5_hash'] = [gcp.extractHash(
    i, "md5") for i in gcp.lsFiles(renamed[onlycol[0]].tolist(), '-L')]
  
  tracker.loc[selected, samplesetname]=1
  tracker.loc[samplesinset, ['low_quality', 'blacklist', 'prioritized']] = 0
  tracker.loc[lowqual,'low_quality']=1
  tracker.loc[failed, 'blacklist'] = 1
  dfToSheet(tracker, sheetname, secret=sheetcreds)
  print("updated the sheet, please reactivate protections")


def loadFromRSEMaggregate(refworkspace, todrop=[], filenames=RSEMFILENAME,
                          sampleset="all", renamingFunc=None):
  """
  #TODO: to document
  Args:
  ----

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
      renaming.update({i:i for i in file.columns[2:] if i not in todrop})
    renaming.update({'transcript_id(s)': 'transcript_id'})
    # we remove the failed samples where we did not found anything else to replace them with
    files[val] = file[file.columns[:2].tolist()+[i for i in file.columns[2:]
                        if i in set(renaming.keys())]].rename(columns=renaming)
  return files, renaming


def subsetGenes(files, gene_rename, filenames=RSEM_TRANSCRIPTS,
                              drop=[], index="transcript_id"):
  """
  # TODO: to document
  """
  print('subsetting transcript columns')
  rename_transcript = {}
  missing = []
  for val in filenames:
    file = files[val].drop(columns=drop).set_index(index)
    file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
    r = [i.split('.')[0] for i in file.index]
    dup = h.dups(r)
    if len(dup) > 0:
      print(dup)
      raise ValueError('duplicate genes')
    file.index = r
    if len(rename_transcript) == 0 and index == "transcript_id":
      for _, v in files[val].iterrows():
        if v.gene_id.split('.')[0] in gene_rename:
          rename_transcript[v.transcript_id] = gene_rename[
              v.gene_id.split('.')[0]].split(' (')[0] + ' (' + v.transcript_id + ')'
        else:
          missing.append(v.gene_id.split('.')[0])
      print('missing: '+str(len(missing))+' genes')
    files[val].rename(index=rename_transcript if len(
      rename_transcript) != 0 else gene_rename).T
  return files


def extractProtCod(files, ensembltohgnc, protcod_rename,
                   filenames=RSEM_TRANSCRIPTS,
                  rep=('genes', 'proteincoding_genes')):
  """
  # TODO: to document
  """
  for val in filenames:
    name = val.replace(rep[0], rep[1])
    files[name] = files[val][files[val].columns.isin(set(ensembltohgnc.ensembl_gene_id))]
    files[name] = files[name].rename(columns=protcod_rename)
    # removing genes that did not match.. pretty unfortunate
    files[name] = files[name][[i for i in files[name].columns if ' (' in i]]
    # check: do we have any duplicates?
    # if we do, managing duplicates
    if len(set(h.dups(files[name].columns.tolist()))) > 0:
      print("we have duplicate gene names!!")
      for dup in h.dups(files[name].columns):
        a = files[name][dup].sum()
        files[name].drop(columns=dup)
        files[name][dup] = a
          
  return files


async def ssGSEA(tpm_genes, pathtogenepy=PATHTOGENEPY,
                 geneset_file=SSGSEAFILEPATH):
  """the way we run ssGSEA on the CCLE dataset

  Args:
      tpm_genes ([type]): [description]
      pathtogenepy ([type], optional): [description]. Defaults to PATHTOGENEPY.
      geneset_file ([type], optional): [description]. Defaults to SSGSEAFILEPATH.

  Returns:
      [type]: [description]
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

  enrichments = (await rna.gsva(tpm_genes.T, pathtogenepy=pathtogenepy,
                                geneset_file=geneset_file, method='ssgsea')).T
  enrichments.index = [i.replace('.', '-') for i in enrichments.index]
  return enrichments


def saveFiles(files, folder=TMP_PATH, rep=('rsem', 'expression')):
  """
  # TODO: to document
  """
  print('storing files in {}'.format(folder))
  for k, val in files.items():
    val.to_csv(os.path.join(folder, k.replace(
      rep[0], rep[1])+'.csv'))
    if 'tpm' in k:
      val.apply(lambda x: np.log2(x + 1)).to_csv(os.path.join(folder,
                                                k.replace(rep[0],rep[1]) + 
                                                '_logp1.csv'))


def postProcess(refworkspace, samplesetname,
              save_output="", doCleanup=False,
              colstoclean=[], ensemblserver=ENSEMBL_SERVER_V,
              todrop=[], samplesetToLoad="all", priority=[],
              geneLevelCols=RSEMFILENAME_GENE,
              trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
              ssGSEAcol="gene_tpm", renamingFunc=None, useCache=False
              ):
  """postprocess a set of aggregated Expression table from RSEM in the CCLE way
  
  (usually using the aggregate_RSEM terra worklow) 

  Args:
      refworkspace ([type]): [description]
      samplesetname ([type]): [description]
      save_output (str, optional): [description]. Defaults to "".
      doCleanup (bool, optional): [description]. Defaults to False.
      colstoclean (list, optional): [description]. Defaults to [].
      ensemblserver ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
      todrop (list, optional): [description]. Defaults to [].
      samplesetToLoad (str, optional): [description]. Defaults to "all".
      priority (list, optional): [description]. Defaults to [].
      geneLevelCols ([type], optional): [description]. Defaults to RSEMFILENAME_GENE.
      trancriptLevelCols ([type], optional): [description]. Defaults to RSEMFILENAME_TRANSCRIPTS.
      ssGSEAcol (str, optional): [description]. Defaults to "gene_tpm".
      renamingFunc ([type], optional): [description]. Defaults to None.
      useCache (bool, optional): [description]. Defaults to False.

  Returns:
      [type]: [description]
  """
  if not samplesetToLoad:
    samplesetToLoad = samplesetname
  refwm = dm.WorkspaceManager(refworkspace)
  if save_output:
    terra.saveConfigs(refworkspace, save_output+'terra/')
  print("load QC and generate QC report")
  samplesinset = [i['entityName'] for i in refwm.get_entities(
    'sample_set').loc[samplesetname].samples]

  _, lowqual, failed = myQC.plot_rnaseqc_results(refworkspace, samplesinset,
                                                  save=bool(save_output),  
                                                  output_path=save_output+"rna_qcs/")

  failed = failed.index.tolist()
  print('those samples completely failed qc: ', failed)
  print("rescuing from whitelist list: ", set(failed) & set(priority) )
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
  gene_rename, protcod_rename, ensembltohgnc = utils.generateGeneNames(
    ensemble_server=ensemblserver, useCache=useCache)
  
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
    files = subsetGenes(files, gene_rename, filenames=geneLevelCols, 
                                drop="transcript_id", index="gene_id")

    files = extractProtCod(files, ensembltohgnc[ensembltohgnc.gene_biotype == 'protein_coding'], 
                            protcod_rename, filenames=geneLevelCols)
  if len(trancriptLevelCols) > 0:
    files = subsetGenes(
      files, gene_rename, filenames=trancriptLevelCols, drop="gene_id", index="transcript_id")
  
  print("doing ssGSEA")
  enrichments = asyncio.run(ssGSEA(files[ssGSEAcol]))
  print("saving files")
  enrichments.to_csv(save_output+'gene_sets_all.csv')
  saveFiles(files, save_output)
  print("done")

  return files, enrichments, failed, samplesinset, renaming, lowqual
  

def CCLEPostProcessing(refworkspace=rnaworkspace, samplesetname=SAMPLESETNAME, refsheet_url=REFSHEET_URL,
                      colstoclean=['fastq1', 'fastq2', 'recalibrated_bam', 'recalibrated_bam_index'], 
                      ensemblserver=ENSEMBL_SERVER_V, doCleanup=True,
                      my_id=MY_ID, mystorage_id=MYSTORAGE_ID, samplesetToLoad="all", 
                      tocompare={"genes_expected_count": "CCLE_RNAseq_reads",
                                "genes_tpm": "CCLE_expression_full",
                                "proteincoding_genes_tpm": "CCLE_expression"},
                      sheetname=SHEETNAME, sheetcreds=SHEETCREDS,
                      prevcounts = tc.get(name=TAIGA_ETERNAL, file='CCLE_RNAseq_reads'),
                      taiga_dataset="expression-d035", minsimi=0.95, 
                      dataset_description=RNAseqreadme, **kwargs):
  """the full CCLE Expression post processing pipeline (used only by CCLE)

  see postprocessing() to reproduce our analysis

  Args:
      refworkspace ([type], optional): [description]. Defaults to rnaworkspace.
      samplesetname ([type], optional): [description]. Defaults to SAMPLESETNAME.
      refsheet_url ([type], optional): [description]. Defaults to REFSHEET_URL.
      colstoclean (list, optional): [description]. Defaults to ['fastq1', 'fastq2', 'recalibrated_bam', 'recalibrated_bam_index'].
      ensemblserver ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
      doCleanup (bool, optional): [description]. Defaults to True.
      my_id ([type], optional): [description]. Defaults to MY_ID.
      mystorage_id ([type], optional): [description]. Defaults to MYSTORAGE_ID.
      samplesetToLoad (str, optional): [description]. Defaults to "all".
      tocompare (dict, optional): [description]. Defaults to {"genes_expected_count": "CCLE_RNAseq_reads", "genes_tpm": "CCLE_expression_full", "proteincoding_genes_tpm": "CCLE_expression"}.
      sheetname ([type], optional): [description]. Defaults to SHEETNAME.
      sheetcreds ([type], optional): [description]. Defaults to SHEETCREDS.
      prevcounts ([type], optional): [description]. Defaults to tc.get(name=TAIGA_ETERNAL, file='CCLE_RNAseq_reads').
      taiga_dataset (str, optional): [description]. Defaults to "expression-d035".
      minsimi (float, optional): [description]. Defaults to 0.95.
      dataset_description ([type], optional): [description]. Defaults to RNAseqreadme.

  Returns:
      [type]: [description]
  """
  
  sheets = Sheets.from_files(my_id, mystorage_id)
  ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)
  
  todrop = ccle_refsamples[ccle_refsamples.blacklist==1].index.tolist()
  priority = ccle_refsamples[ccle_refsamples.prioritized==1].index.tolist()

  def rn(r, todrop):
    renaming = track.removeOlderVersions(
      names=r, refsamples=ccle_refsamples[ccle_refsamples.datatype=="rna"], 
      priority="prioritized")
    # if we have a replaceable failed version in our dataset
    rename = solveQC(ccle_refsamples, todrop)
    for k, _ in renaming.items():
      if k in rename:
        renaming[rename[k]] = renaming.pop(k)
    return renaming

  folder = os.path.join("temp", samplesetname,"")
  files, _, failed, samplesinset, renaming, lowqual = postProcess(refworkspace, samplesetname,
                          save_output=folder, doCleanup=doCleanup, priority=priority,
                          colstoclean=colstoclean, ensemblserver=ensemblserver,
                          todrop=todrop, samplesetToLoad=samplesetToLoad,
                          geneLevelCols=RSEMFILENAME_GENE,
                          trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
                          ssGSEAcol="gene_tpm", renamingFunc=rn, **kwargs)
  
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
    assert omissmatchCols == 0
    assert omissmatchInds == 0
    assert newNAs == 0
    assert new0s == 0
  
  print("updating the tracker")
  updateTracker(refworkspace, set(renaming.keys()) - set(['transcript_id(s)']), failed,
                lowqual[lowqual.sum(1) > 3].index.tolist(),
                ccle_refsamples, samplesinset, samplesetname, 
                sheetname=sheetname, sheetcreds=sheetcreds)
  
  print("uploading to taiga")
  tc.update_dataset(changes_description="new "+samplesetname+" release!",
                    dataset_permaname=taiga_dataset,
                    upload_files=[
                      {
                          "path": "temp/"+samplesetname+"/expression_proteincoding_tpm_logp1.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_transcripts_tpm_logp1.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_genes_tpm_logp1.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_genes_tpm.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_transcripts_tpm.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_proteincoding_tpm.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_transcripts_expectedcount.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_proteincoding_expectedcount.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": "temp/"+samplesetname+"/expression_genes_expectedcount.csv",
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                      {
                          "path": 'temp/'+samplesetname+'/gene_sets_all.csv',
                          "format": "NumericMatrixCSV",
                          "encoding": "utf-8"
                      },
                    ],
                    upload_async=False,
            dataset_description=dataset_description)
  print("done")
