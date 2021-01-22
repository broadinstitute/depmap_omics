# -*- coding: utf-8 -*-
# Jérémie Kalfon
# for BroadInsitute
# in 2019


####
#
# HELPER FUNC  ######################################
#
#
from scipy.stats import pearsonr
import pdb
import pandas as pd
import numpy as np
import dalmatian as dm
from taigapy import TaigaClient
from time import sleep
tc = TaigaClient()
import os
import ipdb
import sys
print("you need to have JKBio in your path:\ne.g. have installed JKBio in the same folder as ccle_processing")
from JKBio import terra
from JKBio import sequencing as seq
from JKBio.utils import helper as h
from JKBio.google import gcp
from collections import Counter
from gsheets import Sheets
import seaborn as sns
from matplotlib import pyplot as plt

#####################
# Const Variables
#####################


sheets = Sheets.from_files('~/.client_secret.json', '~/.storage.json')


CHROMLIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX']


extract_defaults = {
    'name': 'sample_alias',
    'bai': 'crai_or_bai_path',
    'bam': 'cram_or_bam_path',
    'ref_bam': 'internal_bam_filepath',
    'ref_type': 'datatype',
    'ref_bai': "internal_bai_filepath",
    'version': 'version',
    'primary_disease': 'primary_disease',
    'ref_arxspan_id': 'arxspan_id',
    'ref_name': 'stripped_cell_line_name',
    'source': 'source',
    'size': 'size',
    "prev_size":"legacy_size",
    'from_arxspan_id': 'individual_alias',
    'ref_id': 'sample_id',
    'PDO_id':'PDO',
    "update_time":"update_time",
    'from_patient_id': 'individual_alias',
    'patient_id': 'participant_id',
    'ref_date': 'date_sequenced',
    'hs_hs_library_size': 'hs_hs_library_size',
    'hs_het_snp_sensitivity': 'hs_het_snp_sensitivity',
    'hs_mean_bait_coverage': 'hs_mean_bait_coverage',
    'hs_mean_target_coverage': 'hs_mean_target_coverage',
    'hs_on_target_bases': 'hs_on_target_bases',
    'total_reads': 'total_reads',
    'release_date': 'sequencing_date',
    'hash': 'crc32c_hash',
    'mean_depth': 'mean_depth'
}

MINSIZES = {
    'rna': 2000000000,
    'wes': 3000000000,
    'wgs': 50000000000,
}

# found same patient
samepatient = [["ACH-000635", "ACH-000717", "ACH-000864", "ACH-001042", "ACH-001547"],
               ["ACH-002291", "ACH-001672"],
               ["ACH-001706", "ACH-001707"]]

# duplicate ach-id
dup = {"ACH-001620": "ACH-001605",
       "ACH-001621": "ACH-001606"}
# rename ccle_name TODO: ask becky what to do
rename = {"PEDS117": "CCLFPEDS0009T"}


#####################
# Loading Functions
#####################

def GetNewCellLinesFromWorkspaces(wto, wmfroms, sources, stype, maxage, refurl="",
                                  addonly=[], match='ACH', extract={},
                                  extract_defaults=extract_defaults, refsamples=None,
                                  participantslicepos=10, accept_unknowntypes=False,
                                  rename=dict(), recomputehash=False):
  """
  As GP almost always upload their data to a data workspace. we have to merge it to our processing workspace

  Will merge samples from a set of data workspaces to a processing workspace on Terra. Will only
  get a subset of the metadata and rename it.
  Will find out the duplicates based on the file size.
  Can also upload the bam files to a google storage bucket

  Args:
  -----
    wto: str the workspace where you want to create the tsvs
    wfroms: list[str] the workspaces where the Samples to add are stored
    sources: list[str] the corresponding source names
    stype: str sequencing type
    maxage: str earliest date of the bam file upload to be considered new
    refurl: str(url) the reference url for the cell line tracker spreadsheet (only if no refsamples)
    match: list[str]|str the possible values that a sample id need to contain to be considered valid
    refsamples: pdDataFrame with columns matching values is in "extract" for the right keys (see "extract_default")
    participantslicepos: int the length of the sample id string
    accept_unknowntypes: bool whether or not the sample type column for that sample can be different from "Tumor"
    rename: dict(str:str) mapping a wrong arxpand_id to a good arxspan id for known cases of misslabelling
    recomputehash: bool whether or not to recompute the hash of the bam file when loading it
    addonly: list of sample id that you only want to add
    match: list of substring(s) that has to be matched against the id of the samples to add them
    extract: if you want to specify what values should refer to which column names
      dict{
      'name':
      'bai':
      'bam':
      'source':
      'from_arxspan_id':
      ...} (see extract_defaults)
    extract_defaults: the full default dict to specificy what values should refer to which column names

  Returns:
  -------
    samples: a dataframe with the samples that were resolved by the tool (we still need to add some more annotations)
    pairs: the corresponding pair from matching known normals with known tumors
    wrongssamples: a dataframe containing samples that passed most QCs but couldn't be resolved

  Raise:
  -----
    Exception: when no new samples in this matrix
  """
  extract.update(extract_defaults)
  wto = dm.WorkspaceManager(wto)
  if type(match) is str and match:
    match = [match]
  if refurl:
    print('refsamples is overrided by a refurl')
    refsamples = sheets.get(refurl).sheets[0].to_frame().set_index("cds_sample_id")
  if refsamples is None:
    print('we do not have refsamples data. Using the wto workspace sample data instead')
    refsamples = wto.get_samples()
    # TODO: update directly the df if data is not already in here)
    refsamples[extract['ref_arxspan_id']] = [a.split('_')[0] for a in refsamples[extract['ref_arxspan_id']] if type(a) is str]
    if extract['hash'] not in refsamples.columns:
      refsamples[extract['hash']] = [gcp.extractHash(val) for val in gcp.lsFiles(
          [i for i in refsamples[extract["ref_bams"]] if type(i) is str and str(i) != 'NA'], "-L", 200)]
    if extract['size'] not in refsamples.columns:
      refsamples['size'] = [gcp.extractSize(i)[1] for i in gcp.lsFiles(refsamples[extract['bam']].tolist(), '-al', 200)]
    if extract['release_date'] not in refsamples.columns:
      refsamples[extract["ref_bams"]] = seq.getBamDate(refsamples[extract["ref_bams"]])
  refsamples[extract['release_date']] = list(h.datetoint(refsamples[extract["release_date"]].values, '/'))
  if stype not in set(refsamples[extract['ref_type']]):
    h.ask("we have never seen this type: " + stype + ", in the reference, continue?")
  # do NOT make refids a set; we use the num of occurences as way to determine what number to add to the sample id
  # filter refids to only include those that include the strings in the 'match' argument
  refsamples = refsamples[refsamples.index.str.contains('|'.join(match))]
  for match_substring in match:
    refsamples.index = [match_substring + i.split(match_substring)[-1] if match_substring in i else i for i in refsamples.index]
  refsamples.index = [i[:participantslicepos] for i in refsamples.index]
  print("Getting sample infos...")
  if type(sources) is str:
    sources = [sources]
  if type(wmfroms) is str:
    wmfroms = [wmfroms]
  sampless = pd.DataFrame()
  wrongsampless = pd.DataFrame()
  for source, wmfrom in zip(sources, wmfroms):
    broken_bams = []
    wmfrom = dm.WorkspaceManager(wmfrom)
    samples = wmfrom.get_samples().replace(np.nan, '', regex=True).reset_index()
    # keep samples that contain the match requirement (e.g. ACH for DepMap IDs)

    print("\nThe shape of the sample tsv from " + str(wmfrom) + ": " + str(samples.shape))

    # remove true duplicates from consideration
    print("Identifying any true duplicates by checking file hashes (this runs for each data source)...")
    print("This step can take a while as we need to use gsutil to check the size of each potential duplicate...")
    dups_to_remove = []
    # check for broken bam files; if broken, then remove from consideration
    # need to check for broken filepaths before checking if the sample is in Terra so that we don't
    # add a broken file path for a new participant
    foundfiles = gcp.lsFiles(samples[extract['bam']])
    broken_bams = set(samples[extract['bam']]) - set(foundfiles)
    print('These ' + str(len(broken_bams)) + ' bam file path do not exist: ' + str(broken_bams))

    wrongsamples = samples[(~samples[extract['bam']].isin(broken_bams)) & (~samples[extract['from_arxspan_id']].str.contains('|'.join(match)))]
    wrongsamples = extractFromWorkspace(wrongsamples, stype, recomputehash, extract)
    if wrongsamples is not None:
      wrongsamples = mapSamples(wrongsamples, source, extract)
      wrongsampless = pd.concat([wrongsampless, wrongsamples], sort=False)
    samples = samples[(~samples[extract['bam']].isin(broken_bams)) & (samples[extract['from_arxspan_id']].str.contains('|'.join(match)))]
    # getting correct arxspan id
    if samples is None:
      continue
    samples = extractFromWorkspace(samples, stype, recomputehash, extract)
    if samples is None:
      continue
    samples = mapSamples(samples, source, extract)
    samples = resolveFromWorkspace(samples, refsamples[refsamples[extract['ref_type']] == stype], match, participantslicepos,
                                   accept_unknowntypes, addonly, extract)
    if samples is None:
      continue
    sampless = pd.concat([sampless, samples], sort=False)

  if len(sampless) == 0:
    print("no new data available")
    return sampless, pd.DataFrame()

  sampless = assessAllSamples(sampless, refsamples, stype, rename, extract)
  # creating pairs
  pairs = setupPairsFromSamples(sampless, refsamples[refsamples[extract['ref_type']] == stype], extract)
  # I am trying to remove duplicates from samples without arxspan ids to then look more into them
  # and see if I have to get data for them or if I should just throw them out
  toremov = set()
  for k, val in wrongsampless.iterrows():
    withsamesize = wrongsampless[wrongsampless[extract["size"]] == val[extract["size"]]]
    if (val[extract["size"]] in sampless[extract["size"]].tolist()) or (val[extract["size"]] in refsamples[extract["size"]]):
      toremov.add(k)
    if len(withsamesize) > 1:
      for l, _ in withsamesize.iloc[1:].iterrows():
        toremov.add(l)
    #elif len(refsamples[refsamples[extract['size']] == withsamesize[extract["size"]][0]]):
      #toremov.add(k)
  for i in toremov:
    wrongsampless = wrongsampless.drop(i)
  for i, v in wrongsampless.iterrows():
    if not gcp.exists(v[extract['ref_bam']]):
      print(v.ccle_name)
      wrongsampless = wrongsampless.drop(i)
  a = len(sampless)
  sampless = deleteClosest(sampless,refsamples, extract['size'], extract['size'], extract['ref_arxspan_id'])
  sampless = deleteClosest(sampless,refsamples, extract['size'], extract['prev_size'], extract['ref_arxspan_id'])
  print('removed: '+str(a-len(sampless))+" samples from size alone (too similar to a replicate)")
  wrongsampless = wrongsampless[~wrongsampless[extract['size']].isin(set(refsamples[extract['size']]))]
  wrongsampless = wrongsampless[~wrongsampless[extract['size']].isin(set(refsamples[extract['prev_size']]))]
  wrongsampless = deleteClosest(
      wrongsampless, refsamples, extract['size'], extract['size'], extract['ref_arxspan_id'])
  wrongsampless = deleteClosest(
      wrongsampless, refsamples, extract['size'], extract['prev_size'], extract['ref_arxspan_id'])
  #removing duplicate PDOs
  a = len(sampless)
  wrongsampless = wrongsampless[~wrongsampless[extract['PDO_id']].isin(set(refsamples[extract['PDO_id']]))]
  sampless = sampless[~sampless[extract['PDO_id']].isin(
      set(refsamples[extract['PDO_id']]))]
  print('removed: '+str(a-len(sampless)) +
        " samples with duplicat PDO ids ")
  # removing anything too old
  a = len(sampless)
  wrongsampless = wrongsampless[wrongsampless[extract['update_time']] > maxage]
  sampless = sampless[sampless[extract['update_time']]>maxage]
  print('removed: '+str(a-len(sampless))+" samples that have not changed since last time (likely duplicate having been removed)")
  return sampless, pairs, wrongsampless


def deleteClosest(sampless, refsamples, size='size', ref_size='size', arxspid='arxspan_id'):
  """
  for a list of samples and a tracker, will find the index of the sample with the closest size

  if this sample is the same cell line, it will judge it to be a duplicate and remove it

  Args:
  -----
    sampless
    refsamples
    size
    ref_size
    arxspid

  Returns:
  --------
    samples:
  """
  sizes = refsamples[ref_size].tolist()
  for k,v in sampless.iterrows():
    if type(v[size]) is int:
      val = refsamples.iloc[sizes.index(h.closest(sizes, v[size]))]
      if val[arxspid] == v[arxspid]:
        sampless = sampless.drop(v.name)
  return sampless


def extractFromWorkspace(samples, stype, recomputeTime=True, recomputesize=True, recomputehash=True, extract={}):
  """
  # getting the hash
  Args:
  -----
    samples
    stype
    recomputeTime
    recomputesize
    recomputehash
    extract

  """
  extract.update(extract_defaults)
  if extract['hash'] not in samples.columns or recomputehash:
    samples[extract['hash']] = [gcp.extractHash(val) for val in gcp.lsFiles(samples[extract["bam"]].tolist(), "-L", 200)]
  lis = gcp.lsFiles(samples[extract['bam']].tolist(), '-al', 200)
  if extract['size'] not in samples.columns or recomputesize:
    samples[extract['size']] = [gcp.extractSize(i)[1] for i in lis]
  if extract['update_time'] not in samples.columns or recomputeTime:
    samples[extract['update_time']] = [gcp.extractTime(i) for i in lis]
  todrop = []
  for k, val in samples.iterrows():
    if val[extract['size']] < MINSIZES[stype]:
      todrop.append(k)
      print("too small size, removing sample: " + str(val[extract["from_arxspan_id"]]))
  samples = samples.drop(index=todrop)
  # getting the date released
  if len(samples) == 0:
    return None
  if extract['release_date'] not in samples.columns or recomputedate:
    samples[extract["release_date"]] = seq.getBamDate(samples[extract["bam"]])
  samples[extract['release_date']] = list(h.datetoint(samples[extract['release_date']].values))
  return samples


def mapSamples(samples, source, extract={}):
  """
  """
  # creating unique ids
  samples[extract['ref_id']] = ['CDS-' + h.randomString(stringLength=6, stype='all', withdigits=True) for _ in range(len(samples))]
  samples[extract['patient_id']] = ['PT-' + h.randomString(stringLength=8, stype='all', withdigits=True) for _ in range(len(samples))]
  samples.reset_index(drop=True, inplace=True)
  samples[extract['source']] = source

  # renamings
  samples = samples.rename(columns={extract['bam']: extract['ref_bam'],
                                    extract['bai']: extract['ref_bai'],
                                    extract['name']: extract['ref_name'],
                                    extract['from_arxspan_id']: extract['ref_arxspan_id']
                                    }).set_index(extract["ref_id"], drop=True)
  # subsetting
  samples = samples[[extract['ref_bam'], extract['ref_bai'], extract['ref_name'], extract["ref_arxspan_id"], extract["release_date"], extract["patient_id"], extract["hash"], extract['size'], extract['PDO_id'], extract['update_time']]]
  return samples


def resolveFromWorkspace(samples, refsamples, match, participantslicepos=10, accept_unknowntypes=True, addonly=[], extract={}):
  """
  """
  extract.update(extract_defaults)
  prevlen = len(samples)
  for match_substring in match:
    samples[extract['ref_arxspan_id']] = [(match_substring + i.split(match_substring)[-1]) if match_substring in i else i for i in samples[extract['ref_arxspan_id']]]
  samples[extract['ref_arxspan_id']] = [i[:participantslicepos] for i in samples[extract['ref_arxspan_id']]]
  print('we found and removed ' + str(prevlen - len(samples)) + ' samples which did not match our id names: ' + str(match))

  tolookfor = [val[extract['ref_bam']] for _, val in samples.iterrows() if val[extract['ref_arxspan_id']] in set(refsamples[extract['ref_arxspan_id']])]
  print("found " + str(len(tolookfor)) + ' likely replicate')
  sample_hash = {gcp.extractSize(val)[1]: gcp.extractSize(val)[0] for val in gcp.lsFiles(tolookfor, "-la")}
  dups_to_remove = [sample_hash[a] for a in set(sample_hash.keys()) & set(refsamples[extract['size']])]
  dups_to_remove.extend([sample_hash[a] for a in set(sample_hash.keys()) & set(refsamples[extract['prev_size']])])
  # remove the duplicates from consideration
  print("Len of samples before removal: " + str(len(samples)))
  print("Dups from this workspace has len " + str(len(dups_to_remove)) + ":\n " + str(dups_to_remove))
  # remove the samples with broken bam filepaths from consideration
  samples = samples[~samples[extract['ref_bam']].isin(dups_to_remove)]

  print("Len of samples after removal: " + str(len(samples)))
  if len(samples) == 0:
    return None

  # if only add some samples
  if len(addonly) > 0:
    samples = samples[samples[extract['ref_arxspan_id']].isin(addonly)]

  # unknown types
  if 'sample_type' in samples.columns:
    if not accept_unknowntypes:
      samples = samples[samples['sample_type'].isin(['Tumor'])]
  return samples


def assessAllSamples(sampless, refsamples, stype, rename={}, extract={}):
  extract.update(extract_defaults)
  rename.update(dup)
  sample_ids = []
  prevlen = len(sampless)
  sampless[extract['ref_type']] = stype

  # checking no duplicate in the buckets
  for k, val in sampless.iterrows():
    withsamesize = sampless[sampless[extract["size"]] == val[extract["size"]]]
    if len(withsamesize) > 1:
      if len(withsamesize[withsamesize[extract["ref_name"]] == val[extract["ref_name"]]]) < 2:
        raise ValueError('we have duplicate samples with different names!')
      else:
        for l, v in withsamesize.iloc[1:].iterrows():
          sampless = sampless.drop(l)
  print("we had " + str(prevlen - len(sampless)) + " duplicates in the release buckets")
  # check: currently, below lines prevent forcekeep from working on true duplicates
  # (aka same size  file). Need to think about how to bring back the forcekeep functionality
  sampless[extract['ref_arxspan_id']] = [rename[name] if name in rename else name for name in sampless[extract['ref_arxspan_id']]]
  names = []
  # need to keep track of whether we're adding more than one new entry for a given sample id
  subrefsamples = refsamples[refsamples[extract['ref_type']] == stype]
  for k, val in sampless.iterrows():
    val = val[extract['ref_arxspan_id']]
    names.append(val)
    sampless.loc[k, extract['version']] = len(subrefsamples[subrefsamples[extract['ref_arxspan_id']] == val]) + names.count(val)
  sampless[extract['version']] = sampless[extract['version']].astype(int)

  sampless[extract['patient_id']] = [val[extract['patient_id']] if
                                     len(refsamples[refsamples[extract['ref_arxspan_id']] == val[extract['ref_arxspan_id']]]) < 1 else
                                     refsamples[refsamples[extract['ref_arxspan_id']] == val[extract['ref_arxspan_id']]][extract['patient_id']].values[0]
                                     for i, val in sampless.iterrows()]

  return sampless


def setupPairsFromSamples(sampless, refsamples, extract):
  extract.update(extract_defaults)
  pairs = pd.DataFrame()
  normals = refsamples[refsamples[extract['primary_disease']] == 'normal']
  pairs['control_sample'] = ['nan' if len(normals[normals[extract['patient_id']] == val]) < 1 else normals[normals[extract['patient_id']] == val].index.tolist()[0] for val in sampless[extract['patient_id']]]
  pairs['case_sample'] = sampless.index.tolist()
  pairs['participant_id'] = sampless[extract['patient_id']].tolist()
  pairs['pair_id'] = [val['case_sample'] + '_' + val['control_sample'] for i, val in pairs.iterrows()]
  print('found ' + str(len(pairs['control_sample'].unique()) - 1) + ' matched normals')
  return pairs.set_index('pair_id')


def changeCellLineNameInNew(ref, new, datatype, dupdict, toupdate=['stripped_cell_line_name',
                                                                      'arxspan_id', "patient_id",
                                                                      "sex", "primary_disease",
                                                                      "cellosaurus_id", "age",
                                                                      "primary_site", "subtype",
                                                                      "subsubtype"]):
  """
  Args:
  -----
    new: change the cell line name in this dataframe
    dupdict: dict(tochange,newname)
    datatype: str for a ref with many datatype (to get the right version number)

  Returns:
  --------
    the updated dataframe
  """
  for k, v in dupdict.items():
    new.loc[new[new.arxspan_id == k].index, toupdate] = ref[ref.arxspan_id == v][toupdate].values[0]
    new.loc[new[new.arxspan_id == v].index, 'version'] = len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
  return new


def changeCellLineName(ref, datatype, dupdict, toupdate=["stripped_cell_line_name",
                                                         "participant_id",
                                                         "cellosaurus_id",
                                                         "sex",
                                                         "arxspan_id",
                                                         "matched_normal",
                                                         "age",
                                                         "primary_site",
                                                         "primary_disease",
                                                         "subtype",
                                                         "subsubtype",
                                                         "origin"]):
  """
  Args:
  -----
    dupdict: dict(tochange,newname)
    datatype: str for a ref with many datatype (to get the right version number)

  Returns:
  --------
    the updated dataframe
  """
  for k, v in dupdict.items():
    try:
      ref.loc[k, toupdate] = ref[ref.arxspan_id == v][toupdate].values[0]
      ref.loc[k, 'version'] = len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
    except IndexError:
      raise IndexError(str(v)+" not found in tracker")
  return ref

#####################
# CN Functions
#####################

def checkAmountOfSegments(segmentcn, thresh=850, samplecol="DepMap_ID"):
  """
  if there is too many segments, something might be wrong

  will compute the number of segments for each samples from a df of segments from RSEM

  Args:
  ----
    segmentcn: segment dataframe
    thresh: max ok amount
  """
  failed = []
  segmentcn = renameColumns(segmentcn)
  celllines = set(segmentcn[samplecol].tolist())
  amounts = []
  for cellline in celllines:
    val = segmentcn[segmentcn[samplecol] == cellline].shape[0]
    amounts.append(val)
    if val > thresh:
      failed.append(cellline)
      print(cellline, val)
  sns.kdeplot(amounts)
  return failed


def checkGeneChangeAccrossAll(genecn, thresh=0.2):
  """
  used to find poor quality genes in CN data

  compute given a df of gene x sample CN counts, how much change there is accross samples for
  a same gene and returns ones that are below the threshold
  Args:
  -----
    genecn: gene cn data frame
    thresh: threshold in logfold change accross all of them
  """
  return genecn.columns[genecn.var()<thresh].tolist()


def renameColumns(df):
  """
  rename some of the main columns names from RSEM, GATK.. to more readable column names

  Args:
  -----
    df: the df to rename

  Returns:
  ------
    df the renamed df
  """
  return(df.rename(columns={'Sample': 'DepMap_ID', 'CONTIG': 'Chromosome', 'START': 'Start',
                            'END': 'End', 'seqnames': 'Chromosome', 'start': 'Start', 'end': 'End'}))


def checkDifferencesWESWGS(segmentcn_wes, segmentcn_wgs, chromlist=CHROMLIST):
  """
  if the similarity between WES and WGS on same cell line is low there might be problems. Checks for that

  Args:
  ----
    segmentcn_wes: df the RSEM segment data from WES
    segmentcn_wgs: df the RSEM segment data from WGS
    chromlist: list off chromosome to compute on
  """

  # function for overlap between two segments 1D
  def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

  segmentcn_wes = renameColumns(segmentcn_wes)
  segmentcn_wgs = renameColumns(segmentcn_wgs)

  common = set(segmentcn_wgs["DepMap_ID"].tolist()) & set(segmentcn_wes["DepMap_ID"].tolist())
  print('anlaysing on ' + str(common) + " common cell lines")
  segmentcn_wes = segmentcn_wes.loc[segmentcn_wes["DepMap_ID"].isin(common)]
  segmentcn_wgs = segmentcn_wgs.loc[segmentcn_wgs["DepMap_ID"].isin(common)]
  if not segmentcn_wes['Chromosome'].iloc[0] in chromlist:
    print(segmentcn_wes['Chromosome'].iloc[0])
    print("assuming that the values are of type 1,2.. converting to chr1,chr2...")
    segmentcn_wes['Chromosome'] = 'chr' + segmentcn_wes['Chromosome'].astype(str)
  if not segmentcn_wgs['Chromosome'].iloc[0] in chromlist:
    print(segmentcn_wgs['Chromosome'].iloc[0])
    print("assuming that the values are of type 1,2.. converting to chr1,chr2...")
    segmentcn_wgs['Chromosome'] = 'chr' + segmentcn_wgs['Chromosome'].astype(str)

  counts = []
  for name in common:
    count = 0
    for chrom in chromlist:
      i = 0
      cnwes = segmentcn_wes[segmentcn_wes["DepMap_ID"] == name][
          segmentcn_wes["Chromosome"] == chrom][['Start', 'End', 'Segment_Mean']].values
      cnwgs = segmentcn_wgs[segmentcn_wgs["DepMap_ID"] == name][
          segmentcn_wgs["Chromosome"] == chrom][['Start', 'End', 'Segment_Mean']].values
      for val in cnwes:
        if not overlap(val[0], val[1], cnwgs[i][0], cnwgs[i][1]):
          count += 1
          if cnwgs[i][0] > val[1]:
            continue
          elif cnwgs[i][1] < val[0]:
            while cnwgs[i][1] < val[0]:
              i += 1
              count += 1
        else:
          if abs(val[2] - cnwgs[i][2]) > 0.2:
            count += 1
    if count > 500:
      print(name, count)


def manageGapsInSegments(segtocp, Chromosome='Chromosome', End="End", Start="Start", cyto=None):
  """
  """
  prevchr = ''
  prevend = 0
  count = 0
  l = []
  segments = segtocp.copy()
  le = len(segments)
  for k, val in segments.iterrows():
    h.showcount(count, le)
    count += 1
    if val[Chromosome] != prevchr:  # we changed chromosome
      # we extend the previous segment (last of the prev chrom) to.. way enough
      if len(l) > 0:
        l[-1][2] = 1000000000 if cyto is None else cyto[cyto['chrom']
                                                        == prevchr]['end'].values[-1]
      # we extend the first segment to 0
      l.append([val[Chromosome], 0, val[End]])
    else:
      if val[Start] > prevend + 1:  # we have a gap in the same chrom
        sizeofgap = val[Start] - prevend
        # we add to the previous one half of the gap
        l[-1][2] += int(sizeofgap /
                        2) if sizeofgap % 2 == 0 else int(sizeofgap / 2) + 1
        # the rest to the other
        l.append([val[Chromosome], val[Start] - int(sizeofgap / 2), val[End]])
      elif val[Start] < prevend:  # this should never happen
        raise ValueError("start comes after end")
      else:
        l.append([val[Chromosome], val[Start], val[End]])
    prevchr = val[Chromosome]
    prevend = val[End]
  # we extend the last one
  l[-1][2] = 1000000000 if cyto is None else cyto[cyto['chrom']
                                                  == prevchr]['end'].values[-1]
  segments[[Chromosome, Start, End]] = l
  return segments


def toGeneMatrix(segments, gene_mapping, style='weighted', missingchrom=['Y']):
  """
  makes gene matrix from segment level copy number

  Args:
  ----
    style: one of "weighted","mean","closest"
  """
  samples = list(set(segments.DepMap_ID))
  data = np.zeros((len(samples), len(gene_mapping)))
  for i, sample in enumerate(samples):
    segs = segments[segments.DepMap_ID == sample][[
        'Chromosome', 'Start', 'End', "Segment_Mean"]].values
    hasmissing = set(missingchrom) - set(segs[:, 0])
    j = 0
    h.showcount(i, len(samples))
    for k, gene in enumerate(gene_mapping[['Chromosome', 'start', 'end']].values):
        if gene[0] in hasmissing:
          data[i, k] = np.nan
          continue
        while gene[0] != segs[j][0] or gene[1] >= segs[j][2]:
          #print("went beyong",gene, segs[j])
          j += 1
        # some genes are within other genes, we need to go back in the list of segment in that case
        while gene[1] < segs[j][1]:
          j -= 1
          #print("decrease gene",gene)
        # we are entirely within the segment
        c = 1
        if gene[2] <= segs[j][2]:
          data[i, k] = segs[j][3]
        else:
          # how much of the gene is covered by the segment
          coef = (segs[j][2] - gene[1]) / (gene[2] - gene[1])
          # print('coef',coef)
          val = segs[j][3] * coef if style == "weighted" else segs[j][3]
          end = segs[j][2]
          # until the end of a segments goes beyon the end of the gene (say if we have X segments within the gene)
          while end < gene[2]:
            # pdb.set_trace()
            j += 1
            c += 1
            nextend = segs[j][2] if segs[j][2] < gene[2] else gene[2]
            # here, end (of prevsegment) is the next segment's start
            ncoef = (nextend - end) / (gene[2] - gene[1])
            # print('multi',gene, ncoef)
            if style == "closest":
              if ncoef > coef:
                val = segs[j][3]
              else:
                # we switch it back (see line 894)
                ncoef = coef
            else:
              val += segs[j][3] * ncoef if style == "weighted" else segs[j][3]
            end = segs[j][2]
            coef = ncoef
          data[i, k] = val if style == "weighted" else val / c
  return pd.DataFrame(data=data, index=samples, columns=[i['symbol'] + ' (' + str(i['ensembl_id']) + ')' for _, i in gene_mapping.iterrows()])


def plotCNchanges(newgenecn, prevgenecn, newsegments, prevsegments, depmap_id="DepMap_ID", source="Source", prevname='prev', newname="new"):
  """
  makes a Javad Plot on the gene copy number dataset

  Args:
  -----

  """
  grouped = pd.concat(
      [prevgenecn.stack(), newgenecn.stack()], axis=1)
  grouped.columns = [prevname, newname]

  grouped.reset_index(inplace=True)
  grouped.rename(
    columns={'level_0': depmap_id, 'level_1': 'gene'}, inplace=True)
  sources = pd.merge(prevsegments[[depmap_id, source]].drop_duplicates(),
          newsegments[[depmap_id, source]].drop_duplicates(),
          on=depmap_id, suffixes=['_'+prevname, '_'+newname])

  sources['source_change'] = sources.apply(lambda x: '{:s} -> {:s}'.format(x[source+'_'+prevname], x[source+'_'+newname]), axis=1)
  sources['source_has_changed'] = (sources[source+'_'+prevname] != sources[source+"_"+newname])
  grouped = pd.merge(grouped, sources, on=depmap_id)
  plt.figure(figsize=(20,10))
  sns.scatterplot(data=grouped.sample(1000000, random_state=0), x=prevname, y=newname,
                hue='source_change', style='source_has_changed', alpha=0.5, cmap='Tab20')

#####################
# Fusion Functions
#####################

def addToMainFusion(input_filenames, main_filename):
  """
  Given a tsv fusion files from RSEM algorithm, merge it to a tsv set of fusion data

  Args:
  ----
    input_filenames: a set of filepath to input the files should be c|tsv from Terra fusion pipeline
    main_filename: a filepath to input the files should be c|tsv from Terra aggregation pipeline
  """
  maindata = pd.read_csv(main_filename, sep='\t')
  if '.' in maindata["DepMap_ID"][0]:
    maindata["DepMap_ID"] = [i[0] for i in maindata["DepMap_ID"].str.split('.').tolist()]
  samples = set(maindata["DepMap_ID"].tolist())
  with open(main_filename, 'a') as f:
    for input_filename in input_filenames:
      df = pd.read_csv(input_filename, sep='\t')
      input_filename = input_filename.split('/')[-1].split('.')[0]
      if input_filename in samples:
        print(input_filename + " is Already in main fusions")
      df['DepMap_ID'] = pd.Series([input_filename] * len(df.index.tolist()), index=df.index)
      cols = df.columns.tolist()
      cols = cols[-1:] + cols[: -1]
      df = df[cols]
      df.to_csv(f, header=False, sep='\t', index=False)


def renameFusionGene(a):
    return [str(i.split('^')).replace(', ', ' (').replace("'", "")[1:-1]+')' for i in a]


def filterFusions(fusions, maxfreq=0.1, minffpm=0.05, red_herring=['GTEx_recurrent', 'DGD_PARALOGS', 'HGNC_GENEFAM', 'Greger_Normal', 'Babiceanu_Normal', 'ConjoinG', 'NEIGHBORS']):
  fusions = fusions.copy()
  # remove recurrent
  fusions = fusions[fusions['CCLE_count'] <
                    len(set(fusions['DepMap_ID']))*maxfreq]
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


#####################
# Expression Functions
#####################

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
    name = input_filename['rsem_genes'].split('/')[-1].split('.')[0].split('_')[-1]
    rsem_genes = pd.read_csv('temp/' + input_filename['rsem_genes'].split('/')[-1], sep='\t')
    rsem_transcripts = pd.read_csv('temp/' + input_filename['rsem_isoforms'].split('/')[-1], sep='\t')
    genes_count[name] = pd.Series(rsem_genes['expected_count'], index=rsem_genes.index)
    transcripts_tpm[name] = pd.Series(rsem_transcripts['TPM'], index=rsem_transcripts.index)
    genes_tpm[name] = pd.Series(rsem_genes['TPM'], index=rsem_genes.index)

  genes_count.to_csv('temp/' + main_filename['rsem_genes_expected_count'].split('/')[-1], sep='\t',
                     index=False, index_label=False, compression='gzip')
  transcripts_tpm.to_csv('temp/' + main_filename['rsem_transcripts_tpm'].split('/')[-1], sep='\t',
                         index=False, index_label=False, compression='gzip')
  genes_tpm.to_csv('temp/' + main_filename['rsem_genes_tpm'].split('/')[-1], sep='\t',
                   index=False, index_label=False, compression='gzip')


def ExtractStarQualityInfo(samplesetname, workspace, release='temp'):
  """
  put all of the Star Quality results from Terra Star Workflow into one txt file

  Args:
  -----
    samplesetname: the sampleset name for which to grab the samples processed by star.
    wm: the terra workspace
    release: the name of the folder where it will be stored
  """
  a = dm.WorkspaceManager(workspace).get_samples().loc[
      dm.WorkspaceManager(workspace).get_sample_sets().loc[samplesetname].samples].star_logs
  for i, sample in enumerate(a):
    if sample is None:
      print("no log file found for: " + a.index[i])
    for log in sample:
      if 'final.out' in log:
        print("copying " + a.index[i])
        os.system('gsutil cp ' + log + ' temp/')
  os.system("cat data/" + release + "/*.Log.final.out > temp/" + samplesetname + ".txt")
  os.system("rm data/" + release + "/*.Log.final.out")


def findMissAnnotatedReplicates(repprofiles, goodprofile, names, exactMatch=True):
  """
  from a new rnaseq profile on replicate level and a good rnaseq profile on sample level

  will if some replicates are missanotated based on correlation.

  Returns:
  -------
      notindataset: list[str] replicates not in the good dataset
      missannotated: dict(str: tuple(str,str)).  dict containing replicates that are missanotated: for each, gives a tuple (old annotation, right annotation)
  """
  notindataset = []
  missannotated = {}
  unmatched = {}
  if exactMatch:
    res = findClosestMatching(repprofiles, goodprofile)
    for val in repprofiles.index.tolist():
        if val not in res:
            notindataset.append(val)
        elif val not in names:
            unmatched.update({val: res[val]})
        elif res[val] != names[val]:
            missannotated.update({val: (names[val], res[val])})
    return notindataset, missannotated, unmatched
  else:
    corr, closest = findClosestMatching(
        repprofiles, goodprofile, returncorr=True)
    for k, v in corr.iterrows():
      print(k, v.mean())
      try:
          if v[names[k]] < 0.75:
              print(v[[closest[k], names[k]]])
      except:
          a = np.argsort(v.values)[-5:]
          if v.values[a[-1]] > 0.8:
              print(names[k],
                    corr.columns[a], v.values[a])


def findClosestMatching(repprofiles, goodprofile, closest=False, returncorr=False):
  """
  will find what replicate matches best what known profile
  """
  match = {}
  a = set(repprofiles.columns) & set(goodprofile.columns)
  ind = goodprofile.index.tolist()
  corr = []
  for i, (k, v) in enumerate(repprofiles[a].iterrows()):
      h.showcount(i, len(repprofiles))
      res = np.array([np.corrcoef(v, w)[0, 1]
                      for _, w in goodprofile[a].iterrows()])
      if max(res) == 1 or closest:
          match[k] = ind[np.argmax(res)]
      if returncorr:
        corr.append(res)
  if returncorr:
    corr = pd.DataFrame(data=corr, index=repprofiles.index.tolist(
    ), columns=goodprofile.index.tolist())
    return corr, match
  else:
    return match


def rnaseqcorrelation(cn, rna, ax=None, name=None):
  """
  correlates the copy number to the rnaseq in ccle and shows the plot
  """
  a = set(cn.columns) & set(rna.columns)
  ind = set(cn.index) & set(rna.index)
  re = rna.loc[ind]
  ce = cn.loc[ind]
  print(len(ind), len(a))
  corr = np.array([pearsonr(ce[j], re[j])[0] for j in a])
  #corr = pd.DataFrame(data=corr, columns=[name if name is not None else "data"])
  print(np.mean(corr), len(corr))
  sns.kdeplot(corr, ax=ax) if ax is not None else sns.kdeplot(corr)

#####################
# Mutations Functions
#####################

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


def filterAllelicFraction(maf, loc=['CGA_WES_AC'], sep=':', frac=0.1):
    muts = np.zeros((len(maf), 2))
    for val in loc:
        muts += np.array([[v[0], 0] if 'NA' in v else v for v in maf[val].fillna(
            '0'+sep+'0').astype(str).str.split(sep).tolist()]).astype(int)
    muts = muts[:, 0]/(muts[:, 0]+muts[:, 1])
    return maf[muts >= frac]


def filterCoverage(maf, loc=['CGA_WES_AC'], sep=':', cov=4, altloc=0):
    muts = np.zeros((len(maf), 2))
    for val in loc:
        muts += np.array([[v[0], 0] if 'NA' in v else v for v in maf[val].fillna(
            '0'+sep+'0').astype(str).str.split(sep).tolist()]).astype(int)
    return maf[muts[:, altloc] >= cov]


def annotate_likely_immortalized(maf, sample_col="DepMap_ID", genome_change_col="Genome_Change", TCGAlocs=['TCGAhsCnt',
                                                                                                           'COSMIChsCnt'], max_recurrence=0.05, min_tcga_true_cancer=5):
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
    maf['NCBI_Build'] = NCBI_Build
    maf['Strand'] = Strand
    return maf


def mafToMat(maf, boolify=False, freqcol='tumor_f', samplesCol="DepMap_ID", mutNameCol="Hugo_Symbol"):
    maf = maf.sort_values(by=mutNameCol)
    samples = set(maf[samplesCol])
    mut = pd.DataFrame(data=np.zeros((len(set(maf[mutNameCol])), 1)), columns=[
                       'fake'], index=set(maf[mutNameCol])).astype(float)
    for i, val in enumerate(samples):
        h.showcount(i, len(samples))
        mut = mut.join(maf[maf[samplesCol] == val].drop_duplicates(
            mutNameCol).set_index(mutNameCol)[freqcol].rename(val))
    return mut.fillna(0).astype(bool if boolify else float).drop(columns=['fake'])


def mergeAnnotations(firstmaf, additionalmaf, Genome_Change="Genome_Change", Start_position="Start_position", Chromosome="Chromosome", samplename="DepMap_ID", useSecondForConflict=True, dry_run=False):
    mutations = firstmaf.copy()
    mutations['ind'] = mutations[samplename]+"_"+mutations[Genome_Change]
    mutations['loci'] = mutations[samplename] + "_" + \
        mutations[Chromosome] + "_" + mutations[Start_position].astype(str)
    additionalmaf['ind'] = additionalmaf[samplename] + \
        "_"+additionalmaf[Genome_Change]
    additionalmaf['loci'] = additionalmaf[samplename] + "_" + \
        additionalmaf[Chromosome] + "_" + \
        additionalmaf[Start_position].astype(str)
    inboth = set(additionalmaf['loci']) & set(mutations['loci'])
    notineach = set(additionalmaf['ind']) ^ set(mutations['ind'])
    submut = mutations[mutations.loci.isin(
        inboth) & mutations.ind.isin(notineach)]
    subother = additionalmaf[additionalmaf.loci.isin(
        inboth) & additionalmaf.ind.isin(notineach)]
    issues = None
    if len(submut) > 0:
        print("found " + str(len(submut)) + " nonmatching mutations")
        issues = np.vstack([submut.sort_values(by='loci')[
                           Genome_Change].values, subother.sort_values(by='loci')[Genome_Change].values]).T
        if dry_run:
            print(issues)
    if not dry_run:
        if issues is not None:
            if useSecondForConflict:
                mutations = mutations[~mutations.ind.isin(set(submut.ind))]
            else:
                additionalmaf = additionalmaf[~additionalmaf.ind.isin(
                    set(subother.ind))]
            mutations = mutations.append(additionalmaf[additionalmaf['ind'].isin(
                set(additionalmaf['ind']) - set(mutations['ind']))])
        return mutations.drop(columns=['loci', 'ind']).sort_values(by=[samplename, Chromosome, Start_position])
    else:
        return issues

#####################
# Other Helpers
#####################

def AddToVirtual(virtualname, folderfrom = None, files):
  """
  will add some files from a taiga folder to a taiga virtual dataset folder and preserve the previous files

  Args:
  ----
    virtualname: the taiga virtual dataset name
    folderfrom: the taiga folder wher the files are
    files: a list(tuples(newfilename,prevfilename)) can be from the folder
  """
  file_dict = {}
  assert type(files[0]) is tuple
  if folderfrom is not None:
    versiona = max([int(i['name']) for i in tc.get_dataset_metadata(folderfrom)['versions']])
  versionb = max([int(i['name']) for i in tc.get_dataset_metadata(virtualname)['versions']])
  keep = [(i['name'], i['underlying_file_id']) for i in tc.get_dataset_metadata(virtualname, version=versionb)
          ['datasetVersion']['datafiles'] if 'underlying_file_id' in i]

  for i, val in enumerate(files):
    if "/" in val[1]:
      print("assuming "+val[1]+" to be a local file")
      file_dict.update({val[0]:val[1]})
    else:
      files[i] = (val[0], folderfrom + '.' + str(versiona) + '/' + val[1])
  print(files)
  tc.update_dataset(dataset_permaname=virtualname, add_taiga_ids=files, upload_file_path_dict=file_dict, add_all_existing_files=True)

#####################
# DB Functions
#####################

def compareToCuratedGS(url, sample, samplesetname, sample_id='DepMap ID', clientsecret='~/.client_secret.json',
                       storagepath='~/.storage.json', colname='CN New to internal', value='no data yet'):
  """
  from a google spreadsheet, will check that we have all of the samples we should have in our sample
  set name (will parse NAME_additional for sample_id)

  Args:
  -----
    url: str the url of the gsheet
    sample: list(str) the samples to check
    samplesetname: str the name of the sampleset in the googlesheet
    sample_id: str the name of the sample_id column in the google sheet
    clientsecret: str path to your secret google api account file
    storagepath: str path to your secret google api storage file
    colname: str if we need not to include some rows from the spreadsheet that have the value value
    value: str the value for which not to include the rows

  @gmiller
  """
  sheets = Sheets.from_files(clientsecret, storagepath)
  # Cell Line Profiling Status google sheet
  gsheet = sheets.get(url).sheets[0].to_frame()
  gsheet.index = gsheet[sample_id]
  new_cn = gsheet[gsheet[colname] == samplesetname + 'tent']
  if colname and value:
    data_not_ready_cn = gsheet[gsheet[colname] == value]
    print(data_not_ready_cn)
  # these are the "new" samples discovered by our function, createDatasetsFromNewCellLines
  sample_ids = [id.split('_')[0] for id in sample]
  print("We found data for " + str(len(sorted(sample))) + " samples.\n")

  print("Sanity check: Since we have the tacked on number, we should only have 1 each per sample ID:\n")
  Counter(sample)

  in_sheet_not_found = set(new_cn.index.tolist()) - set(sample_ids)
  if len(in_sheet_not_found) > 0:
    print("We have not found " + str(len(in_sheet_not_found)) + " of the samples we're supposed to \
      have this release:\n" + str(sorted(list(in_sheet_not_found))))
  else:
    print("We aren't missing any samples that we're supposed to have this release!")


def removeOlderVersions(names, refsamples, arxspan_id="arxspan_id", version="version"):
  """
  Given a dataframe containing ids, versions, sample_ids and you dataset df indexed by the same ids, will set it to your sample_ids using the latest version available for each sample

  Args:
  -----
    refsamples: df[id, version, arxspan_id,...] the reference metadata
    names: list[id] only do it on this set of samples
    arxspan_id: the name of the id field
    version: the name of the version field

  Returns:
  --------
    the subsetted dataframe

  """
  lennames = len(names)
  res = {}
  refsamples = refsamples[refsamples.index.isin(names)]
  if lennames > len(refsamples):
    ipdb.set_trace()
    raise ValueError('we had some ids in our dataset not registered in this refsample dataframe')
  for arxspan in set(refsamples[arxspan_id]):
    allv = refsamples[refsamples[arxspan_id] == arxspan]
    for k, val in allv.iterrows():
      if val[version] == max(allv.version.values):
        res[k] = arxspan
        break
  print("removed " + str(lennames - len(res)) + " duplicate samples")
  # remove all the reference metadata columns except the arxspan ID
  return res


def getQC(workspace, only=[], qcname=[], match=""):
  """
  Will get from a workspace, the QC data for each samples

  Args:
  -----
    workspace: the workspace name
    only: do it only for this set of samples
    qcname: col name where the QC is in the workspace samples
    match: for example'.Log.final.out' get only that QC if you have a list of QCs in you qcname col

  Returns:
  --------
    a dict(sample_id:list[QC_filepaths])
  """
  if type(qcname) is str:
    qcname = [qcname]
  res = {}
  wm = dm.WorkspaceManager(workspace)
  sam = wm.get_samples()
  if len(only) > 0:
    sam = sam[sam.index.isin(only)]
  for k, val in sam.iterrows():
    res[k] = []
    for i in val[qcname]:
      if type(i) is list:
        if match:
          res[k].extend([e for e in i if match in e])
        else:
          res[k].extend(i)
      else:
        res[k].append(i)
  return res


def updateSamplesSelectedForRelease(refsamples, releaseName, samples):
  refsamples.loc[samples, releaseName] = '1'
  return refsamples


def copyToWorkspace(workspaceID, tracker, columns=["arxspan_id", "version", "sm_id", "datatype", "size",
                                                   "ccle_name", "stripped_cell_line_name", "patient_id", "cellosaurus_id",
"bam_public_sra_path", "internal_bam_filepath", "internal_bai_filepath",
"parent_cell_line", "sex", "matched_normal", "age", "primary_site",
"primary_disease", "subtype", "subsubtype", "origin", "mediatype",
"condition", "sequencing_type", "baits", "source", "legacy_bam_filepath", "legacy_bai_filepath"], rename={'participant_id':'patient_id'}, deleteUnmatched=False):
  """
  will used the current sample tracker to update sample annotation in the workspace

  it can remove samples that are not in the tracker.
  """
  wm = dm.WorkspaceManager(workspaceID).disable_hound()
  sam = wm.get_samples()
  track = tracker[tracker.index.isin(sam.index)].rename(columns=rename)
  miss = set(columns) - set(sam.columns)
  print('found these columns to be missing in workspace: '+str(miss))
  if len(track)==0:
      raise ValueError('wrong tracker or index non matching')
  unmatched = set(sam.index) - (set(tracker.index) | set(['nan']))
  print("found these to be unmatched: "+str(unmatched))
  if deleteUnmatched and len(unmatched)>0:
      terra.removeSamples(workspaceID, unmatched)
  wm.update_sample_attributes(track[columns])


def updatePairs(workspaceID, tracker, removeDataFiles=True, ):
  """
  looks at the current sample tracker and updates the pairs in Terra

  It will add and remove them based on what information of match normal is available in the sample tracker. if an update happens it will remove the data files for the row.
  """

def cleanVersions(tracker, samplecol='arxspan_id', dryrun=False, datatypecol='datatype', versioncol="version"):
  """
  cleans the versions of a sample tracker:

  checks that we get 1,2,3 instead of 2,4,5 when samples are renamed or removed
  """
  tracker = tracker.copy()
  tracker['samval']=tracker[samplecol]+tracker[datatypecol]
  for v in set(tracker['samval']):
      sams = tracker[tracker['samval']==v]
      vs = sams[versioncol].tolist()
      if max(vs)==len(vs):
          continue
      print("found issue")
      if dryrun:
        continue
      vs.sort()
      rn = {}
      for i,v in enumerate(vs):
          rn.update({v:i+1})
      for k, val in sams.iterrows():
          tracker.loc[k, versioncol] = rn[val[versioncol]]
  tracker = tracker.drop(columns='samval')
  return tracker


def setRightName(tracker, name='stripped_cell_line_name', signs=['-','_','.',' ']):
  """
  cell line name needd to be all uppercase and not contain special signs
  """
  new = []
  for val in tracker[name]:
      for s in signs:
          val=val.replace(s,'')
      new.append(val.upper())
  tracker[name]=new
  return tracker


def findLikelyDup(tracker, signs=['-', '_', '.', ' '], name='stripped_cell_line_name', arxspid='arxspan_id', looksub=True):
    """
    find cell lines that are likely to be duplicates

    will return a list of likly duplicate names as tuples (rh13, RH-13), a list of associated arxspan ids as tuples as well,
    and a list of arxspan ids that have multiple names associated
    """
    names = set(tracker[name])
    simi = []
    arxsp = []
    issues = {}
    for i, name1 in enumerate(names):
        h.showcount(i, len(names))
        n1 = name1
        for s in signs:
            name1 = name1.replace(s, '')
        name1 = name1.upper()
        for name2 in names-set([n1]):
            n2 = name2
            for s in signs:
                name2 = name2.replace(s, '')
            name2 = name2.upper()
            if name1 == name2:
                if (looksub and (name1 in name2 or name2 in name1) and abs(len(name1)-len(name2)) < 2) or not looksub:
                    if (n1, n2) not in simi and (n2, n1) not in simi:
                        simi.append((n1, n2))
                        arxsp.append(
                            (tracker[tracker[name] == n1][arxspid][0], tracker[tracker[name] == n2][arxspid][0]))
    for val in set(tracker[name]):
        v = set(tracker[tracker[name] == val][arxspid])
        if len(v) > 1:
            issues.update({val: v})
    return simi, arxsp, issues


def merge(tracker, new, old, arxspid, cols):
  """
  given a tracker a a new and old arxspan id, will merge the two cells lines in the tracker
  """
  #loc = tracker[tracker[arxspid]==old].index


def resolveIssues(tracker, issus, arxspid, cols):
  """
  given a dict of names: [arxp ids] will try to find back the right name for the right
  arxspan id by looking at their rfequncy of occurance along the tracker

  if we have rh12: [ACH-00001,ACH-0002]
  and rh12 is associated 1 time with ach-00002 and 3 with ach-00001
  and rh13 is assocated 2 time with ach-00002, then it associates:
  ach-00001 : rh12
  ach-00002 : rh13
  """
  #for val in issus:


def retrieveFromCellLineName(noarxspan, ccle_refsamples, datatype, extract={},
stripped_cell_line_name="stripped_cell_line_name", arxspan_id="arxspan_id", depmappvlink="https://docs.google.com/spreadsheets/d/1uqCOos-T9EMQU7y2ZUw4Nm84opU5fIT1y7jet1vnScE"):
    # find back from cell line name in ccle ref samples
    noarxspan.arxspan_id = [ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i].arxspan_id[0] if i in ccle_refsamples[stripped_cell_line_name].tolist() else 0 for i in noarxspan[arxspan_id]]
    a = [ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i][arxspan_id][0]
         if i in ccle_refsamples[stripped_cell_line_name].tolist() else 0 for i in noarxspan[stripped_cell_line_name]]
    noarxspan[arxspan_id] = [i if i != 0 else a[e]
                             for e, i in enumerate(noarxspan.arxspan_id)]

    # get depmap pv
    depmap_pv = sheets.get(depmappvlink).sheets[0].to_frame(header=2)
    depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index)

    # find back from depmapPV
    signs = ['-', '_', '.', ' ']
    for k, val in noarxspan[noarxspan[arxspan_id] == 0].iterrows():
        val = val[stripped_cell_line_name].upper()
        for s in signs:
            val = val.replace(s, '')
        a = depmap_pv[depmap_pv['CCLE_name'].str.contains(
            val) | depmap_pv['Stripped Cell Line Name'].str.contains(val) | depmap_pv['Aliases'].str.contains(val)]
        if len(a) > 0:
            noarxspan.loc[k, arxspan_id] = a['DepMap_ID'].values[0]
    noarxspan[arxspan_id] = noarxspan[arxspan_id].astype(str)
    new_noarxspan= resolveFromWorkspace(noarxspan[noarxspan[arxspan_id].str.contains('ACH-')], refsamples=ccle_refsamples[ccle_refsamples['datatype'] == datatype], match=[
                                     'ACH', 'CDS'], participantslicepos=10, accept_unknowntypes=True, extract=extract)
    return pd.concat([new_noarxspan, noarxspan[~noarxspan[arxspan_id].str.contains('ACH-')]])


def updateFromTracker(samples, ccle_refsamples, arxspan_id='arxspan_id', participant_id='participant_id', toupdate={"sex":[],
"primary_disease":[],
"cellosaurus_id":[],
"age":[],
"primary_site":[],
"subtype":[],
"subsubtype":[],
"origin":[],
"parent_cell_line":[],
"matched_normal":[],
"comments":[],
"mediatype":[],
"condition":[],
'stripped_cell_line_name':[],
"participant_id":[]}):
    # If I have a previous samples I can update unknown data directly
    index = []
    notfound = []
    for k, val in samples.iterrows():
        dat = ccle_refsamples[ccle_refsamples[arxspan_id] == val[arxspan_id]]
        if len(dat) > 0:
            index.append(k)
            for k, v in toupdate.items():
                toupdate[k].append(dat[k].tolist()[0])
        else:
            notfound.append(k)
    # doing so..
    for k, v in toupdate.items():
        samples.loc[index, k] = v
    len(samples.loc[notfound][participant_id]
        ), samples.loc[notfound][participant_id].tolist()
    return samples, notfound
