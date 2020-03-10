# Jérémie Kalfon
# for BroadInsitute
# in 2019


####
#
# HELPER FUNC  ######################################
#
#
import pdb
import pandas as pd
import numpy as np
import dalmatian as dm
from taigapy import TaigaClient
tc = TaigaClient()
import os
import sys
print("you need to have installed JKBio in the same folder as ccle_processing")
from JKBio import TerraFunction as terra
from JKBio import GCPFunction as gcp
from collections import Counter
from gsheets import Sheets

CHROMLIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX']


extract_defaults = {
    'name': 'sample_alias',
    'bai': 'crai_or_bai_path',
    'bam': 'cram_or_bam_path',
    'source': 'source',
    'id': 'individual_alias',
    'participant': 'individual_alias',
    'bait_set': 'bait_set',
    'hs_hs_library_size': 'hs_hs_library_size',
    'hs_het_snp_sensitivity': 'hs_het_snp_sensitivity',
    'hs_mean_bait_coverage': 'hs_mean_bait_coverage',
    'hs_mean_target_coverage': 'hs_mean_target_coverage',
    'hs_on_target_bases': 'hs_on_target_bases',
    'total_reads': 'total_reads',
    'release_date': 'release_date',
    'mean_depth': 'mean_depth'
}


def createDatasetWithNewCellLines(wto, samplesetname,
                                  wmfroms, sources,
                                  forcekeep=[], addonly=[], match='ACH', other_to_add=[], extract={},
                                  extract_defaults=extract_defaults,
                                  participantslicepos=10, accept_unknowntypes=False,
                                  gsfolderto=None, dry_run=False, rename=dict()):
  """
  As GP almost always upload their data to a data workspace. we have to merge it to our processing workspace

  Will merge samples from a set of data workspaces to a processing workspace on Terra. Will only 
  get a subset of the metadata and rename it.
  Will find out the duplicates based on the file size. 
  Can also upload the bam files to a google storage bucket

  Args:
  -----
    check: will need to change slicepos based on match;
    may have to take diff approach (regex?) for CCLF samples since CCLF sample IDs are not of consistent length
    wto: dalmatian.workspacemanager the workspace where you want to create the tsvs
    samplesetname: String the name of the sample set created by this new set of samples updated
    wfrom1: dalmatian.workspacemanager the workspace where the Samples to add are stored
    source1: the corresponding source name
    Can add as much as one wants
    dry_run: whether to perform a dry run without actually uploading anything to Terra
    forcekeep: list of sample id that you want to keep even if already in the previous workspace (will
    cause an overwrite) NOTE: currently this functionality is broken due to the implementation of checking for duplicates
    addonly: list of sample id that you only want to add
    match: list of substring(s) that has to be matched against the id of the samples to add them
    other_to_add:
    extract: if you want to specify what values should refer to which column names
      dict{    'name':
      'bai':
      'bam':
      'source':
      'id':
      ...}
    extract_defaults: the full default dict to specificy what values should refer to which column names


  Returns:
  -------
    sample_ids: list(str) the SM-id of the samples that were updated
    refsamples: all of the samples in the data of the workspaces
    CCLE_name: CCLE names of new samples uploaded


  Raise:
  -----
    Exception: when no new samples in this matrix
  """
  wto = dm.WorkspaceManager(wto)
  refsamples = wto.get_samples()
  refids_full = refsamples['participant'].tolist()
  # do NOT make refids a set; we use the num of occurences as way to determine what number to add to the sample id
  # filter refids to only include those that include the strings in the 'match' argument
  refids = []
  if type(match) is str:
    match = [match]
  for match_substring in match:
    refids += [val[val.index(match_substring):] for val in refids_full if match_substring in val]

  extract.update(extract_defaults)

  print("Getting sample infos...")
  if type(sources) is str:
    sources = [sources]
  if type(wmfroms) is str:
    wmfroms = [wmfroms]
  sampless = pd.DataFrame()
  names = [a.split('_')[0] for a in refsamples.index if type(a) is str]

  ref_size = {gcp.extractSize(val): gcp.extractPath(val) for val in gcp.lsFiles(
      [i for i in refsamples[extract["ref_bams"]] if type(i) is str and str(i) != 'NA'], "-al", 200)}

  for source, wmfrom in zip(sources, wmfroms):
    broken_bams = []
    wmfrom = dm.WorkspaceManager(wmfrom)
    samples = wmfrom.get_samples().replace(np.nan, '', regex=True).reset_index()
    # keep samples that contain the match requirement (e.g. ACH for DepMap IDs)
    samples = samples[samples[extract['id']].str.contains('|'.join(match))]
    print("\nThe shape of the sample tsv from " + str(wmfrom) + ": " + str(samples.shape))

    # remove true duplicates from consideration
    print("Identifying any true duplicates by checking file sizes (this runs for each data source)...")
    print("This step can take a while as we need to use gsutil to check the size of each potential duplicate...")
    dups_to_remove = []
    # check for broken bam files; if broken, then remove from consideration
    # need to check for broken filepaths before checking if the sample is in Terra so that we don't
    # add a broken file path for a new participant
    foundfiles = gcp.lsFiles(samples[extract['bam']])
    broken_bams = set(samples[extract['bam']]) - set(foundfiles)
    tolookfor = [val[extract['bam']] for _, val in samples.iterrows() if val[extract['id']][0:participantslicepos] in set(refids)]
    sample_size = {gcp.extractSize(val): gcp.extractPath(val) for val in gcp.lsFiles(tolookfor, "-al")}

    dups_to_remove = [sample_size[a] for a in set(ref_size.keys()) & set(sample_size.keys())]

    # remove the duplicates from consideration
    print("Len of samples before removal: " + str(len(samples)))
    print('These bam file path do not exist: ' + str(broken_bams))
    print("Dups from " + str(wmfrom) + " has len " + str(len(dups_to_remove)) + ":\n " + str(dups_to_remove))

    # remove the samples with broken bam filepaths from consideration
    samples = samples[(~samples[extract['bam']].isin(dups_to_remove)) & (~samples[extract['bam']].isin(broken_bams))]
    samples.reset_index(drop=True, inplace=True)
    # add number to sample ID so runs of same participant have unique sample ID
    new_samples = samples[extract['id']].str.slice(0, participantslicepos)
    # check: currently, below lines prevent forcekeep from working on true duplicates
    # (aka same size  file). Need to think about how to bring back the forcekeep functionality
    print("Len of samples after removal: " + str(len(samples)))

    # need to keep track of whether we're adding more than one new entry for a given sample id
    for i, val in enumerate(new_samples):
      num_in_workspace = names.count(val)
      samples[extract['id']][i] = val + '_' + str(num_in_workspace + 1)
      samples["version"] = int(num_in_workspace + 1)
      names.append(val)
    if len(addonly) > 0:
      samples = samples[samples[extract['id']].isin(addonly)]

    samples[extract['source']] = source
    if 'sample_type' in samples.columns:
      if not accept_unknowntypes:
        samples = samples[samples['sample_type'].isin(['Tumor'])]
    sampless = pd.concat([sampless, samples], sort=False)

  # ## currently, we don't do anything with notfound. Therefore addonly functionality is broken.
  # notfound = set(addonly) - set(sampless.index.tolist())
  # if len(notfound) > 0:
  #   # samples not found in the wto workspace so will be adding them now
  #   print('we did not find:' + str(notfound))
  sample_ids = []
  sampless['version'] = sampless['version'].astype(int)
  # to do the download to the new dataspace
  if not dry_run:
    pdb.set_trace()
    if gsfolderto is not None:
      for i, val in sampless.iterrows():
        if not gcp.exists(gsfolderto + 'v' + str(val["version"]) + '/' + val[extract['bam']].split('/')[-1]):
          cmd = 'gsutil cp ' + val[extract['bam']] + ' ' + val[extract['bai']] + ' ' + gsfolderto + 'v' + str(val["version"]) + '/'
          res = os.system(cmd)
          if res != 0:
            print("error in command '" + cmd + "'")
      sampless[extract['bam']] = [gsfolderto + 'v' + str(a["version"]) + '/' + a[extract['bam']].split('/')[-1] for _, a in sampless.iterrows()]
      sampless[extract['bai']] = [gsfolderto + 'v' + str(a["version"]) + '/' + a[extract['bai']].split('/')[-1] for _, a in sampless.iterrows()]
    print(sampless)
  if len(sampless) == 0:
    raise Exception("no new samples in this matrix")
  for ind, val in sampless.iterrows():
    name = val[extract['id']]
    if name in rename.keys():
      name = rename[name]
    refsamples = refsamples.append(pd.Series(
        {
            "CCLE_name": val[extract['name']],
            extract['ref_bais']: val[extract['bai']],
            extract['ref_bams']: val[extract['bam']],
            "Source": val[extract['source']],
            "participant": val[extract['participant']][:participantslicepos],

        }, name=name))
    sample_ids.append(name)

  if not dry_run:
    print("uploading new samples")
    wto.upload_samples(refsamples)
    print("creating a sample set")
    wto.update_sample_set(sample_set_id=samplesetname, sample_ids=sample_ids)
  return sample_ids, refsamples, refsamples[refsamples.index.isin(sample_ids)].CCLE_name.tolist()

#####################
# VALIDATION
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
  segmentcn = renameColumns(segmentcn)
  celllines = set(segmentcn[samplecol].tolist())
  for cellline in celllines:
    if segmentcn[segmentcn[samplecol] == cellline].shape[0] > thresh:
      print(cellline, segmentcn[segmentcn[samplecol] == cellline].shape[0])


def checkGeneChangeAccrossAll(genecn, thresh=1.5):
  """
  used to find poor quality genes in CN data

  compute given a df of gene x sample CN counts, how much change there is accross samples for
  a same gene and returns ones that are below the threshold
  Args:
  -----
    genecn: gene cn data frame
    thresh: threshold in logfold change accross all of them
  """
  pos = genecn.where((genecn > thresh) & (genecn < 1 / thresh), 0).all(0)
  return pos.loc[pos == True].index.values


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
      cols = cols[-1:] + cols[:-1]
      df = df[cols]
      df.to_csv(f, header=False, sep='\t', index=False)


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


def AddToVirtual(virtualname, folderfrom, files):
  """
  will add some files from a taiga folder to a taiga virtual dataset folder and preserve the previous files

  Args:
  ----
    virtualname: the taiga virtual dataset name
    folderfrom: the taiga folder wher the files are
    files: a list(tuples(newfilename,prevfilename))
  """
  assert type(files[0]) is tuple
  versiona = max([int(i['name']) for i in tc.get_dataset_metadata(folderfrom)['versions']])
  versionb = max([int(i['name']) for i in tc.get_dataset_metadata(virtualname)['versions']])
  keep = [(i['name'], i['underlying_file_id']) for i in tc.get_dataset_metadata(virtualname, version=versionb)
          ['datasetVersion']['datafiles'] if 'underlying_file_id' in i]

  for i, val in enumerate(files):
    files[i] = (val[0], folderfrom + '.' + str(versiona) + '/' + val[1])
    for j, v in enumerate(keep):
      if val[0] == v[0]:
        keep.pop(j)

  files.extend(keep)
  print(files)
  tc.update_dataset(dataset_permaname=virtualname, add_taiga_ids=files, upload_file_path_dict={})


def removeColDuplicates(a, prepended=['dm', 'ibm', 'ccle']):
  """
  This function is used to subset a df to only the columns with the most up to date names

  We consider a naming convention preprended_NAME_version and transform it into NAME with latest NAMES
  

  Args:
  ----
    a:
    prepended:

  Returns:
  ------
    a:
  """
  values = []
  for i in a.columns:
    i = i.split('_')
    if i[0] in prepended:
      values.append('_'.join(i[1:]))
    else:
      values.append('_'.join(i))
  a.columns = values
  a = a[a.columns[np.argsort(values)]]
  todrop = []
  for i in range(len(a.columns) - 1):
    e = a.columns[i + 1].split('_')
    if len(e[-1]) == 1:
      if int(e[-1]) > 1 and e[0] == a.columns[i].split('_')[0]:
        todrop.append(a.columns[i])
        print("removing: " + str(a.columns[i]) + " and replacing by: " + str(e))
  a = a.drop(todrop, 1)
  return a


def removeDuplicates(a, loc, prepended=['dm', 'ibm', 'ccle']):
  """
  This function is used to subset a df to only the columns with the most up to date names

  We consider a naming convention preprended_NAME_version and transform it into NAME with latest NAMES
  
  Args:
  ----

  Returns:
  -------
  """
  values = []
  if len(prepended) > 0:
    print('heyyy')
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
