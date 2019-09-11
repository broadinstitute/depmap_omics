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
from taigapy import TaigaClient
tc = TaigaClient()

CHROMLIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX']


extract_defaults = {
    'name': 'sample_alias',
    'bai_path': 'crai_or_bai_path',
    'bam_path': 'cram_or_bam_path',
    'source': 'source',
    'id': 'individual_alias',
    'participant': 'individual_alias'
}


def read(filename):
  f = open(filename, 'r')
  x = f.read().splitlines()
  f.close()
  return x


def createDatasetWithNewCellLines(wto, samplesetname,
                                  wfrom1, source1, wfrom2=None, source2='U', wfrom3=None, source3='U',
                                  forcekeep=[], addonly=[], match='ACH', other_to_add=[], extract=extract_defaults,
                                  slicepos=10, **kwargs):
  """
  wto: dalmatian.workspacemanager the workspace where you want to create the tsvs
  samplesetname: String the name of the sample set created by this new set of samples updated
  wfrom1: dalmatian.workspacemanager the workspace where the Samples to add are stored
  source1: the corresponding source name
  Can add as much as one wants
  forcekeep: list of sample id that you want to keep even if already in the previous workspace (will
  cause an overwrite)
  addonly: list of sample id that you only want to add
  match: substring that has to be matched against the id of the samples to add them
  other_to_add:
  extract: if you want to specify what values should refer to what
    dict{    'name':
    'bai_path':
    'bam_path':
    'source':
    'id':
    ...}

  """
  refsamples = wto.get_samples()
  refids = refsamples['participant'].tolist()
  refids = [val[val.index('ACH'):] for val in refids if 'ACH' in val]

  samples1 = wfrom1.get_samples().replace(np.nan, '', regex=True).reset_index()
  if len(addonly) > 0:
    samples1 = samples1
  samples1 = samples1[samples1[extract['id']].str.contains(match)][
      (~samples1[extract['id']].str.slice(0, slicepos).isin(refids)) |
      (samples1[extract['id']].isin(forcekeep))]
  samples1[extract['source']] = [source1] * samples1.shape[0]

  samples = []
  for key, val in kwargs.iteritems():
    if key[:-1] == "source":
      sample = val.get_samples().replace(np.nan, '', regex=True).reset_index()
      sample = samples[samples[extract['id']].str.contains(match)][
          (~samples[extract['id']].str.slice(0, slicepos).isin(refids)) |
          (samples[extract['id']].isin(forcekeep))]
        [samples[extract['id']].isin(addonly)]

      sample[extract['source']] = [kwargs["source" + key[-1]]] * samples2.shape[0]
      if 'sample_type' in sample.columns.tolist():
        sample = sample[sample['sample_type'] == 'Tumor']
      samples.append(sample)

  samples.append(samples1)
  samples = pd.concat(samples, sort=False)

  # if wfrom2 is not None:
  #   samples2 = wfrom2.get_samples().replace(np.nan, '', regex=True).reset_index()
  #   samples2 = samples2[samples2[extract['id']].str.contains(match)][
  #       (~samples2[extract['id']].str.slice(0, slicepos).isin(refids)) |
  #       (samples2[extract['id']].isin(forcekeep))][samples2[extract['id']].isin(addonly)]
  #   samples2[extract['source']] = [source2] * samples2.shape[0]
  # else:
  #   samples2 = pd.DataFrame()
  # if wfrom3 is not None:
  #   samples3 = wfrom3.get_samples().replace(np.nan, '', regex=True).reset_index()
  #   samples3 = samples3[samples3[extract['id']].str.contains(match)][
  #       (~samples3[extract['id']].str.slice(0, slicepos).isin(refids)) |
  #       (samples3[extract['id']].isin(forcekeep))][samples3[extract['id']].isin(addonly)]
  #   samples3[extract['source']] = [source3] * samples3.shape[0]
  # else:
  #   samples3 = pd.DataFrame()
  # samples = pd.concat([samples1, samples2, samples3], sort=False)

  notfound = set(samples.index.tolist()) - (set(samples.index.tolist()) & set(addonly))
  if len(notfound) > 0:
    print('we did not found:' + str(notfound))
  sample_ids = []
  for ind, val in samples.iterrows():
    name = val[extract['source']] + '_' + val[extract['id']][:slicepos]
    refsamples = refsamples.append(pd.Series(
        {
            "CCLE_name": val[extract['name']],
            "WES_bai": val[extract['bai']],
            "WES_bam": val[extract['bam']],
            "Source": val[extract['source']],
            "participant": val[extract['participant']][:slicepos],

        }, name=name))
    sample_ids.append(name)

  print("uploading new samples")
  wto.upload_samples(refsamples)
  print("creating a sample set")
  wto.update_sample_set(sample_set_id=samplesetname, sample_ids=sample_ids)
  return sample_ids
#################
#
# VALIDATION
#####################


def checkAmountOfSegments(segmentcn, thresh=850):
  """
  if there is too many segments, something might be wrong
  segmentcn: segment dataframe
  thresh: max ok amount
  """
  segmentcn = renameColumns(segmentcn)
  celllines = set(segmentcn["DepMap_ID"].tolist())
  for cellline in celllines:
    if segmentcn[segmentcn["DepMap_ID"] == cellline].shape[0] > thresh:
      print(cellline, segmentcn[segmentcn["DepMap_ID"] == cellline].shape[0])


def checkGeneChangeAccrossAll(genecn, thresh=1.5):
  """
  genecn: gene cn data frame
  thresh: threshold in logfold change accross all of them
  """
  pos = genecn.where((genecn > thresh) & (genecn < 1 / thresh), 0).all(0)
  return pos.loc[pos == True].index.values


def renameColumns(df):
  return(df.rename(columns={'Sample': 'DepMap_ID', 'CONTIG': 'Chromosome', 'START': 'Start',
                            'END': 'End', 'seqnames': 'Chromosome', 'start': 'Start', 'end': 'End'}))


def checkDifferencesWESWGS(segmentcn_wes, segmentcn_wgs, chromlist=CHROMLIST):
  """
  if the similarity between WES and WGS on same cell line is different there might be problems.
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


def ExtractStarQualityInfo(samplesetname, wm):
  """
  """
  a = wm.get_samples().loc[wm.get_sample_sets().loc[samplesetname].samples].star_logs
  for i, sample in enumerate(a):
    if sample is None:
      print("no log file found for: " + a.index[i])
    for log in sample:
      if 'final.out' in log:
        print("copying " + a.index[i])
        os.system('gsutil cp ' + log + ' temp/')
  os.system("cat temp/*.Log.final.out > temp/" + samplesetname + ".txt")
  os.system("rm temp/*.Log.final.out")
