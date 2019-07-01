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
from taigapy import TaigaClient
tc = TaigaClient()

CHROMLIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX']


def CreateDatasetWithNewCellLines(wfrom, wto, source, samplesetname):
  """

  """
  refsamples = wto.get_samples()
  refids = refsamples['participant'].tolist()
  refids = [val[val.index('ACH'):] for val in refids if 'ACH' in val]
  samples = wfrom.get_samples()
  samples = samples[samples['individual_alias'].str.contains('ACH')][
    ~samples['individual_alias'].str.slice(0, 10).isin(refids)]
  for ind, val in samples.iterrows():
    refsamples = refsamples.append(pd.Series(
        {
            "CCLE_name": val['sample_alias'],
            "WES_bai": val['crai_or_bai_path'],
            "WES_bam": val['cram_or_bam_path'],
            "Source": source,
            "participant": val['individual_alias'][:10],
        }, name=source + '_' + val['individual_alias'][:10]))
  print("uploading new samples")
  wto.upload_samples(refsamples)
  sample_ids = [source + '_' + i for i in samples['individual_alias'].str.slice(0, 10).tolist()]
  print("creating a sample set")
  wto.update_sample_set(sample_set_id=samplesetname, sample_ids=sample_ids)
  return sample_ids
####
#
# VALIDATION
#


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
  pdb.set_trace()
  maindata = pd.read_csv(main_filename, sep='\t')
  if '.' in maindata["DepMap_ID"][0]:
    maindata["DepMap_ID"] = [i[0] for i in maindata["DepMap_ID"].str.split('.').tolist()]
  samples = set(maindata["DepMap_ID"].tolist())
  with open(main_filename, 'a') as f:  
    for input_filename in input_filenames:
      df = pd.read_csv(input_filename,  sep='\t')
      input_filename = input_filename.split('/')[-1].split('.')[0]
      if input_filename in samples:
        print(input_filename + " is Already in main fusions")
      df['DepMap_ID'] = pd.Series([input_filename]*len(df.index.tolist()), index=df.index)
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
    index= False, index_label=False, compression='gzip')
  transcripts_tpm.to_csv('temp/' + main_filename['rsem_transcripts_tpm'].split('/')[-1], sep='\t', 
    index= False, index_label=False, compression='gzip')
  genes_tpm.to_csv('temp/' + main_filename['rsem_genes_tpm'].split('/')[-1], sep='\t', 
    index= False, index_label=False, compression='gzip')



def ExtractStarQualityInfo(samplesetname, wm):
  """
  """
  a = wm.get_samples().loc[wm.get_sample_sets().loc[samplesetname].samples].star_logs
  for i, sample in enumerate(a):
      if sample is None:
          print("no log file found for: "+a.index[i])
      for log in sample:
          if 'final.out' in log:
              print("copying "+a.index[i])
              os.system('gsutil cp '+log+' temp/')
  os.system("cat temp/*.Log.final.out > temp/"+samplesetname+".txt")
  os.system("rm temp/*.Log.final.out")
