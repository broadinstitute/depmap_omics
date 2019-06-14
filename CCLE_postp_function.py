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


def CreateDatasetWithNewCellLines(wfrom, wto, source, samplesetname):
  refsamples = wto.get_samples()
  refids = refsamples['participant'].tolist()
  refids = [val[val.index('ACH'):] for val in refids if 'ACH' in val]
  samples = wfrom.get_samples()
  samples = samples[samples['individual_alias'].str.contains('ACH')][~samples['individual_alias'].str.slice(0, 10).isin(refids)]
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
