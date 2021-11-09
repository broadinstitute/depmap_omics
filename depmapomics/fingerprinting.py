# fingerprinting.py

import pandas as pd
import numpy as np
import dalmatian as dm
from taigapy import TaigaClient
from genepy.google import gcp
import asyncio
from genepy import terra
from genepy.utils import helper as h
from depmapomics import tracker
from depmapomics.config import *
from depmapomics import terra as myterra
import os 


def updateLOD(wm, sampleset, taiga_dataset, taiga_filename, working_dir, upload_to_taiga=True):
  #Here we update the fingerprint LOD matrix on taiga with the new fingerprints
  # Generate matrix with LOD score for new fingerprint vcfs
  new_lod_list = []
  sample_batch_pair_df = wm.get_entities("sample_batch_pair")
  samples_df = sample_batch_pair_df[sample_batch_pair_df.sample_batch_b.apply(lambda x: x['entityName'] == sampleset)]['cross_checks_out'].tolist()
  for batch in samples_df:
    # could be pd concat
    df = pd.read_csv(batch, sep='\t', comment='#')
    lod_mat = df.pivot(index="LEFT_SAMPLE",
                      columns="RIGHT_SAMPLE", values="LOD_SCORE")
    new_lod_list.append(lod_mat)
  new_lod_mat = pd.concat(new_lod_list)
  new_lod_mat.index.name = None
  new_lod_mat = new_lod_mat.T

  # Update LOD matrix ( have to update (A+a)*(B+b) = (AB)+(aB)+(Ab)+(ab))
  tc = TaigaClient()
  prev_lod_mat =  tc.get(name=taiga_dataset,file=taiga_filename)
  new_ids = set(new_lod_mat.index)
  old_ids = set(prev_lod_mat.index) - set(new_ids)
  updated_lod_mat = pd.concat((prev_lod_mat.loc[old_ids,old_ids],
                               new_lod_mat.loc[new_ids,old_ids]), axis=0)
  updated_lod_mat = pd.concat((updated_lod_mat.loc[new_ids.union(old_ids), old_ids], 
                              new_lod_mat.transpose().loc[new_ids.union(old_ids, new_ids)]), axis=1)
  updated_lod_mat.to_csv(working_dir+taiga_filename+'.csv')
  
  # Upload updated LOD matrix to Taiga
  if upload_to_taiga:
    tc.update_dataset(dataset_permaname=taiga_dataset,
                      changes_description="New bam fingerprints added for "+sampleset,
                      upload_files=[
                          {
                              "path": working_dir+taiga_filename+'.csv',
                              "name": taiga_filename,
                              "format": "NumericMatrixCSV",
                              "encoding": "utf-8"
                          }
                      ],
                      add_all_existing_files=True)
  return new_ids, updated_lod_mat

def checkMismatches(v, ref, samples, thr=100):
  should = {}
  print("\n\nsamples that should match but don't:")
  for u in set(samples.arxspan_id):
    res = v.loc[samples[samples.arxspan_id == u].index,
                ref[ref.arxspan_id == u].index.tolist()]
    for i, j in [(res.index[x], res.columns[y]) for x, y in np.argwhere(res.values < thr)]:
      print('__________________________')
      print(res.loc[i, j])
      print(i, ':', tuple(ref.loc[i, ['arxspan_id', 'version', 'datatype', 'participant_id']].values), j, ':', tuple(
          ref.loc[j, ['arxspan_id', 'version', 'datatype', 'participant_id', 'blacklist']]))
      should.update({str(i) + ': ' + str(tuple(ref.loc[i, ['arxspan_id', 'version', 'datatype', 'participant_id']].values)):
                      str(j)} + ': ' + str(tuple(ref.loc[j, ['arxspan_id', 'version', 'datatype', 'participant_id', 'blacklist']])))
  return should

def checkMatches(v, ref, thr=500):
  print("\n\nsamples that shouldn't match but do")
  previ = ''
  shouldnt = {}
  for i, j in [(v.index[x], v.columns[y]) for x, y in np.argwhere(v.values > thr)]:
      if i == j:
          continue
      if ref.loc[i]['participant_id'] == ref.loc[j]['participant_id']:
          continue
      if i != previ:
          if previ != '':
              shouldnt.update({'_'.join(ref.loc[previ, ['arxspan_id', 'version', 'datatype',
                                                'participant_id', 
                                                'stripped_cell_line_name']].astype(str).values.tolist()): n})
          n = [tuple(ref.loc[j, ['arxspan_id', 'version', 'datatype',
                                'participant_id', 'stripped_cell_line_name']].values)]
      else:
          n.append(tuple(ref.loc[j, ['arxspan_id', 'version', 'datatype',
                                    'participant_id', 'stripped_cell_line_name']].values))
      previ = i
  return shouldnt

def add_sample_batch_pairs(wm, working_dir=WORKING_DIR):
  # add and update sample_batch_pairs and sample_batch_pair_set
  all_sample_sets = wm.get_entities("sample_set").index
  sample_set_a_list = []
  sample_set_b_list = []
  pair_ids = []
  for s in all_sample_sets:
      for t in all_sample_sets:
          sample_set_a_list.append(s)
          sample_set_b_list.append(t)
          pair_ids.append(s + '-' + t)

  # pair_df contains all possible pairs between sample_sets         
  pair_df = pd.DataFrame(
      np.array([sample_set_a_list, sample_set_b_list]).T,
      columns=['sample_batch_a', 'sample_batch_b'],
      index=pair_ids
  )
  pair_df.index.name = 'entity:sample_batch_pair_id'
  
  # update sample_batch_pair
  try:
    wm.upload_entities("sample_batch_pair", pair_df)
  except:
    print("still can't update sample_batch_pair")
    #in case it does not work
    pair_df.to_csv(working_dir + "sample_batch_pairs.tsv", sep='\t')
    if h.askif("Please upload the file ../sample_batch_pairs.tsv to the terra workspace as a new data table, and type yes once \
      finished, else write no to quit and they will not be added"):
      print('sample_batch_pair manually updated')
  
  # replace string entries in sample_batch_pairs with references to sample_sets
  myterra.updateReferences(wm, 'sample_batch_pair', pair_df)
  
  unique_pairs = wm.get_entities('sample_batch_pair').index.tolist()
  sample_batch_pair_set_df = pd.DataFrame(np.transpose(unique_pairs), index=['all'] * len(unique_pairs), columns=['sample_batch_pair'])
  sample_batch_pair_set_df.index.name = 'membership:sample_batch_pair_set_id'
  
  # update sample_batch_pair_set
  try:
    wm.upload_entities("sample_batch_pair_set", sample_batch_pair_set_df)
  except:
    print("still can't update sample_batch_pair_set")
    #in case it does not work
    sample_batch_pair_set_df.to_csv(working_dir + "sample_batch_pair_set.tsv", sep='\t')
    if h.askif("Please upload the file ../sample_batch_pair_set.tsv to the terra workspace as a new data table, and type yes once \
      finished, else write no to quit and they will not be added"):
      print('sample_batch_pair_set manually updated')
  

async def _CCLEFingerPrint(rnasamples, wgssamples, sampleset=SAMPLESETNAME, allbatchpairset=FPALLBATCHPAIRSETS, workspace=FPWORKSPACE, 
working_dir=WORKING_DIR, bamcolname=LEGACY_BAM_COLNAMES, taiga_dataset=TAIGA_FP, taiga_filename=TAIGA_FP_FILENAME):
  """1.1  Generate Fingerprint VCFs

  Here we use Dalmatian to run the fingerprint_bam_with_liftover workflow on Terra.
  This workflow calls Picard ExtractFingerprint to generate a fingerprint VCF and 
  then calls Picard LiftoverVcf to covert this vcf to hg38. To fingerprint hg38 bam files 
  just run fingerprint_bam instead.

  Args:
      samples ([type]): [description]
      sampleset ([type]): [description]
      sid ([type]): [description]
      vcf_list_dir ([type]): [description]
      working_dir ([type]): [description]
      crosscheck_batch_size ([type]): [description]
      recreate_batch ([type]): [description]
      bamcolname ([type]): [description]
      workspace ([type], optional): [description]. Defaults to WORKSPACE.
      vcf_list ([type], optional): [description]. Defaults to None.

  Author:
      William Colgan (wcolgan@broadinstitute.org)
  """

  sid = 'id'
  samples = pd.concat([rnasamples, wgssamples])
  bams = samples[bamcolname]
  bams[sid] = bams.index
  print('adding '+str(len(bams))+' new samples to the fingerprint')
  wm = dm.WorkspaceManager(workspace).disable_hound()

  # Upload sample sheet
  samples_df = pd.DataFrame()
  samples_df = pd.DataFrame(bams[bamcolname + [sid, sid]].values, columns = ["bam_filepath", "bai_filepath", "sample_id", "participant_id"])
  samples_df = samples_df.set_index('sample_id')
  wm.upload_samples(samples_df, add_participant_samples=True)
  wm.update_sample_set(sampleset, samples_df.index)
  add_sample_batch_pairs(wm, working_dir=WORKING_DIR)

  # Submit fingerprinting jobs, generate vcf files for all lines
  submission_id = wm.create_submission("fingerprint_bam_with_liftover", sampleset, 
                                       'sample_set', expression='this.samples')
  await terra.waitForSubmission(workspace, submission_id)

  #1.2  Crosscheck Fingerprint VCFs
  # Here we use Dalmation to run the crosscheck_vcfs workflow on Terra. 
  # This workflow calls Picard CrosscheckFingerprints to compare vcfs between every 
  # sample_batch_pair in the sample_batch_pair_set

  # Submit crosscheck jobs
  conf = wm.get_config("crosscheck_vcfs")
  wm.update_config(conf)
  submission_id = wm.create_submission("crosscheck_vcfs", allbatchpairset, 'sample_batch_pair_set', expression='this.sample_batch_pairs')
  await terra.waitForSubmission(workspace, submission_id)

  # save config after jobs finish running
  terra.saveWorkspace(workspace,'data/'+sampleset+'/FPconfig/')

  # Update LOD matrix
  new_ids, updated_lod_mat = updateLOD(wm, sampleset, taiga_dataset, taiga_filename, working_dir, upload_to_taiga=False)

  # finding issues with the dataset
  v = updated_lod_mat.loc[new_ids]
  ref = tracker.getTracker()
  ref = ref.append(samples)

  # find samples that should match but don't
  should = checkMismatches(v, ref, samples)

  # find samples that shouldn't match but do
  shouldnt = checkMatches(v, ref)

  return updated_lod_mat, should, shouldnt