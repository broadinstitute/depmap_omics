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
import os 


def recreateBatch(vcf_list, vcf_list_dir, working_dir, wm, crosscheck_batch_size):
  # Create batch files listing all vcfs in fingerprints dir and upload to bucket
  # (NEW VERSION ONLY) will only needed if need to recreate batches
  if not vcf_list:
    vcf_list = gcp.lsFiles([vcf_list_dir])
  vcf_list = wm.get_samples()["fingerprint_vcf"].tolist()
  batches = []
  for i, l in enumerate(range(0, len(vcf_list), crosscheck_batch_size)):
    f = open(working_dir + "vcf_batch_"+str(i), 'w')
    f.write("\n".join(vcf_list[l:l + crosscheck_batch_size]))
    f.close()
    batches.append(working_dir+"vcf_batch_"+str(i))
  gcp.cpFiles(batches, vcf_list_dir)
  return batches


def updateLOD(wm, sampleset, taiga_dataset, taiga_filename, working_dir, upload_to_taiga=True):
  #Here we update the fingerprint LOD matrix on taiga with the new fingerprints
  # Generate matrix with LOD score for new fingerprint vcfs
  new_lod_list = []
  samples_df = wm.get_entities("sample_group")['cross_checks_out'].tolist()
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

async def _CCLEFingerPrint(samples, sampleset=SAMPLESETNAME, allsampleset=FPALLSAMPLESET, workspace=FPWORKSPACE, vcf_list=None, 
  vcf_list_dir=VCF_LIST_DIR, working_dir=WORKING_DIR, crosscheck_batch_size=CROSSCHECK_BATCH_SIZE, recreate_batch=RECREATE_BATCH, bamcolname=LEGACY_BAM_COLNAMES,
  taiga_dataset=TAIGA_FP, taiga_filename=TAIGA_FP_FILENAME):
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
  bams = samples[bamcolname]
  bams[sid] = bams.index
  print('adding '+str(len(bams))+' new samples to the fingerprint')
  wm = dm.WorkspaceManager(workspace).disable_hound()

  batches = []
  if recreate_batch:
    print('recreating batch')
    batches = recreateBatch(vcf_list, vcf_list_dir, working_dir, wm, crosscheck_batch_size)

  # Upload sample sheet
  samples_df = pd.DataFrame()
  samples_df = pd.DataFrame(bams[bamcolname + [sid, sid]].values, columns = ["bam_filepath", "bai_filepath", "sample_id", "participant_id"])
  samples_df = samples_df.set_index('sample_id')
  wm.upload_samples(samples_df, add_participant_samples=True)
  wm.update_sample_set(sampleset, samples_df.index)

  # Submit jobs 
  submission_id = wm.create_submission("fingerprint_bam_with_liftover", sampleset, 
                                       'sample_set', expression='this.samples')
  await terra.waitForSubmission(workspace, submission_id)

  #1.2  Crosscheck Fingerprint VCFs
  #Here we use Dalmation to run the crosscheck_vcfs workflow on Terra. 
  # This workflow calls Picard CrosscheckFingerprints to compare the new 
  # fingerprint vcfs to batches of existing fingerprint vcfs in fingerprints_dir
  # Create list with new vcfs and upload to bucket
  f = open(working_dir + sampleset, 'w')
  f.write(('\n').join(wm.get_samples().loc[samples_df.index, 'fingerprints'].tolist()))
  f.close()
  gcp.cpFiles(working_dir + sampleset, vcf_list_dir)
  os.system('rm '+working_dir + sampleset)

  # Upload sample sheet
  if recreate_batch:
    sample_group_df = pd.DataFrame(data={"entity:sample_group_id" : batches, "vcf_group" : [vcf_list_dir + x for x in batches]}).set_index('entity:sample_group_id')
  else:
    sample_group_df = pd.DataFrame(data={"entity:sample_group_id" : [sampleset], "vcf_group" : [vcf_list_dir+sampleset]}).set_index('entity:sample_group_id')
  
  print(wm.get_entities('sample_group').index.tolist())

  try:
    wm.upload_entities("sample_group", sample_group_df)
  except:
    print("still can't update entities")
    #in case it does not work
    sample_group_df.to_csv("../temp.tsv", sep='\t')
    if h.askif("Please upload the file ../temp.tsv to the terra workspace as a new sample_group, and type yes once \
      finished, else write no to quit and they will not be added"):
      print('sample_group ' + sampleset + ' manually added as a sample_group')
  
  try:
    wm.update_entity_set("sample_group_set", set_id=allsampleset,
                         entity_ids=wm.get_entities('sample_group').index)
  except:
    print("still can't update entities")
    #in case it does not work
    if h.askif("Please go to the sample_group_set tab on terra and add the sample_group name to the all_samples set, and type yes once \
      finished, else write no to quit and they will not be added"):
      print(allsampleset + ' manually updated')

  # Submit jobs
  conf = wm.get_config("crosscheck_vcfs")
  conf['inputs']['crosscheck.run_crosscheck.vcf_second_input_file'] = '"'+vcf_list_dir+sampleset+'"'
  wm.update_config(conf)
  submission_id = wm.create_submission("crosscheck_vcfs", allsampleset, 
  'sample_group_set', expression='this.sample_groups')
  await terra.waitForSubmission(workspace, submission_id)

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