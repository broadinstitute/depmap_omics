from depmapomics import constants
from depmapomics import env_config

# fingerprinting.py

import pandas as pd
import numpy as np
import dalmatian as dm
from mgenepy import terra
from mgenepy.utils import helper as h
from depmap_omics_upload import tracker as track

from depmapomics import terra as myterra
from depmapomics.fp_snp import updateLOD, checkMismatches, checkMatches


def add_sample_batch_pairs(wm, working_dir=constants.WORKING_DIR):
    """add and update sample_batch_pairs and sample_batch_pair_set in workspace

    Args:
        wm (dm.workspaceManager): dalmatian workspace manager for the terra workspace
        working_dir (str): working directory where we store temp files. optional, defaults to constants.WORKING_DIR

    """
    all_sample_sets = wm.get_entities("sample_set").index
    # only create sample set pairs between merged (hg19 + hg38) sets
    all_sample_sets = [i for i in all_sample_sets if not i.endswith("subset")]
    sample_set_a_list = []
    sample_set_b_list = []
    pair_ids = []
    for s in all_sample_sets:
        for t in all_sample_sets:
            sample_set_a_list.append(s)
            sample_set_b_list.append(t)
            pair_ids.append(s + "-" + t)

    # pair_df contains all possible pairs between all sample_sets
    pair_df = pd.DataFrame(
        np.array([sample_set_a_list, sample_set_b_list]).T,
        columns=["sample_batch_a", "sample_batch_b"],
        index=pair_ids,
    )
    pair_df.index.name = "entity:sample_batch_pair_id"

    # update sample_batch_pair
    try:
        wm.upload_entities("sample_batch_pair", pair_df, model="flexible")
    except:
        print("still can't update sample_batch_pair")
        # in case it does not work
        pair_df.to_csv(working_dir + "sample_batch_pairs.tsv", sep="\t")
        if h.askif(
            "Please upload the file "
            + working_dir
            + "sample_batch_pairs.tsv to the terra workspace as a new data table, and type yes once \
      finished, else write no to quit and they will not be added"
        ):
            print("sample_batch_pair manually updated")
        else:
            raise ValueError("sample_batch_pair not manually updated")

    # replace string entries in sample_batch_pairs with references to sample_sets
    myterra.updateReferences(wm, "sample_batch_pair", pair_df)

    unique_pairs = wm.get_entities("sample_batch_pair").index.tolist()
    existing_pairs = wm.get_entities("sample_batch_pair_set").loc["all", :][0]
    existing_pairs = [x["entityName"] for x in existing_pairs]
    pairs_to_add = list(set(unique_pairs) - set(existing_pairs))

    sample_batch_pair_set_df = pd.DataFrame(
        np.transpose(pairs_to_add),
        index=["all"] * len(pairs_to_add),
        columns=["sample_batch_pair"],
    )
    sample_batch_pair_set_df.index.name = "membership:sample_batch_pair_set_id"

    # update sample_batch_pair_set to include all sample_batch_pairs
    try:
        wm.upload_entities("sample_batch_pair_set", sample_batch_pair_set_df)
    except:
        print("still can't update sample_batch_pair_set")
        # in case it does not work
        sample_batch_pair_set_df.to_csv(
            working_dir + "sample_batch_pair_set.tsv", sep="\t"
        )
        if h.askif(
            "Please upload the file "
            + working_dir
            + "sample_batch_pair_set.tsv to the terra workspace as a new data table, and type yes once \
      finished, else write no to quit and they will not be added"
        ):
            print("sample_batch_pair_set manually updated")
        else:
            raise ValueError("sample_batch_pair_set not manually updated")


async def fingerPrint(
    samples,
    sid="id",
    use_gumbo=False,
    sampleset=constants.SAMPLESETNAME,
    allbatchpairset=constants.FPALLBATCHPAIRSETS,
    workspace=env_config.FPWORKSPACE,
    working_dir=constants.WORKING_DIR,
    bamcolname=constants.LEGACY_BAM_COLNAMES + constants.HG38_CRAM_COLNAMES,
    terrabamcolname=["bam_filepath", "bai_filepath"] + constants.HG38_CRAM_COLNAMES,
    prev_mat_df=None,
    updated_mat_filename="",
):
    """1.1  Generate Fingerprint VCFs

    Here we use Dalmatian to run the fingerprint_bam_with_liftover workflow on Terra.
    This workflow calls Picard ExtractFingerprint to generate a fingerprint VCF and
    then calls Picard LiftoverVcf to covert this vcf to hg38. To fingerprint hg38 bam files
    just run fingerprint_bam instead.

    Args:
        samples (pd.DataFrame): dataframe for samples
        sampleset (str): name for the sample set to be processed
        sid (str): column name for sample id
        working_dir (str): location for temp files to be saved in
        bamcolname (list[str], optional): column names for bam and bai files on terra
        workspace (str, optional): terra workspace for fingerprinting. Defaults to env_config.FPWORKSPACE.

    Returns:
        updated_lod_mat (pd.DataFrame): updated lod matrix
        mismatches (dict): dict representing the mismatches
        matches (dict): dict representing the matches

    Author:
        William Colgan (wcolgan@broadinstitute.org)
    """

    bams = samples[bamcolname]
    bams[sid] = bams.index

    wm = dm.WorkspaceManager(workspace).disable_hound()

    # Upload sample sheet
    samples_df = pd.DataFrame(
        bams[bamcolname + [sid, sid]].values,
        columns=terrabamcolname + ["sample_id", "participant_id"],
    )
    samples_df = samples_df.set_index("sample_id")
    print("adding " + str(len(bams)) + " new samples to the fingerprinting workspace")
    wm.upload_samples(samples_df, add_participant_samples=True)
    # create a sample set containing both rna and wgs data
    wm.update_sample_set(sampleset, samples_df.index)
    # and one for only hg19
    wm.update_sample_set(
        sampleset + "_hg19subset",
        samples_df[~samples_df[terrabamcolname[0]].isna()].index,
    )
    # and one for only hg38
    wm.update_sample_set(
        sampleset + "_hg38subset",
        samples_df[~samples_df[constants.HG38_CRAM_COLNAMES[0]].isna()].index,
    )
    add_sample_batch_pairs(wm, working_dir=constants.WORKING_DIR)

    # Submit fingerprinting jobs, generate vcf files for all lines
    submission_id_hg19 = wm.create_submission(
        "fingerprint_bam_with_liftover",
        sampleset + "_hg19subset",
        "sample_set",
        expression="this.samples",
    )
    submission_id_hg38 = wm.create_submission(
        "fingerprint_cram_hg38",
        sampleset + "_hg38subset",
        "sample_set",
        expression="this.samples",
    )
    submission_id_hipstr = wm.create_submission(
        "hipSTR",
        sampleset + "_hg38subset",
        "sample_set",
        expression="this.samples",
    )
    await terra.waitForSubmission(workspace, submission_id_hg19)
    await terra.waitForSubmission(workspace, submission_id_hg38)
    await terra.waitForSubmission(workspace, submission_id_hipstr)

    # 1.2  Crosscheck Fingerprint VCFs
    # Here we use Dalmation to run the crosscheck_vcfs workflow on Terra.
    # This workflow calls Picard CrosscheckFingerprints to compare vcfs between every
    # sample_batch_pair in the sample_batch_pair_set

    # Submit crosscheck jobs
    conf = wm.get_config("crosscheck_vcfs")
    wm.update_config(conf)
    submission_id = wm.create_submission(
        "crosscheck_vcfs",
        allbatchpairset,
        "sample_batch_pair_set",
        expression="this.sample_batch_pairs",
    )
    await terra.waitForSubmission(workspace, submission_id)

    # save config after jobs finish running
    terra.saveWorkspace(workspace, "data/" + sampleset + "/FPconfig/")

    # Update LOD matrix
    new_ids, updated_lod_mat = updateLOD(
        wm,
        sampleset,
        working_dir,
        save_new_mat=False,
        prev_mat_df=prev_mat_df,
        updated_mat_filename=updated_mat_filename,
    )

    # finding issues with the dataset
    latest_lod_mat = updated_lod_mat.loc[new_ids]

    ref = pd.DataFrame()

    if use_gumbo:
        mytracker = track.SampleTracker()
        ref = mytracker.add_model_cols_to_seqtable(["PatientID", "ModelID"])
        pr_table = mytracker.read_pr_table()
        ref["ModelCondition"] = ref["ProfileID"].apply(
            lambda x: pr_table.loc[x, "ModelCondition"]
        )
        # add model and participant info to new samples that are not in the sequencing table yet
        samples["ModelID"] = samples["ProfileID"].apply(
            lambda x: mytracker.lookup_model_from_pr(x, "ModelID")
        )
        samples["PatientID"] = samples["ProfileID"].apply(
            lambda x: mytracker.lookup_model_from_pr(x, "PatientID")
        )
        samples["ModelCondition"] = samples["ProfileID"].apply(
            lambda x: pr_table.loc[x, "ModelCondition"]
        )

    else:
        ref = ref.append(samples)

    # find samples that should match but don't
    mismatches = checkMismatches(latest_lod_mat, ref, samples)

    # find samples that shouldn't match but do
    matches = checkMatches(latest_lod_mat, ref)

    return updated_lod_mat, mismatches, matches


async def _CCLEFingerPrint(
    rnasamples,
    wgssamples,
    sid="id",
    sampleset=constants.SAMPLESETNAME,
    allbatchpairset=constants.FPALLBATCHPAIRSETS,
    workspace=env_config.FPWORKSPACE,
    working_dir=constants.WORKING_DIR,
    bamcolname=constants.LEGACY_BAM_COLNAMES + constants.HG38_CRAM_COLNAMES,
    taiga_dataset=env_config.TAIGA_FP,
    updated_mat_filename=constants.TAIGA_FP_FILENAME,
    upload_to_taiga=True,
):
    """CCLE fingerprinting function

    Args:
        samples ([str]): list of sample ids to check SNP fingerprints for
        sampleset (str): name of this new sample set to be added to the SNP fingerprinting terra workspce
        sid (str): [description]
        working_dir (str): location where intermediate outputs are stored
        bamcolname ([str]): columns in gumbo where bam file locations are stored
        workspace (str, optional): name of the SNP fingerprinting workspace on terra. Defaults to env_config.FPWORKSPACE.

    Returns:
        updated_lod_mat (pd.DataFrame): updated lod matrix
        mismatches (dict): dict representing the mismatches
        matches (dict): dict representing the matches
    """
    mytracker = track.SampleTracker()
    seq_table = mytracker.read_seq_table()

    samples = seq_table.loc[rnasamples + wgssamples]

    from taigapy import TaigaClient

    tc = TaigaClient()
    prev_lod_mat = tc.get(name=taiga_dataset, file=updated_mat_filename)

    # call generic function
    updated_lod_mat, mismatches, matches = await fingerPrint(
        samples,
        use_gumbo=True,
        sid=sid,
        sampleset=sampleset,
        allbatchpairset=allbatchpairset,
        workspace=workspace,
        working_dir=working_dir,
        bamcolname=bamcolname,
        prev_mat_df=prev_lod_mat,
        updated_mat_filename=updated_mat_filename,
    )

    if upload_to_taiga:
        # Upload updated LOD matrix to Taiga
        print("uploading updated lod matrix to taiga")
        tc.update_dataset(
            dataset_permaname=taiga_dataset,
            changes_description="New bam fingerprints added for " + sampleset,
            upload_files=[
                {
                    "path": working_dir + updated_mat_filename + ".csv",
                    "name": updated_mat_filename,
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                }
            ],
            add_all_existing_files=True,
        )

    return updated_lod_mat, mismatches, matches
