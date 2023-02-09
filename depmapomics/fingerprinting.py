# fingerprinting.py

import pandas as pd
import numpy as np
import dalmatian as dm
from genepy import terra
from genepy.utils import helper as h
from depmap_omics_upload import tracker as track
from depmapomics.config import *
from depmapomics import terra as myterra


def updateLOD(
    wm,
    sampleset,
    working_dir,
    batch_pair_entity="sample_batch_pair",
    save_new_mat=True,
    new_mat_filename="new_fingerprint_lod_matrix.csv",
    prev_mat_df=None,
    updated_mat_filename="updated_lod_matrix.csv",
):
    """Here we update the fingerprint LOD matrix on taiga with the new fingerprints
    and enerate matrix with LOD score for new fingerprint vcfs

    Args:
        wm (wmanager): workspace manager of current FP terra workspace
        sampleset (str): name of current sample set
        working_dir (str): directory where the lod matrix file will be saved
        batch_pair_entity (str): optional. name of entity on terra workspace that represents sample batch pairs. Defaults to "sample_batch_pair".
        save_new_mat (bool): if new lod matrix (not combined with a previous matrix) is saved
        new_mat_filename (str): name of the new lod matrix file if save_new_mat
        prev_mat_df (pd.DatatFrame): a previous lod matrix file that will be merged with the new lod matrix
        updated_mat_filename (str): name of the updated (new merged with prev) lod matrix file

    Returns:
        new_ids (set): set of ids of new samples
        updated_lod_mat (pd.DataFrame): dataframe for the updated lod score matrix
    """
    print("generating lod matrix with new samples")
    new_lod_list = []
    sample_batch_pair_df = wm.get_entities(batch_pair_entity)
    samples_df = sample_batch_pair_df[
        sample_batch_pair_df.sample_batch_b.apply(
            lambda x: x["entityName"] == sampleset
        )
    ]["cross_checks_out"].tolist()
    for batch in samples_df:
        # could be pd concat
        df = pd.read_csv(batch, sep="\t", comment="#")
        lod_mat = df.pivot(
            index="LEFT_SAMPLE", columns="RIGHT_SAMPLE", values="LOD_SCORE"
        )
        new_lod_list.append(lod_mat)
    new_lod_mat = pd.concat(new_lod_list)
    new_lod_mat.index.name = None
    new_lod_mat = new_lod_mat.T

    if save_new_mat:
        print("saving lod matrix with new samples")
        new_lod_mat.to_csv(working_dir + new_mat_filename + ".csv")

    new_ids = set(new_lod_mat.index)
    updated_lod_mat = new_lod_mat

    # Update LOD matrix ( have to update (A+a)*(B+b) = (AB)+(aB)+(Ab)+(ab))
    if prev_mat_df is not None:
        print("merging new lod matrix with the previous one")
        prev_lod_mat = prev_mat_df
        old_ids = set(prev_lod_mat.index) - set(new_ids)
        updated_lod_mat = pd.concat(
            (prev_lod_mat.loc[old_ids, old_ids], new_lod_mat.loc[new_ids, old_ids]),
            axis=0,
        )
        updated_lod_mat = pd.concat(
            (
                updated_lod_mat.loc[new_ids.union(old_ids), old_ids],
                new_lod_mat.transpose().loc[new_ids.union(old_ids, new_ids)],
            ),
            axis=1,
        )
        print("saving merged (updated) lod matrix")
        updated_lod_mat.to_csv(working_dir + updated_mat_filename + ".csv")
    return new_ids, updated_lod_mat


def checkMismatches(lod_mat, ref, samples, thr=100):
    """find samples that should match but don't, by comparing the lod scores

    with thr

    Args:
        lod_mat (pd.DataFrame): dataframe containing lod scores
        ref (pd.DataFrame): dataframe representing the sample tracker
        samples (pd.DataFrame): dataframe representing info for new samples
        thr (int): lod score under which we consider two samples mismatched. optional, defaults to 100

    Returns:
        mismatches (dict): dict representing the mismatches
    """
    mismatches = {}
    print("\n\nsamples that should match but don't:")
    for u in set(samples.ModelID):
        scores = lod_mat.loc[
            samples[samples.ModelID == u].index,
            ref[
                (ref.ModelID == u)
                & (ref.index.isin(lod_mat.columns))
                & (ref.blacklist != 1)
            ].index.tolist(),
        ]
        for i, j in [
            (scores.index[x], scores.columns[y])
            for x, y in np.argwhere(scores.values < thr)
        ]:
            print("__________________________")
            print(scores.loc[i, j])
            print(
                i,
                ":",
                tuple(
                    ref.loc[
                        i, ["ModelID", "version", "expected_type", "PatientID"]
                    ].values
                ),
                j,
                ":",
                tuple(
                    ref.loc[
                        j,
                        [
                            "ModelID",
                            "version",
                            "expected_type",
                            "PatientID",
                            "blacklist",
                        ],
                    ]
                ),
            )
            mismatches.update(
                {
                    str(i)
                    + ": "
                    + str(
                        tuple(
                            ref.loc[
                                i,
                                ["ModelID", "version", "expected_type", "PatientID"],
                            ].values
                        )
                    )
                    + ": "
                    + str(j): str(
                        tuple(
                            ref.loc[
                                j,
                                [
                                    "ModelID",
                                    "version",
                                    "expected_type",
                                    "PatientID",
                                    "blacklist",
                                ],
                            ]
                        )
                    )
                }
            )
    return mismatches


def checkMatches(lod_mat, ref, thr=500):
    """find samples that shouldn't match but do, by comparing the lod scores

    with thr

    Args:
        lod_mat (pd.DataFrame): dataframe containing lod scores
        ref (pd.DataFrame): dataframe representing the sample tracker
        samples (pd.DataFrame): dataframe representing info for new samples. optional, defaults to 500
        thr (int): lod score above which we consider two samples matching

    Returns:
        matches (dict): dict representing the matches
    """

    print("\n\nsamples that shouldn't match but do")
    previ = ""
    matches = {}
    matched_samples = []
    for i, j in [
        (lod_mat.index[x], lod_mat.columns[y])
        for x, y in np.argwhere(lod_mat.values > thr)
    ]:
        if (
            i in ref.index
            and j in ref.index
            and ref.loc[i, "blacklist"] is not True
            and ref.loc[j, "blacklist"] is not True
        ):
            if i == j:
                continue
            if ref.loc[i, "PatientID"] == ref.loc[j, "PatientID"]:
                continue
            if i != previ:
                if previ != "":
                    matches.update(
                        {
                            "_".join(
                                [previ]
                                + ref.loc[
                                    previ,
                                    [
                                        "ModelID",
                                        "ProfileID",
                                        "version",
                                        "expected_type",
                                        "PatientID",
                                    ],
                                ]
                                .astype(str)
                                .values.tolist()
                            ): matched_samples
                        }
                    )
                matched_samples = [
                    tuple(
                        ref.loc[
                            j,
                            [
                                "ModelID",
                                "ProfileID",
                                "version",
                                "expected_type",
                                "PatientID",
                            ],
                        ].values
                    )
                ]
            else:
                matched_samples.append(
                    tuple(
                        ref.loc[
                            j,
                            [
                                "ModelID",
                                "ProfileID",
                                "version",
                                "expected_type",
                                "PatientID",
                            ],
                        ].values
                    )
                )
            previ = i
    return matches


def add_sample_batch_pairs(wm, working_dir=WORKING_DIR):
    """add and update sample_batch_pairs and sample_batch_pair_set in workspace

    Args:
        wm (dm.workspaceManager): dalmatian workspace manager for the terra workspace
        working_dir (str): working directory where we store temp files. optional, defaults to WORKING_DIR

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
        wm.upload_entities("sample_batch_pair", pair_df)
    except:
        print("still can't update sample_batch_pair")
        # in case it does not work
        pair_df.to_csv(working_dir + "sample_batch_pairs.tsv", sep="\t")
        if h.askif(
            "Please upload the file ../sample_batch_pairs.tsv to the terra workspace as a new data table, and type yes once \
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
            "Please upload the file ../sample_batch_pair_set.tsv to the terra workspace as a new data table, and type yes once \
      finished, else write no to quit and they will not be added"
        ):
            print("sample_batch_pair_set manually updated")
        else:
            raise ValueError("sample_batch_pair_set not manually updated")


async def fingerPrint(
    samples,
    sid="id",
    use_gumbo=False,
    sampleset=SAMPLESETNAME,
    allbatchpairset=FPALLBATCHPAIRSETS,
    workspace=FPWORKSPACE,
    working_dir=WORKING_DIR,
    bamcolname=LEGACY_BAM_COLNAMES + HG38_CRAM_COLNAMES,
    terrabamcolname=["bam_filepath", "bai_filepath"] + HG38_CRAM_COLNAMES,
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
        workspace (str, optional): terra workspace for fingerprinting. Defaults to FPWORKSPACE.

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
    # create a sample set containing both hg19 and hg38 samples
    wm.update_sample_set(sampleset, samples_df.index)
    # and one for only hg19
    wm.update_sample_set(
        sampleset + "_hg19subset",
        samples_df[~samples_df[LEGACY_BAM_COLNAMES[0]].isna()].index,
    )
    # and one for only hg38
    wm.update_sample_set(
        sampleset + "_hg38subset",
        samples_df[~samples_df[HG38_CRAM_COLNAMES[0]].isna()].index,
    )
    add_sample_batch_pairs(wm, working_dir=WORKING_DIR)

    # Submit fingerprinting jobs, generate vcf files for all lines
    submission_id = wm.create_submission(
        "fingerprint_bam_with_liftover",
        sampleset,
        "sample_set",
        expression="this.samples",
    )
    await terra.waitForSubmission(workspace, submission_id)

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
        # add model and participant info to new samples that are not in the sequencing table yet
        samples["ModelID"] = samples["ProfileID"].apply(
            lambda x: mytracker.lookup_model_from_pr(x, "ModelID")
        )
        samples["PatientID"] = samples["ProfileID"].apply(
            lambda x: mytracker.lookup_model_from_pr(x, "PatientID")
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
    sampleset=SAMPLESETNAME,
    allbatchpairset=FPALLBATCHPAIRSETS,
    workspace=FPWORKSPACE,
    working_dir=WORKING_DIR,
    bamcolname=LEGACY_BAM_COLNAMES,
    taiga_dataset=TAIGA_FP,
    updated_mat_filename=TAIGA_FP_FILENAME,
    upload_to_taiga=True,
):
    """CCLE fingerprinting function

    Args:
        samples ([type]): [description]
        sampleset ([type]): [description]
        sid ([type]): [description]
        working_dir ([type]): [description]
        bamcolname ([type]): [description]
        workspace ([type], optional): [description]. Defaults to WORKSPACE.

    Returns:
        updated_lod_mat (pd.DataFrame): updated lod matrix
        mismatches (dict): dict representing the mismatches
        matches (dict): dict representing the matches
    """

    samples = pd.concat([rnasamples, wgssamples])

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
