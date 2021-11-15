# fingerprinting.py

import pandas as pd
import numpy as np
import dalmatian as dm
from genepy.google import gcp
from genepy import terra
from genepy.utils import helper as h
from depmapomics import tracker
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
    updated_mat_filename="updated_lod_matrix.csv"
):
    """Here we update the fingerprint LOD matrix on taiga with the new fingerprints
    and enerate matrix with LOD score for new fingerprint vcfs

    Args:
        wm (wmanager): workspace manager of current FP terra workspace
        sampleset (str): name of current sample set
        working_dir (str): directory where the lod matrix file will be saved
        batch_pair_entity (str): name of entity on terra workspace that represents sample batch pairs
        save_new_mat (bool): if new lod matrix (not combined with a previous matrix) is saved
        new_mat_filename (str): name of the new lod matrix file if save_new_mat
        prev_mat_df (pd.DatatFrame): a previous lod matrix file that will be merged with the new lod matrix
        updated_mat_filename (str): name of the updated (new merged with prev) lod matrix file
    """

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
        new_lod_mat.to_csv(working_dir + new_mat_filename + ".csv")

    new_ids = set(new_lod_mat.index)
    updated_lod_mat = new_lod_mat

    # Update LOD matrix ( have to update (A+a)*(B+b) = (AB)+(aB)+(Ab)+(ab))
    if prev_mat_df is not None:
        prev_lod_mat = prev_mat_df  
        old_ids = set(prev_lod_mat.index) - set(new_ids)
        updated_lod_mat = pd.concat(
            (prev_lod_mat.loc[old_ids, old_ids], new_lod_mat.loc[new_ids, old_ids]), axis=0
        )
        updated_lod_mat = pd.concat(
            (
                updated_lod_mat.loc[new_ids.union(old_ids), old_ids],
                new_lod_mat.transpose().loc[new_ids.union(old_ids, new_ids)],
            ),
            axis=1,
        )
        updated_lod_mat.to_csv(working_dir + updated_mat_filename + ".csv")

    return new_ids, updated_lod_mat


def checkMismatches(lod_mat, ref, samples, thr=100):
    mismatches = {}
    print("\n\nsamples that should match but don't:")
    for u in set(samples.arxspan_id):
        scores = lod_mat.loc[
            samples[samples.arxspan_id == u].index,
            ref[ref.arxspan_id == u].index.tolist(),
        ]
        for i, j in [
            (scores.index[x], scores.columns[y]) for x, y in np.argwhere(scores.values < thr)
        ]:
            print("__________________________")
            print(scores.loc[i, j])
            print(
                i,
                ":",
                tuple(
                    ref.loc[
                        i, ["arxspan_id", "version", "datatype", "participant_id"]
                    ].values
                ),
                j,
                ":",
                tuple(
                    ref.loc[
                        j,
                        [
                            "arxspan_id",
                            "version",
                            "datatype",
                            "participant_id",
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
                                ["arxspan_id", "version", "datatype", "participant_id"],
                            ].values
                        )
                    ): str(j)
                }
                + ": "
                + str(
                    tuple(
                        ref.loc[
                            j,
                            [
                                "arxspan_id",
                                "version",
                                "datatype",
                                "participant_id",
                                "blacklist",
                            ],
                        ]
                    )
                )
            )
    return mismatches


def checkMatches(lod_mat, ref, thr=500):
    print("\n\nsamples that shouldn't match but do")
    previ = ""
    matches = {}
    for i, j in [(lod_mat.index[x], lod_mat.columns[y]) for x, y in np.argwhere(lod_mat.values > thr)]:
        if i == j:
            continue
        if ref.loc[i]["participant_id"] == ref.loc[j]["participant_id"]:
            continue
        if i != previ:
            if previ != "":
                matches.update(
                    {
                        "_".join(
                            ref.loc[
                                previ,
                                [
                                    "arxspan_id",
                                    "version",
                                    "datatype",
                                    "participant_id",
                                    "stripped_cell_line_name",
                                ],
                            ]
                            .astype(str)
                            .values.tolist()
                        ): n
                    }
                )
            n = [
                tuple(
                    ref.loc[
                        j,
                        [
                            "arxspan_id",
                            "version",
                            "datatype",
                            "participant_id",
                            "stripped_cell_line_name",
                        ],
                    ].values
                )
            ]
        else:
            n.append(
                tuple(
                    ref.loc[
                        j,
                        [
                            "arxspan_id",
                            "version",
                            "datatype",
                            "participant_id",
                            "stripped_cell_line_name",
                        ],
                    ].values
                )
            )
        previ = i
    return matches


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
            raise ValueError('sample_batch_pair not manually updated')

    # replace string entries in sample_batch_pairs with references to sample_sets
    myterra.updateReferences(wm, "sample_batch_pair", pair_df)

    unique_pairs = wm.get_entities("sample_batch_pair").index.tolist()
    sample_batch_pair_set_df = pd.DataFrame(
        np.transpose(unique_pairs),
        index=["all"] * len(unique_pairs),
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
            raise ValueError('sample_batch_pair_set not manually updated')


async def fingerPrint(
    rnasamples,
    wgssamples,
    sid='id'
    sampleset=SAMPLESETNAME,
    allbatchpairset=FPALLBATCHPAIRSETS,
    workspace=FPWORKSPACE,
    working_dir=WORKING_DIR,
    bamcolname=LEGACY_BAM_COLNAMES,
    prev_mat_df=None,
    updated_mat_filename=""
):
    """1.1  Generate Fingerprint VCFs

    Here we use Dalmatian to run the fingerprint_bam_with_liftover workflow on Terra.
    This workflow calls Picard ExtractFingerprint to generate a fingerprint VCF and 
    then calls Picard LiftoverVcf to covert this vcf to hg38. To fingerprint hg38 bam files 
    just run fingerprint_bam instead.

    Args:
        samples ([type]): [description]
        sampleset ([type]): [description]
        sid ([type]): [description]
        working_dir ([type]): [description]
        bamcolname ([type]): [description]
        workspace ([type], optional): [description]. Defaults to WORKSPACE.

    Author:
        William Colgan (wcolgan@broadinstitute.org)
    """

    samples = pd.concat([rnasamples, wgssamples])
    bams = samples[bamcolname]
    bams[sid] = bams.index
    print("adding " + str(len(bams)) + " new samples to the fingerprint")
    wm = dm.WorkspaceManager(workspace).disable_hound()

    # Upload sample sheet
    samples_df = pd.DataFrame(
        bams[bamcolname + [sid, sid]].values,
        columns=["bam_filepath", "bai_filepath", "sample_id", "participant_id"],
    )
    samples_df = samples_df.set_index("sample_id")
    wm.upload_samples(samples_df, add_participant_samples=True)
    wm.update_sample_set(sampleset, samples_df.index)
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
    new_ids, updated_lod_mat = updateLOD(wm, sampleset, working_dir, save_new_mat=False, prev_mat_df=None, updated_mat_filename=updated_mat_filename)

    # finding issues with the dataset
    v = updated_lod_mat.loc[new_ids]
    ref = tracker.getTracker()
    ref = ref.append(samples)

    # find samples that should match but don't
    mismatches = checkMismatches(v, ref, samples)

    # find samples that shouldn't match but do
    matches = checkMatches(v, ref)

    return updated_lod_mat, mismatches, matches



async def _CCLEFingerPrint(
    rnasamples,
    wgssamples,
    sid='id'
    sampleset=SAMPLESETNAME,
    allbatchpairset=FPALLBATCHPAIRSETS,
    workspace=FPWORKSPACE,
    working_dir=WORKING_DIR,
    bamcolname=LEGACY_BAM_COLNAMES,
    taiga_dataset=TAIGA_FP,
    updated_mat_filename=TAIGA_FP_FILENAME
):
    """ CCLE fingerprinting function
    Args:
        samples ([type]): [description]
        sampleset ([type]): [description]
        sid ([type]): [description]
        working_dir ([type]): [description]
        bamcolname ([type]): [description]
        workspace ([type], optional): [description]. Defaults to WORKSPACE.
    """
    from taigapy import TaigaClient

    tc = TaigaClient()
    prev_lod_mat = tc.get(name=taiga_dataset, file=updated_mat_filename)
    
    # call generic function
    updated_lod_mat, should, shouldnt = fingerPrint(
        rnasamples,
        wgssamples,
        sid=sid,
        sampleset=sampleset,
        allbatchpairset=allbatchpairset,
        workspace=workspace,
        working_dir=working_dir,
        bamcolname=bamcolname,
        prev_mat_df=prev_lod_mat,
        updated_mat_filename=updated_mat_filename
    )

    # Upload updated LOD matrix to Taiga
    tc.update_dataset(
        dataset_permaname=taiga_dataset,
        changes_description="New bam fingerprints added for " + sampleset,
        upload_files=[
            {
                "path": working_dir + taiga_filename,
                "name": taiga_filename,
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            }
        ],
        add_all_existing_files=True,
    )

    return updated_lod_mat, should, shouldnt
