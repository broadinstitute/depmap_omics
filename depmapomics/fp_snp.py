from depmapomics import constants
from depmapomics import env_config

# fingerprinting.py

import pandas as pd
import numpy as np


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
            samples[
                (samples.ModelID == u) & (ref.blacklist != 1) & (ref.profile_blist != 1)
            ].index,
            ref[
                (ref.ModelID == u)
                & (ref.index.isin(lod_mat.columns))
                & (ref.blacklist != 1)
                & (ref.profile_blist != 1)
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
                        i,
                        [
                            "ModelID",
                            "ModelCondition",
                            "version",
                            "expected_type",
                            "PatientID",
                        ],
                    ].values
                ),
                j,
                ":",
                tuple(
                    ref.loc[
                        j,
                        [
                            "ModelID",
                            "ModelCondition",
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
                                [
                                    "ModelID",
                                    "ModelCondition",
                                    "version",
                                    "expected_type",
                                    "PatientID",
                                ],
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
                                    "ModelCondition",
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
            and ref.loc[i, "profile_blist"] is not True
            and ref.loc[j, "profile_blist"] is not True
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
