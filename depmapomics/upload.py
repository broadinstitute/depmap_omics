from __future__ import print_function
import pandas as pd
import os
from datetime import date

from genepy.utils import helper as h
from depmapomics import tracker
from depmapomics.config import *
from taigapy import TaigaClient


def getPRToRelease(trackerobj):
    """generate lists of profiles to release based on date for all portals
    
    Args:
        trackerobj (SampleTracker): tracker object

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    date_col_dict = {
        "internal": "InternalReleaseDate",
        "ibm": "IBMReleaseDate",
        "dmc": "ConsortiumReleaseDate",
        "public": "PublicReleaseDate",
    }
    today = int(str(date.today()).replace("-", ""))
    pr_table = trackerobj.read_pr_table()
    prs = dict()
    for k, v in date_col_dict.items():
        prs_with_date = pr_table[~(pr_table[v] == "")]
        prs[k] = prs_with_date[prs_with_date[v].astype(int) <= today].index.tolist()
    return prs


def makeAchillesChoiceTable(
    trackerobj,
    one_pr_per_type=True,
    source_priority=SOURCE_PRIORITY,
    folder="temp/" + SAMPLESETNAME,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which MC

    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id
        folder (str): location where tables are saved as .csv files

    Returns:
        ach_tables (dict{(str: pd.DataFrame)}): for each portal, a df containing MC-PR mapping
    """
    avail_prs = getPRToRelease(trackerobj)
    pr_table = trackerobj.read_pr_table()
    seq_table = trackerobj.read_seq_table()
    ach_tables = dict()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    for k, v in avail_prs.items():
        rows = []
        subset_pr_table = pr_table.loc[v]
        subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
        mcs = set(subset_pr_table["ModelCondition"])
        # one_pr_per_type assumes we're only picking one PR per datatype (rna/dna) for each MC
        if one_pr_per_type:
            for mc in mcs:
                prs_in_mc = subset_pr_table[(subset_pr_table.ModelCondition == mc)]
                # rna
                if len(prs_in_mc[prs_in_mc.Datatype == "rna"]) == 1:
                    pr = prs_in_mc[prs_in_mc.Datatype == "rna"].index[0]
                    rows.append((mc, pr, "rna"))
                elif len(prs_in_mc[prs_in_mc.Datatype == "rna"]) > 1:
                    cds_ids = prs_in_mc[prs_in_mc.Datatype == "rna"].CDSID.tolist()
                    # at this point it is guaranteed that all cds_ids have different sources
                    subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
                    subset_seq_table.source = subset_seq_table.source.replace(
                        source_priority
                    )
                    latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
                    pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[
                        0
                    ]
                    rows.append((mc, pr, "rna"))
                # dna
                if (
                    len(
                        prs_in_mc[
                            (prs_in_mc.Datatype == "wgs")
                            | (prs_in_mc.Datatype == "wes")
                        ]
                    )
                    == 1
                ):
                    pr = prs_in_mc[
                        (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
                    ].index[0]
                    rows.append((mc, pr, "dna"))
                elif (
                    len(
                        prs_in_mc[
                            (prs_in_mc.Datatype == "wgs")
                            | (prs_in_mc.Datatype == "wes")
                        ]
                    )
                    > 1
                ):
                    cds_ids_wgs = prs_in_mc[prs_in_mc.Datatype == "wgs"].CDSID.tolist()
                    cds_ids_wes = prs_in_mc[prs_in_mc.Datatype == "wes"].CDSID.tolist()
                    pr = ""
                    if len(cds_ids_wgs) == 0:
                        if len(prs_in_mc[prs_in_mc.Datatype == "wes"]) == 1:
                            pr = prs_in_mc[prs_in_mc.Datatype == "wes"].index[0]
                        else:
                            subset_seq_table = seq_table[
                                seq_table.index.isin(cds_ids_wes)
                            ]
                            subset_seq_table.source = subset_seq_table.source.replace(
                                source_priority
                            )
                            latest_cds_id_wes = subset_seq_table.loc[
                                cds_ids_wes, "source"
                            ].idxmin()
                            pr = subset_pr_table[
                                subset_pr_table.CDSID == latest_cds_id_wes
                            ].index[0]
                    else:
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                        subset_seq_table.source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wgs = subset_seq_table.loc[
                            cds_ids_wgs, "source"
                        ].idxmin()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wgs
                        ].index[0]
                    rows.append((mc, pr, "dna"))
        ach_tables[k] = pd.DataFrame(
            rows, columns=["ModelConditionID", "ProfileID", "ProfileType"]
        )
        ach_tables[k].to_csv(folder + "/" + k + "_achilles_choice_table.csv")

    return ach_tables


def makeDefaultModelTable(
    trackerobj,
    one_pr_per_type=True,
    source_priority=SOURCE_PRIORITY,
    folder="temp/" + SAMPLESETNAME,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which modelID

    Args:
        trackerobj (SampleTracker): tracker object
        folder (str): location where tables are saved as .csv files

    Returns:
        ach_tables (dict{(str: pd.DataFrame)}): for each portal, a df containing MC-PR mapping"""
    avail_prs = getPRToRelease(trackerobj)
    pr_table = trackerobj.read_pr_table()
    mc_table = trackerobj.read_mc_table()
    seq_table = trackerobj.read_seq_table()
    default_tables = dict()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    for k, v in avail_prs.items():
        rows = []
        subset_pr_table = pr_table.loc[v]
        subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
        mcs = set(subset_pr_table["ModelCondition"])
        models = set(mc_table.loc[mcs].ModelID)
        # one_pr_per_type assumes we're only picking one PR per datatype (rna/dna) for each MC
        if one_pr_per_type:
            for m in models:
                subset_mc_table = mc_table[mc_table.ModelID == m]
                mcs_in_model = subset_mc_table.index.tolist()
                prs_in_model = subset_pr_table[
                    (subset_pr_table.ModelCondition.isin(mcs_in_model))
                ]
                # rna
                if len(prs_in_model[prs_in_model.Datatype == "rna"]) == 1:
                    pr = prs_in_model[prs_in_model.Datatype == "rna"].index[0]
                    rows.append((m, pr, "rna"))
                elif len(prs_in_model[prs_in_model.Datatype == "rna"]) > 1:
                    cds_ids = prs_in_model[
                        prs_in_model.Datatype == "rna"
                    ].CDSID.tolist()
                    subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
                    subset_seq_table.source = subset_seq_table.source.replace(
                        source_priority
                    )
                    latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
                    pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[
                        0
                    ]
                    rows.append((m, pr, "rna"))
                # dna
                if (
                    len(
                        prs_in_model[
                            (prs_in_model.Datatype == "wgs")
                            | (prs_in_model.Datatype == "wes")
                        ]
                    )
                    == 1
                ):
                    pr = prs_in_model[
                        (prs_in_model.Datatype == "wgs")
                        | (prs_in_model.Datatype == "wes")
                    ].index[0]
                    rows.append((m, pr, "dna"))
                elif (
                    len(
                        prs_in_model[
                            (prs_in_model.Datatype == "wgs")
                            | (prs_in_model.Datatype == "wes")
                        ]
                    )
                    > 1
                ):
                    # assuming SANGER doesn't have wgs
                    cds_ids_wgs = prs_in_model[
                        prs_in_model.Datatype == "wgs"
                    ].CDSID.tolist()
                    cds_ids_wes = prs_in_model[
                        (prs_in_model.Datatype == "wes") & (prs_in_model.CDSID != "")
                    ].CDSID.tolist()  # CDSID is '' when the profile is in legacy
                    pr = ""
                    # if no wgs, look at MC table and select broad wes over sanger wes
                    if len(cds_ids_wgs) == 0:
                        if len(prs_in_model[prs_in_model.Datatype == "wes"]) == 1:
                            pr = prs_in_model[prs_in_model.Datatype == "wes"].index[0]
                        else:
                            subset_seq_table = seq_table[
                                seq_table.index.isin(cds_ids_wes)
                            ]
                            subset_seq_table.source = subset_seq_table.source.replace(
                                source_priority
                            )
                            latest_cds_id_wes = subset_seq_table.loc[
                                cds_ids_wes, "source"
                            ].idxmin()
                            pr = subset_pr_table[
                                subset_pr_table.CDSID == latest_cds_id_wes
                            ].index[0]
                    # if there is wgs, always select wgs
                    else:
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                        subset_seq_table.source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wgs = subset_seq_table.loc[
                            cds_ids_wgs, "source"
                        ].idxmin()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wgs
                        ].index[0]
                    rows.append((m, pr, "dna"))
        default_tables[k] = pd.DataFrame(
            rows, columns=["ModelID", "ProfileID", "ProfileType"]
        )
        default_tables[k].to_csv(folder + "/" + k + "_default_model_table.csv")
    return default_tables


def initVirtualDatasets(
    samplesetname=SAMPLESETNAME, taiga_folder_id=VIRTUAL_FOLDER, portals=DATASETS
):
    """initialize taiga virtual datasets for all portals by uploading an empty dummy file
    """

    with open("temp/dummy.csv", "w") as fp:
        pass
    virtual = dict()
    tc = TaigaClient()
    for p in portals:
        virtual[p] = tc.create_dataset(
            p + "_" + samplesetname,
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "temp/dummy.csv",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
    return virtual


def uploadCNMatrices(renaming_dict, taiga_id="", folder="temp/" + SAMPLESETNAME):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """

    # load cds-id indexed matrices for the current quarter
    genecn = pd.read_csv(folder + "/achilles_gene_cn.csv", index_col=0)
    segmentcn = pd.read_csv(folder + "/achilles_segment.csv")
    wescn = pd.read_csv(folder + "/wes_genecn_latest.csv", index_col=0)
    wessegment = pd.read_csv(folder + "/wes_segments_latest.csv")

    # subset and rename
    segmentcn_renamed = segmentcn[
        segmentcn.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    segmentcn_renamed.to_csv("temp/all_merged_segments.csv", index=False)
    genecn_renamed = genecn[genecn.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    genecn_renamed.to_csv("temp/all_merged_genes_cn.csv")
    wescn_renamed = wescn[wescn.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wescn_renamed.to_csv("temp/wes_genes_cn.csv")
    wessegment_renamed = wessegment[
        wessegment.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    wessegment_renamed.to_csv("temp/wes_segments.csv", index=False)


def makePRLvMatrices(trackerobj, taiga_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(trackerobj)
    pr_table = trackerobj.read_pr_table()
    for portal, prs_to_release in prs_allportals.items():
        subset_pr_table = pr_table[pr_table.index.isin(prs_to_release)]
        renaming_dict = dict(list(zip(subset_pr_table.CDSID, subset_pr_table.index)))
        uploadCNMatrices(renaming_dict, taiga_id=taiga_ids[portal], folder=folder)


def makeModelLvMatrices():

    return True

