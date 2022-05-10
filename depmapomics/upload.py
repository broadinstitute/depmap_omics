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
    taiga_ids=VIRTUAL,
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
    source_priority = dict(
        [(source_priority[i], i) for i in range(len(source_priority))]
    )
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
                    renamed_source = subset_seq_table.source.replace(source_priority)
                    latest_cds_id = renamed_source.loc[cds_ids, "source"].idxmin()
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
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
                        renamed_source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wes = renamed_source.loc[
                            cds_ids_wes, "source"
                        ].idxmin()
                        latest_cds_id_wes = seq_table.loc[
                            cds_ids_wes, "version"
                        ].idxmax()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wes
                        ].index[0]
                    else:
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                        renamed_source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wgs = renamed_source.loc[
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
    taiga_ids=VIRTUAL,
    folder="temp/" + SAMPLESETNAME,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which modelID

    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id
        folder (str): location where tables are saved as .csv files

    Returns:
        ach_tables (dict{(str: pd.DataFrame)}): for each portal, a df containing MC-PR mapping"""
    avail_prs = getPRToRelease(trackerobj)
    pr_table = trackerobj.read_pr_table()
    mc_table = trackerobj.read_mc_table()
    seq_table = trackerobj.read_seq_table()
    default_tables = dict()
    source_priority = dict(
        [(source_priority[i], i) for i in range(len(source_priority))]
    )
    for k, v in avail_prs.items():
        rows = []
        subset_pr_table = pr_table.loc[v]
        subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
        mcs = set(subset_pr_table["ModelCondition"])
        models = set(mc_table.loc[mcs].ModelID)
        # one_pr_per_type assumes we're only picking one PR per datatype (rna/dna) for each MC
        if one_pr_per_type:
            for m in models:
                mcs_in_model = mc_table[mc_table.ModelID == m]
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
                    latest_cds_id = seq_table.loc[cds_ids, "version"].idxmax()
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
                        prs_in_model.Datatype == "wes"
                    ].CDSID.tolist()
                    pr = ""
                    # if no wgs, select broad wes over sanger wes
                    if len(cds_ids_wgs) == 0:
                        latest_cds_id_wes = seq_table.loc[
                            cds_ids_wes, "version"
                        ].idxmax()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wes
                        ].index[0]
                    # if there is wgs, always select wgs
                    else:
                        latest_cds_id_wgs = seq_table.loc[
                            cds_ids_wgs, "version"
                        ].idxmax()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wgs
                        ].index[0]
                    rows.append((mc, pr, "dna"))
        default_tables[k] = pd.DataFrame(
            rows, columns=["ModelID", "ProfileID", "ProfileType"]
        )
        default_tables[k].to_csv(folder + "/" + k + "_default_model_table.csv")
    return True


def makeModelLvMatrices():
    return True


def makePRLvMatrices():
    return True
