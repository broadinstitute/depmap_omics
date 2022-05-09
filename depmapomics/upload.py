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
    trackerobj, one_pr_per_type=True, taiga_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME
):
    """generate and upload a table for each portal that indicates which profiles are released corresponding to which MC

    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id
        folder (str): location where tables are saved as .csv files

    Returns:
        prs (dict{(str: pd.DataFrame)}): for each portal, list of profile IDs
    """
    avail_prs = getPRToRelease(trackerobj)
    pr_table = trackerobj.read_pr_table()
    seq_table = trackerobj.read_seq_table()
    ach_tables = dict()
    for k, v in avail_prs.items():
        rows = []
        subset_pr_table = pr_table.loc[v]
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
                    latest_cds_id = seq_table.loc[cds_ids, "version"].idxmax()
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
                        latest_cds_id_wes = seq_table.loc[
                            cds_ids_wes, "version"
                        ].idxmax()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wes
                        ].index[0]
                    else:
                        latest_cds_id_wgs = seq_table.loc[
                            cds_ids_wgs, "version"
                        ].idxmax()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wgs
                        ].index[0]
                    rows.append((mc, pr, "dna"))
        ach_tables[k] = pd.DataFrame(
            rows, columns=["ModelConditionID", "ProfileID", "ProfileType"]
        )
        ach_tables[k].to_csv(folder + "/" + k + "_achilles_choice_table.csv")

    return ach_tables


def makeDefaultModelTable():
    return True


def makeModelLvMatrices():
    return True


def makePRLvMatrices():
    return True
