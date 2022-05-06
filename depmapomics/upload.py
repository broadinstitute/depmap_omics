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
        prs[k] = pr_table[pr_table[v] <= today].index.tolist()
    return prs


def makeAchillesChoiceTable():
    return True


def makeDefaultModelTable():
    return True


def makeModelLvMatrices():
    return True


def makePRLvMatrices():
    return True
