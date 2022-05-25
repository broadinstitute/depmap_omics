import pandas as pd

from depmapomics.config import *
from depmapomics import tracker as track
from depmapomics import upload


def test_gumbo_formatting():
    # ensures the tables in gsheet are formatted correctly
    trackerobj = track.initTracker()

    pr_table = trackerobj.read_pr_table()
    assert [
        "ModelCondition",
        "CDSID",
        "PublicReleaseDate",
        "ConsortiumReleaseDate",
        "InternalReleaseDate",
        "IBMReleaseDate",
    ] in pr_table.columns, "missing column(s) in omicsProfile table"

    mc_table = trackerobj.read_mc_table()
    assert [
        "Source",
        "ModelID",
    ] in mc_table.columns, "missing column(s) in modelCondition table"


def test_makeAchillesChoiceTable():
    trackerobj = track.initTracker()
    d = upload.makeAchillesChoiceTable(trackerobj)

    print("testing makeAchillesChoiceTable():")
    for k, v in d.items():
        assert ~d[k].isnull().values.any(), "nans found"
        assert d[k].ProfileID.is_unique(), "duplicated profile IDs found"


def test_makeDefaultModelTable():
    trackerobj = track.initTracker()
    d = upload.makeDefaultModelTable(trackerobj)

    print("testing makeDefaultModelTable():")
    for k, v in d.items():
        assert ~d[k].isnull().values.any(), "nans found"
        assert d[k].ProfileID.is_unique(), "duplicated profile IDs found"
