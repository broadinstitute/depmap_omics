from __future__ import print_function
import pandas as pd
import os
from datetime import date

from genepy.utils import helper as h
from depmapomics import tracker
from depmapomics.config import *
from taigapy import TaigaClient


def getPRToRelease(trackerobj, date_col_dict=DATE_COL_DICT):
    """generate lists of profiles to release based on date for all portals
    
    Args:
        trackerobj (SampleTracker): tracker object

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    today = int(str(date.today()).replace("-", ""))
    pr_table = trackerobj.read_pr_table()
    prs = dict()
    for k, v in date_col_dict.items():
        prs_with_date = pr_table[~(pr_table[v] == "")]
        prs[k] = prs_with_date[
            (prs_with_date[v].astype(int) <= today)
            & (prs_with_date.ProfileSource == "bam")
        ].index.tolist()
    return prs


def makeAchillesChoiceTable(
    trackerobj,
    prs,
    one_pr_per_type=True,
    source_priority=SOURCE_PRIORITY,
    colnames=ACH_CHOICE_TABLE_COLS,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which MC

    Args:
        trackerobj (SampleTracker): tracker object
        prs (list): list of profile IDs to be released
        one_pr_per_type (bool, optional): whether to enforce including only one profile type per MC
        source_priority (list, optional): ordered list of how different data sources should be prioritized

    Returns:
        ach_table (pd.DataFrame): a df containing MC-PR mapping
    """
    pr_table = trackerobj.read_pr_table()
    seq_table = trackerobj.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
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
                pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[0]
                rows.append((mc, pr, "rna"))
            # dna
            if (
                len(
                    prs_in_mc[
                        (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
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
                        (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
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
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
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
    ach_table = pd.DataFrame(rows, columns=colnames)

    return ach_table


def makeDefaultModelTable(
    trackerobj,
    prs,
    one_pr_per_type=True,
    source_priority=SOURCE_PRIORITY,
    colnames=DEFAULT_TABLE_COLS,
):
    """generate a table that indicates which profiles are released corresponding to which modelID

    Args:
        trackerobj (SampleTracker): tracker object
        prs (list): list of profile IDs to be released
        one_pr_per_type (bool, optional): whether to enforce including only one profile type per model
        source_priority (list, optional): ordered list of how different data sources should be prioritized

    Returns:
        default_table (pd.DataFrame): a df containing Model-PR mapping
    """
    pr_table = trackerobj.read_pr_table()
    mc_table = trackerobj.read_mc_table()
    seq_table = trackerobj.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
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
                cds_ids = prs_in_model[prs_in_model.Datatype == "rna"].CDSID.tolist()
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[0]
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
                    (prs_in_model.Datatype == "wgs") | (prs_in_model.Datatype == "wes")
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
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
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
    default_table = pd.DataFrame(rows, columns=colnames)
    return default_table


def initVirtualDatasets(
    samplesetname=SAMPLESETNAME, taiga_folder_id=VIRTUAL_FOLDER, portals=DATASETS
):
    """initialize both PR- and Model-level taiga virtual datasets for all 4 portals by uploading an empty dummy file
    """
    virutal_pr = dict()
    virtual_model = dict()
    tc = TaigaClient()
    for p in portals:
        virutal_pr[p] = tc.create_dataset(
            p + "_" + samplesetname + "_profile",
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "/dev/null",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
        virtual_model[p] = tc.create_dataset(
            p + "_" + samplesetname + "_model",
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "/dev/null",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
    return virutal_pr, virtual_model


def uploadPRMatrix(
    prs,
    taiga_latest,
    taiga_virtual,
    latest_fn,
    virtual_fn,
    matrix_format,
    pr_col="index",
    folder="temp",
    change_desc="",
):
    """subset, save and upload to taiga PR-level matrix

    Args:
        prs (list): list of PR-ids to release
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
        latest_fn (str): file name on taiga latest
        virtual_fn (str): file name on taiga virtual
        matrix_format (str): which format this matrix should be uploaded in (NumericMatrixCSV, TableCSV, etc)
        pr_col (str): which column in the matrix contains pr ids
        folder (str): where the file should be stores before uploading to virtual
        change_desc (str): change description on taiga virtual
    """
    print("loading ", latest_fn, " from latest")
    to_subset = tc.get(name=taiga_latest, file=latest_fn)

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(prs)]
        subset_mat.to_csv(folder + virtual_fn + ".csv")
    else:
        subset_mat = subset_mat[subset_mat[pr_col].isin(prs)]
        subset_mat.to_csv(folder + virtual_fn + ".csv", index=False)

    print("uploading ", virtual_fn, " to virtual")
    tc = TaigaClient()
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description=change_desc,
        upload_files=[
            {
                "path": folder + virtual_fn + ".csv",
                "name": virtual_fn,
                "format": matrix_format,
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadModelMatrix(
    pr2model_dict,
    taiga_latest,
    taiga_virtual,
    latest_fn,
    virtual_fn,
    matrix_format,
    pr_col="index",
    folder="temp",
    change_desc="",
):
    """subset, rename, save and upload to taiga model-level matrix

    Args:
        pr2model_dict (dict): dictionary mapping profile ids to model ids
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
        latest_fn (str): file name on taiga latest
        virtual_fn (str): file name on taiga virtual
        matrix_format (str): which format this matrix should be uploaded in (NumericMatrixCSV, TableCSV, etc)
        pr_col (str): which column in the matrix contains pr ids
        folder (str): where the file should be stores before uploading to virtual
        change_desc (str): change description on taiga virtual
    """
    print("loading ", latest_fn, " from latest")
    to_subset = tc.get(name=taiga_latest, file=latest_fn)

    print("subsetting ", latest_fn)
    if pr_col == "index":
        subset_mat = to_subset[to_subset.index.isin(set(pr2model_dict.keys()))].rename(
            index=pr2model_dict
        )
        subset_mat.to_csv(folder + virtual_fn + ".csv")
    else:
        subset_mat = subset_mat[
            subset_mat[pr_col].isin(set(pr2model_dict.keys()))
        ].replace({"DepMap_ID": pr2model_dict})
        subset_mat.to_csv(folder + virtual_fn + ".csv", index=False)

    print("uploading ", virtual_fn, " to virtual")
    tc = TaigaClient()
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description=change_desc,
        upload_files=[
            {
                "path": folder + virtual_fn + ".csv",
                "name": virtual_fn,
                "format": matrix_format,
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadGermlineMatrixModel(
    pr2model_dict, portal, taiga_latest=TAIGA_CN, taiga_virtual="",
):
    """subset, rename, save and upload to taiga germline binary matrix
    
    Args:
        pr2model_dict (dict): renaming scheme mapping from PR-id to model id
        portal (str): which portal the data is being uploaded to
        taiga_latest (str): which dataset the matrices to be subsetted are being read from
        taiga_virtual (str): which dataset the matrices are being uploaded to
    """
    folder = "temp/" + portal + "/model/"
    h.createFoldersFor(folder)
    # load cds-id indexed matrices for the current quarter
    print("Germline matrix: loading from taiga latest")
    tc = TaigaClient()
    germline = tc.get(name=taiga_latest, file="merged_binary_germline")

    # subset and rename
    print("Germline matrix: subsetting and renaming")
    whitelist = [x for x in germline.columns if x in pr2model_dict]
    whitelist_germline = germline[whitelist]
    whitelist_germline = whitelist_germline.rename(columns=pr2model_dict)
    whitelist_germline = whitelist_germline.astype(bool).astype(int)
    sorted_mat = germline.iloc[:, :4].join(whitelist_germline)
    sorted_mat["end"] = sorted_mat["end"].astype(int)
    sorted_mat.to_csv(folder + "merged_binary_germline.csv", index=False)

    # upload to taiga
    print("Germline matrix: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding model-level germline matrix",
        upload_files=[
            {
                "path": folder + "merged_binary_germline.csv",
                "name": "germline_mutation",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadAuxTables(trackerobj, taiga_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME):
    """upload achilles choice and default model table to all portals
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict, optional): dict mapping portal name to taiga virtual dataset id
        folder (str, optional): where the tables are saved
    """
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs in prs_allportals.items():
        achilles_table = makeAchillesChoiceTable(trackerobj, prs)
        default_table = makeDefaultModelTable(trackerobj, prs)
        achilles_table.to_csv(
            folder + portal + "_achilles_choice_table.csv", index=False
        )
        default_table.to_csv(folder + portal + "_default_model_table.csv", index=False)
        tc = TaigaClient()
        tc.update_dataset(
            dataset_id=taiga_ids[portal],
            changes_description="adding mapping tables",
            upload_files=[
                {
                    "path": folder + "/" + portal + "_achilles_choice_table.csv",
                    "name": "Achilles_choice_table",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "/" + portal + "_default_model_table.csv",
                    "name": "default_model_table",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def makePRLvMatrices(
    trackerobj, virtual_ids=VIRTUAL,
):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs_to_release in prs_allportals.items():
        print("uploading profile-level matrices to ", portal)
        for latest_id, fn_dict in LATEST2FN_NUMMAT.items():
            for latest, virtual in fn_dict:
                uploadPRMatrix(
                    prs_to_release,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col="index",
                    folder="temp",
                    change_desc="adding " + virtual,
                )
        for latest_id, fn_dict in LATEST2FN_TABLE.items():
            for latest, virtual in fn_dict:
                uploadPRMatrix(
                    prs_to_release,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "TableCSV",
                    pr_col=SAMPLEID,
                    folder="temp",
                    change_desc="adding " + virtual,
                )
        if portal == "internal":
            for latest, virtual in VIRTUAL_FILENAMES_NUMMAT_EXP_INTERNAL:
                uploadPRMatrix(
                    prs_to_release,
                    TAIGA_EXPRESSION,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col=SAMPLEID,
                    folder="temp",
                    change_desc="adding " + virtual,
                )


def makeModelLvMatrices(
    trackerobj, virtual_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME
):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs_to_release in prs_allportals.items():
        default_table = makeDefaultModelTable(trackerobj, prs_to_release)
        pr2model_dict = dict(list(zip(default_table.ProfileID, default_table.ModelID)))
        h.dictToFile(pr2model_dict, folder + "/" + portal + "_pr2model_renaming.json")
        print("uploading model-level matrices to", portal)
        for latest_id, fn_dict in LATEST2FN_NUMMAT.items():
            for latest, virtual in fn_dict:
                uploadModelMatrix(
                    pr2model_dict,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col="index",
                    folder="temp",
                    change_desc="adding " + virtual,
                )
        for latest_id, fn_dict in LATEST2FN_TABLE.items():
            for latest, virtual in fn_dict:
                uploadModelMatrix(
                    pr2model_dict,
                    latest_id,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "TableCSV",
                    pr_col=SAMPLEID,
                    folder="temp",
                    change_desc="adding " + virtual,
                )
        if portal == "internal":
            for latest, virtual in VIRTUAL_FILENAMES_NUMMAT_EXP_INTERNAL:
                uploadModelMatrix(
                    pr2model_dict,
                    TAIGA_EXPRESSION,
                    virtual_ids[portal],
                    latest,
                    virtual,
                    "NumericMatrixCSV",
                    pr_col=SAMPLEID,
                    folder="temp",
                    change_desc="adding " + virtual,
                )
        uploadGermlineMatrixModel(
            pr2model_dict, portal, taiga_virtual=taiga_ids[portal]
        )


def findLatestVersion(dataset, approved_only=True):
    highest = 0
    latest_version = 0
    tc = TaigaClient()
    data = tc.get_dataset_metadata(dataset)
    for val in data["versions"]:
        if val["state"] == "approved" or not approved_only:
            if int(val["name"]) > highest:
                highest = int(val["name"])
                latest_version = highest
    if latest_version == 0:
        raise ValueError("could not find a version")
    return data["permanames"][0] + "." + str(latest_version)


def updateEternal(
    eternal_id="depmap-a0ab", virtual=VIRTUAL, samplesetname=SAMPLESETNAME
):
    """update taiga eternal dataset by linking to latest virtual internal dataset"""
    latest_version = findLatestVersion(virtual["internal"])

    files = [
        VIRTUAL_FILENAMES_NUMMAT_EXP.values()
        + VIRTUAL_FILENAMES_NUMMAT_EXP_INTERNAL.values()
        + VIRTUAL_FILENAMES_NUMMAT_CN.values()
        + VIRTUAL_FILENAMES_NUMMAT_MUT.values()
        + VIRTUAL_FILENAMES_GERMLINE.values()
        + VIRTUAL_FILENAMES_TABLE_FUSION.values()
        + VIRTUAL_FILENAMES_TABLE_CN.values()
        + VIRTUAL_FILENAMES_TABLE_MUT.values()
    ]

    tc = TaigaClient()
    tc.update_dataset(
        eternal_id,
        changes_description="new " + samplesetname + " omics dataset.",
        add_taiga_ids=[
            {"taiga_id": latest_version + "/" + file, "name": file} for file in files
        ],
        add_all_existing_files=True,
    )


def CCLEupload(trackerobj, taiga_ids=""):
    if taiga_ids == "":
        taiga_ids = initVirtualDatasets()

    makePRLvMatrices(trackerobj, taiga_ids=taiga_ids)
    makeModelLvMatrices(trackerobj, taiga_ids=taiga_ids)
    uploadAuxTables(trackerobj, taiga_ids=taiga_ids)
    updateEternal(virtual=taiga_ids)
