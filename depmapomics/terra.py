# terra.py
from gsheets import Sheets
import pandas as pd
import dalmatian as dm
from genepy import terra
import numpy as np
import traceback
import firecloud.api


def compareToCuratedGS(
    url,
    sample,
    samplesetname,
    sample_id="DepMap ID",
    clientsecret="~/.client_secret.json",
    storagepath="~/.storage.json",
    colname="CN New to internal",
    value="no data yet",
):
    """
    from a google spreadsheet, will check that we have all of the samples we should have in our sample
    set name (will parse NAME_additional for sample_id)

    Args:
    -----
        url: str the url of the gsheet
        sample: list(str) the samples to check
        samplesetname: str the name of the sampleset in the googlesheet
        sample_id: str the name of the sample_id column in the google sheet
        clientsecret: str path to your secret google api account file
        storagepath: str path to your secret google api storage file
        colname: str if we need not to include some rows from the spreadsheet that have the value value
        value: str the value for which not to include the rows

    @gmiller
    """
    sheets = Sheets.from_files(clientsecret, storagepath)
    # Cell Line Profiling Status google sheet
    gsheet = sheets.get(url).sheets[0].to_frame()
    gsheet.index = gsheet[sample_id]
    new_cn = gsheet[gsheet[colname] == samplesetname + "tent"]
    if colname and value:
        data_not_ready_cn = gsheet[gsheet[colname] == value]
        print(data_not_ready_cn)
    # these are the "new" samples discovered by our function, createDatasetsFromNewCellLines
    sample_ids = [id.split("_")[0] for id in sample]
    print("We found data for " + str(len(sorted(sample))) + " samples.\n")

    print(
        "Sanity check: Since we have the tacked on number, we should only have 1 each per sample ID:\n"
    )

    in_sheet_not_found = set(new_cn.index.tolist()) - set(sample_ids)
    if len(in_sheet_not_found) > 0:
        print(
            "We have not found "
            + str(len(in_sheet_not_found))
            + " of the samples we're supposed to \
            have this release:\n"
            + str(sorted(list(in_sheet_not_found)))
        )
    else:
        print("We aren't missing any samples that we're supposed to have this release!")


def getQC(workspace, only=[], qcname=[], match=""):
    """
    Will get from a workspace, the QC data for each samples

    Args:
    -----
        workspace: the workspace name
        only: do it only for this set of samples
        qcname: col name where the QC is in the workspace samples
        match: for example'.Log.final.out' get only that QC if you have a list of QCs in you qcname col

    Returns:
    --------
        a dict(sample_id:list[QC_filepaths])
    """
    if type(qcname) is str:
        qcname = [qcname]
    res = {}
    wm = dm.WorkspaceManager(workspace)
    sam = wm.get_samples()
    if len(only) > 0:
        sam = sam[sam.index.isin(only)]
    for k, val in sam.iterrows():
        res[k] = []
        for i in val[qcname]:
            if type(i) is list:
                if match:
                    res[k].extend([e for e in i if match in e])
                else:
                    res[k].extend(i)
            else:
                res[k].append(i)
    return res


def updatePairs(
    workspaceID, tracker, removeDataFiles=True,
):
    """
    looks at the current sample tracker and updates the pairs in Terra

    It will add and remove them based on what information of match normal is available in the sample tracker. if an update happens it will remove the data files for the row.
    """


def setupPairsFromSamples(sampless, refsamples, extract):
    """
    Given a list of samples, will compute the corresponding pairs (with nan if no matched normals)

    Args:
    -----
        sampless: pd.df samples to compute pairs for
        refsamples: pd.df samples to match to
        extract: str the name of the column to extract from the samples

    Returns:
    -------
        pairs that can be uploaded to the portal team
    """
    pairs = pd.DataFrame()
    normals = refsamples[refsamples[extract["primary_disease"]] == "normal"]
    pairs["control_sample"] = [
        "nan"
        if len(normals[normals[extract["patient_id"]] == val]) < 1
        else normals[normals[extract["patient_id"]] == val].index.tolist()[0]
        for val in sampless[extract["patient_id"]]
    ]
    pairs["case_sample"] = sampless.index.tolist()
    pairs["participant_id"] = sampless[extract["patient_id"]].tolist()
    pairs["pair_id"] = [
        val["case_sample"] + "_" + val["control_sample"] for i, val in pairs.iterrows()
    ]
    print(
        "found " + str(len(pairs["control_sample"].unique()) - 1) + " matched normals"
    )
    return pairs.set_index("pair_id")


def updateAllSampleSet(workspace, Allsample_setname="all"):
    """
    update the previous All Sample sample_set with the new samples that have been added.

    It is especially useful for the aggregate task. Can more generally merge two samplesets together

    Args:
    ----
        workspace: str namespace/workspace from url typically
        newsample_setname: str name of sampleset to add to All_samples
    """
    dm.WorkspaceManager(workspace).update_sample_set(
        Allsample_setname, dm.WorkspaceManager(workspace).get_samples().index.tolist()
    )


def copyToWorkspace(
    workspaceID,
    tracker,
    columns=[
        "arxspan_id",
        "version",
        "sm_id",
        "datatype",
        "size",
        "ccle_name",
        "stripped_cell_line_name",
        "participant_id",
        "cellosaurus_id",
        "bam_public_sra_path",
        "internal_bam_filepath",
        "internal_bai_filepath",
        "parent_cell_line",
        "sex",
        "matched_normal",
        "age",
        "collection_site",
        "primary_disease",
        "subtype",
        "subsubtype",
        "origin",
        "mediatype",
        "condition",
        "sequencing_type",
        "baits",
        "source",
        "legacy_bam_filepath",
        "legacy_bai_filepath",
    ],
    rename={},
    deleteUnmatched=False,
    addMissing=False,
    group=50,
):
    """
    will use the current sample tracker to update samples in the workspace

    it can remove samples that are not in the tracker.

    Args:
    ----
        workspaceID: str the workspace id
        tracker: dataframe the sample tracker
        columns: list[str] the columns to sync
        rename: dict(str:str) columns to rename from sample tracker to workspace
        deleteUnmatched: bool whether or not to delete samples in the workspace and not in the sample tracker
    """
    wm = dm.WorkspaceManager(workspaceID).disable_hound()
    sam = wm.get_samples()
    track = tracker[tracker.index.isin(sam.index)][columns].rename(columns=rename)
    track.index.name = "sample_id"
    if len(track) == 0 and not addMissing:
        raise ValueError("wrong tracker or index non matching")
    unmatched = set(sam.index) - (set(tracker.index) | set(["nan"]))
    if not addMissing:
        print("found these to be unmatched in the tracker: " + str(unmatched))
        if deleteUnmatched and len(unmatched) > 0:
            terra.removeSamples(workspaceID, unmatched)
    unmatched = set(tracker.index) - set(sam.index)
    if len(track) != 0:
        for i in range(0, len(track), group):
            wm.update_sample_attributes(track.iloc[i : i + group])
    if addMissing and len(unmatched) > 0:
        print("found these columns to be missing in workspace: " + str(unmatched))
        track = tracker[tracker.index.isin(unmatched)][columns].rename(columns=rename)
        wm.upload_samples(track)


def updateReferences(wm, etype, attrs):
    """written for FP, where we need to update the sample_batch_pair data table

    where entries are references to sample_sets instead of strings
    
    Args:
    ----
        wm (dm.WorkspaceManager): for the workspace to be updated
        etype (str): entity type to be updated
        attrs (df.DataFrame): updated dataframe for this entity type"""

    reserved_attrs = {}
    if etype == "sample":
        reserved_attrs = {"participant": "participant"}
    elif etype == "pair":
        reserved_attrs = {
            "participant": "participant",
            "case_sample": "sample",
            "control_sample": "sample",
        }
    elif etype == "sample_batch_pair":
        reserved_attrs = {
            "sample_batch_a": "sample_set",
            "sample_batch_b": "sample_set",
        }

    attr_list = []
    for entity, row in attrs.iterrows():
        attr_list.extend(
            [
                {
                    "name": entity,
                    "entityType": etype,
                    "operations": [
                        {
                            "op": "AddUpdateAttribute",
                            "attributeName": i,
                            "addUpdateAttribute": wm._process_attribute_value(
                                i, j, reserved_attrs
                            ),  # pylint: disable=protected-access
                        }
                        for i, j in row.iteritems()
                        if not np.any(pd.isnull(j))
                    ],
                }
            ]
        )

    # try rawls batch call if available
    r = dm.wmanager._batch_update_entities(
        wm.namespace, wm.workspace, attr_list
    )  # pylint: disable=protected-access
    try:
        if r.status_code == 204:
            if isinstance(attrs, pd.DataFrame):
                print(
                    "Successfully updated attributes '{}' for {} {}s.".format(
                        attrs.columns.tolist(), attrs.shape[0], etype
                    )
                )
            elif isinstance(attrs, pd.Series):
                print(
                    "Successfully updated attribute '{}' for {} {}s.".format(
                        attrs.name, len(attrs), etype
                    )
                )
            else:
                print(
                    "Successfully updated attribute '{}' for {} {}s.".format(
                        attrs.name, len(attrs), etype
                    )
                )
        elif r.status_code >= 400:
            raise BlockingIOError("Unable to update entity attributes")
        else:
            print(r.text)
    except:  # revert to public API
        traceback.print_exc()
        print("Failed to use batch update endpoint; switching to slower fallback")
        for update in attr_list:
            r = firecloud.api.update_entity(
                wm.namespace,
                wm.workspace,
                update["entityType"],
                update["name"],
                update["operations"],
            )
            if r.status_code == 200:
                print("Successfully updated {}.".format(update["name"]))
            elif r.status_code >= 400:
                raise BlockingIOError("Unable to update entity attributes")
            else:
                print(r.text)

