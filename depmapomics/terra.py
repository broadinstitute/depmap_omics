# terra.py
import pandas as pd
import dalmatian as dm
from mgenepy import terra
import numpy as np
import traceback
import firecloud.api



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
        "datatype",
        "size",
        "bam_public_sra_path",
        "internal_bam_filepath",
        "internal_bai_filepath",
        "source",
        "hg19_bam_filepath",
        "hg19_bai_filepath",
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
                            "addUpdateAttribute": wm._process_attribute_value(  # pylint: disable=protected-access
                                i, j, reserved_attrs
                            ),
                        }
                        for i, j in row.iteritems()
                        if not np.any(pd.isnull(j))
                    ],
                }
            ]
        )

    # try rawls batch call if available
    r = dm.wmanager._batch_update_entities(  # pylint: disable=protected-access
        wm.namespace, wm.workspace, attr_list
    )
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

