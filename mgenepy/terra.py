import dalmatian as dm
import time
from dalmatian.core import MethodNotFound
from mgenepy.utils import helper as h


def removeSamples(workspace, samples):
    """
    removes a set of samples from a workspace (very usefull when we have linked pairs and pairsets)

    Args:
    -----
      workspace: str workspace name
      samples: list of samples
    """
    wm = dm.WorkspaceManager(workspace).disable_hound()
    try:
        wm.delete_sample(samples)
    except:
        print("we had pairs.")
        pairs = wm.get_pairs()
        if len(pairs) > 0:
            pairid = pairs[pairs.case_sample.isin(samples)].index.tolist()
            for k, val in wm.get_pair_sets().iterrows():
                wm.update_pair_set(k, set(val.tolist()[0]) - set(pairid))
            wm.delete_pair(pairid)
        wm.delete_sample(samples)


def saveWorkspace(workspace, folderpath):
    """
    will save everything about a workspace into a csv and json file

    Args:
    -----
      workspace: str namespace/workspace from url typically
        namespace (str): project to which workspace belongs
        workspace (str): Workspace name
      folderpath: str path to save files
    """
    wm = dm.WorkspaceManager(workspace)
    h.createFoldersFor(folderpath)

    conf = wm.get_configs()
    for k, val in conf.iterrows():
        with open(folderpath + val["name"] + ".wdl", "w") as f:
            if val.sourceRepo == "dockstore":
                name = (
                    "dockstore.org/"
                    + "/".join(val["methodPath"].split("/")[2:4])
                    + "/"
                    + val["methodVersion"]
                )
            else:
                name = "/".join(
                    val[["methodNamespace", "methodName", "methodVersion"]]
                    .astype(str)
                    .tolist()
                )
            try:
                f.write(dm.get_wdl(name))
            except MethodNotFound:
                print(name + " could not be found")
    conf.to_csv(folderpath + "worflow_list.csv")
    params = {}
    params["GENERAL"] = wm.get_workspace_metadata()
    for k, val in conf.iterrows():
        params[k] = wm.get_config(val["name"])
        h.dictToFile(
            params[k]["inputs"], folderpath + "inputs_" + val["name"] + ".json"
        )
        h.dictToFile(params[k], folderpath + "conf_" + val["name"] + ".json")
        h.dictToFile(
            params[k]["outputs"], folderpath + "outputs_" + val["name"] + ".json"
        )
    h.dictToFile(params, folderpath + "all_configs.json")


async def waitForSubmission(workspace, submissions, raise_errors=True):
    """
    wrapper to create many submissions for a workflow

    Args:
    -----
      workspace: str namespace/workspace from url typically
      submissions: list[str] of submission ids
      raise_errors: bool to true if errors should stop your code

    Returns:
    -------
      list of ids of failed submissions
    """
    failed_submission = []
    timing = 0
    wm = dm.WorkspaceManager(workspace).disable_hound()
    assert submissions is not None
    if type(submissions) is type(""):
        submissions = [submissions]
    for scount, submission_id in enumerate(submissions):
        finished = False
        while not finished:
            done = 0
            failed = 0
            finished = True
            submission = wm.get_submission(submission_id)["workflows"]
            for _, i in enumerate(submission):
                if i["status"] not in {"Done", "Aborted", "Failed", "Succeeded"}:
                    finished = False
                elif i["status"] in {"Failed", "Aborted"}:
                    failed += 1
                    if i["workflowEntity"]["entityName"] not in failed_submission:
                        print(i["workflowEntity"]["entityName"])
                        failed_submission.append(i["workflowEntity"]["entityName"])
                elif i["status"] in {"Done", "Succeeded"}:
                    done += 1
            if not finished:
                time.sleep(40)
                print(
                    "status is: Done for "
                    + str(done)
                    + " jobs in submission "
                    + str(scount)
                    + ". "
                    + str(timing)
                    + ",5 mn elapsed.",
                    end="\r",
                )
                timing += 1
                time.sleep(20)
                print(
                    "status is: Failed for "
                    + str(failed)
                    + " jobs in submission "
                    + str(scount)
                    + ". "
                    + str(timing)
                    + " mn elapsed.",
                    end="\r",
                )
            else:
                print(
                    str(done / (done + failed))
                    + " of jobs Succeeded in submission "
                    + str(scount)
                    + "."
                )
    if len(failed_submission) > 0 and raise_errors:
        raise RuntimeError(str(len(failed_submission)) + " failed submission")
    return failed_submission
    # print and return well formated data
