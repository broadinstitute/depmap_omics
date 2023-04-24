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

def changeGSlocation(
    workspacefrom,
    newgs,
    workspaceto=None,
    prevgslist=[],
    index_func=None,
    flag_non_matching=False,
    onlysamples=[],
    bamfilepaths=[],
    entity="samples",
    droplists=True,
    keeppath=True,
    dry_run=True,
    par=20,
):
    """
  Function to move data around from one workspace to a bucket or to another workspace.

  can also work on dataframes containing lists of paths

  Args:
  -----
    workspacefrom: the workspace name where the data is
    newgs: the newgs bucket where to copy the data in
    workspaceto: if we should have these new samples and columns added to another workspace instead \
    of just updating the same one (usefull to copy one workspace to another)
    prevgslist: if providded, will only move files that are in the set of google bucket listed here
    index_func: *WIP* unused
    flag_non_matching: if set to true and prevgslist is set to some value, will return a list of samples that were not
    matched to anything in the prevgslist
    bamfilepaths: do this only on a subset of columns in terra workspace
    entity: the entity in the terra workspace on which to do this
    droplists: if set to true remove all columns containing list of paths (list of path are not uploaded well in terra)
    keeppath: if set to true, will keep the full object path and just change the bucket
    dry_run: if set to true will not update anything on Terra but just return the result
    par: on how many processor do the gs copy commands.

  Returns:
  -------
    torename: the pandas.df containing the new paths
    flaglist: the samples that were non matching (if flag_non_matching is set to true)
  """
    flaglist = []
    wmfrom = dm.WorkspaceManager(workspacefrom)
    a = wmfrom.get_entities(entity)
    if len(onlysamples) > 0:
        a = a[a.index.isin(onlysamples)]
    print("using the data from " + workspacefrom + " " + entity + " list")
    if len(a) == 0:
        raise ValueError("no " + entity)
    if bamfilepaths:
        a = a[bamfilepaths]
    todrop = set()
    torename = {}
    print(
        'this should only contains gs:// paths otherwise precise columns using "bamfilepaths"'
    )
    cmd = []
    for col in a.columns.tolist():
        val = []
        for k, prev in a[col].iteritems():
            if type(prev) is str:
                new = prev
                if newgs not in new:
                    if len(prevgslist) > 0:
                        for prevgs in prevgslist:
                            new = new.replace(prevgs, newgs)
                        if flag_non_matching:
                            if new == prev:
                                flaglist.append(prev)
                    if not keeppath:
                        new = newgs + new.split("/")[-1]
                    else:
                        new = newgs + "/".join(new.split("/")[3:])
                else:
                    print("sample " + str(k) + " was already in the new gs")
                val.append(new)
            # IN CASE WE HAVE A LIST
            if type(prev) is list:
                if droplists:
                    todrop.add(k)
                    continue
                ind = []
                for prevname in prev:
                    newname = prevname
                    if newgs not in newname:
                        if len(prevgslist) > 0:
                            for prevgs in prevgslist:
                                new = new.replace(prevgs, newgs)
                            if flag_non_matching:
                                if new == prev:
                                    flaglist.append(prev)
                        if not keeppath:
                            new = newgs + new.split("/")[-1]
                        else:
                            new = newgs + "/".join(new.split("/")[3:])
                    else:
                        print("sample " + str(k) + " was already in the new gs")
                    ind.append(newname)
                val.append(ind)
        torename.update({col: val})
        if not dry_run:
            if keeppath:
                h.parrun(
                    [
                        "gsutil mv " + a.iloc[i][col] + " " + v
                        for i, v in enumerate(val)
                    ],
                    cores=20,
                )
            else:
                gcp.mvFiles(a[col].tolist(), newgs)
        else:
            if keeppath:
                cmd = [
                    "gsutil mv " + a.iloc[i][col] + " " + v for i, v in enumerate(val)
                ]
            else:
                cmd = "mv " + str(a[col].tolist()) + " " + newgs
    torename = pd.DataFrame(
        data=torename, index=[i for i in a.index.tolist() if i != "nan"]
    )
    if workspaceto is not None:
        wmto = dm.WorkspaceManager(workspaceto)
        if not dry_run:
            wmto.disable_hound().update_entity_attributes(entity, torename)
    return torename, flaglist, cmd