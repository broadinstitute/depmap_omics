# tracker.py
from genepy.utils import helper as h
import pandas as pd
from depmapomics import loading
from gsheets import Sheets
from depmapomics.config import *
from genepy import terra
from genepy.google import gcp


def getTracker():
    """
    get the tracker from google sheets
    """
    return (
        Sheets.from_files(MY_ID, MYSTORAGE_ID)
        .get(REFSHEET_URL)
        .sheets[0]
        .to_frame(index_col=0)
    )


def _getDEPMAPPV(pv_index="arxspan_id", pv_tokeep=[], index="DepMap_ID"):
    """get the DEPMAP master spreadsheet from google sheet

    Args:
        pv_index (str, optional): [description]. Defaults to "arxspan_id".
        pv_tokeep (list, optional): [description]. Defaults to [].
        index (str, optional): [description]. Defaults to "DepMap_ID".

    Returns:
        (pandas.DataFrame): the DEPMAP master spreadsheet
    """
    depmap_pv = (
        Sheets.from_files(MY_ID, MYSTORAGE_ID)
        .get(DEPMAP_PV)
        .sheets[0]
        .to_frame(header=2)
    )
    depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index)
    depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index).set_index(index, drop=True)
    return depmap_pv[pv_tokeep] if pv_tokeep else depmap_pv


def merge(tracker, new, old, arxspid, cols):
    """
    given a tracker a a new and old arxspan id, will merge the two cells lines in the tracker

    Args:
        tracker (pandas.DataFrame): the tracker
        new (str): the new arxspan id
        old (str): the old arxspan id
        arxspid (str): the column name of the arxspan id
        cols (list): the columns to merge on

    Returns:
    """
    # loc = tracker[tracker[arxspid]==old].index
    return False


def findIssue(
    tracker,
    dup=[
        "age",
        "sex",
        "arxspan_id",
        "cellosaurus_id",
        "collection_site",
        "primary_disease",
        "subtype",
        "origin",
        "stripped_cell_line_name",
    ],
):
    """

    """
    print("things that are from the same patient but don't have the same value")
    dup = tracker[dup].set_index("arxspan_id").drop_duplicates()
    print(dup.loc[h.dups(dup.index)])
    print("things that have duplicate versions")
    print(
        h.dups(
            tracker["arxspan_id"]
            + "_"
            + tracker["datatype"]
            + "_"
            + tracker["version"].astype(str)
        )
    )
    print("things that don't have their legacy bam file")
    print(
        tracker[
            tracker["datatype"].isin(["rna", "wes", "wgs"])
            & (
                tracker["legacy_bam_filepath"].isna()
                | tracker["legacy_bai_filepath"].isna()
            )
        ].index
    )
    print("things that don't have their bam file path")
    print(
        tracker[
            (
                tracker["internal_bam_filepath"].isna()
                | tracker["internal_bai_filepath"].isna()
            )
        ].index
    )


def updateFromTracker(
    samples,
    ccle_refsamples,
    arxspan_id="arxspan_id",
    participant_id="participant_id",
    toupdate={},
):
    """update a list of samples' missing information from what is known in the ccle sample tracker

    Args:
        samples (pandas.DataFrame): the samples to update
        ccle_refsamples (pandas.DataFrame): the ccle sample tracker
        arxspan_id (str, optional): the name of the arxspan id column. Defaults to 'arxspan_id'.
        participant_id (str, optional): the name of the participant id column. Defaults to 'participant_id'.
        toupdate (dict(str, []), optional): the columns to update. Defaults to {}.

    Returns:
        (pandas.DataFrame): the updated samples
        (list(str)): the list of samples that were not found in the ccle sample tracker
    """
    # If I have a previous samples I can update unknown data directly
    index = []
    notfound = []
    if len(toupdate) == 0:
        toupdate = {
            "sex": [],
            "primary_disease": [],
            "cellosaurus_id": [],
            "age": [],
            "collection_site": [],
            "subtype": [],
            "subsubtype": [],
            "origin": [],
            "parent_cell_line": [],
            "matched_normal": [],
            "comments": [],
            "mediatype": [],
            "condition": [],
            "stripped_cell_line_name": [],
            "participant_id": [],
        }
    # import pdb;pdb.set_trace()
    for k, val in samples.iterrows():
        dat = ccle_refsamples[ccle_refsamples[arxspan_id] == val[arxspan_id]]
        if len(dat) > 0:
            index.append(k)
            for k, v in toupdate.items():
                toupdate[k].append(dat[k].tolist()[0])
        else:
            notfound.append(k)
    # doing so..
    for k, v in toupdate.items():
        samples.loc[index, k] = v
    len(samples.loc[notfound][participant_id]), samples.loc[notfound][
        participant_id
    ].tolist()
    return samples, notfound


def removeOlderVersions(
    names, refsamples, arxspan_id="arxspan_id", version="version", priority=None
):
    """
    will set it to your sample_ids using the latest version available for each sample

    Given a dataframe containing ids, versions, sample_ids and you dataset df indexed
    by the same ids, will set it to your sample_ids using the latest version
    available for each sample

    Args:
    -----
        refsamples (pd.df): the reference metadata. should contain [id, version, arxspan_id,...]
        names (list[str]): only do it on this set of samples.
        arxspan_id (str, optional): the name of the arxspan id column. Defaults to 'arxspan_id'.
        version (str, optional): the name of the version column. Defaults to 'version'.

    Returns:
    --------
        (dict): the subseted samples

    """
    # pandas throws an error if index is unavailable
    names = [x for x in names if x in refsamples.index.tolist()]

    lennames = len(names)
    res = {}

    refsamples = refsamples.loc[names].copy()
    if lennames > len(refsamples):
        print(set(names) - set(refsamples.index))
        import pdb

        pdb.set_trace()
        raise ValueError(
            "we had some ids in our dataset not registered in this refsample dataframe"
        )
    for arxspan in set(refsamples[arxspan_id]):
        allv = refsamples[refsamples[arxspan_id] == arxspan]
        for k, val in allv.iterrows():
            if priority is None:
                if val[version] == max(allv.version.values):
                    res[k] = arxspan
                    break
            else:
                if val[version] == max(allv.version.values):
                    res[k] = arxspan
                if val[priority] == 1:
                    res[k] = arxspan
                    break

    print("removed " + str(lennames - len(res)) + " duplicate samples")
    # remove all the reference metadata columns except the arxspan ID
    return res


def updateIsogenecity(
    di,
    tracker,
    unset=False,
    toupdate=["participant_id", "age", "sex", "matched_normal"],
):
    """

    Args:
        di ([type]): [description]
        tracker ([type]): [description]
        unset (bool, optional): [description]. Defaults to False.

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    tracker = tracker.copy()
    for k, v in di.items():
        print("________________________________")
        a = tracker[tracker.arxspan_id == k]
        b = tracker[tracker.arxspan_id == v]
        if len(a) == 0:
            print(v, "does not exist")
        if len(b) == 0:
            print(k, "does not exist")
        if len(set(a.participant_id)) > 1 or len(set(b.participant_id)) > 1:
            raise ValueError("not same participant for same cell line")
        if a.participant_id[0] == b.participant_id[0] or (
            unset and a.participant_id[0] != b.participant_id[0]
        ):
            print("already set")
            continue
        print("merging:")
        print(k, v)
        if unset:
            print("changing participant_id of ", v)
            tracker.loc[b.index, "participant_id"] = "PT-" + h.randomString()
        else:
            print("doing:")
            print(a.loc[a.index[0], toupdate].values)
            print("into")
            print(
                tracker.loc[
                    tracker[tracker.participant_id == b.participant_id[0]].index,
                    toupdate,
                ].values
            )
            tracker.loc[
                tracker[tracker.participant_id == b.participant_id[0]].index, toupdate
            ] = a.loc[a.index[0], toupdate].tolist()
    return tracker


def changeCellLineNameInNew(
    ref,
    new,
    datatype,
    dupdict,
    toupdate=[
        "stripped_cell_line_name",
        "arxspan_id",
        "patient_id",
        "sex",
        "primary_disease",
        "cellosaurus_id",
        "age",
        "collection_site",
        "subtype",
        "subsubtype",
    ],
):
    """
    Rename a/some line in a DF and takes care of corresponding metadata and versions from the sample tracker

    Args:
    -----
        new: change the cell line name in this dataframe
        dupdict: dict(tochange,newname)
        datatype: str for a ref with many datatype (to get the right version number)

    Returns:
    --------
        the updated dataframe
    """
    for k, v in dupdict.items():
        new.loc[new[new.arxspan_id == k].index, toupdate] = ref[ref.arxspan_id == v][
            toupdate
        ].values[0]
        new.loc[new[new.arxspan_id == v].index, "version"] = (
            len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
        )
    return new


def changeCellLineName(
    tracker,
    datatype,
    dupdict,
    toupdate=[
        "stripped_cell_line_name",
        "participant_id",
        "cellosaurus_id",
        "sex",
        "arxspan_id",
        "parent_cell_line",
        "matched_normal",
        "age",
        "collection_site",
        "primary_disease",
        "subtype",
        "subsubtype",
        "origin",
    ],
):
    """
    Rename a/some line in our sample tracker and takes care of corresponding metadata and versions from the sample tracker

    Args:
    -----
        dupdict: dict(tochange,newname)
        datatype: str for a tracker with many datatype (to get the right version number)

    Returns:
    --------
        the updated dataframe
    """
    for k, v in dupdict.items():
        try:
            tracker.loc[k, toupdate] = tracker[tracker.arxspan_id == v][
                toupdate
            ].values[0]
            tracker.loc[k, "version"] = len(
                tracker[(tracker.arxspan_id == v) & (tracker.datatype == datatype)]
            )
        except IndexError:
            raise IndexError(str(v) + " not found in tracker")
    return tracker


def updatePairs(
    workspaceID, tracker, removeDataFiles=True,
):
    """
    looks at the current sample tracker and updates the pairs in Terra

    It will add and remove them based on what information of match normal is available in the sample tracker. if an update happens it will remove the data files for the row.
    """
    return False


def cleanVersions(
    tracker,
    samplecol="arxspan_id",
    dryrun=False,
    datatypecol="datatype",
    versioncol="version",
):
    """
    updates and sorts the versions of samples in the sample tracker:

    checks that we get 1,2,3 instead of 2,4,5 when samples are renamed or removed

    Args:
    -----
        tracker: dataframe of the sample tracker
        samplecol: str colname of samples in the trackerr
        dryrun: bool whether or not to apply it or just print if there is issues
        datatypecol: str colname of the datatype values
        versioncol: str colname of the version values

    Returns:
    -------
        an updated sample tracker
    """
    tracker = tracker.copy()
    tracker["samval"] = tracker[samplecol] + tracker[datatypecol]
    for v in set(tracker["samval"]):
        sams = tracker[tracker["samval"] == v]
        vs = sams[versioncol].tolist()
        if max(vs) == len(vs):
            continue
        print("found issue")
        if dryrun:
            continue
        vs.sort()
        rn = {}
        for i, v in enumerate(vs):
            rn.update({v: i + 1})
        for k, val in sams.iterrows():
            tracker.loc[k, versioncol] = rn[val[versioncol]]
    tracker = tracker.drop(columns="samval")
    return tracker


def setRightName(tracker, name="stripped_cell_line_name", signs=["-", "_", ".", " "]):
    """
    cell line name need to be all uppercase and not contain special signs will take carre of ones that are not

    BE CARREFUL, it does not solve issues of duplicate lines (diff sample id, same new name)
    see findLikelyDup for that

    Args:
    -----
        tracker: dataframe of the sample tracker

    Returns:
    -----
        an updated sample tracker
    """
    new = []
    for val in tracker[name]:
        for s in signs:
            val = val.replace(s, "")
        new.append(val.upper())
    tracker[name] = new
    return tracker


def findLikelyDup(
    tracker,
    name="stripped_cell_line_name",
    signs=["-", "_", ".", " "],
    arxspid="arxspan_id",
    looksub=True,
):
    """
    find cell lines that are likely to be duplicates

    will return ,  as well,
    and

    Args:
    -----
        tracker: dataframe of the sample tracker
        looksub: bool, look if a name if within another name (can flag many derivatives)

    Returns:
    --------
        a list[tuples(str,str)] of likly duplicate names as tuples (rh13, RH-13)
        a list[tuples(str,str)] of associated arxspan ids
        a dict[str:set(str)] of arxspan ids that have multiple cell line names associated
    """
    names = set(tracker[name])
    simi = []
    arxsp = []
    issues = {}
    for i, name1 in enumerate(names):
        h.showcount(i, len(names))
        n1 = name1
        for s in signs:
            name1 = name1.replace(s, "")
        name1 = name1.upper()
        for name2 in names - set([n1]):
            n2 = name2
            for s in signs:
                name2 = name2.replace(s, "")
            name2 = name2.upper()
            if name1 == name2:
                if (
                    looksub
                    and (name1 in name2 or name2 in name1)
                    and abs(len(name1) - len(name2)) < 2
                ) or not looksub:
                    if (n1, n2) not in simi and (n2, n1) not in simi:
                        simi.append((n1, n2))
                        arxsp.append(
                            (
                                tracker[tracker[name] == n1][arxspid][0],
                                tracker[tracker[name] == n2][arxspid][0],
                            )
                        )
    for val in set(tracker[name]):
        v = set(tracker[tracker[name] == val][arxspid])
        if len(v) > 1:
            issues.update({val: v})
    return simi, arxsp, issues


def resolveIssues(tracker, issus, arxspid, cols):
    """
    given a dict of names: [arxp ids] will try to find back the right name for the right
    arxspan id by looking at their rfequncy of occurance along the tracker

    if we have rh12: [ACH-00001,ACH-0002]
    and rh12 is associated 1 time with ach-00002 and 3 with ach-00001
    and rh13 is assocated 2 time with ach-00002, then it associates:
    ach-00001 : rh12
    ach-00002 : rh13
    """
    # for val in issus:
    return False


def retrieveFromCellLineName(
    noarxspan,
    ccle_refsamples,
    datatype,
    extract={},
    my_id="~/.client_secret.json",
    stripped_cell_line_name="stripped_cell_line_name",
    arxspan_id="arxspan_id",
    mystorage_id="~/.storage.json",
    depmappvlink=DEPMAP_PV,
):
    """
    Given a List of samples with no arxspan ids, will try to retrieve an arxspan id and associated data from trackers

    For now the sample tracker and paquita's depmap pv are used.

    Args:
    -----
        noarxspan: dataframe of samples with missing arxspan ids
        ccle_refsamples: dataframe of the sample tracker
        datatype: str the datatype we are interested in (wes/rna/..) see the ones in the sample tracker
        extract: dict(str:str) see the extract in the "resolveFromWorkspace" function
        stripped_cell_line_name: str colname where the cell line name is stored
        arxspan_id: str colname wherre the sample id is stored in both noarxspan and ccle_refsamples
        depmappvlink: str the url to the depmap_pv google sheet

    Returns:
    --------
        a new dataframe with filled annotations when possible
    """
    sheets = Sheets.from_files(my_id, mystorage_id)

    # find back from cell line name in ccle ref samples
    noarxspan.arxspan_id = [
        ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i].arxspan_id[0]
        if i in ccle_refsamples[stripped_cell_line_name].tolist()
        else 0
        for i in noarxspan[arxspan_id]
    ]
    a = [
        ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i][arxspan_id][0]
        if i in ccle_refsamples[stripped_cell_line_name].tolist()
        else 0
        for i in noarxspan[stripped_cell_line_name]
    ]
    noarxspan[arxspan_id] = [
        i if i != 0 else a[e] for e, i in enumerate(noarxspan.arxspan_id)
    ]

    # get depmap pv
    depmap_pv = sheets.get(depmappvlink).sheets[0].to_frame(header=2)
    depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index)

    # find back from depmapPV
    signs = ["-", "_", ".", " "]
    for k, val in noarxspan[noarxspan[arxspan_id] == 0].iterrows():
        val = val[stripped_cell_line_name].upper()
        for s in signs:
            val = val.replace(s, "")
        a = depmap_pv[
            depmap_pv["CCLE_name"].str.contains(val)
            | depmap_pv["Stripped Cell Line Name"].str.contains(val)
            | depmap_pv["Aliases"].str.contains(val)
        ]
        if len(a) > 0:
            noarxspan.loc[k, arxspan_id] = a["DepMap_ID"].values[0]
    noarxspan[arxspan_id] = noarxspan[arxspan_id].astype(str)
    new_noarxspan = loading.resolveFromWorkspace(
        noarxspan[noarxspan[arxspan_id].str.contains("ACH-")],
        refsamples=ccle_refsamples[ccle_refsamples["datatype"] == datatype],
        match=["ACH", "CDS"],
        participantslicepos=10,
        accept_unknowntypes=True,
        extract=extract,
    )
    return pd.concat(
        [new_noarxspan, noarxspan[~noarxspan[arxspan_id].str.contains("ACH-")]]
    )


def updateSamplesSelectedForRelease(refsamples, releaseName, samples):
    """
    given a list of samples, a release name and our sample tracker, 
    
    will set these samples as 1 for this releasename and the rest at 0

    Args:
        refsamples (): of the sample tracker
    """
    refsamples[releaseName] = "0"
    refsamples.loc[samples, releaseName] = "1"
    return refsamples


def makeCCLE2(tracker, source="CCLE2"):
    """will turn the ccle sample tracker into the original ccle2 dataset based on the source column

    this means it will return a table with arxspan ids, cell line name, ...[bam file type]

    Args:
        tracker (dataframe): the sample tracker
        source (str, optional): the source column to use. Defaults to 'CCLE2'.

    Returns:
        pd.df: a table with arxspan ids, cell line name, ...[bam file type]
    """
    tracker = tracker[tracker.source == source]
    ccle = pd.DataFrame(
        index=set(tracker.arxspan_id),
        data=tracker["stripped_cell_line_name"],
        columns=["stripped_cell_line_name"],
    )
    ccle = ccle[~ccle.index.duplicated(keep="first")]
    ccle["stripped_cell_line_name"] = tracker[~tracker.index.duplicated(keep="first")][
        "stripped_cell_line_name"
    ]
    for val in ["wes", "wgs", "rna"]:
        a = tracker[tracker.datatype == val].set_index("arxspan_id")
        ccle[val + "_bam"] = a[~a.index.duplicated(keep="first")]["legacy_bam_filepath"]
        ccle[val + "_bai"] = a[~a.index.duplicated(keep="first")]["legacy_bai_filepath"]
    for val in ["hybrid_capture", "targeted", "raindance"]:
        a = tracker[tracker.datatype == val].set_index("arxspan_id")
        ccle[val + "_bam"] = a[~a.index.duplicated(keep="first")][
            "internal_bam_filepath"
        ]
        ccle[val + "_bai"] = a[~a.index.duplicated(keep="first")][
            "internal_bai_filepath"
        ]
    return ccle


def updateParentRelationFromCellosaurus(ref, cellosaurus=None):
    """
    """
    if cellosaurus is None:
        cellosaurus = h.makeCellosaurusExport()
    nol = []
    for val in set(cellosaurus.patient_id):
        mcellosaurus = set(cellosaurus[cellosaurus.patient_id == val].depmap_id) - {""}
        if len(mcellosaurus) > 0:
            pref = set(ref[ref.arxspan_id.isin(mcellosaurus)].participant_id)
            if len(pref) < 1:
                nol.append(mcellosaurus)
                continue
            elif len(pref) > 1:
                print("unmatching lines in the ref:")
                print(mcellosaurus)
            other = set()
            # finding if linked to other
            for val in pref:
                other.update(ref[ref.participant_id == val].arxspan_id.tolist())
            if len(other - mcellosaurus) > 0:
                print("we have some unmatch to the ref:")
                print(other - mcellosaurus)
            ref.loc[
                ref[ref.arxspan_id.isin(other)].index, "participant_id"
            ] = pref.pop()
    print("\n_________________\nno lines found:")
    print(nol)
    print("\n_________________________\nparents:")
    for val in set(ref.arxspan_id):
        pos = cellosaurus[cellosaurus.depmap_id == val]
        if len(pos) > 0:
            pcellosaurus = pos.parent_id[0]
            if pcellosaurus in cellosaurus.index.tolist():
                a = cellosaurus.loc[pcellosaurus, "depmap_id"]
                if a in set(ref.arxspan_id):
                    print(val, "<--", a)
                    ref.loc[ref[ref.arxspan_id == val].index, "parent_id"] = a
    return ref


def update(
    tracker,
    selected,
    samplesetname,
    samplesinset,
    lowqual,
    newgs="",
    sheetcreds=SHEETCREDS,
    sheetname=SHEETNAME,
    refworkspace=None,
    onlycol=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=True,
):
    """updates the sample tracker with the new samples and the QC metrics

    Args:
        tracker (df): [description]
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        onlycol (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to ['internal_bam_filepath', 'internal_bai_filepath'].
    """
    # updating locations of bam files and extracting infos
    if newgs and refworkspace is not None:
        res, _ = terra.changeGSlocation(
            refworkspace,
            newgs=newgs,
            onlycol=onlycol,
            entity="sample",
            keeppath=False,
            dry_run=dry_run,
            onlysamples=samplesinset,
        )
        tracker.loc[res.index.tolist()][
            ["legacy_size", "legacy_crc32c_hash"]
        ] = tracker.loc[res.index.tolist()][["size", "crc32c_hash"]].values
        tracker.loc[res.index.tolist()][HG38BAMCOL] = res[onlycol[:2]].values
        tracker.loc[res.index.tolist(), "size"] = [
            gcp.extractSize(i)[1] for i in gcp.lsFiles(res[onlycol[0]].tolist(), "-l")
        ]
        tracker.loc[res.index.tolist(), "crc32c_hash"] = [
            gcp.extractHash(i) for i in gcp.lsFiles(res[onlycol[0]].tolist(), "-L")
        ]
        tracker.loc[res.index.tolist(), "md5_hash"] = gcp.catFiles(
            dm.WorkspaceManager(refworkspace)
            .get_samples()
            .loc[samplesinset, "analysis_ready_bam_md5"]
            .tolist(),
            cut=32,
        )

    len(selected)
    tracker.loc[selected, samplesetname] = 1
    tracker.loc[samplesinset, ["low_quality", "blacklist", "prioritized"]] = 0
    tracker.loc[lowqual, "low_quality"] = 1
    tracker.loc[lowqual, "blacklist"] = 1
    if dry_run:
        return tracker
    else:
        dfToSheet(tracker, sheetname, secret=sheetcreds)
    print("updated the sheet, please reactivate protections")
