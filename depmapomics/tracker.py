# tracker.py
from genepy.utils import helper as h
import numpy as np
import os
import pandas as pd
from depmapomics import loading
from gsheets import Sheets
import pygsheets
from depmapomics.config import *
from depmapomics import terra as myterra
from genepy import terra
from genepy.google import gcp
from genepy.google.google_sheet import dfToSheet
import dalmatian as dm
import signal


# condense all interactions with tracker (for emeril integration)
class SampleTracker:
    """
    interacts with (read + write) the sample tracker gsheet
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(self):
        self.my_id = MY_ID
        self.mystorage_id = MYSTORAGE_ID
        self.sheetcreds = SHEETCREDS
        self.refsheet_url = REFSHEET_URL
        self.depmap_pv = DEPMAP_PV
        self.samples_not_found_url = SAMPLES_NOT_FOUND_URL
        self.samples_missing_arxspan_url = SAMPLES_MISSING_ARXSPAN_URL
        self.gumbo_url = GUMBO_SHEET
        self.gumbo_sheetname = GUMBO_SHEETNAME
        self.sheetname = SHEETNAME
        self.samples_found_name = SAMPLES_FOUND_NAME
        self.samples_not_found_name = SAMPLES_NOT_FOUND_NAME
        self.samples_missing_arxspan_name = SAMPLES_MISSING_ARXSPAN_NAME
        self.mc_table_name = MC_TABLE_NAME
        self.pr_table_name = PR_TABLE_NAME
        self.seq_table_name = SEQ_TABLE_NAME
        self.sample_table_name = SAMPLE_TABLE_NAME
        self.gc = pygsheets.authorize(service_file=SHEETCREDS)

    def read_tracker(self):
        return (
            Sheets.from_files(self.my_id, self.mystorage_id)
            .get(self.refsheet_url)
            .sheets[0]
            .to_frame(index_col=0)
        )

    def write_tracker(self, df):
        dfToSheet(df, self.sheetname, secret=self.sheetcreds)

    def read_pv(self):
        return (
            Sheets.from_files(self.my_id, self.mystorage_id)
            .get(self.depmap_pv)
            .sheets[0]
            .to_frame(header=2)
        )

    def write_all_samples_found(self, df):
        dfToSheet(df, self.samples_found_name, self.sheetcreds)

    def write_samples_not_found(self, df):
        dfToSheet(df, self.samples_not_found_name, self.sheetcreds)

    def write_samples_missing_arxspan(self, df):
        dfToSheet(df, self.samples_missing_arxspan_name, self.sheetcreds)

    def read_samples_not_found(self):
        return (
            Sheets.from_files(self.my_id, self.mystorage_id)
            .get(self.samples_not_found_url)
            .sheets[0]
            .to_frame()
            .set_index("sample_id")
        )

    def read_samples_missing_arxspan(self):
        return (
            Sheets.from_files(self.my_id, self.mystorage_id)
            .get(self.samples_missing_arxspan_url)
            .sheets[0]
            .to_frame()
            .set_index("sample_id")
        )

    def read_mc_table(self):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.mc_table_name)
        return wksht.get_as_df(index_column=1, start="A3")

    def write_mc_table(self, df):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.mc_table_name)
        wksht.set_dataframe(df, "A3", copy_index=True, nan="")

    def read_pr_table(self):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.pr_table_name)
        return wksht.get_as_df(index_column=1, start="A3")

    def write_pr_table(self, df):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.pr_table_name)
        wksht.set_dataframe(df, "A3", copy_index=True, nan="")

    def read_seq_table(self):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.seq_table_name)
        return wksht.get_as_df(index_column=1, start="A3")

    def write_seq_table(self, df):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.seq_table_name)
        wksht.set_dataframe(df, "A3", copy_index=True, nan="")

    def read_sample_table(self):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.sample_table_name)
        return wksht.get_as_df(index_column=1)

    def write_sample_table(self, df):
        sheet = self.gc.open(self.gumbo_sheetname)
        wksht = sheet.worksheet("title", self.sample_table_name)
        wksht.set_dataframe(df, "A1", copy_index=True, nan="")

    def update_pr_from_seq(
        self,
        cols={
            "bam_public_sra_path": "BamPublicSRAPath",
            "blacklist": "BlacklistOmics",
            "issue": "Issue",
            "prioritized": "Prioritized",
        },
        priority=None,
        dryrun=True,
    ):
        seq_table = self.read_seq_table()
        pr_table = self.read_pr_table()
        prs_in_seq_table = seq_table.ProfileID.unique()

        cds2pr_dict = {}
        for pr in prs_in_seq_table:
            if len(seq_table[seq_table.ProfileID == pr]) == 1:
                pr_table.loc[pr, "CDSID"] = seq_table[seq_table.ProfileID == pr].index
            else:
                allv = seq_table[seq_table["ProfileID"] == pr]
                for k, val in allv.iterrows():
                    if priority is None:
                        if val["version"] == max(allv.version.values):
                            cds2pr_dict[k] = pr
                            break
                    else:
                        if val["version"] == max(allv.version.values):
                            cds2pr_dict[k] = pr
                        if val["priority"] == 1:
                            cds2pr_dict[k] = pr
                            break
        for k, v in cds2pr_dict.items():
            pr_table.loc[v, "CDSID"] = k
        if not dryrun:
            self.write_pr_table(pr_table)
        return pr_table

    def shareCCLEbams(
        self,
        samples,
        users=[],
        groups=[],
        raise_error=True,
        arg_max_length=100000,
        bamcols=["internal_bam_filepath", "internal_bai_filepath"],
        refsheet_url="https://docs.google.com/spreadsheets/d/1Pgb5fIClGnErEqzxpU7qqX6ULpGTDjvzWwDN8XUJKIY",
        privacy_sheeturl="https://docs.google.com/spreadsheets/d/115TUgA1t_mD32SnWAGpW9OKmJ2W5WYAOs3SuSdedpX4",
        unshare=False,
        requesterpays_project="",
    ):
        """
    same as shareTerraBams but is completed to work with CCLE bams from the CCLE sample tracker

    You need to have gsheet installed and you '~/.client_secret.json', '~/.storage.json' set up

    Args:
    ----
        users: list[str] of users' google accounts
        groups: list[str] of groups' google accounts
        samples list[str] of samples cds_ids for which you want to share data
        bamcols: list[str] list of column names where bams/bais are
        raise_error: whether or not to raise an error if we find blacklisted lines
        refsheet_url: the google spreadsheet where the samples are stored
        privacy_sheeturl: the google spreadsheet where the samples are stored
        requesterpays_project: the google project where the requester pays bucket is located

    Returns:
    --------
        a list of the gs path we have been giving access to
    """
        sheets = Sheets.from_files("~/.client_secret.json", "~/.storage.json")
        print(
            "You need to have gsheet installed and you '~/.client_secret.json', '~/.storage.json' set up"
        )
        privacy = sheets.get(privacy_sheeturl).sheets[6].to_frame()
        refdata = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)
        blacklist = [i for i in privacy["blacklist"].values.tolist() if i is not np.nan]
        blacklisted = set(blacklist) & set(samples)
        print("we have " + str(len(blacklist)) + " blacklisted files")
        if len(blacklisted):
            print("these lines are blacklisted " + str(blacklisted))
            if raise_error:
                raise ValueError("blacklistedlines")
        if type(users) is str:
            users = [users]

        togiveaccess = np.ravel(refdata[bamcols].loc[samples].values)
        usrs = ""
        for group in groups:
            usrs += (
                (" -d " if unshare else " -g") + group + (":R" if not unshare else "")
            )
        for user in users:
            usrs += (
                (" -d " if unshare else " -u ") + user + (":R" if not unshare else "")
            )
        requesterpays_cmd = (
            "" if requesterpays_project == "" else "-u " + requesterpays_project
        )
        cmd_prefix = "gsutil {} -m acl ch ".format(requesterpays_cmd) + usrs
        cmd = cmd_prefix
        for n, filename in enumerate(togiveaccess):
            if type(filename) is str and filename:
                oldcmd = cmd
                cmd += " " + filename
                if (len(cmd) > arg_max_length) | (n == len(togiveaccess) - 1):
                    if n < len(togiveaccess) - 1:
                        cmd = oldcmd
                    if unshare:
                        print("preventing access to {:d} files".format(n + 1))
                    else:
                        print("granting access to {:d} files".format(n + 1))
                    with open("/tmp/grantaccess{:d}.sh".format(n), "w") as f:
                        f.write(cmd)
                    code = os.system(cmd)
                    cmd = cmd_prefix + " " + filename
                    if code == signal.SIGINT:
                        print("Awakened")
                        return None

        print("the files are stored here:\n\n" + refsheet_url)
        print("\n\njust install and use gsutil to copy them")
        print("https://cloud.google.com/storage/docs/gsutil_install")
        print("https://cloud.google.com/storage/docs/gsutil/commands/cp")
        return togiveaccess


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
        "lineage",
        "stripped_cell_line_name",
    ],
):
    """
    findIssue looks at a couple metrics:

    'things that are from the same patient but don\'t have the same value'
    'things that have duplicate versions'
    'things that don\'t have their legacy bam file'
    'things that don\'t have their bam file path'

    Args:
        tracker (pandas.DataFrame): the tracker
        dup (list, optional): the list of columns to check for duplicates. 
        Defaults to ['age', 'sex', 'arxspan_id', 'cellosaurus_id', 'primary_site', 'primary_disease',
        'subtype', 'lineage', 'stripped_cell_line_name']
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
            "lineage": [],
            "parent_cell_line": [],
            "matched_normal": [],
            "comments": [],
            # "mediatype": [],
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
    # Commented out the following because it has no effect
    # len(samples.loc[notfound][participant_id]), samples.loc[notfound][
    #    participant_id
    # ].tolist()
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

    !!!! DOES NOT YET WORK !!!! version compute is wrong
    Args:
    -----
        new: change the cell line name in this dataframe 
        dupdict: dict(tochange,newname) (arxspan_id:arxspan_id)
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
        "lineage",
    ],
):
    """
    Rename a/some line in our sample tracker and takes care of corresponding metadata and versions from the sample tracker

    Args:
    -----
        dupdict: dict(tochange,newname): the dict of the new name for the cell line: cds-id: arxspan_id
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


def update(
    table,
    selected,
    samplesetname,
    failed,
    lowqual,
    newgs="",
    refworkspace=None,
    bamfilepaths=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=True,
    samplesinset=[],
    todrop=[],
):
    """updates the sample tracker (or Gumbo omicSequencing table) with the new samples and the QC metrics
    Args:
        table (df): [description]
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to ['internal_bam_filepath', 'internal_bai_filepath'].
    """
    # updating locations of bam files and extracting infos
    if newgs and refworkspace is not None:
        if not samplesinset:
            samplesinset = [
                i["entityName"]
                for i in dm.WorkspaceManager(refworkspace)
                .get_entities("sample_set")
                .loc[samplesetname]
                .samples
            ]
        res, _ = terra.changeGSlocation(
            refworkspace,
            newgs=newgs,
            bamfilepaths=bamfilepaths,
            entity="sample",
            keeppath=False,
            dry_run=dry_run,
            onlysamples=samplesinset,
            workspaceto=refworkspace,
        )
        table.loc[res.index.tolist()][
            ["legacy_size", "legacy_crc32c_hash"]
        ] = table.loc[res.index.tolist()][["size", "crc32c_hash"]].values
        table.loc[res.index.tolist(), HG38BAMCOL] = res[bamfilepaths[:2]].values
        table.loc[res.index.tolist(), "size"] = [
            gcp.extractSize(i)[1]
            for i in gcp.lsFiles(res[bamfilepaths[0]].tolist(), "-l")
        ]
        table.loc[res.index.tolist(), "crc32c_hash"] = [
            gcp.extractHash(i) for i in gcp.lsFiles(res[bamfilepaths[0]].tolist(), "-L")
        ]
        table.loc[res.index.tolist(), "md5_hash"] = [
            gcp.extractHash(i, typ="md5")
            for i in gcp.lsFiles(res[bamfilepaths[0]].tolist(), "-L")
        ]

    table.loc[samplesinset, ["low_quality", "blacklist", "prioritized"]] = 0
    table.loc[lowqual, "low_quality"] = 1
    failed_not_dropped = list(set(failed) - set(todrop))
    # print(todrop)
    table.loc[failed_not_dropped, "blacklist"] = 1
    mytracker = SampleTracker()
    if dry_run:
        return table
    else:
        mytracker.write_seq_table(table)
    print("updated the sheet, please reactivate protections")
    return None


def updateTrackerRNA(
    selected,
    failed,
    lowqual,
    tracker,
    samplesetname,
    refworkspace=RNAWORKSPACE,
    bamfilepaths=STARBAMCOLTERRA,
    newgs=RNA_GCS_PATH_HG38,
    dry_run=False,
    qcname="star_logs",
    match=".Log.final.out",
    samplesinset=[],
    starlogs={},
    todrop=[],
):
    """updates the sample tracker with the new rna samples and the QC metrics

    Args:
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        qcname (str, optional): Terra column containing QC files. Defaults to "star_logs".
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to STARBAMCOLTERRA.
        todrop (list, optional): list of samples to be dropped. Defaults to []
        samplesinset (list[str], optional): list of samples in set when refworkspace is None (bypass interacting with terra)
        starlogs (dict(str:list[str]), optional): dict of samples' star qc log locations when refworkspace is None (bypass interacting with terra)
    """
    refwm = dm.WorkspaceManager(refworkspace)
    samplesinset = [
        i["entityName"]
        for i in refwm.get_entities("sample_set").loc[samplesetname].samples
    ]
    starlogs = myterra.getQC(
        workspace=refworkspace, only=samplesinset, qcname=qcname, match=match
    )
    for k, v in starlogs.items():
        if k == "nan":
            continue
        a = tracker.loc[k, "processing_qc"]
        a = "" if a is np.nan else a
        tracker.loc[k, "processing_qc"] = str(v) + "," + a
        if tracker.loc[k, "bam_qc"] != v[0]:
            tracker.loc[k, "bam_qc"] = v[0]
    tracker.loc[tracker[tracker.datatype.isin(["rna"])].index, samplesetname] = 0
    return update(
        tracker,
        selected,
        samplesetname,
        failed,
        lowqual,
        newgs=newgs,
        refworkspace=refworkspace,
        bamfilepaths=bamfilepaths,
        dry_run=dry_run,
        todrop=todrop,
    )


def updateTrackerWGS(
    tracker,
    selected,
    samplesetname,
    lowqual,
    datatype,
    newgs=WGS_GCS_PATH_HG38,
    samplesinset=[],
    procqc=[],
    bamqc=[],
    refworkspace=None,
    bamfilepaths=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=False,
):
    """updates the sample tracker with the new wgs samples and the QC metrics

    Args:
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        selected (list[str]): which samples were selected in the release of the analysis
        samplesetname (str): the name of the sample set or of the current analysis
        samplesinset (list[str]): list of samples in the analysis.
        lowqual (list[str]): list of samples that failed QC
        newgs (str, optional): google storage path where to move the files. Defaults to ''.
        sheetcreds (str, optional): google sheet service account file path. Defaults to SHEETCREDS.
        sheetname (str, optional): google sheet service account file path. Defaults to SHEETNAME.
        procqc (list, optional): list of Terra columns containing QC files. Defaults to [].
        bamqc (list, optional): list of Terra columns containing bam QC files. Defaults to [].
        refworkspace (str, optional): if provideed will extract workspace values (bam files path, QC,...). Defaults to None.
        bamfilepaths (list, optional): Terra columns containing the bam filepath for which to change the location. Defaults to ['internal_bam_filepath', 'internal_bai_filepath'].
    """
    # computing QC
    print("looking for QC..")
    if refworkspace is not None:
        if not samplesinset:
            samplesinset = [
                i["entityName"]
                for i in dm.WorkspaceManager(refworkspace)
                .get_entities("sample_set")
                .loc[samplesetname]
                .samples
            ]
        dataProc = (
            {}
            if procqc == []
            else myterra.getQC(workspace=refworkspace, only=samplesinset, qcname=procqc)
        )
        dataBam = (
            {}
            if bamqc == []
            else myterra.getQC(workspace=refworkspace, only=samplesinset, qcname=bamqc)
        )
        for k, v in dataProc.items():
            if k == "nan":
                continue
            a = tracker.loc[k, "processing_qc"]
            a = "" if a is np.nan else a
            tracker.loc[k, "processing_qc"] = str(v) + "," + a
        for k, v in dataBam.items():
            if k == "nan":
                continue
            a = tracker.loc[k, "bam_qc"]
            a = "" if a is np.nan else a
            tracker.loc[k, "bam_qc"] = str(v) + "," + a
    if type(datatype) is str:
        datatype = [datatype]
    tracker.loc[tracker[tracker.datatype.isin(datatype)].index, samplesetname] = 0
    update(
        tracker,
        selected,
        samplesetname,
        lowqual,
        lowqual,
        newgs=newgs,
        refworkspace=refworkspace,
        bamfilepaths=bamfilepaths,
        dry_run=dry_run,
        samplesinset=samplesinset,
    )
