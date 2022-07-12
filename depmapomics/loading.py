# -*- coding: utf-8 -*-
# Jérémie Kalfon
# for BroadInsitute
# in 2019

import pandas as pd
import numpy as np
import dalmatian as dm
from depmapomics import tracker as track
from depmapomics.config import *
from depmapomics import terra as myterra
from genepy import terra
from genepy import sequencing as seq
from genepy.utils import helper as h
from genepy.google import gcp

#####################
# Loading Functions
#####################


def loadFromMultipleWorkspaces(
    wsnames,
    wsidcol,
    gumboidcol,
    stype,
    only_load_mapped=False,
    bamcol="cram_or_bam_path",
    load_undefined=False,
    extract=EXTRACT_DEFAULTS,
    accept_unknowntypes=True,
    addonly=[],
    maxage=MAXAGE,
):
    """
    Load new RNAseq samples from multiple terra workspace and attempt to map them to existing profiles in gumbo.

    Args:
    -----
        wsnames (dict): dictionary mapping source to workspace name ({source: ws name})
        wsidcol (str): name of column in wsname's sample table that contains IDs to map to ProfileIDs by
        gumboidcol (str): name of column in gumbo's profile table that contains IDs to map to ProfileIDs by
        stype (str): type of the data (wgs, rna, etc.)
        extract: if you want to specify what values should refer to which column names

    Returns:
    --------
        samples: pd dataframe the filtered sample list
        unmapped_new_lines: list of lines that don't have PR ids assigned in gumbo
    """
    samples = []
    unmapped_samples = []
    for s, wsname in wsnames.items():
        print("loading " + stype + " samples from terra workspace: " + wsname)
        samples_per_ws, unmapped = loadFromTerraWorkspace(
            wsname,
            wsidcol,
            gumboidcol,
            s,
            stype,
            only_load_mapped,
            bamcol,
            load_undefined,
            extract,
            accept_unknowntypes,
            addonly,
        )
        samples.append(samples_per_ws)
        unmapped_samples.extend(unmapped)
    return pd.concat(samples), unmapped_samples


def loadFromTerraWorkspace(
    wsname,
    wsidcol,
    gumboidcol,
    source,
    stype,
    only_load_mapped=False,
    bamcol="cram_or_bam_path",
    load_undefined=False,
    extract=EXTRACT_DEFAULTS,
    accept_unknowntypes=True,
    addonly=[],
    maxage=MAXAGE,
):
    """
    Load new samples from a terra workspace and attempt to map them to existing profiles in gumbo.

    Args:
    -----
        wsname (str): name of terra workspace to load data from
        wsidcol (str): name of column in wsname's sample table that contains IDs to map to ProfileIDs by
        gumboidcol (str): name of column in gumbo's profile table that contains IDs to map to ProfileIDs by
        source (str): source of the data delivered (DEPMAP, IBM, etc.)
        stype (str): type of the data (wgs, rna, etc.)
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':

    Returns:
    --------
        samples: pd dataframe the filtered sample list
        unmapped_new_lines: list of lines that don't have PR ids assigned in gumbo
    """
    mytracker = track.SampleTracker()
    seq_table = mytracker.read_seq_table()
    pr_table = mytracker.read_pr_table()
    wm = dm.WorkspaceManager(wsname).disable_hound()
    samples_in_ws = wm.get_samples().replace(np.nan, "", regex=True).reset_index()

    # check broken bam file paths
    print("checking if there are any broken file paths")
    foundfiles = gcp.lsFiles(samples_in_ws[bamcol])
    broken_bams = set(samples_in_ws[bamcol]) - set(foundfiles)
    print(
        "These "
        + str(len(broken_bams))
        + " bam file path do not exist: "
        + str(broken_bams)
    )
    samples = samples_in_ws[(~samples_in_ws[extract["bam"]].isin(broken_bams))]

    print(
        "extracting information from workspace including hash, size, update time, etc."
    )
    samples = extractFromWorkspace(samples, stype)
    print("generating CDS-ids, annotating source, and renaming columns")
    samples = mapSamples(samples, source, extract=extract)
    print("checking for duplicates in the workspace by comparing file sizes")
    samples = resolveFromWorkspace(
        samples,
        seq_table[seq_table[extract["expected_type"]] == stype],
        wsidcol,
        accept_unknowntypes,
        addonly,
        extract,
    )
    samples = samples[samples[extract["update_time"]] > maxage]
    if len(samples) == 0:
        print("no new data available")
        return None

    # out of the new samples, see if ops has registered them in gumbo
    unmapped_new_lines = []
    # init profile id column
    samples[extract["profile_id"]] = ""
    samples[extract["version"]] = 1
    # FOR NOW, assume the id col to map by is in profile table (might change??)
    for k, v in samples.iterrows():
        if (not pd.isnull(v[wsidcol])) and v[wsidcol] in pr_table[gumboidcol].tolist():
            # different datatypes from the same line might share the same SM-ID
            # so mapping should condition on datatype as well
            pr_id = pr_table[
                (pr_table[gumboidcol] == v[wsidcol]) & (pr_table.Datatype == stype)
            ].index
            if len(pr_id) > 1:
                raise ValueError(
                    "multiple profile ids mapped to the same validation id. check with ops!"
                )
            else:
                samples.loc[k, extract["profile_id"]] = pr_id[0]
        else:
            unmapped_new_lines.append(k)

    return samples, unmapped_new_lines


def loadFromMultipleSources(terraworkspaces=[]):
    """
    Load new samples from multiple sources

    Args:
    -----
        terraworkspaces ([dict]): 

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    return True


def deleteClosest(
    sampless,
    refsamples,
    size="legacy_size",
    ref_size="legacy_size",
    arxspid="arxspan_id",
):
    """
    for a list of samples and a tracker, will find the index of the sample with the closest size

    if this sample is the same cell line, it will judge it to be a duplicate and remove it

    Args:
    -----
        sampless: pd dataframes of samples with at least arxspan ids and sizes
        refsamples: pd dataframe representing a sample tracker
        size: str colname of size in the sample list
        ref_size: str colname of size in the sample tracker
        arxspid: str colnme of sample ids
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    sizes = refsamples[ref_size].tolist()
    print("deleting closest samples:")
    for k, v in sampless.iterrows():
        if type(v[size]) is int:
            val = refsamples.iloc[sizes.index(h.closest(sizes, v[size]))]
            if val[arxspid] == v[arxspid]:
                sampless = sampless.drop(v.name)
                print(v.name)
    return sampless


def extractFromWorkspace(
    samples,
    stype,
    recomputeTime=True,
    recomputesize=True,
    recomputedate=True,
    recompute_hash=True,
    extract={},
):
    """
    Extract more information from a list of samples found on GP workspaces

    Args:
    -----
        samples: pd dataframes of samples with at least arxspan ids and sizes
        stype: str sequencing type
        recomputeTime: bool whether to recompute the date of upload of the bam file
        recomputesize: bool whether to recompute the of the bam file
        recomputehash: bool whether to recompute the of the bam file
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    extract.update(EXTRACT_DEFAULTS)
    if extract["legacy_hash"] not in samples.columns or recompute_hash:
        samples[extract["hash"]] = [
            gcp.extractHash(val)
            for val in gcp.lsFiles(samples[extract["bam"]].tolist(), "-L", 200)
        ]
    lis = gcp.lsFiles(samples[extract["bam"]].tolist(), "-al", 200)
    if extract["legacy_size"] not in samples.columns or recomputesize:
        samples[extract["legacy_size"]] = [gcp.extractSize(i)[1] for i in lis]
    if extract["update_time"] not in samples.columns or recomputeTime:
        samples[extract["update_time"]] = [gcp.extractTime(i) for i in lis]
    todrop = []
    for k, val in samples.iterrows():
        if val[extract["legacy_size"]] < MINSIZES[stype]:
            todrop.append(k)
            print(
                "too small size, removing sample: "
                + str(val[extract["from_arxspan_id"]])
            )
    samples = samples.drop(index=todrop)
    # getting the date released
    if len(samples) == 0:
        return None
    if extract["release_date"] not in samples.columns or recomputedate:
        samples[extract["release_date"]] = seq.getBamDate(samples[extract["bam"]])
    return samples


def mapSamples(samples, source, extract={}):
    """
    Convert samples from a list of GP workspaces to something being able to be merged with the sample tracker

    Args:
    -----
        samples: pd dataframes of samples with at least arxspan ids and sizes
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)
        source:

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    # creating unique ids
    samples[extract["ref_id"]] = [
        "CDS-" + h.randomString(stringLength=6, stype="all", withdigits=True)
        for _ in range(len(samples))
    ]
    samples.reset_index(drop=True, inplace=True)
    samples[extract["source"]] = source

    # renamings
    samples = samples.rename(
        columns={
            extract["bam"]: extract["ref_bam"],
            extract["bai"]: extract["ref_bai"],
            extract["root_sample_id"]: extract["sm_id"],
            extract["PDO_id_terra"]: extract["PDO_id_gumbo"],
        }
    ).set_index(extract["ref_id"], drop=True)
    # subsetting
    print(samples.columns)
    samples = samples[
        list(
            set(
                [
                    extract["ref_bam"],
                    extract["ref_bai"],
                    extract["release_date"],
                    extract["legacy_size"],
                    extract["PDO_id_gumbo"],
                    extract["sm_id"],
                    extract["update_time"],
                    extract["source"],
                ]
            )
        )
    ]
    return samples


def resolveFromWorkspace(
    samples, refsamples, wsidcol, accept_unknowntypes=True, addonly=[], extract={},
):
    """
    Filters our list by trying to find duplicate in our dataset and remove any sample that isn't tumor

    Args:
    -----
        match: list[str]|str the possible values that a sample id need to contain to be considered valid
        participantslicepos: int the length of the sample id string
        accept_unknowntypes: bool whether or not the sample type column for that sample can be different from "Tumor"
        refsamples: pd dataframe representing a sample tracker
        samples: pd dataframes of samples with at least arxspan ids and sizes
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)

    Returns:
    --------
        samples: pd dataframe the filtered sample list
    """
    extract.update(EXTRACT_DEFAULTS)

    sample_size = {
        gcp.extractSize(val)[1]: gcp.extractSize(val)[0]
        for val in gcp.lsFiles(samples[extract["ref_bam"]], "-la")
    }
    dups_to_remove = [
        sample_size[a]
        for a in set(sample_size.keys()) & set(refsamples[extract["legacy_size"]])
    ]
    # remove the duplicates from consideration
    print("Len of samples before removal: " + str(len(samples)))
    print(
        "Dups from this workspace has len "
        + str(len(dups_to_remove))
        + ":\n "
        + str(dups_to_remove)
    )
    # remove the samples with broken bam filepaths from consideration
    samples = samples[~samples[extract["ref_bam"]].isin(dups_to_remove)]

    print("Len of samples after removal: " + str(len(samples)))
    if len(samples) == 0:
        return None

    # if only add some samples
    if len(addonly) > 0:
        samples = samples[samples[wsidcol].isin(addonly)]

    # unknown types
    if "sample_type" in samples.columns:
        if not accept_unknowntypes:
            samples = samples[samples["sample_type"].isin(["Tumor"])]
    return samples


def assessAllSamples(sampless, refsamples, stype, rename={}, extract={}):
    """
    Will look for matching lines and duplicates in our sample tracker and compute version and patient information

    Args:
    -----
        refsamples: pd dataframe representing a sample tracker
        stype: str sequencing type
        rename: dict(str:str) mapping a wrong arxpand_id to a good arxspan id for known cases of misslabelling
        samples: pd dataframes of samples with at least arxspan ids and sizes
        extract: if you want to specify what values should refer to which column names
        dict{
        'name':
        'bai':
        'bam':
        'source':
        'from_arxspan_id':
        ...} (see extract_defaults)
    Returns:
    --------
        samples: pd daataframe the filtered sample list
    """
    extract.update(EXTRACT_DEFAULTS)
    rename.update(DUP_ARXSPANS)
    prevlen = len(sampless)
    sampless[extract["ref_type"]] = stype

    # checking no duplicate in the buckets
    for k, val in sampless.iterrows():
        withsamesize = sampless[
            sampless[extract["legacy_size"]] == val[extract["legacy_size"]]
        ]
        if len(withsamesize) > 1:
            if (
                len(
                    withsamesize[
                        withsamesize[extract["ref_name"]] == val[extract["ref_name"]]
                    ]
                )
                < 2
            ):
                raise ValueError("we have duplicate samples with different names!")
            else:
                for l, v in withsamesize.iloc[1:].iterrows():
                    sampless = sampless.drop(l)
    print(
        "we had " + str(prevlen - len(sampless)) + " duplicates in the release buckets"
    )
    # check: currently, below lines prevent forcekeep from working on true duplicates
    # (aka same size  file). Need to think about how to bring back the forcekeep functionality
    sampless[extract["ref_arxspan_id"]] = [
        rename[name] if name in rename else name
        for name in sampless[extract["ref_arxspan_id"]]
    ]
    names = []
    # need to keep track of whether we're adding more than one new entry for a given sample id
    subrefsamples = refsamples[refsamples[extract["ref_type"]] == stype]
    for k, val in sampless.iterrows():
        val = val[extract["ref_arxspan_id"]]
        names.append(val)
        sampless.loc[k, extract["version"]] = len(
            subrefsamples[subrefsamples[extract["ref_arxspan_id"]] == val]
        ) + names.count(val)
        # TODO: copy the patient id too
    sampless[extract["version"]] = sampless[extract["version"]].astype(int)

    sampless[extract["patient_id"]] = [
        val[extract["patient_id"]]
        if len(
            refsamples[
                refsamples[extract["ref_arxspan_id"]] == val[extract["ref_arxspan_id"]]
            ]
        )
        < 1
        else refsamples[
            refsamples[extract["ref_arxspan_id"]] == val[extract["ref_arxspan_id"]]
        ][extract["patient_id"]].values[0]
        for i, val in sampless.iterrows()
    ]

    return sampless


def addSamplesToGumbo(
    samples,
    stype,
    bucket,
    name_col="index",
    values=["hg19_bam_filepath", "hg19_bai_filepath"],
    filetypes=["bam", "bai"],
):
    """update the samples on gumbo's sequencing sheet

    Args:
        samples ([type]): [description]
        stype ([type]): [description]
        bucket ([type]): [description]
        name_col (str, optional): [description]. Defaults to "index".
        values (list, optional): [description]. Defaults to ['legacy_bam_filepath', 'legacy_bai_filepath'].
        filetypes (list, optional): [description]. Defaults to ['bam', 'bai'].
    """
    # uploading to our bucket (now a new function)
    assert set(values).issubset(set(samples.columns))
    samples = terra.changeToBucket(
        samples,
        bucket,
        name_col=name_col,
        values=values,
        filetypes=filetypes,
        catchdup=True,
        dryrun=False,
    )

    mytracker = track.SampleTracker()
    ccle_refsamples = mytracker.read_seq_table()

    names = []
    subccle_refsamples = ccle_refsamples[ccle_refsamples["ExpectedType"] == stype]
    for k, val in samples.iterrows():
        val = val["ProfileID"]
        if val != "":
            names.append(val)
            samples.loc[k, "version"] = len(
                subccle_refsamples[subccle_refsamples["ProfileID"] == val]
            ) + names.count(val)
    samples["version"] = samples["version"].astype(int)

    ccle_refsamples = ccle_refsamples.append(samples, sort=False)
    mytracker.write_seq_table(ccle_refsamples)


def addSamplesToDepMapWorkspace(
    stype, refworkspace, samplesetname="", add_to_samplesets=[],
):
    """update the samples on a depmapomics terra processing workspace

    Args:
        stype (str): data type
        refworkspace (str): terra processing workspace to import data to
        add_to_samplesets (list, optional): add new samples to additional sample_sets on terra. Defaults to []
    """
    mytracker = track.SampleTracker()
    refwm = dm.WorkspaceManager(refworkspace).disable_hound()

    terra_samples = refwm.get_samples()
    seq_table = mytracker.read_seq_table()
    # check which lines are new and need to be imported to terra
    samples_to_add = seq_table[
        (~seq_table.index.isin(terra_samples.index)) & (seq_table.ExpectedType == stype)
    ]
    print("found " + str(len(samples_to_add)) + " new samples to import!")

    # terra requires a participant id column
    # temp multi-step solution here: get participant id from seq id
    pr_table = mytracker.read_pr_table()
    mc_table = mytracker.read_mc_table()
    model_table = mytracker.read_model_table()

    for i in samples_to_add.index:
        samples_to_add.loc[i, "participant_id"] = mytracker.get_participant_id(
            i, seq_table, pr_table, mc_table, model_table
        )

    # uploading new samples
    samples_to_add.index.name = "sample_id"

    refwm.upload_samples(samples_to_add)

    refwm.update_sample_set(
        sample_set_id="all",
        sample_ids=[i for i in refwm.get_samples().index.tolist() if i != "nan"],
    )

    # creating a sample set, if samplesetname is provided
    if samplesetname != "":
        print("adding new samples to sampleset: " + samplesetname)
        refwm.update_sample_set(
            sample_set_id=samplesetname, sample_ids=samples_to_add.index
        )

    # add new samples to additional existing sample_sets if needed
    for sname in add_to_samplesets:
        samples_in_sname = refwm.get_sample_sets().loc[sname, "samples"]
        new_samples = samples_to_add.index

        refwm.update_sample_set(
            sample_set_id=sname, sample_ids=list(set(samples_in_sname + new_samples))
        )
