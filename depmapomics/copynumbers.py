# cn.py

import numpy as np
from depmapomics import tracker as track
from depmapomics.qc import cn
from depmapomics.config import *
from IPython.display import Image, display
import dalmatian as dm
import pandas as pd
import os
from genepy import mutations as mut
from genepy.utils import helper as h
from genepy import terra
from genepy import rna
import matplotlib.pyplot as plt
from depmapomics import terra as myterra


def renameColumns(df):
    """
    rename some of the main columns names from RSEM, GATK.. to more readable column names

    Args:
    -----
      df: the df to rename

    Returns:
    ------
      df the renamed df
    """
    return df.rename(columns=COLRENAMING)


def loadFromGATKAggregation(
    refworkspace,
    setEntity="sample_set",
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    save_output="",
    doCleanup=True,
    todrop=[],
    showPlots=False,
    colname="combined_seg_file",
    plotColname="modeled_segments_plot_tumor",
    tempFolder="temp/",
    toremove=["readgroup_ubams",],
    sampleset="all",
    colRenaming=COLRENAMING,
):
    """fetching the data from the Terra and loading it into a dataframe

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        doCleanup (bool, optional): if do cleanup of Terra workspace unused output and logs. Defaults to True.
        todrop (list, optional): [description]. Defaults to [].
        showPlots (bool, optional): whether to show plot output from the GATK pipeline. Defaults to False.
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        plotColname (str, optional): the column on terra where the plots exist in sample. Defaults to "modeled_segments_plot_tumor".
        tempFolder (str, optional): where to put temp files. Defaults to "temp/".
        toremove (list, optional): columns in Terra samples to remove (and delete corresponding data). Defaults to ["readgroup_ubams", ].
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to COLRENAMING.

    Returns:
        pd.dataframe: dataframe containing the segments concatenated in a bed like format
    """
    wm = dm.WorkspaceManager(refworkspace)
    if save_output:
        terra.saveConfigs(refworkspace, os.path.join(save_output, "terra/"))

    if doCleanup:
        print("cleaning up")
        res = wm.get_samples()
        for val in toremove:
            if val in res.columns.tolist():
                wm.disable_hound().delete_entity_attributes(
                    "sample", res[val], delete_files=True
                )

    segments = pd.read_csv(
        wm.get_entities(setEntity).loc[sampleset, colname], sep="\t"
    ).rename(columns=colRenaming)

    # removing the duplicates
    segments = segments[~segments[SAMPLEID].isin(todrop)].reset_index(drop=True)
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.Segment_Mean = 2 ** segments.Segment_Mean
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments.loc[
        segments[segments.Chromosome.isin(["X", "Y"])].index, "Segment_Mean"
    ] = (segments[segments.Chromosome.isin(["X", "Y"])]["Segment_Mean"] / 2)
    segments = segments.sort_values(by=sortby)

    print("loading " + str(len(set(segments[SAMPLEID]))) + " rows")
    if showPlots:
        # plotting results of CN calls for this new sample set
        for i, (k, val) in enumerate(
            wm.get_samples().loc[wm.get_sample_sets().loc[sampleset].samples].iterrows()
        ):
            plot = val[plotColname]
            os.system("gsutil cp " + plot + " " + tempFolder)
            print(k)
            print(val["arxspan_id"], val["sex"])
            if i > 30:
                continue
            display(Image(os.path.join(tempFolder, plot.split("/")[-1])))
    return segments


def updateTracker(
    tracker,
    selected,
    samplesetname,
    lowqual,
    datatype,
    newgs=WGS_GCS_PATH_HG38,
    samplesinset=[],
    trackerobj=None,
    procqc=[],
    bamqc=[],
    refworkspace=None,
    bamfilepaths=["internal_bam_filepath", "internal_bai_filepath"],
    dry_run=False,
):
    """updates the sample tracker with the new samples and the QC metrics

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
    track.update(
        tracker,
        selected,
        samplesetname,
        lowqual,
        lowqual,
        newgs,
        refworkspace,
        bamfilepaths,
        dry_run,
        samplesinset,
        trackerobj=trackerobj,
    )


def managingDuplicates(samples, failed, datatype, tracker, newname="arxspan_id"):
    """removes duplicates and solves failed data

    Args:
        failed (list): list of samples that failed QC
        datatype (str): 'wes' or 'wgs'
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        samples (list): list of samples to filter
        newname (str): name of the new sample set

    Returns:
        dict: a dict of duplicates to keep and their newname name
    """
    # selecting the right arxspan id (latest version)
    renaming = track.removeOlderVersions(
        names=samples,
        refsamples=tracker[tracker.datatype == datatype],
        priority="prioritized",
    )

    # reparing QC when we have a better duplicate
    ref = pd.DataFrame(tracker[tracker.datatype == datatype][newname])
    replace = 0
    for val in failed:
        if val in renaming:
            a = ref[ref[newname] == ref.loc[val][newname]].index
            for v in a:
                if v not in failed:
                    renaming[v] = renaming.pop(val)
                    replace += 1
                    break
    print("could replace:")
    print(replace)
    return renaming


def pureCNpostprocess(
    refworkspace,
    setEntity="sample_set",
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    priority=[],
    colname="PureCN_loh_merged",
    sampleset="all",
    colRenaming=PURECN_COLRENAMING,
    terracols=SIGTABLE_TERRACOLS,
    save_output="",
    mappingdf=None,
    genechangethresh=0.025,
    segmentsthresh=2000,
):
    """fetching PureCN data from Terra, generate one matrix for absolute copy number, one matrix for LOH,
    and one cell line signature table

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        todrop (list, optional): [description]. Defaults to [].
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to COLRENAMING.
        genechangethresh (float, optional): above this threshold of variance of gene CN, the sample is considered failed. Defaults to 0.025.
        segmentsthresh (int, optional): above this threshold of number of segments the WGS sample is considered failed. Defaults to 2000.

    Returns:
        pd.dataframe: dataframe containing the segments concatenated in a bed like format
    """

    print("loading PureCN merged LOH file")
    wm = dm.WorkspaceManager(refworkspace)
    segments = pd.read_csv(wm.get_entities(setEntity).loc[sampleset, colname]).rename(
        columns=colRenaming
    )

    # removing the duplicates
    segments = segments[~segments[SAMPLEID].isin(todrop)].reset_index(drop=True)
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.Segment_Mean = 2 ** segments.Segment_Mean
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments.loc[
        segments[segments.Chromosome.isin(["X", "Y"])].index, "Segment_Mean"
    ] = (segments[segments.Chromosome.isin(["X", "Y"])]["Segment_Mean"] / 2)
    segments = segments.sort_values(by=sortby)

    print("loading " + str(len(set(segments[SAMPLEID]))) + " rows")

    # Generate gene-level absolute cn matrix
    absolute_genecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments), mappingdf, value_colname="Segment_Mean"
    )
    print("summary of PureCN absolute cn data:")
    print(
        absolute_genecn.values.min(),
        absolute_genecn.values.mean(),
        absolute_genecn.values.max(),
    )
    mut.checkGeneChangeAccrossAll(absolute_genecn, thresh=genechangethresh)
    failed = mut.checkAmountOfSegments(segments, thresh=segmentsthresh)

    # Generate gene-level LOH status matrix
    segments["type"] = segments["type"].replace(
        ["LOH", "COPY-NEUTRAL LOH", "WHOLE ARM COPY-NEUTRAL LOH"], 1
    )
    segments["type"] = segments["type"].replace("", 0)

    loh_status = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments), mappingdf, value_colname="LOH_status"
    )

    print("PureCN: failed our QC")
    print(failed)

    # Pull additional info directly from terra sample table
    samples = wm.get_samples()
    purecn_table = samples[terracols]

    segments = segments[
        ~segments[SAMPLEID].isin((set(failed) | set(todrop)) - set(priority))
    ].reset_index(drop=True)
    absolute_genecn = absolute_genecn[
        ~absolute_genecn.index.isin((set(failed) | set(todrop)) - set(priority))
    ]
    loh_status = loh_status[
        ~loh_status.index.isin((set(failed) | set(todrop)) - set(priority))
    ]
    purecn_table = purecn_table[
        ~purecn_table.index.isin((set(failed) | set(todrop)) - set(priority))
    ]

    print("PureCN: saving files")
    segments.to_csv(save_output + "purecn_segments_all.csv", index=False)
    absolute_genecn.to_csv(save_output + "purecn_genecn_all.csv")
    loh_status.to_csv(save_output + "purecn_loh_all.csv")
    purecn_table.to_csv(save_output + "purecn_table_all.csv")
    print("done")

    return segments, absolute_genecn, loh_status, purecn_table


def postProcess(
    refworkspace,
    setEntity="sample_set",
    sampleset="all",
    save_output="",
    doCleanup=True,
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    priority=[],
    genechangethresh=0.025,
    segmentsthresh=1500,
    ensemblserver=ENSEMBL_SERVER_V,
    source_rename={},
    useCache=False,
    maxYchrom=150,
):
    """post process an aggregated CN segment file, the CCLE way

    take a CN segment file from the Aggregate_WGS master terra workflow and post process it
    in the CCLE way.

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        genechangethresh (float, optional): above this threshold of variance of gene CN, the sample is considered failed. Defaults to 0.025.
        segmentsthresh (int, optional): above this threshold of number of segments the WGS sample is considered failed. Defaults to 1500.
        ensemblserver (str, optional): ensembl server biomart version . Defaults to ENSEMBL_SERVER_V.
        source_rename (dict, optional): dict to rename the source column if needed. Defaults to {}.
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.

    Returns:
        pd.dataframe: a dataframe with the post processed segment file as a bed file concatenated over all samples.
        pd.dataframe: a dataframe with the post processed gene CN file as a bed file of samples x genes.
        list: list of samples that failed to be processed.
    """
    h.createFoldersFor(save_output)
    print("loading CN from Terra")
    segments = loadFromGATKAggregation(
        refworkspace,
        setEntity=setEntity,
        sampleset=sampleset,
        sortby=sortby,
        todrop=todrop,
        doCleanup=doCleanup,
    )
    print("making gene level copy number")

    mybiomart = h.generateGeneNames(
        ensemble_server=ensemblserver,
        useCache=useCache,
        attributes=["start_position", "end_position", "chromosome_name"],
    )
    mybiomart = mybiomart.rename(
        columns={
            "start_position": "start",
            "end_position": "end",
            "chromosome_name": "Chromosome",
        }
    )
    mybiomart["Chromosome"] = mybiomart["Chromosome"].astype(str)
    mybiomart = mybiomart.sort_values(by=["Chromosome", "start", "end"])
    mybiomart = mybiomart[mybiomart["Chromosome"].isin(set(segments["Chromosome"]))]
    mybiomart = mybiomart.drop_duplicates("hgnc_symbol", keep="first")
    mybiomart["gene_name"] = [
        i["hgnc_symbol"] + " (" + str(i["entrezgene_id"]).split(".")[0] + ")"
        for _, i in mybiomart.iterrows()
    ]
    # drop Ychrom if > maxYchrom
    ychrom = segments[segments.Chromosome.str.contains("Y")]
    countYdrop = [
        i
        for i in set(ychrom[SAMPLEID])
        if len(ychrom[ychrom[SAMPLEID] == i]) > maxYchrom
    ]
    segments = segments[
        ~((segments[SAMPLEID].isin(countYdrop)) & (segments.Chromosome == "Y"))
    ]
    genecn = mut.toGeneMatrix(mut.manageGapsInSegments(segments), mybiomart)
    # validation step
    print("summary of the gene cn data:")
    print(genecn.values.min(), genecn.values.mean(), genecn.values.max())
    mut.checkGeneChangeAccrossAll(genecn, thresh=genechangethresh)
    failed = mut.checkAmountOfSegments(segments, thresh=segmentsthresh)

    print("failed our QC")
    print(failed)

    if source_rename:
        segments = segments.replace({"Source": source_rename})
    if save_output:
        h.listToFile(failed, save_output + "failed.txt")
    # subsetting
    segments = segments[
        ~segments[SAMPLEID].isin((set(failed) | set(todrop)) - set(priority))
    ].reset_index(drop=True)
    genecn = genecn[~genecn.index.isin((set(failed) | set(todrop)) - set(priority))]

    # saving
    print("saving files")
    segments.to_csv(save_output + "segments_all.csv", index=False)
    genecn.to_csv(save_output + "genecn_all.csv")
    print("done")
    purecn_segments, purecn_genecn, loh_status, purecn_table = pureCNpostprocess(
        refworkspace,
        setEntity=setEntity,
        sampleset=sampleset,
        sortby=sortby,
        todrop=todrop,
    return segments, genecn, failed, purecn_segments, purecn_genecn, loh_status, purecn_table


def _CCLEPostProcessing(
    wesrefworkspace=WESCNWORKSPACE,
    wgsrefworkspace=WGSWORKSPACE,
    wessetentity=WESSETENTITY,
    wgssetentity=WGSSETENTITY,
    samplesetname=SAMPLESETNAME,
    AllSamplesetName="all",
    my_id=MY_ID,
    mystorage_id=MYSTORAGE_ID,
    sheetcreds=SHEETCREDS,
    sheetname=SHEETNAME,
    refsheet_url=REFSHEET_URL,
    todrop=KNOWN_DROP,
    prevgenecn="ccle",
    taiga_dataset=TAIGA_CN,
    dataset_description=CNreadme,
    subsetsegs=[
        SAMPLEID,
        "Chromosome",
        "Start",
        "End",
        "Segment_Mean",
        "Num_Probes",
        "Status",
        "Source",
    ],
    bamqc=BAMQC,
    procqc=PROCQC,
    source_rename=SOURCE_RENAME,
    redoWES=False,
    wesfolder="",
    **kwargs
):
    """the full CCLE Copy Number post processing pipeline (used only by CCLE)

    see postprocessing() to reproduce most of our analysis and find out about additional parameters

    Args:
        wesrefworkspace (str): wes terra workspace where the ref data is stored
        wgsrefworkspace (str): wgs terra workspace where the ref data is stored
        samplesetname (str): name of the current release
        AllSamplesetName (str, optional): name of the sample set to get the data from (should contain everything). Defaults to 'all'.
        my_id (str, optional): google sheet user id to access the ref sheets . Defaults to MY_ID.
        mystorage_id (str, optional): google sheet storage id to access the ref sheets . Defaults to MYSTORAGE_ID.
        sheetcreds (str, optional): @see updateTracker. Defaults to SHEETCREDS.
        sheetname (str, optional): @see updateTracker. Defaults to SHEETNAME.
        refsheet_url (str, optional): @see updateTracker. Defaults to REFSHEET_URL.
        prevgenecn (str, optional): where the previous version exists on taiga (for QC purposes). Defaults to tc.get(name=TAIGA_ETERNAL, file='CCLE_gene_cn').
        taiga_dataset (str, optional): where to save the output to on taiga. Defaults to TAIGA_CN.
        dataset_description (str, optional): A long string that will be pushed to taiga to explain the CN dataset. Defaults to CNreadme.
        subsetsegs (list[str], optional): what columns to keep for the segments. Defaults to [SAMPLEID, 'Chromosome', 'Start', 'End', 'Segment_Mean', 'Num_Probes', 'Status', 'Source'].
        bamqc ([type], optional): @see updateTracker. Defaults to BAMQC.
        procqc ([type], optional): @see updateTracker. Defaults to PROCQC.
        source_rename ([type], optional): @see managing duplicates. Defaults to SOURCE_RENAME.
    """
    from taigapy import TaigaClient

    tc = TaigaClient()

    if prevgenecn is "ccle":
        prevgenecn = tc.get(name=TAIGA_ETERNAL, file="CCLE_gene_cn")

    sheets = Sheets.from_files(my_id, mystorage_id)
    tracker = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)

    # doing wes
    folder = os.path.join("temp", samplesetname, "wes_")
    if redoWES:
        print("doing wes")
        priority = tracker[
            (tracker.datatype == "wes") & (tracker.prioritized == 1)
        ].index.tolist()
        todropwes = (
            todrop
            + tracker[
                (tracker.datatype == "wes") & (tracker.blacklist == 1)
            ].index.tolist()
        )
        wessegments, genecn, wesfailed = postProcess(
            wesrefworkspace,
            setEntity=wessetentity,
            sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
            todrop=todropwes,
            save_output=folder,
            priority=priority,
            **kwargs
        )

        wesrenaming = managingDuplicates(
            set(wessegments[SAMPLEID]),
            (set(wesfailed) - set(priority)) | set(todropwes),
            "wes",
            tracker,
        )
        h.dictToFile(wesrenaming, folder + "sample_renaming.json")

        # annotating source
        for v in set(wessegments[SAMPLEID]):
            wessegments.loc[
                wessegments[wessegments[SAMPLEID] == v].index, "Source"
            ] = tracker.loc[v, "source"]

        wessegments.Source = wessegments.Source.replace(source_rename)
        wessegments.Source += " WES"

        print("renaming")
        wespriosegments = (
            wessegments[wessegments[SAMPLEID].isin(set(wesrenaming.keys()))]
            .replace({SAMPLEID: wesrenaming})
            .reset_index(drop=True)
        )
        wespriogenecn = genecn[genecn.index.isin(set(wesrenaming.keys()))].rename(
            index=wesrenaming
        )

        # saving prio
        wespriosegments.to_csv(folder + "segments_latest.csv", index=False)
        wespriogenecn.to_csv(folder + "genecn_latest.csv")
    else:
        print("bypassing WES using folder: " + wesfolder if wesfolder else folder)
        wesfailed = h.fileToList((wesfolder if wesfolder else folder) + "failed.txt")
        wesrenaming = h.fileToDict(
            (wesfolder if wesfolder else folder) + "sample_renaming.json"
        )
        wespriosegments = pd.read_csv(folder + "segments_latest.csv")
        wespriogenecn = pd.read_csv(folder + "genecn_latest.csv", index_col=0)

    # doing wgs
    print("doing wgs")
    folder = os.path.join("temp", samplesetname, "wgs_")
    priority = tracker[
        (tracker.datatype == "wgs") & (tracker.prioritized == 1)
    ].index.tolist()
    todropwgs = (
        todrop
        + tracker[(tracker.datatype == "wgs") & (tracker.blacklist == 1)].index.tolist()
    )
    wgssegments, genecn, wgsfailed = postProcess(
        wgsrefworkspace,
        setEntity=wgssetentity,
        sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
        todrop=todropwgs,
        save_output=folder,
        segmentsthresh=1500,
        priority=priority,
    )

    return (
        segments,
        genecn,
        purecn_segments,
        purecn_genecn,
        loh_status,
        purecn_table,
        failed,
    )

