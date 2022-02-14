import os.path
import dalmatian as dm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from genepy.utils import helper as h
from genepy import rna
from genepy import mutations as mut

from depmapomics import tracker as track
from depmapomics import expressions, mutations
from depmapomics import fusions as fusion
from depmapomics import copynumbers as cn
from depmapomics.config import *


def expressionRenaming(r, todrop):
    trackerobj = track.initTracker()
    ccle_refsamples = trackerobj.read_tracker()

    renaming = track.removeOlderVersions(
        names=r,
        refsamples=ccle_refsamples[ccle_refsamples.datatype == "rna"],
        priority="prioritized",
    )
    # if we have a replaceable failed version in our dataset
    rename = expressions.solveQC(ccle_refsamples, todrop)
    for k, _ in renaming.copy().items():
        if k in rename:
            renaming[rename[k]] = renaming.pop(k)
        elif (k in todrop) and (k not in rename):
            renaming.pop(k)
    return renaming


async def expressionPostProcessing(
    refworkspace=RNAWORKSPACE,
    samplesetname=SAMPLESETNAME,
    colstoclean=["fastq1", "fastq2", "recalibrated_bam", "recalibrated_bam_index"],
    ensemblserver=ENSEMBL_SERVER_V,
    doCleanup=True,
    samplesetToLoad="all",
    tocompare={
        "genes_expected_count": "CCLE_RNAseq_reads",
        "genes_tpm": "CCLE_expression_full",
        "proteincoding_genes_tpm": "CCLE_expression",
    },
    trackerobj=None,
    todrop=KNOWN_DROP,
    prevcounts="ccle",
    taiga_dataset=TAIGA_EXPRESSION,
    minsimi=0.95,
    dropNonMatching=True,
    dataset_description=RNAseqreadme,
    dry_run=False,
    samplesinset=[],
    rsemfilelocs=None,
    rnaqclocs={},
    starlogs={},
    **kwargs,
):
    """the full CCLE Expression post processing pipeline (used only by CCLE)

    @see postprocessing() to reproduce our analysis and for parameters

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        colstoclean (list, optional): the columns to clean in the terra workspace. Defaults to [].
        ensemblserver (str, optional): ensembl server biomart version . Defaults to ENSEMBL_SERVER_V.
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        geneLevelCols (list, optional): the columns that contain the gene level 
        expression data in the workspace. Defaults to RSEMFILENAME_GENE.
        trancriptLevelCols (list, optional): the columns that contain the transcript 
        level expression data in the workspacce. Defaults to RSEMFILENAME_TRANSCRIPTS.
        ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
        renamingFunc (function, optional): the function to use to rename the sample columns
        (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
        compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
        dropNonMatching (bool, optional): whether to drop the non matching genes 
        between entrez and ensembl. Defaults to False.
        recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
        prevcounts (str, optional): the previous counts to use to QC the data for the release. Defaults to 'ccle'.
        taiga_dataset (str, optional): the taiga dataset path to use for uploading results. Defaults to TAIGA_EXPRESSION.
        minsimi (float, optional): the minimum similarity to use for comparison to previous dataset. Defaults to 0.95.
        dataset_description (str, optional): the taiga dataset description to use. Defaults to RNAseqreadme.
        tocompare (dict, optional): the columns to compare. Defaults to {"genes_expected_count": "CCLE_RNAseq_reads", "genes_tpm": "CCLE_expression_full", "proteincoding_genes_tpm": "CCLE_expression"}.
        rsemfilelocs (pd.DataFrame, optional): locations of RSEM output files if refworkspace is not provided (bypass interaction with terra)
        samplesinset (list[str], optional): list of samples in the sampleset if refworkspace is not provided (bypass interaction with terra)
        rnaqclocs (dict(str:list[str]), optional): dict(sample_id:list[QC_filepaths]) of rna qc file locations if refworkspace is not provided (bypass interaction with terra)
    """
    from taigapy import TaigaClient

    tc = TaigaClient()

    if prevcounts == "ccle":
        prevcounts = tc.get(name=TAIGA_ETERNAL, file="CCLE_RNAseq_reads")

    ccle_refsamples = trackerobj.read_tracker()

    todrop += ccle_refsamples[
        (ccle_refsamples.blacklist == 1) & (ccle_refsamples.datatype == "rna")
    ].index.tolist()
    priority = ccle_refsamples[
        (ccle_refsamples.prioritized == 1) & (ccle_refsamples.datatype == "rna")
    ].index.tolist()

    folder = os.path.join("temp", samplesetname, "")

    if dry_run:
        folder = os.path.join("temp", "dryrun", "")

    h.createFoldersFor(folder)
    files, failed, _, renaming, lowqual, _ = await expressions.postProcess(
        refworkspace,
        samplesetname,
        save_output=folder,
        doCleanup=doCleanup,
        priority=priority,
        colstoclean=colstoclean,
        ensemblserver=ensemblserver,
        todrop=todrop,
        samplesetToLoad=samplesetToLoad,
        geneLevelCols=RSEMFILENAME_GENE,
        trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,  # compute_enrichment=False,
        ssGSEAcol="genes_tpm",
        renamingFunc=expressionRenaming,
        dropNonMatching=dropNonMatching,
        dry_run=dry_run,
        samplesinset=samplesinset,
        rsemfilelocs=rsemfilelocs,
        rnaqclocs=rnaqclocs,
        starlogs=starlogs,
        **kwargs,
    )

    print("doing validation")
    nonoverlap = set(prevcounts.columns) ^ set(
        files.get("genes_expected_count").columns
    )
    print("number of non overlaping genes:")
    print(len(nonoverlap))
    # have we lost any samples compared to last release?
    lost = set(prevcounts.index) - set(files.get("genes_expected_count").index)
    print("of which, lost genes:")
    print(lost)
    # do we have samples that are missanotated compared to previous releases (replicate level)
    # notindataset, missannotated, unmatched = findMissAnnotatedReplicates(replevel, prevcounts, renaming)
    # for k,v in unmatched.items():
    #    if ccle_refsamples.loc[k].arxspan_id!=v:
    #        print(k,v)
    # do we have samples that are missanotated compared to previous releases (sample level)
    unmatched = rna.getDifferencesFromCorrelations(
        files.get("genes_expected_count"), prevcounts, minsimi=minsimi
    )
    print("differences in correlations against the previous release")
    print(unmatched)
    # Is it because of  duplicate version?
    print("do we see it as a duplicate in the tracker?")
    rnasamples = ccle_refsamples[ccle_refsamples.datatype == "rna"]
    for i, _ in unmatched:
        print(len(rnasamples[rnasamples.arxspan_id == i]))

    # CCLE_expression, CCLE_expression_full, ,
    print("comparing to previous release")
    # h.compareDfs(files["rsem_transcripts_tpm"], tc.get(name=TAIGA_ETERNAL, file='CCLE_RNAseq_transcripts'))
    # h.compareDfs(files["rsem_transcripts_expected_count"], tc.get(name=TAIGA_ETERNAL, file='CCLE_expression_transcripts_expected_count'))
    # h.compareDfs(enrichments, tc.get(name=TAIGA_ETERNAL, file='CCLE_fusions_unfiltered'))
    for key, val in tocompare.items():
        _, omissmatchCols, _, omissmatchInds, newNAs, new0s = h.compareDfs(
            files[key], tc.get(name=TAIGA_ETERNAL, file=val)
        )
        print(key)
        print(len(omissmatchCols), "columns NOT IN df1")
        print(len(omissmatchInds), "indices NOT IN df1")
        print(newNAs, "new NAs")
        print(new0s, "New 0s")

    print("updating the tracker")

    updated_tracker = expressions.updateTracker(
        set(renaming.keys()) - set(["transcript_id(s)"]),
        failed,
        lowqual[lowqual.sum(1) > 3].index.tolist(),
        ccle_refsamples,
        samplesetname,
        refworkspace,
        samplesinset=samplesinset,
        starlogs=starlogs,
        trackerobj=trackerobj,
        todrop=todrop,
        dry_run=dry_run,
    )

    if not dry_run:
        print("uploading to taiga")
        tc.update_dataset(
            changes_description="new " + samplesetname + " release!",
            dataset_permaname=taiga_dataset,
            upload_files=[
                {
                    "path": folder + "proteincoding_genes_tpm_logp1.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm_logp1.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm_logp1.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_tpm.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_expected_count.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_expected_count.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_expected_count.csv",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                }
                # {
                #     "path": folder+'gene_sets_all.csv',
                #     "format": "NumericMatrixCSV",
                #     "encoding": "utf-8"
                # },
            ],
            upload_async=False,
            dataset_description=dataset_description,
        )
        print("done")

    return updated_tracker


async def fusionPostProcessing(
    refworkspace=RNAWORKSPACE,
    sampleset=SAMPLESETNAME,
    fusionSamplecol=SAMPLEID,
    todrop=KNOWN_DROP,
    taiga_dataset=TAIGA_FUSION,
    dataset_description=FUSIONreadme,
    prevdataset="ccle",
    trackerobj=None,
    **kwargs,
):
    """the full CCLE Fusion post processing pipeline (used only by CCLE)

    see postprocessing() to reproduce our analysis

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether and where to save our data. Defaults to "".
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        fusionSamplecol ([type], optional): [description]. Defaults to SAMPLEID.
        taiga_dataset (str, optional): the taiga dataset path to use for uploading results. Defaults to TAIGA_EXPRESSION.
        dataset_description (str, optional): the taiga dataset description to use. Defaults to RNAseqreadme.
        sheetcreds (str, optional): path to the google sheet credentials file to use. Defaults to SHEETCREDS.
        refsheet_url (str, optional): the url of the google sheet containing the data. Defaults to REFSHEET_URL.
        my_id (str, optional): path to the id containing file for google sheet. Defaults to MY_ID.
        mystorage_id (str, optional): path to the id containing file for google storage. Defaults to MYSTORAGE_ID.
        prevdataset (str, optional): the previous dataset to use for the taiga upload. Defaults to 'ccle'.
    
    Returns:
        (pd.df): fusion dataframe
        (pd.df): filtered fusion dataframe
    """
    from taigapy import TaigaClient

    tc = TaigaClient()

    if prevdataset is "ccle":
        prevdataset = tc.get(name=TAIGA_ETERNAL, file="CCLE_fusions_unfiltered")

    ccle_refsamples = trackerobj.read_tracker()

    previousQCfail = ccle_refsamples[ccle_refsamples.low_quality == 1].index.tolist()

    folder = os.path.join("temp", sampleset, "")
    renaming = h.fileToDict(folder + "rna_sample_renaming.json")
    # TODO: include in rna_sample_renaming.json instead
    # lower priority versions of these lines were used

    fusions, _ = fusion.postProcess(
        refworkspace,
        todrop=previousQCfail,
        renaming=renaming,
        save_output=folder,
        **kwargs,
    )

    print("comparing to previous version")
    print("new")
    print(set(fusions[fusionSamplecol]) - set(prevdataset[fusionSamplecol]))

    print("removed")
    print(set(prevdataset[fusionSamplecol]) - set(fusions[fusionSamplecol]))

    print("changes in fusion names")
    pf = prevdataset.copy()
    pf["id"] = pf[fusionSamplecol] + "_" + pf["FusionName"]
    f = fusions.copy()
    f["id"] = f[fusionSamplecol] + "_" + f["FusionName"]
    print(len(set(pf[~pf.id.isin(f.id.tolist())][fusionSamplecol])))

    print("changes in junction read counts")
    f["sid"] = (
        f[fusionSamplecol]
        + "_"
        + f["FusionName"]
        + "_"
        + f["JunctionReadCount"].astype(str)
    )
    pf["sid"] = (
        pf[fusionSamplecol]
        + "_"
        + pf["FusionName"]
        + "_"
        + pf["JunctionReadCount"].astype(str)
    )
    print(len(set(pf[~pf.sid.isin(f.sid.tolist())][fusionSamplecol])))

    # taiga
    print("uploading to taiga")
    tc.update_dataset(
        dataset_permaname=taiga_dataset,
        changes_description="new " + sampleset + " release!",
        upload_files=[
            {
                "path": "temp/" + sampleset + "/fusions_latest.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/" + sampleset + "/filteredfusions_latest.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/" + sampleset + "/fusions_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        dataset_description=dataset_description,
    )
    print("done")
    return fusions


def cnPostProcessing(
    wesrefworkspace=WESCNWORKSPACE,
    wgsrefworkspace=WGSWORKSPACE,
    wessetentity=WESSETENTITY,
    wgssetentity=WGSSETENTITY,
    samplesetname=SAMPLESETNAME,
    trackerobj=None,
    AllSamplesetName="all",
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
    **kwargs,
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

    tracker = pd.DataFrame()
    if trackerobj is not None:
        tracker = trackerobj.read_tracker()

    assert len(tracker) != 0, "broken source for sample tracker"

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
        wessegments, genecn, wesfailed = cn.postProcess(
            wesrefworkspace,
            setEntity=wessetentity,
            sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
            todrop=todropwes,
            save_output=folder,
            priority=priority,
            **kwargs,
        )

        wesrenaming = cn.managingDuplicates(
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
    wgssegments, genecn, wgsfailed = cn.postProcess(
        wgsrefworkspace,
        setEntity=wgssetentity,
        sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
        todrop=todropwgs,
        save_output=folder,
        segmentsthresh=2000,
        priority=priority,
        **kwargs,
    )

    wgsrenaming = cn.managingDuplicates(
        set(wgssegments[SAMPLEID]),
        (set(wgsfailed) - set(priority)) | set(todropwgs),
        "wgs",
        tracker,
    )

    h.dictToFile(wgsrenaming, folder + "sample_renaming.json")

    # annotating source
    for v in set(wgssegments[SAMPLEID]):
        wgssegments.loc[
            wgssegments[wgssegments[SAMPLEID] == v].index, "Source"
        ] = tracker.loc[v, "source"]

    wgssegments.Source = wgssegments.Source.replace(source_rename)
    wgssegments.Source += " WGS"

    print("renaming")
    wgspriosegments = (
        wgssegments[wgssegments[SAMPLEID].isin(set(wgsrenaming.keys()))]
        .replace({SAMPLEID: wgsrenaming})
        .reset_index(drop=True)
    )
    wgspriogenecn = genecn[genecn.index.isin(set(wgsrenaming.keys()))].rename(
        index=wgsrenaming
    )
    # saving prio
    wgspriosegments.to_csv(folder + "segments_latest.csv", index=False)
    wgspriogenecn.to_csv(folder + "genecn_latest.csv")

    print("comparing to previous version")
    h.compareDfs(wespriogenecn, prevgenecn)

    # adding to the sample tracker the sequencing that were selected and the ones that failed QC
    selected = {i for i, j in wgsrenaming.items()}
    selected.update({i for i, j in wesrenaming.items()})

    try:
        cn.updateTracker(
            tracker,
            selected,
            samplesetname,
            wgsfailed,
            datatype=["wgs", "wes"],
            trackerobj=trackerobj,
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wgsrefworkspace,
        )
    except:
        print("no wgs for this sampleset")

    try:
        cn.updateTracker(
            tracker,
            selected,
            samplesetname,
            list(wesfailed),
            datatype=["wes", "wgs"],
            trackerobj=trackerobj,
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wesrefworkspace,
        )
    except:
        print("no wes for this sampleset")

    # merging WES/WGS
    folder = os.path.join("temp", samplesetname, "")
    mergedsegments = wgspriosegments.append(
        wespriosegments[~wespriosegments[SAMPLEID].isin(set(wgspriosegments[SAMPLEID]))]
    )[subsetsegs]
    mergedgenecn = wgspriogenecn.append(
        wespriogenecn[~wespriogenecn.index.isin(set(wgspriogenecn.index))]
    )

    mergedgenecn.to_csv(folder + "merged_genecn_all.csv")
    mergedsegments.to_csv(folder + "merged_segments_all.csv", index=False)

    # uploading to taiga
    print("uploading to taiga")
    tc.update_dataset(
        changes_description="new "
        + samplesetname
        + " release! (removed misslabellings, see changelog)",
        dataset_permaname=taiga_dataset,
        upload_files=[
            {
                "path": folder + "wes_segments_latest.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_genecn_latest.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_segments_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_genecn_all.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_genecn_all.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_segments_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wgs_segments_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wgs_genecn_all.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wgs_segments_latest.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wgs_genecn_latest.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        dataset_description=dataset_description,
    )
    print("done")
    return wespriosegments, wgspriosegments


def cnProcessForAchilles(
    wespriosegs,
    wgspriosegs,
    gene_mapping,
    cytobandloc=CYTOBANDLOC,
    samplesetname=SAMPLESETNAME,
    bad=[],
    taiga_legacy_loc=TAIGA_LEGACY_CN,
    taiga_legacy_filename="legacy_segments",
    taiga_dataset=TAIGA_CN_ACHILLES,
    dataset_description=Achillesreadme,
    prevsegments="ccle",
    prevgenecn="ccle",
    gene_expected_count="ccle",
):
    """runs the Achilles specific part of the pipeline

    Args:
        wespriosegs (pd.dataframe): set of wes to priorise, output of _CCLEPostProcessing
        wgspriosegs ([type]): set of wgs to priorise, output of _CCLEPostProcessing
        gene_mapping (str): 'biomart' or local bed file containing gene names
        samplesetname ([type], optional): . Defaults to SAMPLESETNAME.
        bad (list, optional): list of known samples that should not be inclusded.
        taiga_legacy_loc (str, optional): where the legacy segments file lies. Defaults to 'depmap-wes-cn-data--08f3'.
        taiga_legacy_filename (str, optional): what the legacy segments file is named. Defaults to 'legacy_segments'.
        taiga_dataset (str, optional): where we should upload the output on taiga to. Defaults to "cn-wes-achilles-4dcd".
        dataset_description ([type], optional): README to add to the taiga dataset. Defaults to Achillesreadme.
        cytobandloc (str, optional): bed file containing genomic chromosomal regions to extend segments (needed for Achilles). Defaults to 'data/hg38_cytoband.gz'.
        prevsegments (str, optional): if ccle, gets taiga, else can provide your own file instead of taiga's (used for QC). Defaults to "ccle".
        prevgenecn (str, optional): if ccle, gets taiga, else can provide your own file instead of taiga's (used for QC). Defaults to "ccle".
        gene_expected_count (str, optional): if ccle, gets taiga's expression file, else can provide your own file instead of taiga's (used for QC). Defaults to "ccle".
    """
    # load legacy_segments
    from taigapy import TaigaClient

    tc = TaigaClient()

    if prevsegments == "ccle":
        prevsegments = tc.get(name=TAIGA_ETERNAL, file="CCLE_segment_cn")
    if prevgenecn == "ccle":
        prevgenecn = (2 ** tc.get(name=TAIGA_ETERNAL, file="CCLE_gene_cn")) - 1
    if gene_expected_count == "ccle":
        gene_expected_count = tc.get(
            name=TAIGA_ETERNAL,
            file="CCLE_expression_proteincoding_genes_expected_count",
        )

    # load the legacy taiga dataset
    legacy_segments = tc.get(name=taiga_legacy_loc, file=taiga_legacy_filename).drop(
        columns="Unnamed: 0"
    )
    legacy_segments["Status"] = "U"
    legacy_segments["Chromosome"] = legacy_segments["Chromosome"].map(
        lambda x: x[3:] if x.startswith("chr") else x
    )

    data_sources = pd.crosstab(
        index=legacy_segments["DepMap_ID"], columns=legacy_segments["Source"]
    )
    data_source_duplicates = data_sources[(data_sources > 0).sum(axis=1) > 1]
    duplicate_arxspans = data_source_duplicates.index.tolist()

    # TODO: replace with assert if the data is updated on taiga
    assert (
        data_source_duplicates.empty
    ), "Duplicate data sources found in the taiga dataset"
    onlyinleg = set(legacy_segments[SAMPLEID]) - (
        set(wespriosegs[SAMPLEID]) | (set(wgspriosegs[SAMPLEID]))
    )
    # samegenes = set(prevgenecn.columns) & set(priogenecn.columns)

    onlyinleg = onlyinleg - set(bad)
    print("found samples that are only in the legacy datasets")
    print(onlyinleg)
    # merging
    print("merging wes/wgs/legacy")
    mergedsegments = (
        wespriosegs[~wespriosegs[SAMPLEID].isin(list(onlyinleg))]
        .append(legacy_segments[legacy_segments[SAMPLEID].isin(list(onlyinleg))])
        .reset_index(drop=True)
    )
    # adding wgs to wes
    mergedsegments = wgspriosegs.append(
        mergedsegments[~mergedsegments[SAMPLEID].isin(set(wgspriosegs[SAMPLEID]))]
    )

    mergedsegments = (
        mergedsegments[
            [
                SAMPLEID,
                "Chromosome",
                "Start",
                "End",
                "Segment_Mean",
                "Num_Probes",
                "Status",
                "Source",
            ]
        ]
        .sort_values(by=[SAMPLEID, "Chromosome", "Start", "End"])
        .reset_index(drop=True)
    )

    # setting amplification status to U for X chromosome as it is artificially
    # amplified in female samples:
    mergedsegments.loc[
        mergedsegments[mergedsegments.Chromosome == "X"].index, "Status"
    ] = "U"
    # making the gene cn matrix
    print("making the gene cn")
    cyto = pd.read_csv(
        cytobandloc, sep="\t", names=["chrom", "start", "end", "loc", "stains"]
    ).iloc[:-1]
    cyto["chrom"] = [i[3:] for i in cyto["chrom"]]

    gene_mapping_df = pd.DataFrame()
    if gene_mapping == "biomart":
        mybiomart = h.generateGeneNames(
            ensemble_server=ENSEMBL_SERVER_V,
            useCache=False,
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
        mybiomart = mybiomart[
            mybiomart["Chromosome"].isin(set(mergedsegments["Chromosome"]))
        ]
        mybiomart = mybiomart[
            ~mybiomart.entrezgene_id.isna()
        ]  # dropping all nan entrez id cols
        gene_mapping_df = mybiomart.drop_duplicates("hgnc_symbol", keep="first")
        gene_mapping_df["gene_name"] = [
            i["hgnc_symbol"] + " (" + str(i["entrezgene_id"]).split(".")[0] + ")"
            for _, i in gene_mapping_df.iterrows()
        ]

    else:
        gene_mapping["Chromosome"] = gene_mapping["Chromosome"].astype(str)
        gene_mapping = gene_mapping.sort_values(by=["Chromosome", "start", "end"])
        gene_mapping = gene_mapping[
            gene_mapping["Chromosome"].isin(set(mergedsegments["Chromosome"]))
        ]
        gene_mapping["gene_name"] = [
            i["symbol"] + " (" + str(i["ensembl_id"]).split(".")[0] + ")"
            for _, i in gene_mapping.iterrows()
        ]
        gene_mapping_df = gene_mapping.rename(columns={"ensembl_id": "entrezgene_id"})

    mergedsegments = mut.manageGapsInSegments(mergedsegments, cyto=cyto)
    mergedgenecn = mut.toGeneMatrix(mergedsegments, gene_mapping_df,).apply(
        lambda x: np.log2(1 + x)
    )
    wesgenecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(wespriosegs), gene_mapping_df
    ).apply(lambda x: np.log2(1 + x))
    wesgenecn.to_csv("temp/" + samplesetname + "/wes_genecn_latest.csv")

    # some QC
    print("copy number change with previous release")
    cn.plotCNchanges(
        mergedgenecn,
        prevgenecn.apply(lambda x: np.log2(1 + x)),
        mergedsegments,
        prevsegments,
    )

    if mergedgenecn.values.max() > 100:
        print("\n\n\nTOO HIGH, not LOG2 transformed!")
    if len(mergedgenecn.index) > len(set(mergedgenecn.index)):
        print("Duplicate CL, not reprioritized well!")

    # computing relationship with RNAseq
    print("correlation with RNAseq:")
    _, ax = plt.subplots()
    rna.rnaseqcorrelation(
        mergedgenecn.fillna(0), gene_expected_count.fillna(0), ax, name="current"
    )
    rna.rnaseqcorrelation(
        prevgenecn[prevgenecn.index.isin(mergedgenecn.index.tolist())].fillna(0),
        gene_expected_count.fillna(0),
        ax,
        name="prev",
    )

    # TODO: compute sample specific correlation (see James' function)
    h.compareDfs(mergedgenecn, prevgenecn)
    # h.compareDfs(mergedsegments, tc.get(name='depmap-a0ab', file='CCLE_segment_cn'))

    # saving
    print("saving")
    mergedgenecn.to_csv("temp/" + samplesetname + "/achilles_gene_cn.csv")
    mergedsegments.to_csv(
        "temp/" + samplesetname + "/achilles_segment.csv", index=False
    )

    # saving to taiga
    print("uploading to taiga")

    tc.update_dataset(
        changes_description="updated to new "
        + samplesetname
        + " release! (updated from relabelling see google drive file for more info)",
        dataset_permaname=taiga_dataset,
        upload_files=[
            {
                "path": "temp/" + samplesetname + "/achilles_segment.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/" + samplesetname + "/achilles_gene_cn.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        dataset_description=dataset_description,
    )
    print("done")


async def mutationPostProcessing(
    wesrefworkspace=WESMUTWORKSPACE,
    wgsrefworkspace=WGSWORKSPACE,
    samplesetname=SAMPLESETNAME,
    todrop=KNOWN_DROP,
    AllSamplesetName="all",
    doCleanup=False,
    taiga_description=Mutationsreadme,
    taiga_dataset=TAIGA_MUTATION,
    mutation_groups=MUTATION_GROUPS,
    tokeep_wes=RESCUE_FOR_MUTATION_WES,
    tokeep_wgs=RESCUE_FOR_MUTATION_WGS,
    prev="ccle",
    minfreqtocall=0.25,
    **kwargs,
):
    """the full CCLE mutations post processing pipeline (used only by CCLE)

    see postprocess() to reproduce our analysis

    Args:
        wesrefworkspace (str, optional): the reference workspace for WES. Defaults to WESMUTWORKSPACE.
        wgsrefworkspace (str, optional): the reference workspace for WGS. Defaults to WGSWORKSPACE.
        samplesetname (str, optional): the sample set name to use (for the release). Defaults to SAMPLESETNAME.
        AllSamplesetName (str, optional): the sample set to use for all samples. Defaults to 'all'.
        doCleanup (bool, optional): whether to clean up the workspace. Defaults to False.
        my_id (str, optional): filepath for google cloud access id file. Defaults to MY_ID.
        mystorage_id (str, optional): filepath to a google cloud storage access file. Defaults to MYSTORAGE_ID.
        refsheet_url (str, optional): path to the sample tracker google sheet. Defaults to REFSHEET_URL.
        taiga_description (str, optional): description of the dataset on taiga. Defaults to Mutationsreadme.
        taiga_dataset (str, optional): taiga folder location. Defaults to TAIGA_MUTATION.
        mutation_groups (dict, optional): a dict to group mutations annotations into bigger groups. Defaults to MUTATION_GROUPS.
        tokeep_wes (dict, optional): a dict of wes lines that are blacklisted on the tracker due to CN qc but we want to keep their mutation data. Defaults to RESCUE_FOR_MUTATION_WES.
        tokeep_wgs (dict, optional): a dict of wgs lines that are blacklisted on the tracker due to CN qc but we want to keep their mutation data. Defaults to RESCUE_FOR_MUTATION_WGS.
        prev (pd.df, optional): the previous release dataset (to do QC). 
            Defaults to ccle =>(tc.get(name=TAIGA_ETERNAL, file='CCLE_mutations')).
        minfreqtocall (float, optional): the minimum frequency to call a mutation. Defaults to 0.25.
    """
    from taigapy import TaigaClient

    tc = TaigaClient()
    if prev == "ccle":
        prev = tc.get(name=TAIGA_ETERNAL, file="CCLE_mutations")
    if doCleanup:
        # TODO:
        val = ""
        # gcp.rmFiles('gs://fc-secure-012d088c-f039-4d36-bde5-ee9b1b76b912/$val/**/call-tumorMM_Task/*.cleaned.bam')
    # sometimes it does not work so better check again

    # doing wes
    print("doing wes")
    folder = os.path.join("temp", samplesetname, "wes_")

    wesmutations = mutations.postProcess(
        wesrefworkspace,
        AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        doCleanup=True,
        mutCol="CGA_WES_AC",
        **kwargs,
    )

    # renaming
    print("renaming")
    wesrefwm = dm.WorkspaceManager(wesrefworkspace)
    wesrenaming = track.removeOlderVersions(
        names=set(wesmutations[SAMPLEID]),
        refsamples=wesrefwm.get_samples(),
        arxspan_id="arxspan_id",
        version="version",
    )

    wesrenaming = h.fileToDict(folder + "sample_renaming.json")
    wesrenaming.update(tokeep_wes)

    wesmutations = wesmutations[
        wesmutations[SAMPLEID].isin(wesrenaming.keys())
    ].replace({SAMPLEID: wesrenaming})
    wesmutations.to_csv(folder + "wes_somatic_mutations_all.csv", index=None)

    # doing wgs
    print("doing wgs")
    folder = os.path.join("temp", samplesetname, "wgs_")

    wgsmutations = mutations.postProcess(
        wgsrefworkspace,
        sampleset="allcurrent",  # AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        doCleanup=True,
        mutCol="CGA_WES_AC",
        **kwargs,
    )

    # renaming
    print("renaming")
    wgsrefwm = dm.WorkspaceManager(wgsrefworkspace)
    wgsrenaming = track.removeOlderVersions(
        names=set(wesmutations[SAMPLEID]),
        refsamples=wgsrefwm.get_samples(),
        arxspan_id="arxspan_id",
        version="version",
    )

    wgsrenaming = h.fileToDict(folder + "sample_renaming.json")
    wgsrenaming.update(tokeep_wgs)

    wgsmutations = wgsmutations[
        wgsmutations[SAMPLEID].isin(wgsrenaming.keys())
    ].replace({SAMPLEID: wgsrenaming})
    wgsmutations.to_csv(folder + "wgs_somatic_mutations_all.csv", index=None)

    # merge
    print("merging")
    folder = os.path.join("temp", samplesetname, "merged_")
    toadd = set(wgsmutations[SAMPLEID]) - set(wesmutations[SAMPLEID])
    priomutations = wesmutations.append(
        wgsmutations[wgsmutations[SAMPLEID].isin(toadd)]
    ).reset_index(drop=True)
    # normals = set(ccle_refsamples[ccle_refsamples.primary_disease=="normal"].arxspan_id)
    # mutations = mutations[~mutations[SAMPLEID].isin(normals)]
    priomutations.to_csv(folder + "somatic_mutations.csv", index=False)

    # making binary mutation matrices
    print("creating mutation matrices")
    # binary mutations matrices
    mut.mafToMat(
        priomutations[(priomutations.isDeleterious)], minfreqtocall=minfreqtocall
    ).astype(int).T.to_csv(folder + "somatic_mutations_boolmatrix_deleterious.csv")
    mut.mafToMat(
        priomutations[
            ~(
                priomutations.isDeleterious
                | priomutations.isCOSMIChotspot
                | priomutations.isTCGAhotspot
                | priomutations["Variant_Classification"]
                == "Silent"
            )
        ],
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(folder + "somatic_mutations_boolmatrix_other.csv")
    mut.mafToMat(
        priomutations[(priomutations.isCOSMIChotspot | priomutations.isTCGAhotspot)],
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(folder + "somatic_mutations_boolmatrix_hotspot.csv")

    # genotyped mutations matrices
    mut.mafToMat(
        priomutations[(priomutations.isDeleterious)], mode="genotype",
    ).T.to_csv(folder + "somatic_mutations_matrix_deleterious.csv")
    mut.mafToMat(
        priomutations[
            ~(
                priomutations.isDeleterious
                | priomutations.isCOSMIChotspot
                | priomutations.isTCGAhotspot
                | priomutations["Variant_Classification"]
                == "Silent"
            )
        ],
        mode="genotype",
    ).T.to_csv(folder + "somatic_mutations_matrix_other.csv")
    mut.mafToMat(
        priomutations[(priomutations.isCOSMIChotspot | priomutations.isTCGAhotspot)],
        mode="genotype",
    ).T.to_csv(folder + "somatic_mutations_matrix_hotspot.csv")

    # adding lgacy datasetss
    print("add legacy datasets")
    legacy_hybridcapture = tc.get(name="mutations-da6a", file="legacy_hybridcapture")
    legacy_raindance = tc.get(name="mutations-da6a", file="legacy_raindance")
    legacy_rna = tc.get(name="mutations-da6a", file="legacy_rna")
    legacy_wes_sanger = tc.get(name="mutations-da6a", file="legacy_wes_sanger")
    legacy_wgs_exoniconly = tc.get(name="mutations-da6a", file="legacy_wgs_exoniconly")

    merged = mut.mergeAnnotations(
        priomutations,
        legacy_hybridcapture,
        "HC_AC",
        useSecondForConflict=True,
        dry_run=False,
    )
    merged = mut.mergeAnnotations(
        merged, legacy_raindance, "RD_AC", useSecondForConflict=True, dry_run=False
    )
    merged = mut.mergeAnnotations(
        merged,
        legacy_wgs_exoniconly,
        "WGS_AC",
        useSecondForConflict=False,
        dry_run=False,
    )
    merged = mut.mergeAnnotations(
        merged,
        legacy_wes_sanger,
        "SangerWES_AC",
        useSecondForConflict=False,
        dry_run=False,
    )
    merged = mut.mergeAnnotations(
        merged, legacy_rna, "RNAseq_AC", useSecondForConflict=False, dry_run=False
    )

    merged = merged[merged["tumor_f"] > 0.05]
    merged = mutations.annotateLikelyImmortalized(
        merged,
        TCGAlocs=["TCGAhsCnt", "COSMIChsCnt"],
        max_recurrence=0.05,
        min_tcga_true_cancer=5,
    )
    print("changing variant annotations")
    rename = {}
    for k, v in mutation_groups.items():
        for e in v:
            rename[e] = k
    merged["Variant_annotation"] = [
        rename[i] for i in merged["Variant_Classification"].tolist()
    ]

    print("compare to previous release")
    a = set(merged[SAMPLEID])
    b = set(prev[SAMPLEID])
    print("new lines:")
    print(a - b)
    print("lost lines:")
    print(b - a)

    # making a depmap version
    # removing immortalized ffor now
    merged = merged[~merged["is_likely_immortalization"]]
    # reverting to previous versions
    merged_maf = merged[MUTCOL_DEPMAP].rename(
        columns={"Tumor_Allele": "Alternate_Allele"}
    )
    merged_maf.to_csv(folder + "somatic_mutations_withlegacy.csv", index=False)

    # making binary matrices
    merged = merged[merged["Entrez_Gene_Id"] != 0]
    merged["mutname"] = (
        merged["Hugo_Symbol"] + " (" + merged["Entrez_Gene_Id"].astype(str) + ")"
    )
    mut.mafToMat(
        merged[(merged.Variant_annotation == "damaging")],
        mode="bool",
        mutNameCol="mutname",
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(
        folder + "somatic_mutations_boolmatrix_fordepmap_damaging.csv"
    )
    mut.mafToMat(
        merged[(merged.Variant_annotation == "other conserving")],
        mode="bool",
        mutNameCol="mutname",
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(
        folder + "somatic_mutations_boolmatrix_fordepmap_othercons.csv"
    )
    mut.mafToMat(
        merged[(merged.Variant_annotation == "other non-conserving")],
        mode="bool",
        mutNameCol="mutname",
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(
        folder + "somatic_mutations_boolmatrix_fordepmap_othernoncons.csv"
    )
    mut.mafToMat(
        merged[(merged.isCOSMIChotspot | merged.isTCGAhotspot)],
        mode="bool",
        mutNameCol="mutname",
        minfreqtocall=minfreqtocall,
    ).astype(int).T.to_csv(
        folder + "somatic_mutations_boolmatrix_fordepmap_hotspot.csv"
    )

    # uploading to taiga
    tc.update_dataset(
        changes_description="new " + samplesetname + " release!",
        dataset_permaname=taiga_dataset,
        upload_files=[
            # for depmap
            # {
            #     "path": folder+"somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8"
            # },
            # {
            #     "path": folder+"somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8"
            # },
            # {
            #     "path": folder+"somatic_mutations_boolmatrix_fordepmap_damaging.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8"
            # },
            # genotyped
            {
                "path": folder + "somatic_mutations_matrix_hotspot.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_matrix_other.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_matrix_deleterious.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            # new
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder
                + "somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_othercons.csv",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_withlegacy.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/" + samplesetname + "/wes_somatic_mutations_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/" + samplesetname + "/wgs_somatic_mutations_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        upload_async=False,
        dataset_description=taiga_description,
    )


async def mutationAnalyzeUnfiltered(
    workspace=WGSWORKSPACE,
    allsampleset="all",
    folder="temp/",
    subsetcol=[
        SAMPLEID,
        "Hugo_Symbol",
        "Entrez_Gene_Id",
        "Chromosome",
        "Start_position",
        "End_position",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Allele",
        "dbSNP_RS",
        "dbSNP_Val_Status",
        "Genome_Change",
        "Annotation_Transcript",
        "cDNA_Change",
        "Codon_Change",
        "HGVS_protein_change",
        "Protein_Change",
        "t_alt_count",
        "t_ref_count",
        "tumor_f",
        "CGA_WES_AC",
    ],
    taiga_dataset=TAIGA_MUTATION,
):
    """_CCLEAnalyzeUnfiltered function to subset and filter the CGA unfiltered maf file.

    This will output a much bigger maf file without CGA filters (usefull for QC and more).
    Will take a lot of memory (expect ~32GB minimum). if you don't have that amount of RAM, don't use.

    Args:
        workspace (str): workspace name. Default is WGSWORKSPACE.
        allsampleset (str, optional): sampleset name. Default is 'all'.
        folder (str, optional): folder name. Default is 'temp/'.
        subsetcol (list, optional): list of columns to subset the maf file on. 
            will also output the unfiltered version of themaf file.
            Defaults to [SAMPLEID, 'Hugo_Symbol', 'Entrez_Gene_Id', 
                        'Chromosome', 'Start_position', 'End_position', 
                        'Variant_Classification', 'Variant_Type', 'Reference_Allele', 
                        'Tumor_Allele', 'dbSNP_RS', 'dbSNP_Val_Status', 'Genome_Change', 
                        'Annotation_Transcript', 'cDNA_Change', 'Codon_Change', 
                        'HGVS_protein_change',  'Protein_Change', 't_alt_count', 
                        't_ref_count', 'tumor_f', 'CGA_WES_AC'].
        taiga_dataset (str, optional): taiga dataset path. Default is TAIGA_MUTATION.
    """
    print("retrieving unfiltered")
    ####### WES
    from taigapy import TaigaClient

    tc = TaigaClient()
    res = dm.WorkspaceManager(workspace).get_sample_sets()
    unfiltered = pd.read_csv(
        res.loc[allsampleset, "unfiltered_CGA_MAF_aggregated"],
        sep="\t",
        encoding="L6",
        na_values=["__UNKNOWN__", "."],
        engine="c",
        dtype=str,
    )
    unfiltered["somatic"] = unfiltered["somatic"].replace("nan", "False")
    unfiltered["HGNC_Status"] = unfiltered["HGNC_Status"].replace("nan", "Unapproved")
    unfiltered["judgement"] = unfiltered["judgement"].replace("nan", "REMOVE")
    unfiltered = unfiltered.rename(
        columns={
            "i_ExAC_AF": "ExAC_AF",
            "Tumor_Sample_Barcode": SAMPLEID,
            "Tumor_Seq_Allele2": "Tumor_Allele",
        }
    ).drop(columns=["Tumor_Seq_Allele1"])
    unfiltered["CGA_WES_AC"] = [
        str(i[0]) + ":" + str(i[1])
        for i in np.nan_to_num(
            unfiltered[["t_alt_count", "t_ref_count"]].values.astype(float), 0
        ).astype(int)
    ]
    toremove = []
    subunfilt = unfiltered.iloc[:10000]
    for i, val in enumerate(unfiltered.columns):
        h.showcount(i, len(unfiltered.columns))
        if len(set(subunfilt[val]) - set(["nan"])) == 1:
            if len(set(unfiltered[val]) - set(["nan"])) == 1:
                toremove.append(val)
    unfiltered = unfiltered.drop(columns=set(toremove))
    toint = ["Start_position", "End_position"]
    for val in toint:
        unfiltered[val] = unfiltered[val].astype(int)
    unfiltered.to_csv(
        folder + "mutation_somatic_unfiltered_withreplicates.csv.gz", index=False
    )
    unfiltered = unfiltered[subsetcol]
    unfiltered.to_csv(
        folder + "mutation_somatic_unfiltered_withreplicates_subseted.csv.gz",
        index=False,
    )
    os.system("gunzip " + folder + "mutation_somatic_unfiltered_withreplicates.csv.gz")
    del unfiltered
    tc.update_dataset(
        changes_description="adding unfiltered mutations",
        dataset_permaname=taiga_dataset,
        upload_files=[
            {
                "path": folder
                + "mutation_somatic_unfiltered_withreplicates_subseted.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "mutation_somatic_unfiltered_withreplicates.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
        upload_async=False,
    )

