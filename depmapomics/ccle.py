import os.path
import dalmatian as dm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat

from genepy.utils import helper as h
from genepy import rna
from genepy import mutations as mut

from depmapomics import tracker as track
from depmapomics import expressions, mutations
from depmapomics import fusions as fusion
from depmapomics import copynumbers as cn
from depmapomics.config import *

import multiprocessing
from genepy.epigenetics import chipseq as chip


async def expressionPostProcessing(
    trackerobj,
    gumbo,
    refworkspace=RNAWORKSPACE,
    samplesetname=SAMPLESETNAME,
    ensemblserver=ENSEMBL_SERVER_V,
    geneLevelCols=RSEMFILENAME_GENE,
    trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
    rsemfilename=RSEMFILENAME,
    ssgseafilepath=SSGSEAFILEPATH,
    colstoclean=["fastq1", "fastq2", "recalibrated_bam", "recalibrated_bam_index"],
    doCleanup=True,
    samplesetToLoad="all",
    tocompare={
        "genes_expected_count": "CCLE_RNAseq_reads",
        "genes_tpm": "CCLE_expression_full",
        "proteincoding_genes_tpm": "CCLE_expression",
    },
    todrop=KNOWN_DROP,
    prevcounts="ccle",
    taiga_dataset=TAIGA_EXPRESSION,
    minsimi=RNAMINSIMI,
    dropNonMatching=True,
    dataset_description=RNAseqreadme,
    dry_run=False,
    samplesinset=[],
    rsemfilelocs=None,
    rnaqclocs={},
    starlogs={},
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
    if gumbo:
        ccle_refsamples = trackerobj.read_seq_table()

    todrop += ccle_refsamples[
        (ccle_refsamples.blacklist == 1) & (ccle_refsamples.datatype == "rna")
    ].index.tolist()
    priority = ccle_refsamples[
        (ccle_refsamples.prioritized == 1) & (ccle_refsamples.datatype == "rna")
    ].index.tolist()

    folder = os.path.join("output", samplesetname, "")

    if dry_run:
        folder = os.path.join("output", "dryrun", "")

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
        rsemfilename=rsemfilename,
        ssgseafilepath=ssgseafilepath,
        samplesetToLoad=samplesetToLoad,
        geneLevelCols=geneLevelCols,
        trancriptLevelCols=trancriptLevelCols,  # compute_enrichment=False,
        ssGSEAcol="genes_tpm",
        dropNonMatching=dropNonMatching,
        dry_run=dry_run,
        samplesinset=samplesinset,
        rsemfilelocs=rsemfilelocs,
        rnaqclocs=rnaqclocs,
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
        refworkspace=refworkspace,
        onlycol=STARBAMCOLTERRA,
        newgs=RNA_GCS_PATH_HG38,
        samplesinset=samplesinset,
        starlogs=starlogs,
        trackerobj=trackerobj,
        todrop=todrop,
        dry_run=True,
    )

    # subset and rename, include all PRs that have associated CDS-ids
    pr_table = track.update_pr_from_seq(trackerobj)
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
    h.dictToFile(renaming_dict, folder + "rna_seq2pr_renaming.json")
    pr_files = dict()
    for k, v in files.items():
        pr_files[k + "_profile"] = v[v.index.isin(set(renaming_dict.keys()))].rename(
            index=renaming_dict
        )
    expressions.saveFiles(pr_files, folder)

    if not dry_run:
        print("uploading to taiga")
        tc.update_dataset(
            changes_description="new " + samplesetname + " release!",
            dataset_permaname=taiga_dataset,
            upload_files=[
                {
                    "path": folder + "proteincoding_genes_tpm_logp1.csv",
                    "name": "proteincoding_genes-tpm_logp1",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm_logp1.csv",
                    "name": "transcripts-tpm_logp1",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm_logp1.csv",
                    "name": "genes-tpm_logp1",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm.csv",
                    "name": "genes-tpm_logp1",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm.csv",
                    "name": "transcripts-tpm",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_tpm.csv",
                    "name": "proteincoding_genes-tpm",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_expected_count.csv",
                    "name": "transcripts-expected_count",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_expected_count.csv",
                    "name": "proteincoding_genes-expected_count",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_expected_count.csv",
                    "name": "genes-expected_count",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_tpm_profile_logp1.csv",
                    "name": "proteincoding_genes-tpm_logp1-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm_profile_logp1.csv",
                    "name": "transcripts-tpm_logp1-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm_profile_logp1.csv",
                    "name": "genes-tpm_logp1-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm_profile.csv",
                    "name": "genes-tpm-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm_profile.csv",
                    "name": "transcripts-tpm-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_tpm_profile.csv",
                    "name": "proteincoding_genes-tpm-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_expected_count_profile.csv",
                    "name": "transcripts-expected_count-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_expected_count_profile.csv",
                    "name": "proteincoding_genes-expected_count-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_expected_count_profile.csv",
                    "name": "genes-expected_count-profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                # {
                #     "path": folder+'gene_sets_all.csv',
                #     "name": "gene_sets",
                #     "format": "NumericMatrixCSV",
                #     "encoding": "utf-8"
                # },
                # {
                #     "path": folder+'gene_sets_profile.csv',
                #     "name": "gene_sets-profile",
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
    trackerobj,
    refworkspace=RNAWORKSPACE,
    sampleset=SAMPLESETNAME,
    fusionSamplecol=SAMPLEID,
    todrop=KNOWN_DROP,
    taiga_dataset=TAIGA_FUSION,
    dataset_description=FUSIONreadme,
    prevdataset="ccle",
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

    ccle_refsamples = trackerobj.read_tracker()

    previousQCfail = ccle_refsamples[ccle_refsamples.low_quality == 1].index.tolist()

    folder = os.path.join("output", sampleset, "")
    # TODO: include in rna_sample_renaming.json instead
    # lower priority versions of these lines were used

    fusions, fusions_filtered = fusion.postProcess(
        refworkspace, todrop=previousQCfail, save_output=folder, **kwargs,
    )

    # subset, rename from seqid to prid, and save pr-indexed matrices
    pr_table = trackerobj.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
    fusions_pr = fusions[
        fusions[fusionSamplecol].isin(set(renaming_dict.keys()))
    ].replace({fusionSamplecol: renaming_dict})
    fusions_filtered_pr = fusions_filtered[
        fusions_filtered[fusionSamplecol].isin(set(renaming_dict.keys()))
    ].replace({fusionSamplecol: renaming_dict})

    fusions_pr.to_csv(os.path.join(folder, "fusions_all_profile.csv"), index=False)
    fusions_filtered_pr.to_csv(
        os.path.join(folder, "filteredfusions_latest_profile.csv"), index=False
    )

    if prevdataset == "ccle":
        prevdataset = tc.get(name=TAIGA_ETERNAL, file="CCLE_fusions_unfiltered")

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
                "path": "output/" + sampleset + "/filteredfusions_latest.csv",
                "name": "filtered_fusions",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/" + sampleset + "/fusions_all.csv",
                "name": "fusions-all",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/" + sampleset + "/filteredfusions_latest_profile.csv",
                "name": "filtered_fusion-profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/" + sampleset + "/fusions_all_profile.csv",
                "name": "fusion-profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        dataset_description=dataset_description,
    )
    print("done")
    return fusions


def cnPostProcessing(
    trackerobj,
    gumbo,
    wesrefworkspace=WESCNWORKSPACE,
    wgsrefworkspace=WGSWORKSPACE,
    wessetentity=WESSETENTITY,
    wgssetentity=WGSSETENTITY,
    samplesetname=SAMPLESETNAME,
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
    genechangethresh=GENECHANGETHR,
    segmentsthresh=SEGMENTSTHR,
    maxYchrom=MAXYCHROM,
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

    if prevgenecn == "ccle":
        prevgenecn = tc.get(name=TAIGA_ETERNAL, file="CCLE_gene_cn")

    tracker = pd.DataFrame()
    if trackerobj is not None:
        tracker = trackerobj.read_seq_table()

    assert len(tracker) != 0, "broken source for sample tracker"

    # doing wes
    folder = os.path.join("output", samplesetname, "wes_")
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
        (
            wessegments,
            genecn,
            wesfailed,
            wes_purecn_segments,
            wes_purecn_genecn,
            wes_loh,
            wes_purecn_table,
            _,
        ) = cn.postProcess(
            wesrefworkspace,
            setEntity=wessetentity,
            sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
            todrop=todropwes,
            save_output=folder,
            priority=priority,
            genechangethresh=genechangethresh,
            segmentsthresh=segmentsthresh,
            maxYchrom=maxYchrom,
            **kwargs,
        )

        # annotating source
        for v in set(wessegments[SAMPLEID]):
            wessegments.loc[
                wessegments[wessegments[SAMPLEID] == v].index, "Source"
            ] = tracker.loc[v, "source"]

        wessegments.Source = wessegments.Source.replace(source_rename)
        wessegments.Source += " WES"
        # wes_purecn_segments.to_csv(folder + "purecn_segments_latest.csv", index=False)
        # wes_purecn_genecn.to_csv(folder + "purecn_genecn_latest.csv")
        # wes_loh.to_csv(folder + "purecn_loh_latest.csv")
        # wes_purecn_table.to_csv(folder + "purecn_table_latest.csv")
        # subset and rename to PR-indexed matrices
        # assuming no wgs
        pr_table = trackerobj.read_pr_table()
        renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
        wessegments_pr = (
            wessegments[wessegments[SAMPLEID].isin(set(renaming_dict.keys()))]
            .replace({SAMPLEID: renaming_dict})
            .reset_index(drop=True)
        )
        # wessegments_pr.to_csv(folder + "segments_all_profile.csv", index=False)
        # wescn_pr.to_csv(folder + "genecn_all_profile.csv")
    else:
        print("bypassing WES using folder: " + (wesfolder if wesfolder else folder))
        wesfailed = h.fileToList((wesfolder if wesfolder else folder) + "failed.txt")
        wessegments = pd.read_csv(folder + "segments_all.csv")
        genecn = pd.read_csv(folder + "genecn_all.csv", index_col=0)
        wessegments_pr = (
            wessegments[wessegments[SAMPLEID].isin(set(renaming_dict.keys()))]
            .replace({SAMPLEID: renaming_dict})
            .reset_index(drop=True)
        )

    # doing wgs
    print("doing wgs")
    folder = os.path.join("output", samplesetname, "wgs_")
    priority = tracker[
        (tracker.datatype == "wgs") & (tracker.prioritized == 1)
    ].index.tolist()
    todropwgs = (
        todrop
        + tracker[(tracker.datatype == "wgs") & (tracker.blacklist == 1)].index.tolist()
    )
    (
        wgssegments,
        genecn,
        wgsfailed,
        wgs_purecn_segments,
        wgs_purecn_genecn,
        wgs_loh,
        wgs_purecn_table,
        mybiomart,
    ) = cn.postProcess(
        wgsrefworkspace,
        setEntity=wgssetentity,
        sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
        todrop=todropwgs,
        save_output=folder,
        priority=priority,
        genechangethresh=genechangethresh,
        segmentsthresh=segmentsthresh,
        maxYchrom=maxYchrom,
        **kwargs,
    )

    # annotating source
    for v in set(wgssegments[SAMPLEID]):
        wgssegments.loc[
            wgssegments[wgssegments[SAMPLEID] == v].index, "Source"
        ] = tracker.loc[v, "source"]

    wgssegments.Source = wgssegments.Source.replace(source_rename)
    wgssegments.Source += " WGS"

    # wgs_purecn_segments.to_csv(folder + "purecn_segments_latest.csv", index=False)
    # wgs_purecn_genecn.to_csv(folder + "purecn_genecn_latest.csv")
    # wgs_loh.to_csv(folder + "purecn_loh_latest.csv")
    # wgs_purecn_table.to_csv(folder + "purecn_table_latest.csv")

    # print("comparing to previous version")
    # h.compareDfs(wespriogenecn, prevgenecn)

    # with gumbo, no need to mark this selected field
    selected = []

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
            dry_run=True,
            gumbo=gumbo,
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
            dry_run=True,
            gumbo=gumbo,
        )
    except:
        print("no wes for this sampleset")

    pr_table = track.update_pr_from_seq(trackerobj)
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
    wgssegments_pr = (
        wgssegments[wgssegments[SAMPLEID].isin(set(renaming_dict.keys()))]
        .replace({SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wgscn_pr = genecn[genecn.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wgssegments_pr.to_csv(folder + "segments_all_profile.csv", index=False)
    wgscn_pr.to_csv(folder + "genecn_all_profile.csv")

    print("merging PR-level seg file")
    mergedsegments_pr = wgssegments_pr.append(wessegments_pr).reset_index(drop=True)
    mergedsegments_pr = (
        mergedsegments_pr[
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
    mergedsegments_pr.loc[
        mergedsegments_pr[mergedsegments_pr.Chromosome == "X"].index, "Status"
    ] = "U"
    # merging wes and wgs
    folder = os.path.join("output", samplesetname, "")
    mergedsegments_pr = mut.manageGapsInSegments(mergedsegments_pr)
    mergedsegments_pr.to_csv(folder + "merged_segments_all_profile.csv", index=False)

    mergedgenecn_pr = mut.toGeneMatrix(mergedsegments_pr, mybiomart).apply(
        lambda x: np.log2(1 + x)
    )
    mergedgenecn_pr.to_csv(folder + "merged_genecn_all_profile.csv")

    # uploading to taiga
    print("uploading to taiga")
    tc.update_dataset(
        changes_description="new "
        + samplesetname
        + " release! (removed misslabellings, see changelog)",
        dataset_permaname=taiga_dataset,
        upload_files=[
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
                "path": folder + "merged_segments_all_profile.csv",
                "name": "merged-segments-profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_genecn_all_profile.csv",
                "name": "merged-genecn-profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            # Pure CN outputs
            # {
            #     "path": folder + "wes_purecn_segments_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_genecn_latest.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_segments_all.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_genecn_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_loh_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_loh_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_table_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wes_purecn_table_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_segments_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_genecn_latest.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_segments_all.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_genecn_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_loh_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_loh_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_table_latest.csv",
            #     "format": "TableCSV",
            #     "encoding": "utf-8",
            # },
            # {
            #     "path": folder + "wgs_purecn_table_all.csv",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
        ],
        dataset_description=dataset_description,
    )
    print("done")
    return wessegments, wgssegments


async def mutationPostProcessing(
    trackerobj,
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
    folder = os.path.join("output", samplesetname, "wes_")

    wesmutations = mutations.postProcess(
        wesrefworkspace,
        AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        doCleanup=True,
        mutCol="CGA_WES_AC",
        **kwargs,
    )

    pr_table = trackerobj.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
    wesmutations_pr = wesmutations[
        wesmutations[SAMPLEID].isin(renaming_dict.keys())
    ].replace({SAMPLEID: renaming_dict})
    wesmutations_pr.to_csv(folder + "somatic_mutations_all_profile.csv", index=None)

    # doing wgs
    print("doing wgs")
    folder = os.path.join("output", samplesetname, "wgs_")

    wgsmutations = mutations.postProcess(
        wgsrefworkspace,
        sampleset="allcurrent",  # AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        doCleanup=True,
        mutCol="CGA_WES_AC",
        **kwargs,
    )

    wgsmutations_pr = wgsmutations[
        wgsmutations[SAMPLEID].isin(renaming_dict.keys())
    ].replace({SAMPLEID: renaming_dict})

    wgsmutations_pr.to_csv(folder + "somatic_mutations_all_profile.csv", index=None)

    # merge
    print("merging")
    folder = os.path.join("output", samplesetname, "merged_")
    priomutations = wgsmutations_pr.append(wesmutations_pr).reset_index(drop=True)
    # normals = set(ccle_refsamples[ccle_refsamples.primary_disease=="normal"].arxspan_id)
    # mutations = mutations[~mutations[SAMPLEID].isin(normals)]
    # priomutations.to_csv(folder + "somatic_mutations_profile.csv", index=False)

    # making binary mutation matrices
    print("creating mutation matrices")

    merged = priomutations[priomutations["tumor_f"] > 0.05]
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
    merged_maf.to_csv(folder + "somatic_mutations_fordepmap_profile.csv", index=False)

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
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                "name": "somatic_mutations_boolmatrix-hotspot-profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder
                + "somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                "name": "somatic_mutations_boolmatrix-othernoncons-profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                "name": "somatic_mutations_boolmatrix-damaging-profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            # {
            #     "path": folder + "somatic_mutations_boolmatrix_fordepmap_othercons.csv",
            #     "name": "somatic_mutations_boolmatrix-othercons-profile",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            {
                "path": folder + "somatic_mutations_fordepmap_profile.csv",
                "name": "somatic_mutations-profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/" + samplesetname + "/wes_somatic_mutations_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/" + samplesetname + "/wgs_somatic_mutations_all.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/"
                + samplesetname
                + "/wes_somatic_mutations_all_profile.csv",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "output/"
                + samplesetname
                + "/wgs_somatic_mutations_all_profile.csv",
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
    folder="output/",
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
        folder (str, optional): folder name. Default is 'output/'.
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
    # WES
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
