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
    todrop=KNOWN_DROP,
    taiga_dataset=TAIGA_EXPRESSION,
    save_output=WORKING_DIR,
    minsimi=RNAMINSIMI,
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

    mytracker = track.SampleTracker()

    ccle_refsamples = mytracker.read_seq_table()

    todrop += ccle_refsamples[
        (ccle_refsamples.blacklist == 1) & (ccle_refsamples.ExpectedType == "rna")
    ].index.tolist()
    priority = ccle_refsamples[
        (ccle_refsamples.prioritized == 1) & (ccle_refsamples.ExpectedType == "rna")
    ].index.tolist()

    folder = save_output + samplesetname + "/"

    if dry_run:
        folder = save_output + "dryrun/"

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
        dropNonMatching=dropNonMatching,
        dry_run=dry_run,
        samplesinset=samplesinset,
        rsemfilelocs=rsemfilelocs,
        rnaqclocs=rnaqclocs,
        **kwargs,
    )

    print("updating the tracker")

    track.updateTrackerRNA(
        set(renaming.keys()) - set(["transcript_id(s)"]),
        failed,
        lowqual[lowqual.sum(1) > 3].index.tolist(),
        ccle_refsamples,
        samplesetname,
        refworkspace,
        samplesinset=samplesinset,
        starlogs=starlogs,
        todrop=todrop,
        dry_run=True,
        newgs=None,
    )

    # subset and rename, include all PRs that have associated CDS-ids
    pr_table = track.update_pr_from_seq()
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
                    "path": folder + "transcripts_expected_count.csv",
                    "name": "transcripts_expectedCount_withReplicates",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_expected_count.csv",
                    "name": "genes_expectedCount_withReplicates",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_tpm_profile_logp1.csv",
                    "name": "proteinCoding_genes_tpm_logp1_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_tpm_profile_logp1.csv",
                    "name": "transcripts_tpm_logp1_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_tpm_profile_logp1.csv",
                    "name": "genes_tpm_logp1_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "transcripts_expected_count_profile.csv",
                    "name": "transcripts_expectedCount_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "proteincoding_genes_expected_count_profile.csv",
                    "name": "proteinCoding_genes_expectedCount_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_expected_count_profile.csv",
                    "name": "genes_expectedCount_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                # {
                #     "path": folder+'gene_sets_all.csv',
                #     "name": "gene_set_enrichment_withReplicates",
                #     "format": "NumericMatrixCSV",
                #     "encoding": "utf-8"
                # },
                # {
                #     "path": folder+'gene_sets_profile.csv',
                #     "name": "gene_set_enrichment_profile",
                #     "format": "NumericMatrixCSV",
                #     "encoding": "utf-8"
                # },
            ],
            upload_async=False,
            dataset_description=dataset_description,
        )
        print("done")


async def fusionPostProcessing(
    refworkspace=RNAWORKSPACE,
    sampleset=SAMPLESETNAME,
    fusionSamplecol=SAMPLEID,
    todrop=KNOWN_DROP,
    taiga_dataset=TAIGA_FUSION,
    dataset_description=FUSIONreadme,
    folder=WORKING_DIR,
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
    
    Returns:
        (pd.df): fusion dataframe
        (pd.df): filtered fusion dataframe
    """
    from taigapy import TaigaClient

    tc = TaigaClient()

    mytracker = track.SampleTracker()
    ccle_refsamples = mytracker.read_seq_table()

    previousQCfail = ccle_refsamples[ccle_refsamples.low_quality == 1].index.tolist()

    # TODO: include in rna_sample_renaming.json instead
    # lower priority versions of these lines were used

    fusions, fusions_filtered = fusion.postProcess(
        refworkspace, todrop=previousQCfail, save_output=folder, **kwargs,
    )

    # subset, rename from seqid to prid, and save pr-indexed matrices
    pr_table = mytracker.read_pr_table()
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

    # taiga
    print("uploading to taiga")
    tc.update_dataset(
        dataset_permaname=taiga_dataset,
        changes_description="new " + sampleset + " release!",
        upload_files=[
            {
                "path": folder + sampleset + "/fusions_all.csv",
                "name": "fusions_unfiltered_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + sampleset + "/filteredfusions_latest_profile.csv",
                "name": "fusions_filtered_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + sampleset + "/fusions_all_profile.csv",
                "name": "fusions_unfiltered_profile",
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
    AllSamplesetName="all",
    todrop=KNOWN_DROP,
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
    save_dir=WORKING_DIR,
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
        taiga_dataset (str, optional): where to save the output to on taiga. Defaults to TAIGA_CN.
        dataset_description (str, optional): A long string that will be pushed to taiga to explain the CN dataset. Defaults to CNreadme.
        subsetsegs (list[str], optional): what columns to keep for the segments. Defaults to [SAMPLEID, 'Chromosome', 'Start', 'End', 'Segment_Mean', 'Num_Probes', 'Status', 'Source'].
        bamqc ([type], optional): @see updateTracker. Defaults to BAMQC.
        procqc ([type], optional): @see updateTracker. Defaults to PROCQC.
        source_rename ([type], optional): @see managing duplicates. Defaults to SOURCE_RENAME.
    """
    from taigapy import TaigaClient

    tc = TaigaClient()

    mytracker = track.SampleTracker()
    tracker = mytracker.read_seq_table()

    assert len(tracker) != 0, "broken source for sample tracker"
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))

    # doing wes
    folder = save_dir + "wes_"
    if wesfolder == "":
        print("doing wes")
        todropwes = (
            todrop
            + tracker[
                (tracker.ExpectedType == "wes") & (tracker.blacklist == 1)
            ].index.tolist()
        )
        (
            wessegments,
            wesgenecn,
            wesfailed,
            wes_purecn_segments,
            wes_purecn_genecn,
            wes_loh,
            wes_feature_table,
        ) = cn.postProcess(
            wesrefworkspace,
            setEntity=wessetentity,
            sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
            todrop=todropwes,
            save_output=folder,
            genechangethresh=genechangethresh,
            segmentsthresh=segmentsthresh,
            maxYchrom=maxYchrom,
            **kwargs,
        )

        # wessegments_pr.to_csv(folder + "segments_all_profile.csv", index=False)
        # wescn_pr.to_csv(folder + "genecn_all_profile.csv")
    else:
        print("bypassing WES using folder: " + wesfolder)
        wesfailed = h.fileToList(wesfolder + "failed.txt")
        wessegments = pd.read_csv(wesfolder + "segments_all.csv")
        wesgenecn = pd.read_csv(wesfolder + "genecn_all.csv", index_col=0)
        wes_purecn_segments = pd.read_csv(wesfolder + "purecn_segments_all.csv")
        wes_purecn_genecn = pd.read_csv(
            wesfolder + "purecn_genecn_all.csv", index_col=0
        )
        wes_loh = pd.read_csv(wesfolder + "purecn_loh_all.csv", index_col=0)
        wes_feature_table = pd.read_csv(
            wesfolder + "globalGenomicFeatures_all.csv", index_col=0
        )
    # subset and rename to PR-indexed matrices
    wessegments_pr = (
        wessegments[wessegments[SAMPLEID].isin(set(renaming_dict.keys()))]
        .replace({SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wes_purecn_segments_pr = (
        wes_purecn_segments[
            wes_purecn_segments[SAMPLEID].isin(set(renaming_dict.keys()))
        ]
        .replace({SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wes_genecn_pr = wesgenecn[wesgenecn.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wes_purecn_genecn_pr = wes_purecn_genecn[
        wes_purecn_genecn.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    wes_loh_pr = wes_loh[wes_loh.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wes_feature_table_pr = wes_feature_table[
        wes_feature_table.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)

    # doing wgs
    print("doing wgs")
    folder = save_dir + "wgs_"
    (
        wgssegments,
        wgsgenecn,
        wgsfailed,
        wgs_purecn_segments,
        wgs_purecn_genecn,
        wgs_loh,
        wgs_feature_table,
    ) = cn.postProcess(
        wgsrefworkspace,
        setEntity=wgssetentity,
        sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        genechangethresh=genechangethresh,
        segmentsthresh=segmentsthresh,
        maxYchrom=maxYchrom,
        **kwargs,
    )

    # with gumbo, no need to mark this selected field
    selected = []

    try:
        track.updateTrackerWGS(
            tracker,
            selected,
            samplesetname,
            wgsfailed,
            datatype=["wgs", "wes"],
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wgsrefworkspace,
            dry_run=True,
        )
    except:
        print("no wgs for this sampleset")

    try:
        track.updateTrackerWGS(
            tracker,
            selected,
            samplesetname,
            list(wesfailed),
            datatype=["wes", "wgs"],
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wesrefworkspace,
            dry_run=True,
        )
    except:
        print("no wes for this sampleset")

    wgssegments_pr = (
        wgssegments[wgssegments[SAMPLEID].isin(set(renaming_dict.keys()))]
        .replace({SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wgs_purecn_segments_pr = (
        wgs_purecn_segments[
            wgs_purecn_segments[SAMPLEID].isin(set(renaming_dict.keys()))
        ]
        .replace({SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wgs_genecn_pr = wgsgenecn[wgsgenecn.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wgs_purecn_genecn_pr = wgs_purecn_genecn[
        wgs_purecn_genecn.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    wgs_loh_pr = wgs_loh[wgs_loh.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    wgs_feature_table_pr = wgs_feature_table[
        wgs_feature_table.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)

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

    print("merging PR-level absolute seg file")
    merged_purecn_segments_pr = wgs_purecn_segments_pr.append(
        wes_purecn_segments_pr
    ).reset_index(drop=True)
    merged_purecn_segments_pr = (
        merged_purecn_segments_pr[
            [
                SAMPLEID,
                "Chromosome",
                "Start",
                "End",
                "Absolute_CN",
                "Minor_allele_absolute_CN",
                "LOH_status",
            ]
        ]
        .sort_values(by=[SAMPLEID, "Chromosome", "Start", "End"])
        .reset_index(drop=True)
    )

    # merging wes and wgs
    # CDS-ID level
    print("saving merged files")
    folder = save_dir
    mergedsegments = wgssegments.append(wessegments).reset_index(drop=True)
    mergedsegments.to_csv(folder + "merged_segments.csv", index=False)
    mergedcn = wgsgenecn.append(wesgenecn)
    mergedcn.to_csv(folder + "merged_genecn.csv")
    merged_purecn_segments = wgs_purecn_segments.append(
        wes_purecn_segments
    ).reset_index(drop=True)
    merged_purecn_segments.to_csv(folder + "merged_absolute_segments.csv", index=False)
    merged_purecn_genecn = wgs_purecn_genecn.append(wes_purecn_genecn)
    merged_purecn_genecn.to_csv(folder + "merged_absolute_genecn.csv")
    merged_loh = wgs_loh.append(wgs_loh)
    merged_loh.to_csv(folder + "merged_loh.csv")
    merged_feature_table = wgs_feature_table.append(wes_feature_table)
    merged_feature_table.to_csv(folder + "merged_feature_table.csv")

    # profile-ID level
    mergedsegments_pr.to_csv(folder + "merged_segments_profile.csv", index=False)
    mergedgenecn_pr = wgs_genecn_pr.append(wes_genecn_pr)
    mergedgenecn_pr.to_csv(folder + "merged_genecn_profile.csv")
    merged_purecn_segments_pr.to_csv(
        folder + "merged_absolute_segments_profile.csv", index=False
    )
    merged_purecn_genecn_pr = wgs_purecn_genecn_pr.append(wes_purecn_genecn_pr)
    merged_purecn_genecn_pr.to_csv(folder + "merged_absolute_genecn_profile.csv")
    merged_loh_pr = wgs_loh_pr.append(wes_loh_pr)
    merged_loh_pr.to_csv(folder + "merged_loh_profile.csv")
    merged_feature_table_pr = wgs_feature_table_pr.append(wes_feature_table_pr)
    merged_feature_table_pr.to_csv(folder + "merged_feature_table_profile.csv")

    # uploading to taiga
    print("uploading to taiga")
    tc.update_dataset(
        changes_description="new "
        + samplesetname
        + " release! (removed misslabellings, see changelog)",
        dataset_permaname=taiga_dataset,
        upload_files=[
            {
                "path": folder + "merged_segments.csv",
                "name": "merged_segments_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_genecn.csv",
                "name": "merged_gene_cn_withReplicates",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_segments_all_profile.csv",
                "name": "merged_segments_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_genecn_all_profile.csv",
                "name": "merged_gene_cn_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            # Pure CN outputs
            {
                "path": folder + "merged_absolute_segments.csv",
                "name": "merged_absolute_segments_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_absolute_genecn.csv",
                "name": "merged_absolute_gene_cn_withReplicates",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_loh.csv",
                "name": "merged_loh_withReplicates",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_feature_table.csv",
                "name": "globalGenomicFeatures_withReplicates",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_absolute_segments_profile.csv",
                "name": "merged_absolute_segments_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_absolute_genecn_profile.csv",
                "name": "merged_absolute_gene_cn_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_loh_profile.csv",
                "name": "merged_loh_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_feature_table_profile.csv",
                "name": "globalGenomicFeatures_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        dataset_description=dataset_description,
    )
    print("done")
    return wessegments, wgssegments


async def mutationPostProcessing(
    wesrefworkspace=WESMUTWORKSPACE,
    wescnworkspace=WESCNWORKSPACE,
    wgsrefworkspace=WGSWORKSPACE,
    vcfdir=VCFDIR,
    vcf_colname=VCFCOLNAME,
    samplesetname=SAMPLESETNAME,
    AllSamplesetName="all",
    doCleanup=False,
    taiga_description=Mutationsreadme,
    taiga_dataset=TAIGA_MUTATION,
    mutation_groups=MUTATION_GROUPS,
    bed_location=GUIDESBED,
    minfreqtocall=MINFREQTOCALL,
    sv_col=SV_COLNAME,
    sv_filename=SV_FILENAME,
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

    # doing wes
    print("doing wes")
    folder = WORKING_DIR + samplesetname + "/wes_"

    wesmutations, wessvs = mutations.postProcess(
        wesrefworkspace,
        AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        doCleanup=True,
        mutCol="CGA_WES_AC",
        **kwargs,
    )

    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))
    wesmutations_pr = wesmutations[
        wesmutations[SAMPLEID].isin(renaming_dict.keys())
    ].replace({SAMPLEID: renaming_dict})

    # doing wgs
    print("doing wgs")
    folder = WORKING_DIR + samplesetname + "/wgs_"

    wgsmutations, wgssvs = mutations.postProcess(
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

    # merge
    print("merging")
    folder = WORKING_DIR + samplesetname + "/merged_"
    mergedmutations = wgsmutations.append(wesmutations).reset_index(drop=True)
    mergedmutations.to_csv(folder + "somatic_mutations.csv", index=False)

    mergedsvs = wgssvs.append(wessvs).reset_index(drop=True)
    mergedsvs.to_csv(folder + "svs.bedpe", index=False)
    mergedsvs_pr = mergedsvs[mergedsvs[SAMPLEID].isin(renaming_dict.keys())].replace(
        {SAMPLEID: renaming_dict}
    )
    mergedsvs_pr.to_csv(folder + "svs_profile.bedpe", index=False)

    merged = wgsmutations_pr.append(wesmutations_pr).reset_index(drop=True)

    # making binary mutation matrices
    print("creating mutation matrices")
    print("changing variant annotations")
    rename = {}
    for k, v in mutation_groups.items():
        for e in v:
            rename[e] = k
    merged["Variant_annotation"] = [
        rename[i] for i in merged["Variant_Classification"].tolist()
    ]

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

    # generate germline binary matrix
    wgswm = dm.WorkspaceManager(wgsrefworkspace)
    weswm = dm.WorkspaceManager(wescnworkspace)
    wgs_samples = wgswm.get_samples()
    wes_samples = weswm.get_samples()
    wgs_vcfs = wgs_samples[vcf_colname]
    wes_vcfs = wes_samples[vcf_colname]
    vcflist = wgs_vcfs[~wgs_vcfs.isna()].tolist() + wes_vcfs[~wes_vcfs.isna()].tolist()

    print("generate germline binary matrix")
    germline_mat = mut.generateGermlineMatrix(
        vcflist,
        vcfdir=vcfdir,
        savedir=WORKING_DIR + samplesetname + "/",
        filename="binary_mutguides.tsv.gz",
        bed_location=bed_location,
    )
    # merging wes and wgs
    print("merging and renaming wes and wgs germline matrices")
    germline_mat_noguides = germline_mat.iloc[:, 4:]

    pr_table = mytracker.read_pr_table()
    # transform from CDSID-level to PR-level
    renaming_dict = dict(list(zip(pr_table.CDSID, pr_table.index)))

    whitelist_cols = [x for x in germline_mat_noguides.columns if x in renaming_dict]
    whitelist_germline_mat = germline_mat_noguides[whitelist_cols]
    mergedmat = whitelist_germline_mat.rename(columns=renaming_dict)

    mergedmat = mergedmat.astype(bool).astype(int)
    sorted_mat = germline_mat.iloc[:, :4].join(mergedmat)
    sorted_mat["end"] = sorted_mat["end"].astype(int)
    print("saving wes and wgs germline matrices")
    sorted_mat.to_csv(folder + "binary_germline.csv", index=False)
    # uploading to taiga
    tc.update_dataset(
        changes_description="new " + samplesetname + " release!",
        dataset_permaname=taiga_dataset,
        upload_files=[
            # for depmap
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                "name": "somaticMutations_boolMatrix_hotspot_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder
                + "somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                "name": "somaticMutations_boolMatrix_otherNonconserving_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                "name": "somaticMutations_boolMatrix_damaging_profile",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            # {
            #     "path": folder + "somatic_mutations_boolmatrix_fordepmap_othercons.csv",
            #     "name": "somaticMutations_boolMatrix_otherConserving_profile",
            #     "format": "NumericMatrixCSV",
            #     "encoding": "utf-8",
            # },
            {
                "path": folder + "somatic_mutations_fordepmap_profile.csv",
                "name": "somaticMutations_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "somatic_mutations.csv",
                "name": "somaticMutations_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_binary_germline.csv",
                "name": "binary_germline_mutation",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "svs.bedpe",
                "name": "structuralVariants_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "svs_profile.bedpe",
                "name": "structuralVariants_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        upload_async=False,
        dataset_description=taiga_description,
    )

