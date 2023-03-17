from depmapomics import constants
from depmapomics import env_config
import os.path
import dalmatian as dm
import pandas as pd
import numpy as np
from taigapy import TaigaClient

from mgenepy.utils import helper as h

from depmap_omics_upload import tracker as track

from depmapomics import expressions
from depmapomics import mutations
from depmapomics import fusions as fusion
from depmapomics import copynumbers as cn


async def expressionPostProcessing(
    refworkspace=env_config.RNAWORKSPACE,
    samplesetname=constants.SAMPLESETNAME,
    colstoclean=["fastq1", "fastq2", "recalibrated_bam", "recalibrated_bam_index"],
    ensemblserver=constants.ENSEMBL_SERVER_V,
    doCleanup=True,
    samplesetToLoad="all",
    taiga_dataset=env_config.TAIGA_EXPRESSION,
    save_output=constants.WORKING_DIR,
    minsimi=constants.RNAMINSIMI,
    dropNonMatching=True,
    dataset_description=constants.RNAseqreadme,
    dry_run=False,
    samplesinset=[],
    rsemfilelocs=None,
    rnaqclocs={},
    starlogs={},
    compute_enrichment=False,
    billing_proj=constants.GCS_PAYER_PROJECT,
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
        ensemblserver (str, optional): ensembl server biomart version . Defaults to constants.ENSEMBL_SERVER_V.
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        geneLevelCols (list, optional): the columns that contain the gene level
        expression data in the workspace. Defaults to constants.RSEMFILENAME_GENE.
        trancriptLevelCols (list, optional): the columns that contain the transcript
        level expression data in the workspacce. Defaults to constants.RSEMFILENAME_TRANSCRIPTS.
        ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
        renamingFunc (function, optional): the function to use to rename the sample columns
        (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
        compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
        dropNonMatching (bool, optional): whether to drop the non matching genes
        between entrez and ensembl. Defaults to False.
        recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
        taiga_dataset (str, optional): the taiga dataset path to use for uploading results. Defaults to env_config.TAIGA_EXPRESSION.
        minsimi (float, optional): the minimum similarity to use for comparison to previous dataset. Defaults to 0.95.
        dataset_description (str, optional): the taiga dataset description to use. Defaults to constants.RNAseqreadme.
        tocompare (dict, optional): the columns to compare. Defaults to {"genes_expected_count": "CCLE_RNAseq_reads", "genes_tpm": "CCLE_expression_full", "proteincoding_genes_tpm": "CCLE_expression"}.
        rsemfilelocs (pd.DataFrame, optional): locations of RSEM output files if refworkspace is not provided (bypass interaction with terra)
        samplesinset (list[str], optional): list of samples in the sampleset if refworkspace is not provided (bypass interaction with terra)
        rnaqclocs (dict(str:list[str]), optional): dict(sample_id:list[QC_filepaths]) of rna qc file locations if refworkspace is not provided (bypass interaction with terra)
    """
    tc = TaigaClient()

    mytracker = track.SampleTracker()

    ccle_refsamples = mytracker.read_seq_table()

    folder = save_output + samplesetname + "/"

    if dry_run:
        folder = save_output + "dryrun/"

    h.createFoldersFor(folder)
    files, failed, _, renaming, lowqual, enrichments = await expressions.postProcess(
        refworkspace,
        samplesetname,
        save_output=folder,
        doCleanup=doCleanup,
        colstoclean=colstoclean,
        ensemblserver=ensemblserver,
        samplesetToLoad=samplesetToLoad,
        geneLevelCols=constants.RSEMFILENAME_GENE,
        trancriptLevelCols=constants.RSEMFILENAME_TRANSCRIPTS,
        ssGSEAcol="genes_tpm",
        dropNonMatching=dropNonMatching,
        dry_run=dry_run,
        samplesinset=samplesinset,
        rsemfilelocs=rsemfilelocs,
        rnaqclocs=rnaqclocs,
        compute_enrichment=compute_enrichment,
        **kwargs,
    )

    print("updating the tracker")

    track.updateTrackerRNA(
        failed,
        lowqual[lowqual.sum(1) > 3].index.tolist(),
        ccle_refsamples,
        samplesetname,
        refworkspace,
        samplesinset=samplesinset,
        starlogs=starlogs,
        dry_run=dry_run,
        newgs=None,
        billing_proj=billing_proj,
    )

    # subset and rename, include all PRs that have associated CDS-ids
    pr_table = mytracker.update_pr_from_seq(["rna"])

    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))
    h.dictToFile(renaming_dict, folder + "rna_seq2pr_renaming.json")
    pr_files = dict()
    for k, v in files.items():
        pr_files[k + "_profile"] = v[v.index.isin(set(renaming_dict.keys()))].rename(
            index=renaming_dict
        )
    if enrichments != None:
        enrichments = enrichments[
            enrichments.index.isin(set(renaming_dict.keys()))
        ].rename(index=renaming_dict)
        enrichments.to_csv(folder + "gene_sets_profile.csv")
    expressions.saveFiles(pr_files, folder)
    mytracker.close_gumbo_client()

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
                    "path": folder + "proteincoding_genes_tpm_logp1.csv",
                    "name": "proteinCoding_genes_tpm_logp1_withReplicates",
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
            ],
            upload_async=False,
            dataset_description=dataset_description,
        )
        if enrichments != None:
            tc.update_dataset(
                changes_description="adding enrichments for new "
                + samplesetname
                + " release!",
                dataset_permaname=taiga_dataset,
                upload_files=[
                    {
                        "path": folder + "gene_sets_all.csv",
                        "name": "gene_set_enrichment_withReplicates",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": folder + "gene_sets_profile.csv",
                        "name": "gene_set_enrichment_profile",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                ],
                upload_async=False,
                dataset_description=dataset_description,
            )
        print("done")


async def fusionPostProcessing(
    refworkspace=env_config.RNAWORKSPACE,
    sampleset=constants.SAMPLESETNAME,
    fusionSamplecol=constants.SAMPLEID,
    taiga_dataset=env_config.TAIGA_FUSION,
    dataset_description=constants.FUSIONreadme,
    folder=constants.WORKING_DIR,
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
        fusionSamplecol ([type], optional): [description]. Defaults to constants.SAMPLEID.
        taiga_dataset (str, optional): the taiga dataset path to use for uploading results. Defaults to env_config.TAIGA_EXPRESSION.
        dataset_description (str, optional): the taiga dataset description to use. Defaults to constants.RNAseqreadme.

    Returns:
        (pd.df): fusion dataframe
        (pd.df): filtered fusion dataframe
    """
    tc = TaigaClient()

    mytracker = track.SampleTracker()
    ccle_refsamples = mytracker.read_seq_table()

    previousQCfail = ccle_refsamples[ccle_refsamples.low_quality == 1].index.tolist()

    # TODO: include in rna_sample_renaming.json instead
    # lower priority versions of these lines were used

    folder = folder + sampleset + "/"

    fusions, fusions_filtered = fusion.postProcess(
        refworkspace,
        todrop=previousQCfail,
        save_output=folder,
        **kwargs,
    )

    # subset, rename from seqid to prid, and save pr-indexed matrices
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))
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

    mytracker.close_gumbo_client()

    # taiga
    print("uploading to taiga")
    tc.update_dataset(
        dataset_permaname=taiga_dataset,
        changes_description="new " + sampleset + " release!",
        upload_files=[
            {
                "path": folder + "/fusions_all.csv",
                "name": "fusions_unfiltered_withReplicates",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "/filteredfusions_latest_profile.csv",
                "name": "fusions_filtered_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "/fusions_all_profile.csv",
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
    wesrefworkspace=env_config.WESCNWORKSPACE,
    wgsrefworkspace=env_config.WGSWORKSPACE,
    wessetentity=constants.WESSETENTITY,
    wgssetentity=constants.WGSSETENTITY,
    samplesetname=constants.SAMPLESETNAME,
    purecnsampleset=constants.PURECN_SAMPLESET,
    AllSamplesetName="all",
    taiga_dataset=env_config.TAIGA_CN,
    dataset_description=constants.CNreadme,
    subsetsegs=[
        constants.SAMPLEID,
        "Chromosome",
        "Start",
        "End",
        "Segment_Mean",
        "Num_Probes",
        "Status",
        "Source",
    ],
    bamqc=constants.BAMQC,
    procqc=constants.PROCQC,
    save_dir=constants.WORKING_DIR,
    wesfolder="",
    segmentsthresh=constants.SEGMENTSTHR,
    maxYchrom=constants.MAXYCHROM,
    billing_proj=constants.GCS_PAYER_PROJECT,
    dryrun=False,
    **kwargs,
):
    """the full CCLE Copy Number post processing pipeline (used only by CCLE)
    see postprocessing() to reproduce most of our analysis and find out about additional parameters
    Args:
        wesrefworkspace (str): wes terra workspace where the ref data is stored
        wgsrefworkspace (str): wgs terra workspace where the ref data is stored
        samplesetname (str): name of the current release
        AllSamplesetName (str, optional): name of the sample set to get the data from (should contain everything). Defaults to 'all'.
        taiga_dataset (str, optional): where to save the output to on taiga. Defaults to env_config.TAIGA_CN.
        dataset_description (str, optional): A long string that will be pushed to taiga to explain the CN dataset. Defaults to constants.CNreadme.
        subsetsegs (list[str], optional): what columns to keep for the segments. Defaults to [constants.SAMPLEID, 'Chromosome', 'Start', 'End', 'Segment_Mean', 'Num_Probes', 'Status', 'Source'].
        bamqc ([type], optional): @see updateTracker. Defaults to constants.BAMQC.
        procqc ([type], optional): @see updateTracker. Defaults to constants.PROCQC.
        source_rename ([type], optional): @see managing duplicates. Defaults to constants.SOURCE_RENAME.
    """
    tc = TaigaClient()

    mytracker = track.SampleTracker()
    tracker = mytracker.read_seq_table()

    assert len(tracker) != 0, "broken source for sample tracker"
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))

    mytracker.close_gumbo_client()

    save_dir = save_dir + samplesetname + "/"
    # doing wes
    folder = save_dir + "wes_"
    if wesfolder == "":
        print("doing wes")
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
            save_output=folder,
            segmentsthresh=segmentsthresh,
            maxYchrom=maxYchrom,
            purecnsampleset=purecnsampleset,
            **kwargs,
        )

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
        wessegments[wessegments[constants.SAMPLEID].isin(set(renaming_dict.keys()))]
        .replace({constants.SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wes_purecn_segments_pr = (
        wes_purecn_segments[
            wes_purecn_segments[constants.SAMPLEID].isin(set(renaming_dict.keys()))
        ]
        .replace({constants.SAMPLEID: renaming_dict})
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
        segmentsthresh=segmentsthresh,
        maxYchrom=maxYchrom,
        purecnsampleset=purecnsampleset,
        **kwargs,
    )

    try:
        track.updateTrackerWGS(
            tracker,
            samplesetname,
            wgsfailed,
            datatype=["wgs", "wes"],
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wgsrefworkspace,
            dry_run=dryrun,
            billing_proj=billing_proj,
        )
    except:
        print("no wgs for this sampleset")

    try:
        track.updateTrackerWGS(
            tracker,
            samplesetname,
            list(wesfailed),
            datatype=["wes", "wgs"],
            bamqc=bamqc,
            procqc=procqc,
            refworkspace=wesrefworkspace,
            dry_run=dryrun,
            billing_proj=billing_proj,
        )
    except:
        print("no wes for this sampleset")

    wgssegments_pr = (
        wgssegments[wgssegments[constants.SAMPLEID].isin(set(renaming_dict.keys()))]
        .replace({constants.SAMPLEID: renaming_dict})
        .reset_index(drop=True)
    )
    wgs_purecn_segments_pr = (
        wgs_purecn_segments[
            wgs_purecn_segments[constants.SAMPLEID].isin(set(renaming_dict.keys()))
        ]
        .replace({constants.SAMPLEID: renaming_dict})
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
                constants.SAMPLEID,
                "Chromosome",
                "Start",
                "End",
                "SegmentMean",
                "NumProbes",
                "Status",
            ]
        ]
        .sort_values(by=[constants.SAMPLEID, "Chromosome", "Start", "End"])
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
                constants.SAMPLEID,
                "Chromosome",
                "Start",
                "End",
                "MajorAlleleAbsoluteCN",
                "MinorAlleleAbsoluteCN",
                "LoHStatus",
            ]
        ]
        .sort_values(by=[constants.SAMPLEID, "Chromosome", "Start", "End"])
        .reset_index(drop=True)
    )

    # merging wes and wgs
    # CDS-ID level
    print("saving merged files")
    folder = save_dir
    mergedsegments = wgssegments.append(wessegments).reset_index(drop=True)
    mergedsegments.to_csv(folder + "merged_segments.csv", index=False)
    mergedcn = (wgsgenecn.append(wesgenecn)).apply(lambda x: np.log2(1 + x))
    mergedcn.to_csv(folder + "merged_genecn.csv")
    merged_purecn_segments = wgs_purecn_segments.append(
        wes_purecn_segments
    ).reset_index(drop=True)
    merged_purecn_segments.to_csv(folder + "merged_absolute_segments.csv", index=False)
    merged_purecn_genecn = wgs_purecn_genecn.append(wes_purecn_genecn)
    merged_purecn_genecn.to_csv(folder + "merged_absolute_genecn.csv")
    merged_loh = wgs_loh.append(wes_loh)
    merged_loh.to_csv(folder + "merged_loh.csv")
    merged_feature_table = wgs_feature_table.append(wes_feature_table)
    merged_feature_table.to_csv(folder + "merged_feature_table.csv")

    # profile-ID level
    mergedsegments_pr.to_csv(folder + "merged_segments_profile.csv", index=False)
    mergedgenecn_pr = wgs_genecn_pr.append(wes_genecn_pr).apply(
        lambda x: np.log2(1 + x)
    )
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
                "path": folder + "merged_segments_profile.csv",
                "name": "merged_segments_profile",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "merged_genecn_profile.csv",
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
        upload_async=False,
    )
    print("done")
    return wessegments, wgssegments


async def mutationPostProcessing(
    wesrefworkspace: str = env_config.WESCNWORKSPACE,
    wgsrefworkspace: str = env_config.WGSWORKSPACE,
    vcfdir: str = constants.VCFDIR,
    vcf_colname: str = constants.VCFCOLNAME,
    samplesetname: str = constants.SAMPLESETNAME,
    AllSamplesetName: str = "all",
    taiga_description: str = constants.Mutationsreadme,
    taiga_dataset: str = env_config.TAIGA_MUTATION,
    bed_locations: dict[str, str] = constants.GUIDESBED,
    sv_col: str = constants.SV_COLNAME,
    sv_filename: str = constants.SV_FILENAME,
    mutcol: dict[str, str] = constants.MUTCOL_DEPMAP,
    mafcol: str = constants.MAF_COL,
    doCleanup: bool = False,
    run_sv: bool = False,
    run_guidemat: bool = False,
    upload_taiga: bool = False,
    **kwargs,
):
    """The full CCLE mutations post processing pipeline (used only by CCLE)

    see postprocess() to reproduce our analysis

    Args:
        wesrefworkspace (str, optional): the reference workspace for WES. Defaults to env_config.WESCNWORKSPACE.
        wgsrefworkspace (str, optional): the reference workspace for WGS. Defaults to env_config.WGSWORKSPACE.
        samplesetname (str, optional): the sample set name to use (for the release). Defaults to constants.SAMPLESETNAME.
        AllSamplesetName (str, optional): the sample set to use for all samples. Defaults to 'all'.
        doCleanup (bool, optional): whether to clean up the workspace. Defaults to False.
        upload_taiga (bool, optional): whether to upload to taiga. Defaults to False.
    """
    tc = TaigaClient()

    wes_wm = dm.WorkspaceManager(wesrefworkspace)
    wgs_wm = dm.WorkspaceManager(wgsrefworkspace)

    # doing wes
    print("DOING WES")
    folder = constants.WORKING_DIR + samplesetname + "/wes_"

    wesmutations, wessvs = mutations.postProcess(
        wes_wm,
        AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        sv_col=sv_col,
        sv_filename=sv_filename,
        mafcol=mafcol,
        run_sv=run_sv,
        **kwargs,
    )

    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))
    mytracker.close_gumbo_client()

    wesmutations_pr = wesmutations[
        wesmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace({constants.SAMPLEID: renaming_dict})

    # doing wgs
    print("DOING WGS")
    folder = constants.WORKING_DIR + samplesetname + "/wgs_"

    wgsmutations, wgssvs = mutations.postProcess(
        wgs_wm,
        sampleset="all",  # AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        sv_col=sv_col,
        sv_filename=sv_filename,
        mafcol=mafcol,
        run_sv=run_sv,
        **kwargs,
    )

    wgsmutations_pr = wgsmutations[
        wgsmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace({constants.SAMPLEID: renaming_dict})

    # merge
    print("merging WES and WGS")
    folder = constants.WORKING_DIR + samplesetname + "/merged_"
    mergedmutations = pd.concat([wgsmutations, wesmutations], axis=0).reset_index(
        drop=True
    )

    # some hgnc symbols in the maf are outdated, we are renaming them here and then dropping ones that aren't in biomart
    print("replacing outdated hugo symbols and dropping ones that aren't in biomart")
    hugo_mapping = pd.read_csv(constants.HGNC_MAPPING, sep="\t")
    hugo_mapping = {
        b: a for a, b in hugo_mapping[~hugo_mapping["Previous symbol"].isna()].values
    }

    mybiomart = h.generateGeneNames()
    mybiomart = mybiomart.drop_duplicates("hgnc_symbol", keep="first")

    genes_in_maf = set(mergedmutations.hugo_symbol)
    genes_not_in_biomart = genes_in_maf - set(mybiomart.hgnc_symbol)
    maf_gene_renaming = dict()
    maf_genes_to_drop = []
    for gene in genes_not_in_biomart:
        # if the hugo symbol in maf is outdated, and the new name is in biomart,
        # we will rename it to the new name in the maf
        if gene in hugo_mapping and hugo_mapping[gene] in set(mybiomart.hgnc_symbol):
            maf_gene_renaming[gene] = hugo_mapping[gene]
        # if the hugo symbol can't be found in biomart with or without hugo_mapping,
        # we will drop that gene from the maf
        else:
            maf_genes_to_drop.append(gene)
    mergedmutations = mergedmutations[
        ~mergedmutations.hugo_symbol.isin(maf_genes_to_drop)
    ]
    mergedmutations = mergedmutations.replace({"hugo_symbol": maf_gene_renaming})

    # add entrez id column
    symbol_to_entrez_dict = dict(zip(mybiomart.hgnc_symbol, mybiomart.entrezgene_id))
    mergedmutations["EntrezGeneID"] = mergedmutations["hugo_symbol"].map(
        symbol_to_entrez_dict
    )
    mergedmutations["EntrezGeneID"] = mergedmutations["EntrezGeneID"].fillna("Unknown")
    mergedmutations = mergedmutations.drop(columns=["achilles_top_genes"])
    mergedmutations = mergedmutations.rename(columns=mutcol)

    # https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#somatic-maf-file-generation
    # For all columns, convert "Y" to True/False
    for col in mergedmutations.columns:
        if "Y" in mergedmutations[col].values:
            mergedmutations.loc[:, col] = np.where(
                mergedmutations[col].values == "Y", True, False
            )

    mergedmutations.to_csv(folder + "somatic_mutations.csv", index=False)

    if run_sv:
        if wgssvs is not None:
            mergedsvs = wgssvs.append(wessvs).reset_index(drop=True)
            mergedsvs.to_csv(folder + "svs.csv", index=False)
            mergedsvs_pr = mergedsvs[
                mergedsvs[constants.SAMPLEID].isin(renaming_dict.keys())
            ].replace({constants.SAMPLEID: renaming_dict})
            print("saving somatic svs")
            mergedsvs_pr.to_csv(folder + "svs_profile.csv", index=False)

    merged = pd.concat([wgsmutations_pr, wesmutations_pr], axis=0).reset_index(
        drop=True
    )
    merged["EntrezGeneID"] = merged["hugo_symbol"].map(symbol_to_entrez_dict)
    merged["EntrezGeneID"] = merged["EntrezGeneID"].fillna("Unknown")
    merged = merged.drop(columns=["achilles_top_genes"])
    merged = merged.rename(columns=mutcol)
    merged.to_csv(folder + "somatic_mutations_profile.csv", index=False)

    # making genotyped mutation matrices
    print("creating mutation matrices")
    hotspot_mat, lof_mat, driver_mat = mutations.makeMatrices(merged)
    # add entrez ids to column names
    mybiomart["gene_name"] = [
        i["hgnc_symbol"] + " (" + str(i["entrezgene_id"]).split(".")[0] + ")"
        if not pd.isna(i["entrezgene_id"])
        else i["hgnc_symbol"] + " (Unknown)"
        for _, i in mybiomart.iterrows()
    ]
    symbol_to_symbolentrez_dict = dict(zip(mybiomart.hgnc_symbol, mybiomart.gene_name))
    hotspot_mat = hotspot_mat.rename(columns=symbol_to_symbolentrez_dict)
    lof_mat = lof_mat.rename(columns=symbol_to_symbolentrez_dict)
    driver_mat = driver_mat.rename(columns=symbol_to_symbolentrez_dict)

    hotspot_mat.to_csv(folder + "somatic_mutations_genotyped_hotspot_profile.csv")
    lof_mat.to_csv(folder + "somatic_mutations_genotyped_damaging_profile.csv")
    driver_mat.to_csv(folder + "somatic_mutations_genotyped_driver_profile.csv")

    merged.rename(
        columns={
            "HugoSymbol": "Hugo_Symbol",
            "Chrom": "Chromosome",
            "Pos": "Start_Position",
            "VariantType": "Variant_Type",
            "Ref": "Reference_Allele",
            "Alt": "Alternate_Allele",
            "DepMap_ID": "Tumor_Sample_Barcode",
            "VariantInfo": "Variant_Classification",
            "ProteinChange": "Protein_Change",
        },
        inplace=True,
    )

    merged.loc[:, "Chromosome"] = merged.loc[:, "Chromosome"].str.replace("chr", "")

    def assign_end_pos(
        *,
        Start_Position: int,
        Variant_Type: str,
        Reference_Allele: str,
        Alternate_Allele: str,
        **kwargs,
    ) -> int:
        """Assign End_Position to different Variant_Type

        NOTE: MAF is 1-based coordinate, https://www.biostars.org/p/84686/
        If SNP/DNP, end_pos = start_pos + len(ref)
        If INS, end_pos = start_pos + 1
        If DEL, end_pos = start_pos + len(alt) - 1

        Ignore multiple alleles now

        Parameter
        ----------
        Start_Position: int,
        Variant_Type: str,
        Reference_Allele: str,
        Alternate_Allele: str,

        Return
        ---------
        End_Position: int
        """
        end_pos = Start_Position
        if Variant_Type in ["SNP", "DNP", "TNP"]:
            end_pos = Start_Position + len(Reference_Allele)
        if Variant_Type == "INS":
            end_pos = Start_Position + 1
        if Variant_Type == "DEL":
            end_pos = Start_Position + len(Alternate_Allele) - 1
        # TODO add SV types
        return end_pos

    merged.loc[:, "End_Position"] = merged.apply(
        lambda row: assign_end_pos(**row), axis=1
    )

    merged.loc[:, "NCBI_Build"] = "GRCh38"  # or 38?
    merged.loc[:, "Strand"] = "+"  # TODO: need to check vcf2maf later
    merged.loc[:, "Tumor_Seq_Allele1"] = merged.loc[
        :, "Reference_Allele"
    ]  # TODO: need to check vcf2maf later
    merged.loc[:, "Tumor_Seq_Allele2"] = merged.loc[:, "Alternate_Allele"]

    merged = merged.loc[
        :,
        [
            "Hugo_Symbol",
            "NCBI_Build",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Variant_Type",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
            "Variant_Classification",
            "Protein_Change",
        ],
    ]

    # DepMap 20Q1 mutation issues.. https://github.com/PoisonAlien/maftools/issues/644
    # Mapping of variant classification https://github.com/PoisonAlien/maftools/blob/41ddaabbd824f24a99f66439366be98775edebb2/R/icgc_to_maf.R#L52
    merged.loc[:, "Variant_Classification"] = merged.loc[
        :, "Variant_Classification"
    ].map(
        {
            "MISSENSE": "Missense_Mutation",
            "SILENT": "Silent",
            "IN_FRAME_INS": "In_Frame_Ins",
            "IN_FRAME_DEL": "In_Frame_Del",
            "SPLICE_SITE": "Splice_Site",
            "NONSENSE": "Nonsense_Mutation",
            "FRAME_SHIFT_DEL": "Frame_Shift_Del",
            "FRAME_SHIFT_INS": "Frame_Shift_Ins",
            "NONSTOP": "Nonstop_Mutation",
            "START_CODON_SNP": "Silent",
            "START_CODON_INS": "Silent",
        }
    )

    # sort by chrom, start, and end columns for IGV import
    merged["Chromosome"] = merged["Chromosome"].replace({"X": 23, "Y": 24, "M": 25})
    merged = merged.sort_values(by=["Chromosome", "Start_Position", "End_Position"])
    merged["Chromosome"] = merged["Chromosome"].replace({23: "X", 24: "Y", 25: "M"})

    # TODO: add pandera type validation

    merged.to_csv(folder + "somatic_mutations_profile.maf.txt", index=False, sep="\t")

    if run_guidemat:
        # generate germline binary matrix
        wgs_samples = dm.WorkspaceManager(wgsrefworkspace).get_samples()
        wes_samples = dm.WorkspaceManager(wesrefworkspace).get_samples()
        wgs_vcfs = wgs_samples[vcf_colname]
        wes_vcfs = wes_samples[vcf_colname]
        vcflist = (
            wgs_vcfs[~wgs_vcfs.isna()].tolist() + wes_vcfs[~wes_vcfs.isna()].tolist()
        )
        vcflist = [v for v in vcflist if v.startswith("gs://")]

        print("generating germline binary matrix")
        germline_mats = mutations.generateGermlineMatrix(
            vcflist,
            vcfdir=vcfdir,
            savedir=constants.WORKING_DIR + samplesetname + "/",
            filename="binary_mutguides.tsv.gz",
            bed_locations=bed_locations,
        )
        for lib, mat in germline_mats.items():
            # merging wes and wgs
            print("renaming merged wes and wgs germline matrix for library: ", lib)
            germline_mat_noguides = mat.iloc[:, 4:]

            # transform from CDSID-level to PR-level
            whitelist_cols = [
                x for x in germline_mat_noguides.columns if x in renaming_dict
            ]
            whitelist_germline_mat = germline_mat_noguides[whitelist_cols]
            mergedmat = whitelist_germline_mat.rename(columns=renaming_dict)

            mergedmat = mergedmat.astype(bool).astype(int)
            sorted_mat = mat.iloc[:, :4].join(mergedmat)
            sorted_mat["end"] = sorted_mat["end"].astype(int)
            print("saving merged binary matrix for library: ", lib)
            sorted_mat.to_csv(
                folder + "binary_germline" + "_" + lib + ".csv", index=False
            )
    # uploading to taiga
    if upload_taiga:
        tc.update_dataset(
            changes_description="new " + samplesetname + " release!",
            dataset_permaname=taiga_dataset,
            upload_files=[
                {
                    "path": folder + "somatic_mutations_genotyped_driver_profile.csv",
                    "name": "somaticMutations_genotypedMatrix_driver_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations_genotyped_hotspot_profile.csv",
                    "name": "somaticMutations_genotypedMatrix_hotspot_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations_genotyped_damaging_profile.csv",
                    "name": "somaticMutations_genotypedMatrix_damaging_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations_profile.csv",
                    "name": "somaticMutations_profile",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations_profile.maf.csv",
                    "name": "somaticMutations_profile_maf",
                    "format": "MAF",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations.csv",
                    "name": "somaticMutations_withReplicates",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "merged_binary_germline_avana.csv",
                    "name": "binary_mutation_avana",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "merged_binary_germline_ky.csv",
                    "name": "binary_mutation_ky",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "merged_binary_germline_humagne.csv",
                    "name": "binary_mutation_humagne",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "svs.csv",
                    "name": "structuralVariants_withReplicates",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "svs_profile.csv",
                    "name": "structuralVariants_profile",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            upload_async=False,
            dataset_description=taiga_description,
        )
