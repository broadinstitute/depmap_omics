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

from .mutations import postprocess_main_steps


async def expressionPostProcessing(
    refworkspace=env_config.RNAWORKSPACE,
    samplesetname=constants.SAMPLESETNAME,
    samplesetname_stranded=constants.SAMPLESETNAME_STRANDED,
    colstoclean=["fastq1", "fastq2", "recalibrated_bam", "recalibrated_bam_index"],
    ensemblserver=constants.ENSEMBL_SERVER_V,
    doCleanup=True,
    samplesetToLoad="all",
    strandedSamplesetToLoad="all_stranded",
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
    generate_count_matrix=True,
    run_stranded=True,
    rnaseqc2_gene_count_col=constants.RNASEQC2_GENE_COUNT_COL,
    rnaseqc2_gene_count_col_stranded=constants.RNASEQC2_GENE_COUNT_COL_STRANDED,
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

    files_stranded = dict()
    if run_stranded:
        files_stranded = await expressions.postProcessStranded(
            refworkspace,
            samplesetname_stranded,
            failed,
            save_output=folder,
            ensemblserver=ensemblserver,
            samplesetToLoad=samplesetname_stranded,
            geneLevelCols=constants.RSEMFILENAME_GENE_STRANDED,
            trancriptLevelCols=constants.RSEMFILENAME_TRANSCRIPTS_STRANDED,
        )

    if not dry_run:
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
            billing_proj=billing_proj,
        )

    pr_table = mytracker.read_pr_table()

    if not dry_run:
        # subset and rename, include all PRs that have associated CDS-ids
        pr_table = mytracker.update_pr_from_seq(["rna"])

    mytracker.close_gumbo_client()

    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))
    h.dictToFile(renaming_dict, folder + "rna_seq2pr_renaming.json")
    pr_files = dict()
    for k, v in files.items():
        pr_files[k + "_profile"] = v[v.index.isin(set(renaming_dict.keys()))].rename(
            index=renaming_dict
        )
    if enrichments is not None:
        enrichments = enrichments[
            enrichments.index.isin(set(renaming_dict.keys()))
        ].rename(index=renaming_dict)
        enrichments.to_csv(folder + "gene_sets_profile.csv")
    expressions.saveFiles(pr_files, folder)
    if run_stranded:
        pr_files_stranded = dict()
        tpm_mat = files_stranded["proteincoding_genes_tpm_stranded"]
        pr_files_stranded["proteincoding_genes_tpm_profile_stranded"] = tpm_mat[tpm_mat.index.isin(set(renaming_dict.keys()))].rename(
            index=renaming_dict
        )
        expressions.saveFiles(pr_files_stranded, folder)

    if generate_count_matrix:
        print("generating rnaseqc gene count matrix")
        rnaseqc_count_dfs = expressions.parse_rnaseqc_counts(refworkspace, samplesetToLoad, rnaseqc2_gene_count_col)
        rnaseqc_count_mat = pd.concat(rnaseqc_count_dfs, axis=1)
        rnaseqc_count_mat = rnaseqc_count_mat.T
        rnaseqc_count_mat.to_csv(folder + "rnaseqc_count_mat.csv")
        rnaseqc_count_mat_pr = rnaseqc_count_mat[rnaseqc_count_mat.index.isin(set(renaming_dict.keys()))].rename(index=renaming_dict)
        rnaseqc_count_mat_pr.to_csv(folder + "rnaseqc_count_mat_pr.csv")
        if run_stranded:
            print("generating rnaseqc gene count matrix for stranded subset")
            rnaseqc_count_dfs = expressions.parse_rnaseqc_counts(refworkspace, strandedSamplesetToLoad, rnaseqc2_gene_count_col_stranded)
            rnaseqc_count_mat = pd.concat(rnaseqc_count_dfs, axis=1)
            rnaseqc_count_mat = rnaseqc_count_mat.T
            rnaseqc_count_mat.to_csv(folder + "stranded_rnaseqc_count_mat.csv")
            rnaseqc_count_mat_pr = rnaseqc_count_mat[rnaseqc_count_mat.index.isin(set(renaming_dict.keys()))].rename(index=renaming_dict)
            rnaseqc_count_mat_pr.to_csv(folder + "stranded_rnaseqc_count_mat_pr.csv")

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
                {
                    "path": folder + "genes_effective_length_profile.csv",
                    "name": "genes_effectiveLength_profile",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "genes_effective_length.csv",
                    "name": "genes_effectiveLength_withReplicates",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "rna_qcs/all_qc.csv",
                    "name": "all_samples_qc",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
            ],
            upload_async=False,
            dataset_description=dataset_description,
        )
        if enrichments is not None:
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
                add_all_existing_files=True,
                upload_async=False,
                dataset_description=dataset_description,
            )
        if generate_count_matrix:
            tc.update_dataset(
                changes_description="adding rnaseqc2 gene counts for new "
                + samplesetname
                + " release!",
                dataset_permaname=taiga_dataset,
                upload_files=[
                    {
                        "path": folder + "rnaseqc_count_mat.csv",
                        "name": "rnaseqc_count_mat_withReplicates",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": folder + "rnaseqc_count_mat_pr.csv",
                        "name": "rnaseqc_count_mat_profile",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                ],
                add_all_existing_files=True,
                upload_async=False,
                dataset_description=dataset_description,
            )
        if run_stranded:
            tc.update_dataset(
                changes_description="adding stranded data sets for "
                + samplesetname
                + " release!",
                dataset_permaname=taiga_dataset,
                upload_files=[
                    {
                        "path": folder + "stranded_rnaseqc_count_mat_pr.csv",
                        "name": "stranded_rnaseqc_count_mat_profile",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": folder + "proteincoding_genes_tpm_profile_stranded_logp1.csv",
                        "name": "stranded_proteinCoding_genes_tpm_logp1_profile",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    
                ],
                add_all_existing_files=True,
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
            wes_arm_cna,
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
        wes_arm_cna = pd.read_csv(wesfolder + "arm_cna_all.csv", index_col=0)
        wes_feature_table = pd.read_csv(
            wesfolder + "globalGenomicFeaturesWithAneuploidy_all.csv", index_col=0
        )

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
        wgs_arm_cna,
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
            wesfailed,
            datatype=["wes", "wgs"],
            dry_run=dryrun,
            billing_proj=billing_proj,
        )
    except:
        print("no wes for this sampleset")

    pr_table = mytracker.update_pr_from_seq(["wgs"])
    pr_table = mytracker.update_pr_from_seq(["wes"])

    mytracker.close_gumbo_client()

    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))

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
    wes_arm_cna_pr = wes_arm_cna[
        wes_arm_cna.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    wes_feature_table_pr = wes_feature_table[
        wes_feature_table.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)

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
    wgs_arm_cna_pr = wgs_arm_cna[
        wgs_arm_cna.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
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
                "SegmentAbsoluteCN",
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
    mergedcn = (wgsgenecn.append(wesgenecn))
    mergedcn.to_csv(folder + "merged_genecn.csv")
    merged_purecn_segments = wgs_purecn_segments.append(
        wes_purecn_segments
    ).reset_index(drop=True)
    merged_purecn_segments.to_csv(folder + "merged_absolute_segments.csv", index=False)
    merged_purecn_genecn = wgs_purecn_genecn.append(wes_purecn_genecn)
    merged_purecn_genecn.to_csv(folder + "merged_absolute_genecn.csv")
    merged_loh = wgs_loh.append(wes_loh)
    merged_loh.to_csv(folder + "merged_loh.csv")
    merged_arm_cna = wes_arm_cna.append(wgs_arm_cna)
    merged_arm_cna.to_csv(folder + "merged_arm_cna.csv")
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
    merged_arm_cna_pr = wes_arm_cna_pr.append(wgs_arm_cna_pr)
    merged_arm_cna_pr.to_csv(folder + "merged_arm_cna_profile.csv")
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
                "path": folder + "merged_arm_cna.csv",
                "name": "armLevelCNA_withReplicates",
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
            {
                "path": folder + "merged_arm_cna_profile.csv",
                "name": "armLevelCNA_profile",
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
    bed_locations: dict = constants.GUIDESBED,
    sv_col: str = constants.SV_COLNAME,
    sv_filename: str = constants.SV_FILENAME,
    mutcol: dict = constants.MUTCOL_DEPMAP,
    standardmafcol: dict = constants.MUTCOL_STANDARDMAF,
    mafcol: str = constants.MAF_COL,
    run_sv: bool = True,
    run_guidemat: bool = True,
    upload_taiga: bool = True,
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

    # TODO: replace with multiprocessing 
    # ./sandbox/dna_eval/combine_mafs.py
    wesmutations, _ = mutations.postProcess(
        wes_wm,
        AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        sv_col=sv_col,
        sv_filename=sv_filename,
        mafcol=mafcol,
        run_sv=False,
        debug=False,
        **kwargs,
    )

    mytracker = track.SampleTracker()
    pr_table = mytracker.read_pr_table()
    renaming_dict = dict(list(zip(pr_table.MainSequencingID, pr_table.index)))
    mytracker.close_gumbo_client()

    wesmutations_pr = wesmutations[
        wesmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace({constants.SAMPLEID: renaming_dict, "Tumor_Sample_Barcode": renaming_dict})

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
        debug=False,
        **kwargs,
    )

    wgsmutations_pr = wgsmutations[
        wgsmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace({constants.SAMPLEID: renaming_dict, "Tumor_Sample_Barcode": renaming_dict})

    # merge
    print("merging WES and WGS")
    folder = constants.WORKING_DIR + samplesetname + "/merged_"
    # if not os.path.exists(constants.WORKING_DIR + samplesetname):
    #     os.mkdir(constants.WORKING_DIR + samplesetname)

    mergedmutations = pd.concat([wgsmutations, wesmutations], axis=0).reset_index(
        drop=True
    )

    mutcol.update(constants.MUTCOL_ADDITIONAL)
    mergedmutations = mergedmutations.rename(columns=mutcol)

    mergedmutations = mutations.addEntrez(mergedmutations, ensembl_col="EnsemblGeneID", entrez_col="EntrezGeneID")

    # https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#somatic-maf-file-generation
    # For all columns, convert "Y" to True/False
    for col in mergedmutations.columns:
        if "Y" in mergedmutations[col].values:
            mergedmutations.loc[:, col] = np.where(
                mergedmutations[col].values == "Y", True, False
            )

    mergedmutations[list(mutcol.values()) + ["EntrezGeneID"]].to_csv(folder + "somatic_mutations.csv", index=False)

    if run_sv:
        if wgssvs is not None:
            wgssvs.to_csv(folder + "svs.csv", index=False)
            wgssvs_pr = wgssvs[
                wgssvs[constants.SAMPLEID].isin(renaming_dict.keys())
            ].replace({constants.SAMPLEID: renaming_dict})
            print("saving somatic svs")
            wgssvs_pr.to_csv(folder + "svs_profile.csv", index=False)

    merged = pd.concat([wgsmutations_pr, wesmutations_pr], axis=0).reset_index(
        drop=True
    )

    # For all columns, convert "Y" to True/False
    for col in merged.columns:
        if "Y" in merged[col].values:
            merged.loc[:, col] = np.where(merged[col].values == "Y", True, False)

    merged = merged.rename(columns=mutcol)
    merged = mutations.addEntrez(merged, ensembl_col="EnsemblGeneID", entrez_col="EntrezGeneID")
    merged.to_csv(folder + "somatic_mutations_all_cols_profile.csv", index=False)
    merged[list(mutcol.values()) + ["EntrezGeneID"]].to_csv(folder + "somatic_mutations_profile.csv", index=False)
    merged[standardmafcol.keys()].to_csv(folder + "somatic_mutations_profile.maf.csv", index=False)

    # making genotyped mutation matrices
    print("creating mutation matrices")
    hotspot_mat, lof_mat = mutations.makeMatrices(merged)
    # add entrez ids to column names
    merged["gene_name"] = [
        i["HugoSymbol"] + " (" + str(i["EntrezGeneID"]).split(".")[0] + ")"
        if i["EntrezGeneID"] != ""
        else i["HugoSymbol"] + " (Unknown)"
        for _, i in merged.iterrows()
    ]
    symbol_to_symbolentrez_dict = dict(zip(merged.HugoSymbol, merged.gene_name))
    hotspot_mat = hotspot_mat.rename(columns=symbol_to_symbolentrez_dict)
    lof_mat = lof_mat.rename(columns=symbol_to_symbolentrez_dict)

    hotspot_mat.to_csv(folder + "somatic_mutations_genotyped_hotspot_profile.csv")
    lof_mat.to_csv(folder + "somatic_mutations_genotyped_damaging_profile.csv")

    # TODO: add pandera type validation

    if run_guidemat:
        # aggregate germline binary matrix
        print("aggregating binary guide mutation matrices")
        print("aggregating wes")
        wes_germline_mats = mutations.aggregateGermlineMatrix(wes_wm, AllSamplesetName, save_output=folder)
        print("aggregating wgs")
        wgs_germline_mats = mutations.aggregateGermlineMatrix(wgs_wm, AllSamplesetName, save_output=folder)

        for lib, _ in bed_locations.items():
            assert lib in wes_germline_mats, "library missing in wes"
            assert lib in wgs_germline_mats, "library missing in wgs"
            # merging wes and wgs
            print("renaming merged wes and wgs germline matrix for library: ", lib)
            germline_mat_merged = pd.concat([wes_germline_mats[lib], wgs_germline_mats[lib].iloc[:, 4:]], axis=1)
            germline_mat_merged_noguides = germline_mat_merged.iloc[:, 4:]

            # transform from CDSID-level to PR-level
            whitelist_cols = [
                x for x in germline_mat_merged_noguides.columns if x in renaming_dict
            ]
            whitelist_germline_mat = germline_mat_merged_noguides[whitelist_cols]
            mergedmat = whitelist_germline_mat.rename(columns=renaming_dict)

            mergedmat = mergedmat.astype(bool).astype(int)
            sorted_mat = germline_mat_merged.iloc[:, :4].join(mergedmat)
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
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations_all_cols_profile.csv",
                    "name": "somaticMutations_profile_all_cols",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "somatic_mutations.csv",
                    "name": "somaticMutations_withReplicates",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ], upload_async=False,
            dataset_description=taiga_description,)
        if run_guidemat:
            tc.update_dataset(
                changes_description="new " + samplesetname + " release!",
                dataset_permaname=taiga_dataset,
                upload_files=[
                    {
                        "path": folder + "binary_germline_avana.csv",
                        "name": "binary_mutation_avana",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": folder + "binary_germline_ky.csv",
                        "name": "binary_mutation_ky",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": folder + "binary_germline_humagne.csv",
                        "name": "binary_mutation_humagne",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                ], 
                add_all_existing_files=True,
                upload_async=False,
                dataset_description=taiga_description,)
        if run_sv:
            tc.update_dataset(
                changes_description="new " + samplesetname + " release!",
                dataset_permaname=taiga_dataset,
                upload_files=[           
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
                add_all_existing_files=True,
                upload_async=False,
                dataset_description=taiga_description,)
