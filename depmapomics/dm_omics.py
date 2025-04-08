import os.path

import dalmatian as dm
import numpy as np
import pandas as pd
from depmap_omics_upload import tracker as track
from taigapy import TaigaClient, create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat

from depmapomics import constants, env_config, expressions, mutations
from depmapomics import copynumbers as cn
from depmapomics import fusions as fusion
from mgenepy.utils import helper as h

from .mutations import postprocess_main_steps


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
    save_dir=constants.WORKING_DIR,
    wesfolder="",
    segmentsthresh=constants.SEGMENTSTHR,
    maxYchrom=constants.MAXYCHROM,
    hgnc_mapping_taiga=constants.HGNC_MAPPING_TABLE_TAIGAID,
    hgnc_mapping_table_name=constants.HGNC_MAPPING_TABLE_NAME,
    hgnc_mapping_table_version=constants.HGNC_MAPPING_TABLE_VERSION,
    dryrun=False,
    masked_gene_list=constants.MASKED_GENE_LIST,
    omics_id_mapping_table_taigaid=constants.OMICS_ID_MAPPING_TABLE_TAIGAID,
    omics_id_mapping_table_name=constants.OMICS_ID_MAPPING_TABLE_NAME,
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
    client = create_taiga_client_v3()

    with open(masked_gene_list, "r") as f:
        genes_to_mask = f.read().splitlines()

    # read cds->pr mapping table and construct renaming dictionary
    # always read latest version
    print("reading omics ID mapping table from taiga")
    omics_id_mapping_table = client.get(
        name=omics_id_mapping_table_taigaid, file=omics_id_mapping_table_name
    )
    renaming_dict = dict(
        list(
            zip(
                omics_id_mapping_table["sequencing_id"],
                omics_id_mapping_table["profile_id"],
            )
        )
    )

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
            wes_ms_df,
        ) = cn.postProcess(
            wesrefworkspace,
            run_gatk_relative=True,
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
        wes_ms_df = pd.read_csv(wesfolder + "ms_repeats_all.csv")

    print("masking wes")
    cols_to_drop = [col for col in genes_to_mask if col in wesgenecn.columns]
    print("dropping " + str(len(cols_to_drop)) + " genes from WES relative CN")
    wesgenecn = wesgenecn.drop(columns=cols_to_drop)
    cols_to_drop = [col for col in genes_to_mask if col in wes_purecn_genecn.columns]
    wes_purecn_genecn = wes_purecn_genecn.drop(columns=cols_to_drop)
    cols_to_drop = [col for col in genes_to_mask if col in wes_loh.columns]
    wes_loh = wes_loh.drop(columns=cols_to_drop)

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
        wgs_ms_df,
    ) = cn.postProcess(
        wgsrefworkspace,
        run_gatk_relative=False,
        setEntity=wgssetentity,
        sampleset=AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        segmentsthresh=segmentsthresh,
        maxYchrom=maxYchrom,
        purecnsampleset=purecnsampleset,
        **kwargs,
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

    print("map hugo symbols and entrez ids to gene CN columns")
    # pull the gene id mapping table from taiga dataset maintained by the portal
    hgnc_table = cn.make_hgnc_table(
        taiga_id=hgnc_mapping_taiga,
        dataset_version=hgnc_mapping_table_version,
        dataset_file=hgnc_mapping_table_name,
    ).drop(columns="par")
    ensg2hugo_entrez_dict = dict(
        zip(hgnc_table["ensembl_gene_id"], hgnc_table["hugo_entrez"])
    )

    # merging wes and wgs
    # CDS-ID level
    print("saving merged files")
    folder = save_dir
    mergedsegments = wgssegments.append(wessegments).reset_index(drop=True)
    mergedsegments.to_csv(folder + "merged_segments.csv", index=False)
    mergedcn = wgsgenecn.append(wesgenecn)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(mergedcn.columns) & set(hgnc_table["ensembl_gene_id"])
    mergedcn = mergedcn[cols_in_portal_table].rename(columns=ensg2hugo_entrez_dict)
    mergedcn.to_csv(folder + "merged_genecn.csv")
    merged_purecn_segments = wgs_purecn_segments.append(
        wes_purecn_segments
    ).reset_index(drop=True)
    merged_purecn_segments.to_csv(folder + "merged_absolute_segments.csv", index=False)
    merged_purecn_genecn = wgs_purecn_genecn.append(wes_purecn_genecn)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(merged_purecn_genecn.columns) & set(
        hgnc_table["ensembl_gene_id"]
    )
    merged_purecn_genecn = merged_purecn_genecn[cols_in_portal_table].rename(
        columns=ensg2hugo_entrez_dict
    )
    merged_purecn_genecn.to_csv(folder + "merged_absolute_genecn.csv")
    merged_loh = wgs_loh.append(wes_loh)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(merged_loh.columns) & set(hgnc_table["ensembl_gene_id"])
    merged_loh = merged_loh[cols_in_portal_table].rename(columns=ensg2hugo_entrez_dict)
    merged_loh.to_csv(folder + "merged_loh.csv")
    merged_arm_cna = wes_arm_cna.append(wgs_arm_cna)
    merged_arm_cna.to_csv(folder + "merged_arm_cna.csv")
    merged_feature_table = wgs_feature_table.append(wes_feature_table)
    merged_feature_table.to_csv(folder + "merged_feature_table.csv")

    # profile-ID level
    mergedsegments_pr.to_csv(folder + "merged_segments_profile.csv", index=False)
    mergedgenecn_pr = wgs_genecn_pr.append(wes_genecn_pr)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(mergedgenecn_pr.columns) & set(
        hgnc_table["ensembl_gene_id"]
    )
    mergedgenecn_pr = mergedgenecn_pr[cols_in_portal_table].rename(
        columns=ensg2hugo_entrez_dict
    )
    mergedgenecn_pr.to_csv(folder + "merged_genecn_profile.csv")
    merged_purecn_segments_pr.to_csv(
        folder + "merged_absolute_segments_profile.csv", index=False
    )
    merged_purecn_genecn_pr = wgs_purecn_genecn_pr.append(wes_purecn_genecn_pr)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(merged_purecn_genecn_pr.columns) & set(
        hgnc_table["ensembl_gene_id"]
    )
    merged_purecn_genecn_pr = merged_purecn_genecn_pr[cols_in_portal_table].rename(
        columns=ensg2hugo_entrez_dict
    )
    merged_purecn_genecn_pr.to_csv(folder + "merged_absolute_genecn_profile.csv")
    merged_loh_pr = wgs_loh_pr.append(wes_loh_pr)
    # rename ensg -> hugo (entrez)
    cols_in_portal_table = set(merged_loh_pr.columns) & set(
        hgnc_table["ensembl_gene_id"]
    )
    merged_loh_pr = merged_loh_pr[cols_in_portal_table].rename(
        columns=ensg2hugo_entrez_dict
    )
    merged_loh_pr.to_csv(folder + "merged_loh_profile.csv")
    merged_arm_cna_pr = wes_arm_cna_pr.append(wgs_arm_cna_pr)
    merged_arm_cna_pr.to_csv(folder + "merged_arm_cna_profile.csv")
    merged_feature_table_pr = wgs_feature_table_pr.append(wes_feature_table_pr)
    merged_feature_table_pr.to_csv(folder + "merged_feature_table_profile.csv")

    # merging microsatellite repeats
    pd.testing.assert_frame_equal(wes_ms_df.loc[:, :5], wgs_ms_df.loc[:, :5])
    ms_mat_merged = pd.concat([wes_ms_df, wgs_ms_df.iloc[:, 5:]], axis=1)
    ms_mat_merged_no_coords = ms_mat_merged.iloc[:, 5:]

    # transform from CDSID-level to PR-level
    whitelist_cols = [
        x for x in ms_mat_merged_no_coords.columns if x in set(renaming_dict.keys())
    ]
    whitelist_ms_mat = ms_mat_merged_no_coords[whitelist_cols]
    mergedmat = whitelist_ms_mat.rename(columns=renaming_dict)

    ms_mat = ms_mat_merged.iloc[:, :5].join(mergedmat)
    print("saving microsatellite matrix")
    ms_mat.to_csv(folder + "ms_repeat_profile.csv", index=False)

    # uploading to taiga
    print("uploading to taiga")

    client.update_dataset(
        reason="new "
        + samplesetname
        + " release! (removed misslabellings, see changelog)",
        permaname=taiga_dataset,
        additions=[
            UploadedFile(
                local_path=folder + "merged_segments.csv",
                name="merged_segments_withReplicates",
                format=LocalFormat.CSV_TABLE,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_genecn.csv",
                name="merged_gene_cn_withReplicates",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_segments_profile.csv",
                name="merged_segments_profile",
                format=LocalFormat.CSV_TABLE,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_genecn_profile.csv",
                name="merged_gene_cn_profile_for_achilles",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_absolute_segments.csv",
                name="merged_absolute_segments_withReplicates",
                format=LocalFormat.CSV_TABLE,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_absolute_genecn.csv",
                name="merged_absolute_gene_cn_withReplicates",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_loh.csv",
                name="merged_loh_withReplicates",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_feature_table.csv",
                name="globalGenomicFeatures_withReplicates",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_arm_cna.csv",
                name="armLevelCNA_withReplicates",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_absolute_segments_profile.csv",
                name="merged_absolute_segments_profile",
                format=LocalFormat.CSV_TABLE,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_absolute_genecn_profile.csv",
                name="merged_absolute_gene_cn_profile",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_loh_profile.csv",
                name="merged_loh_profile",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_feature_table_profile.csv",
                name="globalGenomicFeatures_profile",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "merged_arm_cna_profile.csv",
                name="armLevelCNA_profile",
                format=LocalFormat.CSV_MATRIX,
                encoding="utf8",
            ),
            UploadedFile(
                local_path=folder + "ms_repeat_profile.csv",
                name="ms_repeat_profile",
                format=LocalFormat.CSV_TABLE,
                encoding="utf8",
            ),
        ],
    )
    print("done")
    return wessegments, wgssegments


async def mutationPostProcessing(
    wesrefworkspace: str = env_config.WESCNWORKSPACE,
    wgsrefworkspace: str = env_config.WGSWORKSPACE,
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
    sv_af_cutoff: float = constants.SV_INTERNAL_AF_CUTOFF,
    run_guidemat: bool = True,
    upload_taiga: bool = True,
    hgnc_mapping_taiga: str = constants.HGNC_MAPPING_TABLE_TAIGAID,
    hgnc_mapping_table_name: str = constants.HGNC_MAPPING_TABLE_NAME,
    hgnc_mapping_table_version: int = constants.HGNC_MAPPING_TABLE_VERSION,
    omics_id_mapping_table_taigaid=constants.OMICS_ID_MAPPING_TABLE_TAIGAID,
    omics_id_mapping_table_name=constants.OMICS_ID_MAPPING_TABLE_NAME,
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
    client = create_taiga_client_v3()

    # read cds->pr mapping table and construct renaming dictionary
    # always read latest version
    print("reading omics ID mapping table from taiga")
    omics_id_mapping_table = client.get(
        name=omics_id_mapping_table_taigaid, file=omics_id_mapping_table_name
    )
    renaming_dict = dict(
        list(
            zip(
                omics_id_mapping_table["sequencing_id"],
                omics_id_mapping_table["profile_id"],
            )
        )
    )

    wes_wm = dm.WorkspaceManager(wesrefworkspace)
    wgs_wm = dm.WorkspaceManager(wgsrefworkspace)

    # doing wes
    print("DOING WES")
    folder = constants.WORKING_DIR + samplesetname + "/wes_"

    # TODO: replace with multiprocessing
    # ./sandbox/dna_eval/combine_mafs.py
    wesmutations, _, _ = mutations.postProcess(
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

    wesmutations_pr = wesmutations[
        wesmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace(
        {constants.SAMPLEID: renaming_dict, "Tumor_Sample_Barcode": renaming_dict}
    )

    # doing wgs
    print("DOING WGS")
    folder = constants.WORKING_DIR + samplesetname + "/wgs_"

    wgsmutations, wgssvs, wgs_sv_mat = mutations.postProcess(
        wgs_wm,
        sampleset="all",  # AllSamplesetName if AllSamplesetName else samplesetname,
        save_output=folder,
        sv_col=sv_col,
        sv_filename=sv_filename,
        mafcol=mafcol,
        run_sv=run_sv,
        sv_af_cutoff=sv_af_cutoff,
        debug=False,
        **kwargs,
    )

    wgsmutations_pr = wgsmutations[
        wgsmutations[constants.SAMPLEID].isin(renaming_dict.keys())
    ].replace(
        {constants.SAMPLEID: renaming_dict, "Tumor_Sample_Barcode": renaming_dict}
    )

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

    mergedmutations = mutations.addEntrez(
        mergedmutations, ensembl_col="EnsemblGeneID", entrez_col="EntrezGeneID"
    )

    # https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#somatic-maf-file-generation
    # For all columns, convert "Y" to True/False
    for col in mergedmutations.columns:
        if "Y" in mergedmutations[col].values:
            mergedmutations.loc[:, col] = np.where(
                mergedmutations[col].values == "Y", True, False
            )

    mergedmutations[list(mutcol.values()) + ["EntrezGeneID"]].to_csv(
        folder + "somatic_mutations.csv", index=False
    )

    if run_sv:
        if wgssvs is not None and wgs_sv_mat is not None:
            print("saving WGS svs")
            wgssvs.to_csv(folder + "svs.csv", index=False)
            wgssvs_pr = wgssvs[
                wgssvs[constants.SAMPLEID].isin(renaming_dict.keys())
            ].replace({constants.SAMPLEID: renaming_dict})
            wgssvs_pr.to_csv(folder + "svs_profile.csv", index=False)

            print("map entrez ids to SV matrix columns")
            # pull the gene id mapping table from taiga dataset maintained by the portal
            hgnc_table = client.get(
                name=hgnc_mapping_taiga,
                version=hgnc_mapping_table_version,
                file=hgnc_mapping_table_name,
            )
            # some rows in the table are missing entrez ids, replace them with "Unknown"
            hgnc_table["entrez_id"] = hgnc_table["entrez_id"].fillna("Unknown")
            hugo_entrez_pairs = list(zip(hgnc_table.symbol, hgnc_table.entrez_id))
            # generate a dictionary, key: hugo symbol, value: hugo symbol (entrez id)
            gene_renaming_dict = dict(
                [
                    (
                        (e[0], e[0] + " (" + str(int(e[1])) + ")")
                        if e[1] != "Unknown"
                        else (e[0], e[0] + " (Unknown)")
                    )
                    for e in hugo_entrez_pairs
                ]
            )
            wgs_sv_mat = wgs_sv_mat[
                list(set(wgs_sv_mat.columns) & set(gene_renaming_dict.keys()))
            ].rename(columns=gene_renaming_dict)
            wgs_sv_mat.to_csv(folder + "sv_mat_with_entrez.csv")
            wgs_sv_mat_pr = wgs_sv_mat[
                wgs_sv_mat.index.isin(renaming_dict.keys())
            ].rename(index=renaming_dict)
            wgs_sv_mat_pr.to_csv(folder + "sv_mat_with_entrez_profile.csv")
        else:
            print("no WGS SVs processed")

    merged = pd.concat([wgsmutations_pr, wesmutations_pr], axis=0).reset_index(
        drop=True
    )

    # For all columns, convert "Y" to True/False
    for col in merged.columns:
        if "Y" in merged[col].values:
            merged.loc[:, col] = np.where(merged[col].values == "Y", True, False)

    merged = merged.rename(columns=mutcol)
    merged = mutations.addEntrez(
        merged, ensembl_col="EnsemblGeneID", entrez_col="EntrezGeneID"
    )
    merged.to_csv(folder + "somatic_mutations_all_cols_profile.csv", index=False)
    merged[list(mutcol.values()) + ["EntrezGeneID"]].to_csv(
        folder + "somatic_mutations_profile.csv", index=False
    )
    merged[standardmafcol.keys()].to_csv(
        folder + "somatic_mutations_profile.maf.csv", index=False
    )

    # making genotyped mutation matrices
    print("creating mutation matrices")
    hotspot_mat, lof_mat = mutations.makeMatrices(merged)
    # add entrez ids to column names
    merged["gene_name"] = [
        (
            i["HugoSymbol"] + " (" + str(i["EntrezGeneID"]).split(".")[0] + ")"
            if i["EntrezGeneID"] != ""
            else i["HugoSymbol"] + " (Unknown)"
        )
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
        wes_germline_mats = mutations.aggregateGermlineMatrix(
            wes_wm, AllSamplesetName, save_output=folder
        )
        print("aggregating wgs")
        wgs_germline_mats = mutations.aggregateGermlineMatrix(
            wgs_wm, AllSamplesetName, save_output=folder
        )

        for lib, _ in bed_locations.items():
            assert lib in wes_germline_mats, "library missing in wes"
            assert lib in wgs_germline_mats, "library missing in wgs"
            # merging wes and wgs
            print("renaming merged wes and wgs germline matrix for library: ", lib)
            germline_mat_merged = pd.concat(
                [wes_germline_mats[lib], wgs_germline_mats[lib].iloc[:, 4:]], axis=1
            )
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
        client.update_dataset(
            reason="new " + samplesetname + " release!",
            permaname=taiga_dataset,
            additions=[
                UploadedFile(
                    local_path=folder
                    + "somatic_mutations_genotyped_hotspot_profile.csv",
                    name="somaticMutations_genotypedMatrix_hotspot_profile",
                    format=LocalFormat.CSV_MATRIX,
                    encoding="utf8",
                ),
                UploadedFile(
                    local_path=folder
                    + "somatic_mutations_genotyped_damaging_profile.csv",
                    name="somaticMutations_genotypedMatrix_damaging_profile",
                    format=LocalFormat.CSV_MATRIX,
                    encoding="utf8",
                ),
                UploadedFile(
                    local_path=folder + "somatic_mutations_profile.csv",
                    name="somaticMutations_profile",
                    format=LocalFormat.CSV_TABLE,
                    encoding="utf8",
                ),
                UploadedFile(
                    local_path=folder + "somatic_mutations_profile.maf.csv",
                    name="somaticMutations_profile_maf",
                    format=LocalFormat.CSV_TABLE,
                    encoding="utf8",
                ),
                UploadedFile(
                    local_path=folder + "somatic_mutations_all_cols_profile.csv",
                    name="somaticMutations_profile_all_cols",
                    format=LocalFormat.CSV_TABLE,
                    encoding="utf8",
                ),
                UploadedFile(
                    local_path=folder + "somatic_mutations.csv",
                    name="somaticMutations_withReplicates",
                    format=LocalFormat.CSV_TABLE,
                    encoding="utf8",
                ),
            ],
        )
        if run_guidemat:
            client.update_dataset(
                reason="new " + samplesetname + " release!",
                permaname=taiga_dataset,
                additions=[
                    UploadedFile(
                        local_path=folder + "binary_germline_avana.csv",
                        name="binary_mutation_avana",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "binary_germline_ky.csv",
                        name="binary_mutation_ky",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "binary_germline_humagne.csv",
                        name="binary_mutation_humagne",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "binary_germline_brunello.csv",
                        name="binary_mutation_brunello",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "binary_germline_tkov3.csv",
                        name="binary_mutation_tkov3",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                ],
            )
        if run_sv:
            client.update_dataset(
                reason="new " + samplesetname + " release!",
                permaname=taiga_dataset,
                additions=[
                    UploadedFile(
                        local_path=folder + "svs.csv",
                        name="structuralVariants_withReplicates",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "svs_profile.csv",
                        name="structuralVariants_profile",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "sv_mat_with_entrez.csv",
                        name="structuralVariants_geneLevelMatrix_withReplicates",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                    UploadedFile(
                        local_path=folder + "sv_mat_with_entrez_profile.csv",
                        name="structuralVariants_geneLevelMatrix_profile",
                        format=LocalFormat.CSV_TABLE,
                        encoding="utf8",
                    ),
                ],
            )
