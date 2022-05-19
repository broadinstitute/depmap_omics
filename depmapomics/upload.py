from __future__ import print_function
import pandas as pd
import os
from datetime import date

from genepy.utils import helper as h
from depmapomics import tracker
from depmapomics.config import *
from taigapy import TaigaClient


def getPRToRelease(trackerobj):
    """generate lists of profiles to release based on date for all portals
    
    Args:
        trackerobj (SampleTracker): tracker object

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    date_col_dict = {
        "internal": "InternalReleaseDate",
        "ibm": "IBMReleaseDate",
        "dmc": "ConsortiumReleaseDate",
        "public": "PublicReleaseDate",
    }
    today = int(str(date.today()).replace("-", ""))
    pr_table = trackerobj.read_pr_table()
    prs = dict()
    for k, v in date_col_dict.items():
        prs_with_date = pr_table[~(pr_table[v] == "")]
        prs[k] = prs_with_date[
            (prs_with_date[v].astype(int) <= today)
            & (prs_with_date.ProfileSource == "bam")
        ].index.tolist()
    return prs


def makeAchillesChoiceTable(
    trackerobj, prs, one_pr_per_type=True, source_priority=SOURCE_PRIORITY,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which MC

    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id
        folder (str): location where tables are saved as .csv files

    Returns:
        ach_tables (dict{(str: pd.DataFrame)}): for each portal, a df containing MC-PR mapping
    """
    pr_table = trackerobj.read_pr_table()
    seq_table = trackerobj.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
    subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
    mcs = set(subset_pr_table["ModelCondition"])
    # one_pr_per_type assumes we're only picking one PR per datatype (rna/dna) for each MC
    if one_pr_per_type:
        for mc in mcs:
            prs_in_mc = subset_pr_table[(subset_pr_table.ModelCondition == mc)]
            # rna
            if len(prs_in_mc[prs_in_mc.Datatype == "rna"]) == 1:
                pr = prs_in_mc[prs_in_mc.Datatype == "rna"].index[0]
                rows.append((mc, pr, "rna"))
            elif len(prs_in_mc[prs_in_mc.Datatype == "rna"]) > 1:
                cds_ids = prs_in_mc[prs_in_mc.Datatype == "rna"].CDSID.tolist()
                # at this point it is guaranteed that all cds_ids have different sources
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[0]
                rows.append((mc, pr, "rna"))
            # dna
            if (
                len(
                    prs_in_mc[
                        (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
                    ]
                )
                == 1
            ):
                pr = prs_in_mc[
                    (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
                ].index[0]
                rows.append((mc, pr, "dna"))
            elif (
                len(
                    prs_in_mc[
                        (prs_in_mc.Datatype == "wgs") | (prs_in_mc.Datatype == "wes")
                    ]
                )
                > 1
            ):
                cds_ids_wgs = prs_in_mc[prs_in_mc.Datatype == "wgs"].CDSID.tolist()
                cds_ids_wes = prs_in_mc[prs_in_mc.Datatype == "wes"].CDSID.tolist()
                pr = ""
                if len(cds_ids_wgs) == 0:
                    if len(prs_in_mc[prs_in_mc.Datatype == "wes"]) == 1:
                        pr = prs_in_mc[prs_in_mc.Datatype == "wes"].index[0]
                    else:
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
                        subset_seq_table.source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wes = subset_seq_table.loc[
                            cds_ids_wes, "source"
                        ].idxmin()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wes
                        ].index[0]
                else:
                    subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                    subset_seq_table.source = subset_seq_table.source.replace(
                        source_priority
                    )
                    latest_cds_id_wgs = subset_seq_table.loc[
                        cds_ids_wgs, "source"
                    ].idxmin()
                    pr = subset_pr_table[
                        subset_pr_table.CDSID == latest_cds_id_wgs
                    ].index[0]
                rows.append((mc, pr, "dna"))
    ach_table = pd.DataFrame(
        rows, columns=["ModelConditionID", "ProfileID", "ProfileType"]
    )

    return ach_table


def makeDefaultModelTable(
    trackerobj, prs, one_pr_per_type=True, source_priority=SOURCE_PRIORITY,
):
    """generate a table for each portal that indicates which profiles are released corresponding to which modelID

    Args:
        trackerobj (SampleTracker): tracker object
        folder (str): location where tables are saved as .csv files

    Returns:
        ach_tables (dict{(str: pd.DataFrame)}): for each portal, a df containing MC-PR mapping"""
    pr_table = trackerobj.read_pr_table()
    mc_table = trackerobj.read_mc_table()
    seq_table = trackerobj.read_seq_table()
    source_priority = {source_priority[i]: i for i in range(len(source_priority))}
    rows = []
    subset_pr_table = pr_table.loc[prs]
    subset_pr_table = subset_pr_table[subset_pr_table.BlacklistOmics != 1]
    mcs = set(subset_pr_table["ModelCondition"])
    models = set(mc_table.loc[mcs].ModelID)
    # one_pr_per_type assumes we're only picking one PR per datatype (rna/dna) for each MC
    if one_pr_per_type:
        for m in models:
            subset_mc_table = mc_table[mc_table.ModelID == m]
            mcs_in_model = subset_mc_table.index.tolist()
            prs_in_model = subset_pr_table[
                (subset_pr_table.ModelCondition.isin(mcs_in_model))
            ]
            # rna
            if len(prs_in_model[prs_in_model.Datatype == "rna"]) == 1:
                pr = prs_in_model[prs_in_model.Datatype == "rna"].index[0]
                rows.append((m, pr, "rna"))
            elif len(prs_in_model[prs_in_model.Datatype == "rna"]) > 1:
                cds_ids = prs_in_model[prs_in_model.Datatype == "rna"].CDSID.tolist()
                subset_seq_table = seq_table[seq_table.index.isin(cds_ids)]
                subset_seq_table.source = subset_seq_table.source.replace(
                    source_priority
                )
                latest_cds_id = subset_seq_table.loc[cds_ids, "source"].idxmin()
                pr = subset_pr_table[subset_pr_table.CDSID == latest_cds_id].index[0]
                rows.append((m, pr, "rna"))
            # dna
            if (
                len(
                    prs_in_model[
                        (prs_in_model.Datatype == "wgs")
                        | (prs_in_model.Datatype == "wes")
                    ]
                )
                == 1
            ):
                pr = prs_in_model[
                    (prs_in_model.Datatype == "wgs") | (prs_in_model.Datatype == "wes")
                ].index[0]
                rows.append((m, pr, "dna"))
            elif (
                len(
                    prs_in_model[
                        (prs_in_model.Datatype == "wgs")
                        | (prs_in_model.Datatype == "wes")
                    ]
                )
                > 1
            ):
                # assuming SANGER doesn't have wgs
                cds_ids_wgs = prs_in_model[
                    prs_in_model.Datatype == "wgs"
                ].CDSID.tolist()
                cds_ids_wes = prs_in_model[
                    (prs_in_model.Datatype == "wes") & (prs_in_model.CDSID != "")
                ].CDSID.tolist()  # CDSID is '' when the profile is in legacy
                pr = ""
                # if no wgs, look at MC table and select broad wes over sanger wes
                if len(cds_ids_wgs) == 0:
                    if len(prs_in_model[prs_in_model.Datatype == "wes"]) == 1:
                        pr = prs_in_model[prs_in_model.Datatype == "wes"].index[0]
                    else:
                        subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wes)]
                        subset_seq_table.source = subset_seq_table.source.replace(
                            source_priority
                        )
                        latest_cds_id_wes = subset_seq_table.loc[
                            cds_ids_wes, "source"
                        ].idxmin()
                        pr = subset_pr_table[
                            subset_pr_table.CDSID == latest_cds_id_wes
                        ].index[0]
                # if there is wgs, always select wgs
                else:
                    subset_seq_table = seq_table[seq_table.index.isin(cds_ids_wgs)]
                    subset_seq_table.source = subset_seq_table.source.replace(
                        source_priority
                    )
                    latest_cds_id_wgs = subset_seq_table.loc[
                        cds_ids_wgs, "source"
                    ].idxmin()
                    pr = subset_pr_table[
                        subset_pr_table.CDSID == latest_cds_id_wgs
                    ].index[0]
                rows.append((m, pr, "dna"))
    default_table = pd.DataFrame(rows, columns=["ModelID", "ProfileID", "ProfileType"])
    return default_table


def initVirtualDatasets(
    samplesetname=SAMPLESETNAME, taiga_folder_id=VIRTUAL_FOLDER, portals=DATASETS
):
    """initialize taiga virtual datasets for all 4 portals by uploading an empty dummy file
    """

    with open("temp/dummy.csv", "w") as fp:
        pass
    virutal_pr = dict()
    virtual_model = dict()
    tc = TaigaClient()
    for p in portals:
        virutal_pr[p] = tc.create_dataset(
            p + "_" + samplesetname + "_profile",
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "temp/dummy.csv",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
        virtual_model[p] = tc.create_dataset(
            p + "_" + samplesetname + "_model",
            dataset_description=samplesetname
            + " release of the DepMap dataset for the DepMap Portal. Please look at the README file for additional information about this dataset. ",
            upload_files=[
                {
                    "path": "temp/dummy.csv",
                    "name": "init",
                    "format": "Raw",
                    "encoding": "utf-8",
                }
            ],
            folder_id=taiga_folder_id,
        )
    return virutal_pr, virtual_model


def uploadCNMatricesPR(prs, portal, taiga_latest=TAIGA_CN, taiga_virtual=""):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """
    folder = "temp/" + portal + "/profile/"
    h.createFoldersFor(folder)
    # load cds-id indexed matrices for the current quarter
    print("CN: loading from taiga latest")
    tc = TaigaClient()
    genecn = tc.get(name=taiga_latest, file="merged_genecn_all_profile")
    segmentcn = tc.get(name=taiga_latest, file="merged_segments_all_profile")
    wescn = tc.get(name=taiga_latest, file="wes_genecn_all_profile")
    wessegment = tc.get(name=taiga_latest, file="wes_segments_all_profile")

    # subset and rename
    print("CN: subsetting and renaming")
    segmentcn_renamed = segmentcn[segmentcn.DepMap_ID.isin(prs)]
    segmentcn_renamed.to_csv(folder + "all_merged_segments.csv", index=False)
    genecn_renamed = genecn[genecn.index.isin(prs)]
    genecn_renamed.to_csv(folder + "all_merged_genes_cn.csv")
    wescn_renamed = wescn[wescn.index.isin(prs)]
    wescn_renamed.to_csv(folder + "wes_genes_cn.csv")
    wessegment_renamed = wessegment[wessegment.DepMap_ID.isin(prs)]
    wessegment_renamed.to_csv(folder + "wes_segments.csv", index=False)

    # upload to taiga
    print("CN: uploading to taiga")
    tc = TaigaClient()
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding CN data",
        upload_files=[
            {
                "path": folder + "all_merged_genes_cn.csv",
                "name": "CCLE_gene_cn",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "all_merged_segments.csv",
                "name": "CCLE_segment_cn",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_genes_cn.csv",
                "name": "CCLE_wes_gene_cn",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_segments.csv",
                "name": "CCLE_wes_segment_cn",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadCNMatricesModel(
    pr2model_dict, portal, taiga_latest=TAIGA_CN, taiga_virtual=""
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """
    folder = "temp/" + portal + "/model/"
    h.createFoldersFor(folder)
    # load cds-id indexed matrices for the current quarter
    print("CN: loading from taiga latest")
    tc = TaigaClient()
    genecn = tc.get(name=taiga_latest, file="merged_genecn_all_profile")
    segmentcn = tc.get(name=taiga_latest, file="merged_segments_all_profile")
    wescn = tc.get(name=taiga_latest, file="wes_genecn_all_profile")
    wessegment = tc.get(name=taiga_latest, file="wes_segments_all_profile")

    # subset and rename
    print("CN: subsetting and renaming")
    segmentcn_renamed = segmentcn[
        segmentcn.DepMap_ID.isin(set(pr2model_dict.keys()))
    ].replace({"DepMap_ID": pr2model_dict})
    segmentcn_renamed.to_csv(folder + "all_merged_segments.csv", index=False)
    genecn_renamed = genecn[genecn.index.isin(set(pr2model_dict.keys()))].rename(
        index=pr2model_dict
    )
    genecn_renamed.to_csv(folder + "all_merged_genes_cn.csv")
    wescn_renamed = wescn[wescn.index.isin(set(pr2model_dict.keys()))].rename(
        index=pr2model_dict
    )
    wescn_renamed.to_csv(folder + "wes_genes_cn.csv")
    wessegment_renamed = wessegment[
        wessegment.DepMap_ID.isin(set(pr2model_dict.keys()))
    ].replace({"DepMap_ID": pr2model_dict})
    wessegment_renamed.to_csv(folder + "wes_segments.csv", index=False)

    # upload to taiga
    print("CN: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding CN data",
        upload_files=[
            {
                "path": folder + "all_merged_genes_cn.csv",
                "name": "CCLE_gene_cn",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "all_merged_segments.csv",
                "name": "CCLE_segment_cn",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_genes_cn.csv",
                "name": "CCLE_wes_gene_cn",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "wes_segments.csv",
                "name": "CCLE_wes_segment_cn",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadMutationMatrices(
    renaming_dict, taiga_virtual="", folder="temp/" + SAMPLESETNAME
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """

    # load cds-id indexed matrices for the current quarter
    print("Mutation: loading")
    mutations = pd.read_csv(folder + "/merged_somatic_mutations_withlegacy.csv")
    damaging = pd.read_csv(
        folder + "/merged_somatic_mutations_boolmatrix_fordepmap_damaging.csv",
        index_col=0,
    )
    othercons = pd.read_csv(
        folder + "/merged_somatic_mutations_boolmatrix_fordepmap_othercons.csv",
        index_col=0,
    )
    othernoncons = pd.read_csv(
        folder + "/merged_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
        index_col=0,
    )
    hotspot = pd.read_csv(
        folder + "/merged_somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
        index_col=0,
    )

    # subset and rename
    print("Mutation: subsetting and renaming")
    mutations_renamed = mutations[
        mutations.DepMap_ID.isin(set(renaming_dict.keys()))
    ].replace({"DepMap_ID": renaming_dict})
    mutations_renamed.to_csv("temp/all_somatic_mutations_withlegacy.csv", index=False)
    damaging_renamed = damaging[damaging.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    damaging_renamed.to_csv(
        "temp/all_somatic_mutations_boolmatrix_fordepmap_damaging.csv"
    )
    othercons_renamed = othercons[
        othercons.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    othercons_renamed.to_csv(
        "temp/all_somatic_mutations_boolmatrix_fordepmap_othercons.csv"
    )
    othernoncons_renamed = othernoncons[
        othernoncons.index.isin(set(renaming_dict.keys()))
    ].rename(index=renaming_dict)
    othernoncons_renamed.to_csv(
        "temp/all_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv"
    )
    hotspot_renamed = hotspot[hotspot.index.isin(set(renaming_dict.keys()))].rename(
        index=renaming_dict
    )
    hotspot_renamed.to_csv(
        "temp/all_somatic_mutations_boolmatrix_fordepmap_hotspot.csv"
    )

    # upload to taiga
    print("Mutation: uploading to taiga")
    tc = TaigaClient()
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding mutations",
        upload_files=[
            {
                "path": "temp/all_somatic_mutations_withlegacy.csv",
                "name": "CCLE_mutations",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                "name": "CCLE_mutations_bool_damaging",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                "name": "CCLE_mutations_bool_nonconserving",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_othercons.csv",
                "name": "CCLE_mutations_bool_otherconserving",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                "name": "CCLE_mutations_bool_hotspot",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadExpressionMatricesPR(
    prs, portal, include_ssgsea=False, taiga_latest=TAIGA_EXPRESSION, taiga_virtual="",
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
        include_ssgsea (bool): whether or not to upload ssGSEA data to this portal
    """
    folder = "temp/" + portal + "/profile/"
    h.createFoldersFor(folder)
    # load pr-id indexed matrices for the current quarter from latest
    print("Expression: loading from taiga latest")
    tc = TaigaClient()
    genes_expected_count = tc.get(
        name=taiga_latest, file="genes_expected_count_profile"
    )
    genes_tpm = tc.get(name=taiga_latest, file="genes_tpm_profile_logp1")
    proteincoding_genes_expected_count = tc.get(
        name=taiga_latest, file="proteincoding_genes_expected_count_profile"
    )
    proteincoding_genes_tpm = tc.get(
        name=taiga_latest, file="proteincoding_genes_tpm_profile_logp1"
    )
    transcripts_expected_count = tc.get(
        name=taiga_latest, file="transcripts_expected_count_profile"
    )
    transcripts_tpm = tc.get(name=taiga_latest, file="transcripts_tpm_profile_logp1")

    # subset
    print("Expression: subsetting and renaming")
    genes_expected_count_renamed = genes_expected_count[
        genes_expected_count.index.isin(prs)
    ]
    genes_expected_count_renamed = genes_expected_count_renamed[
        (genes_expected_count_renamed.sum(1) != 0)
        & (genes_expected_count_renamed.var(1) != 0)
    ]
    genes_expected_count_renamed.to_csv(folder + "expression_genes_expected_count.csv")

    genes_tpm_renamed = genes_tpm[genes_tpm.index.isin(prs)]
    genes_tpm_renamed = genes_tpm_renamed[
        (genes_tpm_renamed.sum(1) != 0) & (genes_tpm_renamed.var(1) != 0)
    ]
    genes_tpm_renamed.to_csv(folder + "expression_genes_tpm.csv")

    proteincoding_genes_tpm_renamed = proteincoding_genes_tpm[
        proteincoding_genes_tpm.index.isin(prs)
    ]
    proteincoding_genes_tpm_renamed = proteincoding_genes_tpm_renamed[
        (proteincoding_genes_tpm_renamed.sum(1) != 0)
        & (proteincoding_genes_tpm_renamed.var(1) != 0)
    ]
    proteincoding_genes_tpm_renamed.to_csv(
        folder + "expression_proteincoding_genes_tpm.csv"
    )
    transcripts_tpm_renamed = transcripts_tpm[transcripts_tpm.index.isin(prs)]
    transcripts_tpm_renamed = transcripts_tpm_renamed[
        (transcripts_tpm_renamed.sum(1) != 0) & (transcripts_tpm_renamed.var(1) != 0)
    ]
    transcripts_tpm_renamed.to_csv(folder + "expression_transcripts_tpm.csv")
    proteincoding_genes_expected_count_renamed = proteincoding_genes_expected_count[
        proteincoding_genes_expected_count.index.isin(prs)
    ]
    proteincoding_genes_expected_count_renamed = proteincoding_genes_expected_count_renamed[
        (proteincoding_genes_expected_count_renamed.sum(1) != 0)
        & (proteincoding_genes_expected_count_renamed.var(1) != 0)
    ]
    proteincoding_genes_expected_count_renamed.to_csv(
        folder + "expression_proteincoding_genes_expected_count.csv"
    )
    transcripts_expected_count_renamed = transcripts_expected_count[
        transcripts_expected_count.index.isin(prs)
    ]
    transcripts_expected_count_renamed = transcripts_expected_count_renamed[
        (transcripts_expected_count_renamed.sum(1) != 0)
        & (transcripts_expected_count_renamed.var(1) != 0)
    ]
    transcripts_expected_count_renamed.to_csv(
        folder + "expression_transcripts_expected_count.csv"
    )

    # upload to taiga
    print("Expression: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding expression",
        upload_files=[
            {
                "path": folder + "expression_genes_expected_count.csv",
                "name": "CCLE_RNAseq_reads",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_transcripts_tpm.csv",
                "name": "CCLE_RNAseq_transcripts",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_genes_tpm.csv",
                "name": "CCLE_expression_full",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_proteincoding_genes_tpm.csv",
                "name": "CCLE_expression",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_proteincoding_genes_expected_count.csv",
                "name": "CCLE_expression_proteincoding_genes_expected_count",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_transcripts_expected_count.csv",
                "name": "CCLE_expression_transcripts_expected_count",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        upload_async=False,
        add_all_existing_files=True,
    )
    if include_ssgsea:
        print("ssGSEA: loading")
        enrichments = pd.read_csv(folder + "/gene_sets_all.csv", index_col=0)
        print("ssGSEA: subsetting and renaming")
        enrichments_renamed = enrichments[enrichments.index.isin(prs)]
        enrichments_renamed.to_csv("temp/enrichments_ssGSEA.csv")
        print("ssGSEA: uploading to taiga")
        tc.update_dataset(
            dataset_id=taiga_virtual,
            changes_description="adding enrichments",
            upload_files=[
                {
                    "path": "temp/enrichments_ssGSEA.csv",
                    "name": "CCLE_ssGSEA",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
            ],
            upload_async=False,
            add_all_existing_files=True,
        )


def uploadExpressionMatricesModel(
    pr2model_dict,
    portal,
    include_ssgsea=False,
    taiga_latest=TAIGA_EXPRESSION,
    taiga_virtual="",
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
        include_ssgsea (bool): whether or not to upload ssGSEA data to this portal
    """
    folder = "temp/" + portal + "/model/"
    h.createFoldersFor(folder)
    # load pr-id indexed matrices for the current quarter from latest
    print("Expression: loading from taiga latest")
    tc = TaigaClient()
    genes_expected_count = tc.get(
        name=taiga_latest, file="genes_expected_count_profile"
    )
    genes_tpm = tc.get(name=taiga_latest, file="genes_tpm_profile_logp1")
    proteincoding_genes_expected_count = tc.get(
        name=taiga_latest, file="proteincoding_genes_expected_count_profile"
    )
    proteincoding_genes_tpm = tc.get(
        name=taiga_latest, file="proteincoding_genes_tpm_profile_logp1"
    )
    transcripts_expected_count = tc.get(
        name=taiga_latest, file="transcripts_expected_count_profile"
    )
    transcripts_tpm = tc.get(name=taiga_latest, file="transcripts_tpm_profile_logp1")

    # subset and rename
    print("Expression: subsetting and renaming")
    genes_expected_count_renamed = genes_expected_count[
        genes_expected_count.index.isin(pr2model_dict.keys())
    ].rename(index=pr2model_dict)
    genes_expected_count_renamed = genes_expected_count_renamed[
        (genes_expected_count_renamed.sum(1) != 0)
        & (genes_expected_count_renamed.var(1) != 0)
    ]
    genes_expected_count_renamed.to_csv(folder + "expression_genes_expected_count.csv")

    genes_tpm_renamed = genes_tpm[genes_tpm.index.isin(pr2model_dict.keys())].rename(
        index=pr2model_dict
    )
    genes_tpm_renamed = genes_tpm_renamed[
        (genes_tpm_renamed.sum(1) != 0) & (genes_tpm_renamed.var(1) != 0)
    ]
    genes_tpm_renamed.to_csv(folder + "expression_genes_tpm.csv")

    proteincoding_genes_tpm_renamed = proteincoding_genes_tpm[
        proteincoding_genes_tpm.index.isin(pr2model_dict.keys())
    ].rename(index=pr2model_dict)
    proteincoding_genes_tpm_renamed = proteincoding_genes_tpm_renamed[
        (proteincoding_genes_tpm_renamed.sum(1) != 0)
        & (proteincoding_genes_tpm_renamed.var(1) != 0)
    ]
    proteincoding_genes_tpm_renamed.to_csv(
        folder + "expression_proteincoding_genes_tpm.csv"
    )

    transcripts_tpm_renamed = transcripts_tpm[
        transcripts_tpm.index.isin(pr2model_dict.keys())
    ].rename(index=pr2model_dict)
    transcripts_tpm_renamed = transcripts_tpm_renamed[
        (transcripts_tpm_renamed.sum(1) != 0) & (transcripts_tpm_renamed.var(1) != 0)
    ]
    transcripts_tpm_renamed.to_csv(folder + "expression_transcripts_tpm.csv")

    proteincoding_genes_expected_count_renamed = proteincoding_genes_expected_count[
        proteincoding_genes_expected_count.index.isin(pr2model_dict.keys())
    ].rename(index=pr2model_dict)
    proteincoding_genes_expected_count_renamed = proteincoding_genes_expected_count_renamed[
        (proteincoding_genes_expected_count_renamed.sum(1) != 0)
        & (proteincoding_genes_expected_count_renamed.var(1) != 0)
    ]
    proteincoding_genes_expected_count_renamed.to_csv(
        folder + "expression_proteincoding_genes_expected_count.csv"
    )

    transcripts_expected_count_renamed = transcripts_expected_count[
        transcripts_expected_count.index.isin(pr2model_dict.keys())
    ].rename(index=pr2model_dict)
    transcripts_expected_count_renamed = transcripts_expected_count_renamed[
        (transcripts_expected_count_renamed.sum(1) != 0)
        & (transcripts_expected_count_renamed.var(1) != 0)
    ]
    transcripts_expected_count_renamed.to_csv(
        folder + "expression_transcripts_expected_count.csv"
    )

    # upload to taiga
    print("Expression: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding expression",
        upload_files=[
            {
                "path": folder + "expression_genes_expected_count.csv",
                "name": "CCLE_RNAseq_reads",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_transcripts_tpm.csv",
                "name": "CCLE_RNAseq_transcripts",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_genes_tpm.csv",
                "name": "CCLE_expression_full",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_proteincoding_genes_tpm.csv",
                "name": "CCLE_expression",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_proteincoding_genes_expected_count.csv",
                "name": "CCLE_expression_proteincoding_genes_expected_count",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "expression_transcripts_expected_count.csv",
                "name": "CCLE_expression_transcripts_expected_count",
                "format": "NumericMatrixCSV",
                "encoding": "utf-8",
            },
        ],
        upload_async=False,
        add_all_existing_files=True,
    )
    if include_ssgsea:
        print("ssGSEA: loading")
        enrichments = pd.read_csv(folder + "/gene_sets_all.csv", index_col=0)
        print("ssGSEA: subsetting and renaming")
        enrichments_renamed = enrichments[enrichments.index.isin(pr2model_dict.keys())]
        enrichments_renamed.to_csv("temp/enrichments_ssGSEA.csv")
        print("ssGSEA: uploading to taiga")
        tc.update_dataset(
            dataset_id=taiga_virtual,
            changes_description="adding enrichments",
            upload_files=[
                {
                    "path": "temp/enrichments_ssGSEA.csv",
                    "name": "CCLE_ssGSEA",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
            ],
            upload_async=False,
            add_all_existing_files=True,
        )


def uploadFusionMatricesPR(
    prs,
    portal,
    taiga_latest=TAIGA_FUSION,
    taiga_virtual="",
    countCol="CCLE_count",
    maxfreq=FUSION_MAXFREQ,
    sampleCol=SAMPLEID,
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """
    folder = "temp/" + portal + "/profile/"
    h.createFoldersFor(folder)
    # load cds-id indexed matrices for the current quarter
    print("Fusion: loading from taiga latest")
    tc = TaigaClient()
    fusions = tc.get(name=taiga_latest, file="fusions_all_profile")
    filtered = tc.get(name=taiga_latest, file="filteredfusions_latest_profile")

    # subset and rename
    print("Fusion: subsetting and renaming")
    fusions_renamed = fusions[fusions.DepMap_ID.isin(prs)]
    fusions_renamed.to_csv(folder + "fusions.csv", index=False)

    filtered_renamed = filtered[filtered.DepMap_ID.isin(prs)]
    filtered_renamed = filtered_renamed[
        filtered_renamed[countCol] < len(set(filtered_renamed[sampleCol])) * maxfreq
    ]
    filtered_renamed.to_csv(folder + "filtered_fusions.csv", index=False)

    # upload to taiga
    print("Fusion: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding fusions",
        upload_files=[
            {
                "path": folder + "fusions.csv",
                "name": "CCLE_fusions_unfiltered",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "filtered_fusions.csv",
                "name": "CCLE_fusions",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadFusionMatricesModel(
    pr2model_dict,
    portal,
    taiga_latest=TAIGA_FUSION,
    taiga_virtual="",
    countCol="CCLE_count",
    maxfreq=FUSION_MAXFREQ,
    sampleCol=SAMPLEID,
):
    """subset, rename, save and upload to taiga CN matrices
    
    Args:
        renaming_dict (dict): renaming scheme mapping from CDS-id to either model or MC id
        taiga_id (str): which dataset the matrices are being uploaded to
        folder (str): where the files to be subsetted are saved
    """
    folder = "temp/" + portal + "/model/"
    h.createFoldersFor(folder)
    # load cds-id indexed matrices for the current quarter
    print("Fusion: loading from taiga latest")
    tc = TaigaClient()
    fusions = tc.get(name=taiga_latest, file="fusions_all_profile")
    filtered = tc.get(name=taiga_latest, file="filteredfusions_latest_profile")

    # subset and rename
    print("Fusion: subsetting and renaming")
    fusions_renamed = fusions[
        fusions.DepMap_ID.isin(set(pr2model_dict.keys()))
    ].replace({sampleCol: pr2model_dict})
    fusions_renamed.to_csv(folder + "fusions.csv", index=False)

    filtered_renamed = filtered[
        filtered.DepMap_ID.isin(set(pr2model_dict.keys()))
    ].replace({sampleCol: pr2model_dict})
    filtered_renamed = filtered_renamed[
        filtered_renamed[countCol] < len(set(filtered_renamed[sampleCol])) * maxfreq
    ]
    filtered_renamed.to_csv(folder + "filtered_fusions.csv", index=False)

    # upload to taiga
    print("Fusion: uploading to taiga")
    tc.update_dataset(
        dataset_id=taiga_virtual,
        changes_description="adding fusions",
        upload_files=[
            {
                "path": folder + "fusions.csv",
                "name": "CCLE_fusions_unfiltered",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
            {
                "path": folder + "filtered_fusions.csv",
                "name": "CCLE_fusions",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )


def uploadAuxTables(trackerobj, taiga_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME):
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs in prs_allportals.items():
        achilles_table = makeAchillesChoiceTable(trackerobj, prs)
        default_table = makeDefaultModelTable(trackerobj, prs)
        achilles_table.to_csv(
            folder + portal + "_achilles_choice_table.csv", index=False
        )
        default_table.to_csv(folder + portal + "_default_model_table.csv", index=False)
        tc = TaigaClient()
        tc.update_dataset(
            dataset_id=taiga_ids[portal],
            changes_description="adding mapping tables",
            upload_files=[
                {
                    "path": folder + "/" + portal + "_achilles_choice_table.csv",
                    "name": "Achilles_choice_table",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": folder + "/" + portal + "_default_model_table.csv",
                    "name": "default_model_table",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ],
            add_all_existing_files=True,
        )


def makePRLvMatrices(trackerobj, taiga_ids=VIRTUAL):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs_to_release in prs_allportals.items():
        print("uploading profile-level matrices to ", portal)
        uploadCNMatricesPR(prs_to_release, portal, taiga_virtual=taiga_ids[portal])
        # uploadMutationMatrices(prs_to_release, portal, taiga_virtual=taiga_ids[portal])
        uploadFusionMatricesPR(prs_to_release, portal, taiga_virtual=taiga_ids[portal])
        # TODO: don't have ssgsea data for test run yet
        if portal == "":
            uploadExpressionMatricesPR(
                prs_to_release,
                portal,
                include_ssgsea=True,
                taiga_virtual=taiga_ids[portal],
            )
        else:
            uploadExpressionMatricesPR(
                prs_to_release,
                portal,
                include_ssgsea=False,
                taiga_virtual=taiga_ids[portal],
            )


def makeModelLvMatrices(trackerobj, taiga_ids=VIRTUAL, folder="temp/" + SAMPLESETNAME):
    """for each portal, save and upload profile-indexed data matrices
    
    Args:
        trackerobj (SampleTracker): tracker object
        taiga_ids (dict): dictionary that maps portal name to virtual taiga dataset id

    Returns:
        prs (dict{(portal: list of PRs)}): for each portal, list of profile IDs
    """
    prs_allportals = getPRToRelease(trackerobj)
    for portal, prs_to_release in prs_allportals.items():
        default_table = makeDefaultModelTable(trackerobj, prs_to_release)
        pr2model_dict = dict(list(zip(default_table.ProfileID, default_table.ModelID)))
        h.dictToFile(pr2model_dict, folder + "/" + portal + "_pr2model_renaming.json")
        print("uploading model-level matrices to", portal)
        uploadCNMatricesModel(pr2model_dict, taiga_virtual=taiga_ids[portal])
        # uploadMutationMatricesModel(pr2model_dict, taiga_virtual=taiga_ids[portal])
        uploadFusionMatricesModel(
            pr2model_dict, portal, taiga_virtual=taiga_ids[portal]
        )
        # TODO: don't have ssgsea data for test run yet
        if portal == "":
            uploadExpressionMatricesModel(
                pr2model_dict,
                portal,
                include_ssgsea=True,
                taiga_virtual=taiga_ids[portal],
            )
        else:
            uploadExpressionMatricesModel(
                pr2model_dict,
                portal,
                include_ssgsea=False,
                taiga_virtual=taiga_ids[portal],
            )


def findLatestVersion(dataset, approved_only=True):
    highest = 0
    latest_version = 0
    tc = TaigaClient()
    data = tc.get_dataset_metadata(dataset)
    for val in data["versions"]:
        if val["state"] == "approved" or not approved_only:
            if int(val["name"]) > highest:
                highest = int(val["name"])
                latest_version = highest
    if latest_version == 0:
        raise ValueError("could not find a version")
    return data["permanames"][0] + "." + str(latest_version)


def updateEternal(
    eternal_id="depmap-a0ab", virtual=VIRTUAL, samplesetname=SAMPLESETNAME
):
    """update taiga eternal dataset by linking to latest virtual internal dataset"""
    latest_version = findLatestVersion(virtual["internal"])

    files = [
        "CCLE_gene_cn",
        "CCLE_segment_cn",
        "CCLE_mutations",
        "CCLE_mutations_bool_damaging",
        "CCLE_mutations_bool_nonconserving",
        "CCLE_mutations_bool_otherconserving",
        "CCLE_mutations_bool_hotspot",
        "CCLE_expression_full",
        "CCLE_RNAseq_transcripts",
        "CCLE_RNAseq_reads",
        "CCLE_expression",
        "CCLE_expression_proteincoding_genes_expected_count",
        "CCLE_expression_transcripts_expected_count",
        "CCLE_fusions_unfiltered",
        "CCLE_fusions",
        "CCLE_ssGSEA",
        "CCLE_wes_gene_cn",
        "CCLE_wes_segment_cn",
    ]

    tc = TaigaClient()
    tc.update_dataset(
        eternal_id,
        changes_description="new " + samplesetname + " omics dataset.",
        add_taiga_ids=[
            {"taiga_id": latest_version + "/" + file, "name": file} for file in files
        ],
        add_all_existing_files=True,
    )


def CCLEupload(trackerobj, taiga_ids=""):
    if taiga_ids == "":
        taiga_ids = initVirtualDatasets()

    makePRLvMatrices(trackerobj, taiga_ids=taiga_ids)
    makeModelLvMatrices(trackerobj, taiga_ids=taiga_ids)
    uploadAuxTables(trackerobj, taiga_ids=taiga_ids)
    updateEternal(virtual=taiga_ids)
