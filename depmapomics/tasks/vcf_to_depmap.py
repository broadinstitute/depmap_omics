from mgenepy import mutations
import pyarrow as pa
from mgenepy.utils import helper as h
import os
import pandas as pd
import pandas.api.types as ptypes
import argparse
import pyarrow.parquet as pq


import re
import numpy as np


def to_bool(x):
    "Like bool(x) but more paranoid about input values. ie: python's bool('false') is True"
    x = x.lower()
    if x == "true":
        return True
    assert x == "false"
    return False


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_filename")
    parser.add_argument("sample_name", nargs="?", default=None)
    parser.add_argument("--n_rows", default=500_000, type=int)
    parser.add_argument("--use_multi", default=False, type=to_bool)
    parser.add_argument("--force_keep", default=[], type=lambda x: x.split(","))
    parser.add_argument("--whitelist", default=False, type=to_bool)
    parser.add_argument("--drop_clustered_events", default=True, type=to_bool)
    parser.add_argument("--version", default="", type=str)
    args = parser.parse_args()

    vcf_filename = args.vcf_filename

    sample_name = args.sample_name
    if not sample_name:
        sample_name = vcf_filename.split("/")[-1].split(".")[0]

    n_rows = args.n_rows
    use_multi = args.use_multi
    force_keep = args.force_keep
    whitelist = args.whitelist
    drop_clustered_events = args.drop_clustered_events
    version = args.version

    prev_cols = []

    print(
        "inputs: vcf_filename:",
        vcf_filename,
        ", sample_name:",
        sample_name,
        ", n_rows:",
        n_rows,
        ", use_multi:",
        use_multi,
        ", force_keep:",
        force_keep,
        ", version:",
        version,
    )

    tobreak = False

    loc = os.path.dirname(os.path.abspath(__file__))
    oncogene = h.fileToList(loc + "/oncokb_dm/data/oncogene_oncokb_20250205.txt")
    tumor_suppressor_list = h.fileToList(
        loc + "/oncokb_dm/data/tumor_suppressor_oncokb_20250205.txt"
    )
    civic_df = pd.read_csv(loc + "/civic_export_09212022.csv").drop(
        columns=["chromosome_37", "start_37"]
    )

    """
    we are running through these likely very large files by loading a chunk at a time

    the issue is to make sure that each chunk is the same as the previous chunk (we don't remove differ
    set of columns etc..) 

    """

    processed_count = 0
    for i in range(10_000):
        # boolean annotations in vcf
        bool_cols = [
            "PON",
            "ASP",
            "R3",
            "R5",
            "VLD",
            "INT",
            "SLO",
            "KGPhase1",
            "KGPhase3",
            "GNO",
            "HD",
            "RV",
            "G5",
            "G5A",
            "NSM",
            "REF",
            "NOV",
            "SYN",
            "S3D",
            "LSD",
            "PM",
            "PMC",
            "U3",
            "U5",
            "NSF",
            "NSN",
            "NOC",
            "MTP",
            "DSS",
            "TPA",
            "OTH",
            "ASS",
            "CFL",
            "OM",
        ]
        # read in vcf as a df
        vcf_file, _, _ = mutations.vcf_to_df(
            vcf_filename,
            additional_cols=bool_cols,
            parse_filter=True,
            force_keep=force_keep
            + list(TO_RENAME_OC.keys())
            + list(TO_RENAME_HGVS.keys())
            + list(TO_RENAME_BASE.keys()),
            drop_null=False,
            cols_to_drop=[
                "clinvar_vcf_mc",
                "oreganno_build",
                "gt",
                "ad",
                "af",
                "dp",
                "f1r2",
                "f2r1",
                "fad",
                "sb",
                "pid",
                "pl",
                "ps",
                "gq",
                "pgt",
                "gencode_34_chromosome",
            ],
            nrows=n_rows,
            skiprows=n_rows * i,
        )
        if "PID" not in vcf_file.columns.tolist():
            vcf_file["PID"] = ""
        filen = len(vcf_file)
        processed_count += filen
        if filen < n_rows:
            # we have reached the end:
            tobreak = True

        # improve
        vcf_file = improve(
            vcf_file,
            force_list=["oc_genehancer__feature_name"],
            split_multiallelic=use_multi,
            oncogene_list=oncogene,
            tumor_suppressor_list=tumor_suppressor_list,
            civic_df=civic_df,
        )

        # checking we have the same set of columns
        cols = vcf_file.columns.tolist()
        if i == 0:
            prev_cols = cols
        elif set(cols) != set(prev_cols):
            raise ValueError(
                "we are removing different sets of columns",
                cols,
                list(set(cols) ^ set(prev_cols)),
            )
        elif len(cols) != len(prev_cols):
            raise ValueError("some columns have duplicate values", prev_cols, cols)
        elif cols != prev_cols:
            vcf_file = vcf_file[prev_cols]

        # save full
        # need pyarrows
        print("saving " + str(len(vcf_file)) + " variants in parquet")
        pq.write_to_dataset(
            pa.Table.from_pandas(vcf_file), root_path=sample_name + "-maf-full.parquet"
        )

        # save maf
        print("saving maf")
        if i == 0:
            to_maf(
                vcf_file,
                sample_name,
                only_somatic=True,
                only_coding=True,
                whitelist=whitelist,
                drop_multi=True,
                drop_clustered_events=drop_clustered_events,
                tokeep={**TOKEEP_BASE, **TOKEEP_ADD},
                index=True,
                version=version,
            )
        else:
            to_maf(
                vcf_file,
                sample_name,
                only_somatic=True,
                only_coding=True,
                whitelist=whitelist,
                drop_multi=True,
                drop_clustered_events=drop_clustered_events,
                mode="a",
                header=False,
                tokeep={**TOKEEP_BASE, **TOKEEP_ADD},
                index=True,
            )
        del vcf_file
        if tobreak:
            break

    print(f"finished, processed {processed_count} rows")


REPLACE_EMPTY = {
    np.nan: "",
    None: "",
    "Unknown": "",
    "unknown": "",
    "UNKNOWN": "",
    ",": "",
    ",,": "",
    ",,,": "",
    ",%3B,": "",
    ",;,": "",
    ",%2C,": "",
    ",,,,": "",
    ",,,,,": "",
    ",,,,,,": "",
    ",;,,;,": "",
    ",%3B,,%3B,": "",
    ",%2C,,%2C,": "",
    ",,,,,,,": "",
    ",,,,,,,,": "",
    ",,,,,,,,,": "",
    ",;,,;,,;,": "",
    ",%3B,,%3B,,%3B,": "",
    ",%2C,,%2C,,%2C,": "",
    "|": "",
    "||": "",
    "|||": "",
    "||||": "",
    "|||||": "",
    "||||||": "",
    "Approved": "Y",
    "yes": "Y",
    "Yes": "Y",
    "YES": "Y",
    "y": "Y",
    "True": "Y",
    "False": "",
    True: "Y",
    False: "",
    "true": "Y",
    "false": "",
}


REPLACE_SPECIAL_CHAR = {
    "%2C": ",",
    "%3B": ";",
    "%3A": ":",
    "%3D": "=",
    "%25": "%",
}

TO_RENAME_BASE = {
    "AF": "af",
    "AD": "ad",
    "DP": "dp",
    "GT": "gt",
    "PS": "ps",
    "vep_consequence": "variant_info",
    "vep_variant_class": "variant_type",
    "vep_hgvsc": "dna_change",
    "vep_symbol": "hugo_symbol",
    "vep_hgvsp": "protein_change",
    "vep_gene": "ensembl_gene_id",
    "vep_exon": "exon",
    "vep_intron": "intron",
    "vep_uniprot_isoform": "uniprot_id",
    "hgnc_approved_name": "hgnc_name",
    "hgnc_gene_family_name": "hgnc_family",
    # "gencode_34_transcriptexon": "transcript_exon",
    # "gencode_34_transcriptstrand": "transcript_strand",
    # "simple_uniprot_go_biological_process": "go_biological_process",
    # "simple_uniprot_go_cellular_component": "go_cellular_component",
    # "simple_uniprot_go_molecular_function": "go_molecular_function",
    # "cgc_tissue_type": "lineage_association",
    # "cgc_cancer_molecular_genetics": "cancer_molecular_genetics",
    "gencode_34_gccontent": "gc_content",
}

TO_RENAME_HGVS = {
    "vep_sift": "sift",
    "vep_polyphen": "polyphen",
    "vep_gnomade_af": "gnomade_af",
    "vep_gnomadg_af": "gnomadg_af",
    "vep_am_class": "am_class",
    "vep_am_pathogenicity": "am_pathogenicity",
    "vep_feature": "ensembl_feature_id",
    "rs": "dbsnp_rs_id",
    "mc": "molecular_consequence",
}


TO_RENAME_OC = {
    "oc_oncokb__mutation_effect": "oncokb_effect",
    "oc_oncokb__hotspot": "oncokb_hotspot",
    "oc_oncokb__oncogenic": "oncokb_oncogenic",
    ## same
    "oc_civic__description": "civic_description",
    "oc_civic__clinical_a_score": "civic_score",
    "oc_civic__id": "civic_id",
    "oc_pharmgkb__id": "pharmgkb_id",
    "oc_pharmgkb__chemicals": "pharmgkb_chemicals",
    "oc_pharmgkb__pheno_cat": "pharmgkb_pheno_cat",
    "oc_dida__id": "dida_id",
    "oc_dida__name": "dida_name",
    "oc_dida__effect": "dida_effect",
    "oc_dida__relation": "dida_relation",
    "oc_dida__fam": "dida_fam",
    "oc_dida__funct": "dida_funct",
    "oc_dida__dist": "dida_dist",
    "oc_dida__pub": "dida_pub",
    "oc_gwas_catalog__disease": "gwas_disease",
    "oc_gwas_catalog__pmid": "gwas_pmid",
    "oc_ccre_screen___group": "encode_group",
    "oc_ccre_screen__bound": "encode_bound",
    "oc_brca1_func_assay__score": "brca1_func_score",
    "oc_revel__score": "revel_score",
    "oc_spliceai__ds_ag": "spliceai_ds_ag",
    "oc_spliceai__ds_al": "spliceai_ds_al",
    "oc_spliceai__ds_dg": "spliceai_ds_dg",
    "oc_spliceai__ds_dl": "spliceai_ds_dl",
    "oc_spliceai__dp_ag": "spliceai_dp_ag",
    "oc_spliceai__dp_al": "spliceai_dp_al",
    "oc_spliceai__dp_dg": "spliceai_dp_dg",
    "oc_spliceai__dp_dl": "spliceai_dp_dl",
    "oc_gtex__gtex_gene": "gtex_gene",
    "oc_hess_drivers__is_driver": "hess_driver",
    "oc_hess_drivers__signature": "hess_signture",
    "oc_cosmic_sig__cosmic_tier": "cosmic_tier",
    "oc_provean__prediction": "provean_prediction",
}

TOKEEP_BASE = {
    "chrom": "str",
    "pos": "int",
    "ref": "str",
    "alt": "str",
    "af": "float",
    "dp": "int",
    "ref_count": "int",
    "alt_count": "int",
    "gt": "str",
    "ps": "str",
    "variant_type": "str",
    "variant_info": "str",
    "dna_change": "str",
    "protein_change": "str",
    "hugo_symbol": "str",
    "exon": "str",
    "intron": "str",
    "ensembl_gene_id": "str",
    "ensembl_feature_id": "str",
    "hgnc_name": "str",
    "hgnc_family": "str",
    "uniprot_id": "str",
    "dbsnp_rs_id": "str",
    "gc_content": "float",
}

TOKEEP_ADD = {
    "lof_gene_name": "str",
    "lof_gene_id": "str",
    "lof_number_of_transcripts_in_gene": "str",
    "lof_percent_of_transcripts_affected": "str",
    "nmd": "str",
    "clnsig": "str",
    "molecular_consequence": "str",
    "af_exac": "str",
    "af_esp": "str",
    "af_tgp": "str",
    "vep_impact": "str",
    "vep_biotype": "str",
    "vep_hgnc_id": "str",
    "vep_existing_variation": "str",
    "vep_mane_select": "str",
    "vep_ensp": "str",
    "vep_swissprot": "str",
    "sift": "str",
    "polyphen": "str",
    "gnomade_af": "float",
    "gnomadg_af": "float",
    "am_class": "str",
    "am_pathogenicity": "str",
    "vep_clin_sig": "str",
    "vep_somatic": "str",
    "vep_pli_gene_value": "str",
    "vep_loftool": "str",
    "oncogene_high_impact": "str",
    "tumor_suppressor_high_impact": "str",
    ###############
    "achilles_top_genes": "str",
    "structural_relation": "str",
    "associated_with": "str",
    "lof": "str",
    "driver": "str",
    "likely_driver": "str",
    "transcript_likely_lof": "str",
    "brca1_func_score": "str",
    "civic_id": "str",
    "civic_description": "str",
    "civic_score": "str",
    "popaf": "str",
    "likely_gof": "str",
    "likely_lof": "str",
    "hess_driver": "str",
    "hess_signture": "str",
    "revel_score": "str",
    "pharmgkb_id": "str",
    "dida_id": "str",
    "dida_name": "str",
    "gwas_disease": "str",
    "gwas_pmid": "str",
    "gtex_gene": "str",
    "cosmic_tier": "str",
    "oncokb_effect": "str",
    "oncokb_hotspot": "str",
    "oncokb_oncogenic": "str",
    "provean_prediction": "str",
    "spliceai_ds_ag": "str",
    "spliceai_ds_al": "str",
    "spliceai_ds_dg": "str",
    "spliceai_ds_dl": "str",
    "spliceai_dp_ag": "str",
    "spliceai_dp_al": "str",
    "spliceai_dp_dg": "str",
    "spliceai_dp_dl": "str",
    "segdup": "str",
    "rm": "str",
}


def improve(
    vcf,
    force_list=["oc_genehancer__feature_name"],
    torename={**TO_RENAME_BASE, **TO_RENAME_HGVS, **TO_RENAME_OC},
    special_rep=REPLACE_SPECIAL_CHAR,
    replace_empty=REPLACE_EMPTY,
    split_multiallelic=False,
    oncogene_list=[],
    tumor_suppressor_list=[],
    civic_df=None,
):
    """
    given a dataframe representing vcf annotated with opencravat, will improve it.


    Args:
    -----
        revel, spliceai, gtex, funseq2, pharmgkb, dida, gwas_catalog]
        vcf: a df from geney.mutations.vcf_to_df(): the input vcf annotated with opencravat
        force_list: list of elements we know have ',' separated values.
        torename: dict(str: str) renaming dict for columns
        special_rep: dict(str: str) special characters to replace with something else
        replace_empty: dict(str: str) values representing empty fields to replace with something else
        split_multiallelic: bool if True, will split multiallelic variants into separate rows
        min_count_hotspot: int minimum number of mutations in cosmic to consider the loci a hotspot
        civic_df (pd.DataFrame): dataframe containing civic annotations

    Returns:
        the imrpoved vcf
    """

    # mutect2 generates 2 DP columns, one in INFO and the other in FORMAT
    # according to mutect2, INFO fields are for the variant as a whole (over all samples),
    # while FORMAT fields are for individual samples.
    # https://github.com/broadinstitute/gatk/issues/6067
    # here we are dropping the one in INFO
    if vcf[["DP"]].shape[1] > 1:
        vcf["DP_keep"] = vcf["DP"].iloc[:, 1]
        vcf = vcf.drop(columns=["DP"]).rename(columns={"DP_keep": "DP"})

    # solving the special characters
    for val in vcf.columns[9:]:
        try:
            loc = vcf[val].str.contains("%")
        except AttributeError:
            continue
        if loc.sum() > 0:
            loc[loc.isna()] = False
            li = vcf.loc[loc, val]
            # print("replacing: ", val)
            for k, v in special_rep.items():  # "%23": "#",
                li = li.str.replace(k, v)
            vcf.loc[loc, val] = li

    # solving multi allelic sites
    if split_multiallelic:
        for k, row in vcf[vcf.multiallelic].iterrows():
            row.split

    else:
        for val in vcf.columns[9:]:
            try:
                loc = vcf[val].str.contains(",")
            except AttributeError:
                continue
            if loc.sum() > 0:
                # to distinguish between "," from multiallelic and regular ","
                if (
                    vcf[val].str.contains(", ") | vcf[val].str.contains(" ,")
                ).sum() > 0:
                    vcf[val] = vcf[val].str.replace(
                        r"([^\s-]),([^\s-])", r"\1%3B\2", regex=True
                    )
                    loc = vcf[val].str.contains("%3B")
                    spliton = "%3B"
                else:
                    spliton = ","
                loc[loc.isna()] = False
                if loc.sum() == 0:
                    continue
                a, b = np.array(
                    vcf.loc[loc, val].iloc[:1000].str.split(spliton).str[:2].tolist()
                ).T
                if (a == b).sum() == len(a) or val in force_list:
                    print("dropping from: ", val)
                    vcf[val] = vcf[val].str.split(spliton).str[0]

    # replace empty characters:
    print("replacing empty characters:")
    vcf = vcf.replace(replace_empty)
    vcf["oc_brca1_func_assay__score"] = vcf["oc_brca1_func_assay__score"].replace(
        "", np.nan
    )
    vcf["oc_brca1_func_assay__score"] = vcf["oc_brca1_func_assay__score"].astype(
        "float64"
    )

    print("re-annotating CIVIC using static dataframe:")
    if civic_df is not None:
        vcf = civic_df.merge(vcf, on=["chrom", "pos", "ref", "alt"], how="right")

    vcf = vcf.rename(
        columns={
            "description": "oc_civic__description",
            "civic_actionability_score": "oc_civic__clinical_a_score",
            "civic_id": "oc_civic__id",
        }
    )

    print("making new annotations")
    # creating merged annotations

    # ccle_deleterious
    vcf["likely_lof"] = ""
    loc = vcf["vep_impact"] == "HIGH"
    vcf.loc[loc, "likely_lof"] = "Y"

    # structural_relation
    vcf["structural_relation"] = vcf["cgc_translocation_partner"]
    loc = (vcf["cgc_translocation_partner"] == "") & ~(
        vcf["cosmicfusion_fusion_genes"] == ""
    )
    vcf.loc[loc, "structural_relation"] = (
        vcf.loc[loc, "cosmicfusion_fusion_genes"]
        .str.split("_")
        .str[2]
        .str.split("{EN")
        .str[0]
    )

    # DNArepair
    vcf["dna_repair"] = ""
    loc = ~(vcf["dnarepairgenes_activity_linked_to_omim"] == "")
    vcf["dna_repair"] = (
        vcf.loc[loc, "dnarepairgenes_accession_number_linked_to_ncbi_entrez"]
        + ": "
        + vcf.loc[loc, "dnarepairgenes_activity_linked_to_omim"]
    )
    vcf["dna_repair"] = vcf["dna_repair"].replace(
        {": ,%3B,,": "", ": ,%3B,": "", np.nan: ""}
    )
    ############################ ADDITIONNAL ANNOTATOR #########################

    vcf["associated_with"] = ""

    if "oc_oncokb__oncogenic" in vcf.columns.tolist():
        # lof
        vcf["lof"] = ""
        loc = vcf["oc_oncokb__mutation_effect"] == "Loss-of-function"
        vcf.loc[loc, "lof"] = "Y"

        # likely lof
        vcf.loc[loc, "likely_lof"] = "Y"
        loc = vcf["oc_oncokb__mutation_effect"] == "Likely Loss-of-function"  # |
        vcf.loc[loc, "likely_lof"] = "Y"

        # gof
        vcf["gof"] = ""
        loc = vcf["oc_oncokb__mutation_effect"] == "Gain-of-function"
        vcf.loc[loc, "gof"] = "Y"

        # likely_gof
        vcf["likely_gof"] = ""
        vcf.loc[loc, "likely_gof"] = "Y"
        loc = vcf["oc_oncokb__mutation_effect"] == "Likely Gain-of-function"  # |
        vcf.loc[loc, "likely_gof"] = "Y"

    # clinical evidence:
    # https://civic.readthedocs.io/en/latest/model/variants/evidence_score.html
    if "oc_civic__clinical_a_score" in vcf.columns.tolist():
        if "driver" not in vcf.columns.tolist():
            vcf["driver"] = ""
        loc = (~vcf["oc_civic__clinical_a_score"].isnull()) & (
            vcf["multiallelic"] != "Y"
        )
        subvcf = vcf.loc[loc][["oc_civic__clinical_a_score"]]
        vcf.loc[
            subvcf[subvcf["oc_civic__clinical_a_score"].astype(float) >= 8].index,
            "driver",
        ] = "Y"

        if "likely_driver" not in vcf.columns.tolist():
            vcf["likely_driver"] = ""
        loc = ~vcf["oc_civic__clinical_a_score"].isnull()
        vcf.loc[loc, "likely_driver"] = "Y"

    # lof more
    if "oc_brca1_func_assay__score" in vcf.columns.tolist():
        if "lof" not in vcf.columns.tolist():
            vcf["lof"] = ""
        loc = (vcf["oc_brca1_func_assay__score"] != "") & (vcf["multiallelic"] != "Y")
        loc = vcf[loc][
            vcf[loc]["oc_brca1_func_assay__score"].astype(float) <= -1.328
        ].index
        vcf.loc[loc, "lof"] = "Y"

    # lof revel
    if "oc_revel__score" in vcf.columns.tolist():
        # trancript lof
        loc = (vcf["oc_revel__score"] != "") & (vcf["multiallelic"] != "Y")
        vcf["transcript_likely_lof"] = ""
        trscs = []
        for k, val in vcf[loc][["oc_revel__all"]].iterrows():
            trsc = ""
            for v in [i.split(",") for i in val.oc_revel__all[2:-2].split("],[")]:
                if float(v[1]) >= 0.7:
                    trsc += v[0][1:-1] + ";"
            trscs.append(trsc)

        vcf.loc[loc, "transcript_likely_lof"] = trscs

    # generic annotation
    if "oc_spliceai__ds_ag" in vcf.columns.tolist():
        subvcf = vcf[(vcf["oc_spliceai__ds_ag"] != "") & (vcf["multiallelic"] != "Y")]
        loc = subvcf[
            (subvcf["oc_spliceai__ds_ag"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_al"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_dg"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_dl"].astype(float) >= 0.5)
        ].index
        vcf.loc[loc, "associated_with"] += "splicing;"

    loc_e = []
    if "oc_gtex__gtex_gene" in vcf.columns.tolist():
        loc_e += vcf[vcf["oc_gtex__gtex_gene"] != ""].index.tolist()

    vcf.loc[list(set(loc_e)), "associated_with"] += "expression;"

    loc = []
    if "pharmoc_pharmgkb__chemicalsgkb" in vcf.columns.tolist():
        loc += vcf[vcf["oc_pharmgkb__chemicals"] != ""].index.tolist()
    if "oc_oncokb_dm__highestsensitivelevel" in vcf.columns.tolist():
        loc += vcf[vcf["oc_oncokb_dm__highestsensitivelevel"] != ""].index.tolist()
    vcf.loc[list(set(loc)), "associated_with"] += "chemicals;"

    if "likely_gof" not in vcf.columns.tolist():
        vcf["likely_gof"] = ""

    if "oc_hess_drivers__is_driver" in vcf.columns.tolist():
        name = "oc_hess_drivers__is_driver"
        vcf.loc[(vcf[name] == "Y"), "likely_driver"] = "Y"

    if "likely_driver" in vcf.columns.tolist():
        vcf.loc[
            (vcf["likely_driver"] == "Y")
            & vcf["gencode_34_hugosymbol"].isin(tumor_suppressor_list),
            "likely_lof",
        ] = "Y"
        vcf.loc[
            (vcf["likely_driver"] == "Y")
            & vcf["gencode_34_hugosymbol"].isin(oncogene_list),
            "likely_gof",
        ] = "Y"
        vcf.loc[
            (vcf["likely_driver"] == "Y"),  # | (vcf["clinically_significant"] == "Y"),
            "associated_with",
        ] += "cancer;"

    if "likely_gof" in vcf.columns.tolist():
        loc = vcf[vcf["likely_gof"] == "Y"].index.tolist()
        vcf.loc[loc, "associated_with"] += "gene_function_gain;"

    if "likely_lof" in vcf.columns.tolist():
        loc = vcf[vcf["likely_lof"] == "Y"].index.tolist()
    vcf.loc[loc, "associated_with"] += "gene_function_loss;"

    if "oc_dida__id" in vcf.columns.tolist():
        loc_dida = vcf["oc_dida__id"] != ""
        vcf.loc[loc_dida, "associated_with"] += "digenic_disease;"

    # gwas
    if "oc_gwas_catalog__disease" in vcf.columns.tolist():
        loc_g = vcf["oc_gwas_catalog__disease"] != ""
        vcf.loc[loc_g, "associated_with"] += "gwas;"

    vcf.loc[vcf["dna_repair"] != "", "associated_with"] += "dna_repair;"
    vcf.loc[
        vcf["structural_relation"] != "", "associated_with"
    ] += "structural_relation;"

    # high impact oncogenes and tumor suppressor
    if "vep_impact" in vcf.columns.tolist() and "vep_symbol" in vcf.columns.tolist():
        vcf.loc[:, "oncogene_high_impact"] = False
        vcf.loc[:, "tumor_suppressor_high_impact"] = False
        onco_loc = (vcf["vep_impact"] == "HIGH") & (
            vcf["vep_symbol"].isin(oncogene_list)
        )
        ts_loc = (vcf["vep_impact"] == "HIGH") & (
            vcf["vep_symbol"].isin(tumor_suppressor_list)
        )
        vcf.loc[onco_loc, "oncogene_high_impact"] = True
        vcf.loc[ts_loc, "tumor_suppressor_high_impact"] = True

    vcf["vep_gnomade_af"] = vcf["vep_gnomade_af"].astype(str)
    vcf["vep_gnomadg_af"] = vcf["vep_gnomadg_af"].astype(str)

    # rename columns
    vcf = vcf.rename(columns=torename)

    return vcf


def drop_lowqual(
    vcf,
    min_freq=0.15,
    min_depth=2,
):
    loc = (
        # drops 30% of the variants
        (vcf["af"].astype(float) >= min_freq)
        & (vcf["dp"].astype(int) >= min_depth)
        & (~(vcf["chrom"].str.contains("_")))
        # drops 90% of the variants
        & ~(
            (vcf["map_qual"] == "Y")
            | (vcf["slippage"] == "Y")
            | (vcf["strand_bias"] == "Y")
            | (vcf["weak_evidence"] == "Y")
            | (vcf["base_qual"] == "Y")
        )
    )

    vcf = vcf[loc]

    return vcf


def to_maf(
    vcf,
    sample_id,
    tokeep=TOKEEP_BASE,
    whitelist=False,
    drop_multi=True,
    max_pop_af=0.00001,
    only_coding=True,
    only_somatic=True,
    mask_segdup_and_rm=True,
    drop_clustered_events=True,
    version="",
    **kwargs,
):
    """to_maf

    Args:
        vcf (_type_): _description_
        tokeep (_type_, optional): _description_. Defaults to TOKEEP_SMALL.
        whitelist (bool): set it to true to whitelist some mutatios and prevent them
            from being dropped being germline. needs output from vcf.improve and annotators
        sample_id (str): sample id to use for the maf file
        drop_multi (bool): set it to true to drop multi-allelic variants
        min_freq (float): minimum allelic frequency to keep a variant
        min_depth (int): minimum read depth to keep a variant
        max_log_pop_af (float): maximum -log10 population allele frequency to keep a variant
        only_coding (bool): set it to true to keep only coding variants
        only_somatic (bool): set it to true to keep only somatic variants

    """
    # dropping
    initsize = len(vcf)
    if drop_multi:
        #  drops 2% of the variants
        vcf = vcf[(vcf["multiallelic"] is not True) & (vcf["multiallelic"] != "Y")]

    # drop low quality and low coverage
    vcf = drop_lowqual(vcf)

    vcf["gnomade_af"] = vcf["gnomade_af"].replace("", np.nan)
    vcf["gnomadg_af"] = vcf["gnomadg_af"].replace("", np.nan)

    assert ptypes.is_string_dtype(vcf["cosmic_tier"])

    important = pd.Series(False, index=vcf.index)
    loc = pd.Series(True, index=vcf.index)

    if whitelist:
        # if missing columns print issue
        if (
            len(
                set(
                    [
                        "oncokb_effect",
                        "oncokb_oncogenic",
                        "oncokb_hotspot",
                        "cosmic_tier",
                        "brca1_func_score",
                        "clnsig",
                        "civic_score",
                        "hugo_symbol",
                        "hess_driver",
                        "oncogene_high_impact",
                        "tumor_suppressor_high_impact",
                    ]
                )
                - set(vcf.columns)
            )
            > 0
        ):
            print(
                "missing columns to perform whitelisting",
                set(
                    [
                        "oncokb_effect",
                        "oncokb_oncogenic",
                        "oncokb_hotspot",
                        "cosmic_tier",
                        "brca1_func_score",
                        "hugo_symbol",
                        "hess_driver",
                        "oncogene_high_impact",
                        "tumor_suppressor_high_impact",
                    ]
                )
                - set(vcf.columns),
            )
        print("performing whitelisting")
        important = (
            (vcf["oncokb_effect"].isin(["Loss-of-function", "Gain-of-function"]))
            | (vcf["oncokb_oncogenic"] == "Oncogenic")
            | (vcf["oncokb_hotspot"] == "Y")
            | (vcf["cosmic_tier"] == "1")
            | (vcf["brca1_func_score"].astype(float) <= -1.328)
            | (vcf["oncogene_high_impact"])
            | (vcf["tumor_suppressor_high_impact"])
            | (vcf["hess_driver"] == "Y")
            # rescue TERT intronic mutations
            | (
                (vcf["hugo_symbol"] == "TERT")
                & (vcf["pos"] >= 1295054)
                & (vcf["pos"] <= 1295365)
            )
            # rescue MET intron13 PPT mutations
            | (
                (vcf["hugo_symbol"] == "MET")
                & (vcf["pos"] >= 116771825)
                & (vcf["pos"] <= 116771840)
            )
        )
    if only_coding:
        print("only keeping coding mutations")
        vcf["protein_change"] = vcf["protein_change"].fillna("")
        loc = (
            (
                (vcf["variant_info"].str.contains("splice"))
                & (vcf["vep_impact"].isin(["HIGH", "MODERATE"]))
            )
            | (
                (vcf["protein_change"] != "")
                & (vcf["protein_change"].str.endswith("=") == False)
            )
            | important
        )
    if drop_clustered_events:
        print("dropping clustered_events variants")
        loc = ((vcf["clustered_events"] != "Y") | important) & loc
    if mask_segdup_and_rm:
        print("removing variants in segmental duplication and repeatmasker regions")
        loc = (((vcf["segdup"] != "Y") & (vcf["rm"] != "Y")) | important) & loc

    # redefine somatic (not germline or pon and a log pop. af of > max_log_pop_af)
    # drops 80% of the variants
    if only_somatic:
        print("only keeping somatic mutations")
        loc = (
            (
                ~(vcf["pon"] == "Y")
                & (~(vcf["gnomade_af"].astype(float) > max_pop_af))
                & (~(vcf["gnomadg_af"].astype(float) > max_pop_af))
            )
            | important
        ) & loc

    vcf = vcf[loc]

    print(
        "new size: "
        + str(len(vcf))
        + ". removed: {:2f}%".format((1 - (len(vcf) / initsize)) * 100)
    )

    vcf["ref_count"] = None
    vcf["alt_count"] = None

    if len(vcf) > 0:
        # creating count columns
        vcf[["ref_count", "alt_count"]] = vcf["ad"].str.split(",", expand=True)

    # subsetting
    vcf = vcf[list(tokeep.keys())]

    # setting the right type
    for k, v in tokeep.items():
        vcf[k] = vcf[k].astype(v)
        if v == "str":
            vcf[k] = vcf[k].replace(",", "%2C")

    vcf["rescue"] = False
    if len(vcf) > 0:
        vcf.loc[important, "rescue"] = True

    if version != "":
        vcf.index.name = version

    vcf.to_csv(
        sample_id + "-maf-coding_somatic-subset.csv.gz", index_label=version, **kwargs
    )


if __name__ == "__main__":
    main()
