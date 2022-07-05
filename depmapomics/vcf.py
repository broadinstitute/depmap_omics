import re
import numpy as np

DROPCOMMA = [
    "gencode_34_hugosymbol",
    "gencode_34_chromosome",
    "gencode_34_secondaryvariantclassification",
    "gencode_34_annotationtranscript",
    "gencode_34_transcriptstrand",
    "gencode_34_transcriptexon",
    "gencode_34_transcriptpos",
    "achilles_top_genes",
    "cgc_geneid",
    "cgc_chr",
    "cgc_chr_band",
    "cgc_cancer_somatic_mut",
    "cgc_cancer_germline_mut",
    "cgc_cancer_syndrome",
    "cgc_tissue_type",
    "cgc_cancer_molecular_genetics",
    "cgc_translocation_partner",
    "cgc_other_germline_mut",
    "cgc_other_syndrome/disease",
    "clinvar_vcf_af_esp",
    "clinvar_vcf_af_exac",
    "clinvar_vcf_af_tgp",
    "clinvar_vcf_alleleid",
    "clinvar_vcf_clndisdbincl",
    "clinvar_vcf_clndnincl",
    "clinvar_vcf_clnhgvs",
    "clinvar_vcf_clnsig",
    "clinvar_vcf_clnsigincl",
    "clinvar_vcf_clnvc",
    "clinvar_vcf_clnvcso",
    "clinvar_vcf_geneinfo",
    "clinvar_vcf_origin",
    "clinvar_vcf_rs",
    "clinvar_vcf_ssr",
    "clinvar_vcf_id",
    "cosmic_overlapping_mutations",
    "cosmicfusion_fusion_genes",
    "cosmictissue_total_alterations_in_gene",
    "cosmictissue_tissue_types_affected",
    "dnarepairgenes_chromosome_location_linked_to_ncbi_mapview",
    "dnarepairgenes_accession_number_linked_to_ncbi_entrez",
    "familial_cancer_genes_syndrome",
    "familial_cancer_genes_reference",
    "gencode_xhgnc_hgnc_id",
    "gencode_xrefseq_mrna_id",
    "gencode_xrefseq_prot_acc",
    "hgnc_hgnc_id",
    "hgnc_status",
    "hgnc_locus_group",
    "hgnc_chromosome",
    "hgnc_date_modified",
    "hgnc_date_symbol_changed",
    "hgnc_date_name_changed",
    "hgnc_enzyme_ids",
    "hgnc_entrez_gene_id",
    "hgnc_ensembl_gene_id",
    "hgnc_refseq_ids",
    "hgnc_gene_family_id",
    "hgnc_vega_id",
    "hgnc_entrez_gene_id(supplied_by_ncbi)",
    "hgnc_omim_id(supplied_by_omim)",
    "hgnc_refseq(supplied_by_ncbi)",
    "hgnc_uniprot_id(supplied_by_uniprot)",
    "hgnc_ensembl_id(supplied_by_ensembl)",
    "hgnc_ucsc_id(supplied_by_ucsc)",
    "oreganno_build",
    "oreganno_id",
    "oreganno_values",
    "simple_uniprot_uniprot_entry_name",
    "simple_uniprot_drugbank",
    "simple_uniprot_alt_uniprot_accessions",
    "simple_uniprot_uniprot_accession",
    "oc_oncokb__oncogenic",
    "oc_original_input__chrom",
    "oc_oncokb__knowneffect",
]

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
    "GT": "gt",
    "PS": "ps",
    "gencode_34_variantclassification": "variant_info",
    "gencode_34_cdnachange": "dna_change",
    "gencode_34_hugosymbol": "hugo_symbol",
    "gencode_34_proteinchange": "protein_change",
    "hgnc_uniprot_id(supplied_by_uniprot)": "uniprot_id",
    "gencode_34_annotationtranscript": "transcript",
    "hgnc_approved_name": "hgnc_name",
    "hgnc_gene_family_name": "hgnc_family",
    "gencode_34_transcriptexon": "transcript_exon",
    "gencode_34_transcriptstrand": "transcript_strand",
    "simple_uniprot_go_biological_process": "go_biological_process",
    "simple_uniprot_go_cellular_component": "go_cellular_component",
    "simple_uniprot_go_molecular_function": "go_molecular_function",
    "cgc_tissue_type": "lineage_association",
    "cgc_cancer_molecular_genetics": "cancer_molecular_genetics",
    "gencode_34_gccontent": "gc_content",
}


TO_RENAME_OC = {
    "oc_oncokb_dm__highestsensitivelevel": "oncokb_sensitivelevel",
    "oc_oncokb_dm__highestresistancelevel": "oncokb_resistancelevel",
    "oc_oncokb_dm__highestdiagnosticimplicationlevel": "oncokb_diagnosticlevel",
    "oc_oncokb_dm__highestprognosticimplicationlevel": "oncokb_prognosticlevel",
    "oc_oncokb_dm__variantsummary": "oncokb_variant_summary",
    "oc_oncokb_dm__pmids": "oncokb_pmids",
    "oc_civic__description": "civic_description",
    "oc_civic__clinical_a_score": "civic_score",
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
    "oc_cscape__score": "cscape_score",
    "oc_brca1_func_assay__score": "brca1_func_score",
    "oc_dann__score": "dann_score",
    "oc_revel__score": "revel_score",
    "oc_spliceai__ds_ag": "spliceai_ds_ag",
    "oc_spliceai__ds_al": "spliceai_ds_al",
    "oc_spliceai__ds_dg": "spliceai_ds_dg",
    "oc_spliceai__ds_dl": "spliceai_ds_dl",
    "oc_gtex__gtex_gene": "gtex_gene",
    "oc_funseq2__score": "funseq2_score",
    "oc_funseq2__motif": "funseq2_motif",
    "oc_funseq2__hot": "funseq2_hot",
}

TOKEEP_BASE = {
    "chrom": "str",
    "pos": "int",
    "ref": "str",
    "alt": "str",
    "af": "float",
    "ref_count": "int",
    "alt_count": "int",
    "gt": "str",
    "ps": "str",
    "variant_info": "str",
    "dna_change": "str",
    "protein_change": "str",
    "hugo_symbol": "str",
    "hgnc_name": "str",
    "hgnc_family": "str",
    "transcript": "str",
    "transcript_exon": "str",
    "transcript_strand": "str",
    "uniprot_id": "str",
    "str": "str",
    "dbsnp_id": "str",
    "dbsnp_filter": "str",
    "issues": "str",
    "gc_content": "float",
}

TOKEEP_ADD = {
    ###############
    "lineage_association": "str",
    "achilles_top_genes": "str",
    "cancer_molecular_genetics": "str",
    "hotspot": "str",
    "ccle_deleterious": "str",
    "structural_relation": "str",
    "cosmic_hotspot": "str",
    "associated_with": "str",
    "clinically_significant": "str",
    "gof": "str",
    "lof": "str",
    "driver": "str",
    "likely_driver": "str",
    "pathogenic": "str",
    "trancript_likely_pathogenic": "str",
    "transcript_likely_lof": "str",
    "oncokb_sensitivelevel": "str",
    "oncokb_resistancelevel": "str",
    "oncokb_diagnosticlevel": "str",
    "oncokb_prognosticlevel": "str",
    "oncokb_pmids": "str",
}


TOKEEP_LARGE_ADD = {
    "go_biological_process": "str",
    "go_cellular_component": "str",
    "go_molecular_function": "str",
    ###############
    # "quality_score": "float",
    # "somatic_score": "float",
    "oncokb_variant_summary": "str",
    "civic_description": "str",
    "civic_score": "str",
    "pharmgkb_id": "str",
    "pharmgkb_chemicals": "str",
    "pharmgkb_pheno_cat": "str",
    "dida_id": "str",
    "dida_name": "str",
    "dida_effect": "str",
    "dida_relation": "str",
    "dida_fam": "str",
    "dida_funct": "str",
    "dida_dist": "str",
    "dida_pub": "str",
    "gwas_disease": "str",
    "gwas_pmid": "str",
    "encode_group": "str",
    "encode_bound": "str",
    "cscape_score": "str",
    "brca1_func_score": "str",
    "dann_score": "str",
    "revel_score": "str",
    "spliceai_ds_ag": "str",
    "spliceai_ds_al": "str",
    "spliceai_ds_dg": "str",
    "spliceai_ds_dl": "str",
    "gtex_gene": "str",
    "funseq2_score": "str",
    "funseq2_motif": "str",
    "funseq2_hot": "str",
    "likely_gof": "str",
    "likely_lof": "str",
    "likely_pathogenic": "str",
}


ALL_ANNOTATORS = [
    "oncokb",
    "cscape",
    "civic",
    "brca1_func",
    "sift",
    "provean",
    "dann",
    "revel",
    "spliceai",
    "gtex",
    "funseq2",
    "pharmgkb",
    "dida",
    "gwas_catalog",
]


def improve(
    vcf,
    force_list=["oc_genehancer__feature_name"],
    torename={**TO_RENAME_BASE, **TO_RENAME_OC},
    special_rep=REPLACE_SPECIAL_CHAR,
    replace_empty=REPLACE_EMPTY,
    annotators=[],
    split_multiallelic=False,
    min_count_hotspot=5,
    **kwargs
):
    """
    annotators: list(str) many of: [oncokb, cscape, civic, brca1_func, sift, provean, dann, 
    revel, spliceai, gtex, funseq2, pharmgkb, dida, gwas_catalog]
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
            print("replacing: ", val)
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

    print("making new annotations")
    # creating merged annotations

    # issue
    to_add = []
    for val in vcf[["dbsnp_cfl", "dbsnp_asp"]].values:
        if val[0] == "Y":
            to_add.append("assembly_conflict")
        elif val[1] == "Y":
            to_add.append("as_specific")
        else:
            to_add.append("")
    vcf["issues"] = to_add
    todrop = ["dbsnp_asp", "dbsnp_cfl"]

    # defining hotspot
    vcf["cosmic_hotspot"] = ""

    hotspot_l = []
    for cosmic_over in list(
        set(
            vcf[vcf["cosmic_overlapping_mutations"] != ""][
                "cosmic_overlapping_mutations"
            ]
        )
    ):
        # finding the number in " p.I517T(13)", " p.I517T(13), p.G202G(15), p.?(56)"
        res = sum(
            [
                int(val.group(0)[1:-1])
                for val in re.finditer(r"([(])\d+([)])", cosmic_over)
            ]
        )
        if res > min_count_hotspot:
            hotspot_l.append(cosmic_over)
    loc = vcf["cosmic_overlapping_mutations"].isin(hotspot_l)
    vcf["cosmic_hotspot"] = "Y"

    # ccle_deleterious
    vcf["ccle_deleterious"] = ""
    vcf.loc[
        vcf["gencode_34_variantclassification"].isin(
            [
                "DE_NOVO_START_OUT_FRAME",
                "DE_NOVO_START_IN_FRAME",
                "FRAME_SHIFT_DEL",
                "FRAME_SHIFT_INS",
                "START_CODON_INS",
                "START_CODON_DEL",
                "NONSTOP",
                "NONSENSE",
            ]
        ),
        "ccle_deleterious",
    ] = "Y"

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

    todrop.extend(
        [
            "dnarepairgenes_accession_number_linked_to_ncbi_entrez",
            "dnarepairgenes_activity_linked_to_omim",
        ]
    )
    ############################ ADDITIONNAL ANNOTATOR #########################

    vcf["associated_with"] = ""

    if "oncokb" in annotators:
        vcf["hotspot"] = ""
        loc = vcf["oc_oncokb_dm__hotspot"] == "Y"
        vcf.loc[loc, "hotspot"] = "Y"
        todrop.append("oc_oncokb_dm__hotspot")

        # defining drivers
        vcf["driver"] = ""
        loc = vcf["oc_oncokb_dm__oncogenic"] == "Oncogenic"
        vcf.loc[loc, "driver"] = "Y"

        vcf["likely_driver"] = ""
        vcf.loc[loc, "likely_driver"] = "Y"
        loc = vcf["oc_oncokb_dm__oncogenic"] == "Likely Oncogenic"
        vcf.loc[loc, "likely_driver"] = "Y"

        todrop.append("oc_oncokb_dm__oncogenic")

        # clinical significance
        vcf["clinically_significant"] = ""
        loc = ~(
            (vcf["oc_oncokb_dm__highestprognosticimplicationlevel"] == "")
            & (vcf["oc_oncokb_dm__highestdiagnosticimplicationlevel"] == "")
            & (vcf["oc_oncokb_dm__highestsensitivelevel"] == "")
            & (vcf["oc_oncokb_dm__highestresistancelevel"] == "")
        )
        vcf.loc[loc, "clinically_significant"] = "Y"

        # lof
        vcf["lof"] = ""
        loc = vcf["oc_oncokb_dm__knowneffect"] == "Loss-of-function"
        vcf.loc[loc, "lof"] = "Y"

        # likely lof
        vcf["likely_lof"] = ""
        vcf.loc[loc, "likely_lof"] = "Y"
        loc = vcf["oc_oncokb_dm__knowneffect"] == "Likely Loss-of-function"  # |
        vcf.loc[loc, "likely_lof"] = "Y"

        # gof
        vcf["gof"] = ""
        loc = vcf["oc_oncokb_dm__knowneffect"] == "Gain-of-function"
        vcf["gof"] = "Y"

        # likely_gof
        vcf["likely_gof"] = ""
        vcf.loc[loc, "likely_gof"] = "Y"
        loc = vcf["oc_oncokb_dm__knowneffect"] == "Likely Gain-of-function"  # |
        vcf.loc[loc, "likely_gof"] = "Y"

    # additional defining drivers
    if "cscape" in annotators:
        if "likely_driver" not in vcf.columns.tolist():
            vcf["likely_driver"] = ""
        subvcf = vcf[(vcf["oc_cscape__score"] != "") & (vcf["multiallelic"] != "Y")][
            ["oc_base__coding", "oc_cscape__score"]
        ]
        vcf.loc[
            subvcf[
                (
                    (subvcf["oc_base__coding"] != "")
                    & (subvcf["oc_cscape__score"].astype(float) >= 0.89)
                )
                | (
                    (subvcf["oc_base__coding"] == "")
                    & (subvcf["oc_cscape__score"].astype(float) >= 0.80)
                )
            ].index,
            "likely_driver",
        ] = "Y"

    # clinical significance
    if "civic" in annotators:
        if "clinically_significant" not in vcf.columns.tolist():
            vcf["clinically_significant"] = ""
        loc = vcf["oc_civic__clinical_a_score"] != ""
        vcf.loc[loc, "clinically_significant"] = "Y"

    # lof more
    if "brca1_func" in annotators:
        if "lof" not in vcf.columns.tolist():
            vcf["lof"] = ""
        loc = (vcf["oc_brca1_func_assay__score"] != "") & (vcf["multiallelic"] != "Y")
        loc = vcf[loc][
            vcf[loc]["oc_brca1_func_assay__score"].astype(float) <= -1.328
        ].index
        vcf.loc[loc, "lof"] = "Y"

    # transcript lof
    if "sift" in annotators or "provean" in annotators:
        if "provean":
            name = "oc_provean__all"
            score = -2.5
        else:  # 'sift' in annotators:
            name = "oc_sift__all"
            score = 0.05
        vcf["transcript_likely_lof"] = ""
        trscs = []
        loc = vcf[name] != ""
        for k, val in vcf[loc][[name]].iterrows():
            trsc = ""
            for trsc in [i.split(",") for i in val[name][2:-2].split("],[")]:
                if float(trsc[2]) <= score:
                    trsc += trsc[0][1:-1] + ";"
            trscs.append(trsc)
        vcf.loc[loc, "transcript_likely_lof"] = trscs

    # likely lof in DANN
    if "dann" in annotators:
        if "likely_lof" not in vcf.columns.tolist():
            vcf["likely_lof"] = ""
        # http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
        loc = (vcf["oc_dann__score"] != "") & (vcf["multiallelic"] != "Y")
        loc = vcf[loc][vcf[loc]["oc_dann__score"].astype(float) >= 0.96].index
        vcf.loc[loc, "likely_lof"] = "Y"

    # pathogenic
    if "revel" in annotators:
        vcf["pathogenic"] = ""
        loc = (vcf["oc_revel__score"] != "") & (vcf["multiallelic"] != "Y")
        vcf.loc[
            vcf[loc][vcf[loc]["oc_revel__score"].astype(float) >= 0.9].index,
            "pathogenic",
        ] = "Y"

        # likely pathogenic
        vcf["likely_pathogenic"] = ""
        vcf.loc[
            vcf[loc][vcf[loc]["oc_revel__score"].astype(float) >= 0.7].index,
            "likely_pathogenic",
        ] = "Y"
        # trancript pathogenic
        vcf["trancript_likely_pathogenic"] = ""
        trscs = []
        loc = vcf["oc_revel__all"] != ""
        for k, val in vcf[loc][["oc_revel__all"]].iterrows():
            trsc = ""
            for trsc in [i.split(",") for i in val.oc_revel__all[2:-2].split("],[")]:
                if float(trsc[1]) >= 0.7:
                    trsc += trsc[0][1:-1] + ";"
            trscs.append(trsc)

        vcf.loc[loc, "trancript_likely_pathogenic"] = trscs

        loc = vcf["likely_pathogenic"] == "Y"
        vcf.loc[loc, "associated_with"] += "disease;"

    # generic annotation
    if "spliceai" in annotators:
        subvcf = vcf[(vcf["oc_spliceai__ds_ag"] != "") & (vcf["multiallelic"] != "Y")]
        loc = subvcf[
            (subvcf["oc_spliceai__ds_ag"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_al"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_dg"].astype(float) >= 0.5)
            | (subvcf["oc_spliceai__ds_dl"].astype(float) >= 0.5)
        ].index
        vcf.loc[loc, "associated_with"] += "splicing;"

    loc_e = []
    if "gtex" in annotators:
        loc_e += vcf[vcf["oc_gtex__gtex_gene"] != ""].index.tolist()

    if "funseq2" in annotators:
        loc = (vcf["oc_funseq2__score"] != "") & (vcf["multiallelic"] != "Y")
        loc_e += vcf[loc][
            (vcf[loc].oc_funseq2__score.astype(float) >= 0.5)
            & (
                (vcf[loc]["oc_funseq2__motif"] != "")
                | (vcf[loc]["oc_funseq2__hot"] != "")
            )
        ].index.tolist()

    vcf.loc[list(set(loc_e)), "associated_with"] += "expression;"

    loc = []
    if "pharmgkb" in annotators:
        loc += vcf[vcf["oc_pharmgkb__chemicals"] != ""].index.tolist()
    if "oncokb" in annotators:
        loc += vcf[vcf["oc_oncokb_dm__highestsensitivelevel"] != ""].index.tolist()
    vcf.loc[list(set(loc)), "associated_with"] += "chemicals;"

    if "oncokb" in annotators:
        loc = vcf[vcf["likely_gof"] == "Y"].index.tolist()
        vcf.loc[loc, "associated_with"] += "gene_function_gain;"

    if "likely_lof" in vcf.columns.tolist():
        loc = vcf[vcf["likely_lof"] == "Y"].index.tolist()
    vcf.loc[loc, "associated_with"] += "gene_function_loss;"

    loc = []
    if "likely_driver" in vcf.columns.tolist():
        loc += vcf[(vcf["likely_driver"] == "Y")].index.tolist()
    if "hotspot" in vcf.columns.tolist():
        loc += vcf[
            (vcf["hotspot"] == "Y") | (vcf["clinically_significant"] == "Y")
        ].index.tolist()
    vcf.loc[list(set(loc)), "associated_with"] += "cancer;"

    if "dida" in annotators:
        loc_dida = vcf["oc_dida__id"] != ""
        vcf.loc[loc_dida, "associated_with"] += "digenic_disease;"

    # gwas
    if "gwas_catalog" in annotators:
        loc_g = vcf["oc_gwas_catalog__disease"] != ""
        vcf.loc[loc_g, "associated_with"] += "gwas;"

    vcf.loc[vcf["dna_repair"] != "", "associated_with"] += "dna_repair;"
    vcf.loc[
        vcf["structural_relation"] != "", "associated_with"
    ] += "structural_relation;"

    vcf["is_coding"] = ""
    vcf.loc[vcf["gencode_34_proteinchange"] != "", "is_coding"] = "Y"
    # else it is an oncokb value

    for k, df in kwargs.items():
        missing_col = set([""]) - set(df.columns)
        if len(missing_col) > 0:
            raise ValueError("missing columns", missing_col)

    # somatic_score

    # quality_score

    # rename columns
    vcf = vcf.drop(columns=todrop).rename(columns=torename)

    # setting the right types
    # vcf["is_coding"].replace({"Yes": True, "": False})
    # ...

    return vcf


def to_maf(
    vcf,
    sample_id,
    tokeep=TOKEEP_BASE,
    drop_multi=True,
    min_freq=0.15,
    min_depth=2,
    max_log_pop_af=3,
    only_coding=True,
    only_somatic=True,
    oncogenic_list=[],
    tumor_suppressor_list=[],
    **kwargs
):
    """to_maf 

    Args:
        vcf (_type_): _description_
        tokeep (_type_, optional): _description_. Defaults to TOKEEP_SMALL.
    """
    # dropping
    initsize = len(vcf)
    if drop_multi:
        #  drops 2% of the variants
        vcf = vcf[vcf["multiallelic"] != "Y"]

    loc = (
        # drops 30% of the variants
        (vcf["af"].astype(float) >= min_freq)
        & (vcf["ad"].str.split(",").str[1].astype(int) >= min_depth)
        # drops 90% of the variants
        & ~(
            (vcf["map_qual"] == "Y")
            | (vcf["slippage"] == "Y")
            | (vcf["strand_bias"] == "Y")
            | (vcf["weak_evidence"] == "Y")
            | (vcf["clustered_events"] == "Y")
            | (vcf["base_qual"] == "Y")
            # | vcf["fragments"]
        )
    )

    important = (
        (vcf["driver"] == "Y")
        | (vcf["pathogenic"] == "Y")
        | (vcf["likely_gof"] == "Y")
        | (vcf["clinically_significant"] == "Y")
        | (vcf["lof"] == "Y")
        | (
            (vcf["likely_lof"] == "Y")
            & (vcf["hugo_symbol"].isin(tumor_suppressor_list))
        )
        | (
            (vcf["likely_driver"] == "Y")
            & vcf["hugo_symbol"].isin(oncogenic_list + tumor_suppressor_list)
        )
    )
    if only_coding:
        # drops 99.5% of the variants
        loc = loc & (
            (vcf["is_coding"] == "Y")
            | (vcf["variant_info"] == "SPLICE_SITE")
            | important
        )

    # we will drop 99.993% of the variants and 90% of the columns
    vcf = vcf[loc]

    # redefine somatic (not germline or pon and a log pop. af of > max_log_pop_af)
    # drops 80% of the variants
    if only_somatic:
        loc = (
            ~((vcf["germline"] == "Y") | (vcf["pon"] == "Y"))
            & (vcf["popaf"].astype(float) > max_log_pop_af)
        ) | important
        vcf = vcf[loc]
    print(
        "new size: "
        + str(len(vcf))
        + ". removed: {:2f}%".format((1 - (len(vcf) / initsize)) * 100)
    )
    # creating count columns
    vcf[["ref_count", "alt_count"]] = np.array(vcf["ad"].str.split(",").to_list())

    # subsetting
    vcf = vcf[list(tokeep.keys())]

    # setting the right type
    for k, v in tokeep.items():
        vcf[k] = vcf[k].astype(v)
        if v == "str":
            vcf[k] = vcf[k].replace(",", "%2C")
    vcf.to_csv(sample_id + "-maf-coding_somatic-subset.csv.gz", **kwargs)

