import pdb
import numpy as np

BOOLIFY = [
    "hgnc_status",
    "str",
    "cgc_cancer_germline_mut",
    "cgc_cancer_somatic_mut",
    "pon",
]

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
    ",": "",
    "Unknown": "",
    "UNKNOWN": "",
    ",,": "",
    ",,,": "",
    ",,,,": "",
    ",,,,,": "",
    ",,,,,,": "",
    "|": "",
    "||": "",
    "|||": "",
    "||||": "",
    "|||||": "",
    "||||||": "",
}

TOKEEP_ONCO = {
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
    "is_coding": "bool",
    "hugo_symbol": "str",
    "hgnc_name": "str",
    "hgnc_family": "str",
    "transcript": "str",
    "transcript_exon": "str",
    "protein_change": "str",
    "transcript_strand": "str",
    "uniprot_id": "str",
    "go_biological_process": "str",
    "go_cellular_component": "str",
    "go_molecular_function": "str",
    "feature_name": "str",
    "str": "bool",
    "lineage_association": "str",
    "disease_association": "str",
    "chemical_association": "str",
    "achilles_top_genes": "str",
    "cancer_molecular_genetics": "str",
    "hotspot": "bool",
    "dbsnp_id": "str",
    "dbsnp_filter": "str",
    "issues": "str",
    "quality_score": "float",
    "somatic_score": "float",
    "gc_content": "float",
    "dna_repair": "str",
    "splice_ai_info": "",
    "ccle_deleterious": "bool",
    "deleteriousness_score": "float",
    "functionality_score": "float",
    "effect_score": "float",
    "transcript_specific_pathogenecity": "",
    "structural_relation": "str",
}


TOKEEP_SMALL = {
    "chrom": "str",
    "pos": "int",
    "ref": "str",
    "alt": "str",
    "af": "float",
    "ref_count": "int",
    "alt_count": "int",
    "gt": "str",
    "ps": "str",
    "dna_change": "str",
    "protein_change": "str",
    "variant_info": "str",
    "hugo_symbol": "str",
    "hgnc_name": "str",
    "hgnc_family": "str",
    "transcript_exon": "str",
    "transcript_strand": "str",
    "uniprot_id": "str",
    "str": "bool",
    "lineage_association": "str",
    "achilles_top_genes": "str",
    "cancer_molecular_genetics": "str",
    "hotspot": "bool",
    "dbsnp_id": "str",
    "dbsnp_filter": "str",
    "issues": "str",
    "gc_content": "float",
    "dna_repair": "str",
    "ccle_deleterious": "bool",
    "structural_relation": "str",
}

TO_RENAME_OC = {
    "AF": "af",
    "AD": "ad",
    "GT": "gt",
    "PS": "ps",
    "oc_base__so": "variant_info",
    "oc_base__cchange": "dna_change",
    "oc_base__coding": "is_coding",
    "oc_base__hugo": "hugo_symbol",
    "oc_base__achange": "protein_change",
    "oc_provean__uniprot": "uniprot_id",
    "oc_base__transcript": "transcript",
    "hgnc_approved_name": "hgnc_name",
    "hgnc_gene_family_name": "hgnc_family",
    "gencode_34_transcriptexon": "transcript_exon",
    "gencode_34_transcriptstrand": "transcript_strand",
    "simple_uniprot_go_biological_process": "go_biological_process",
    "simple_uniprot_go_cellular_component": "go_cellular_component",
    "simple_uniprot_go_molecular_function": "go_molecular_function",
    "cgc_tissue_type": "lineage_association",
    "cgc_cancer_molecular_genetics": "cancer_molecular_genetics",
    "somatic_score": "somatic_score",
    "gencode_34_gccontent": "gc_content",
}


TO_RENAME = {
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
    "somatic_score": "somatic_score",
    "gencode_34_gccontent": "gc_content",
}

REPLACE_SPECIAL_CHAR = {
    "%2C": ";",
    "%3B": ",",
    "%3A": ":",
    "%3D": "=",
    "%25": "%",
}


def improve(
    vcf,
    force_list=["oc_genehancer__feature_name"],
    torename=TO_RENAME,
    special_rep=REPLACE_SPECIAL_CHAR,
    replace_empty=REPLACE_EMPTY,
    boolify=BOOLIFY,
    with_onco_kb=False,
    split_multiallelic=False,
):
    """
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

    # creating merged annotations

    # issue
    to_add = []
    for val in vcf[["dbsnp_cfl", "dbsnp_asp"]].values:
        if val[0] == "true":
            to_add.append("assembly_conflict")
        elif val[1] == "true":
            to_add.append("as_specific")
        else:
            to_add.append("")
    vcf["issues"] = to_add
    vcf.drop(columns=["dbsnp_asp", "dbsnp_cfl"], inplace=True)

    # defining hotspot
    vcf["hotspot"] = False
    loc = ~(
        vcf["cgc_cancer_somatic_mut"] == ""
    )  # | ~(vcf["cgc_cancer_germline_mut"] == "")
    if with_onco_kb:
        loc = loc | (vcf["oc_oncokb__oncogenic"] == "Likely Oncogenic")
    vcf.loc[loc, "hotspot"] = True

    # ccle_deleterious
    vcf["ccle_deleterious"] = False
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
    ] = True

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

    vcf.drop(
        columns=[
            "dnarepairgenes_accession_number_linked_to_ncbi_entrez",
            "dnarepairgenes_activity_linked_to_omim",
        ],
        inplace=True,
    )

    if not with_onco_kb:
        vcf["is_coding"] = vcf["gencode_34_proteinchange"] != ""

    for val in boolify:
        vcf[val] = (
            vcf[val]
            .astype(str)
            .replace({"": False, "Approved": True, "yes": True})
            .astype(bool)
        )

    # somatic_score

    # quality_score

    # deleteriousness_score

    # transcript_specific_pathogenecity

    # functionality_score

    # effect_score

    # disease_association

    # chemical_association

    # splice_AI_info

    # rename columns
    vcf = vcf.rename(columns=torename)

    # setting the right types
    # vcf["is_coding"].replace({"Yes": True, "": False})
    # ...

    return vcf


def to_maf(
    vcf,
    sample_id,
    tokeep=TOKEEP_SMALL,
    drop_multi=True,
    min_freq=0.15,
    min_depth=2,
    max_log_pop_af=3,
    only_coding=True,
    only_somatic=True,
    **kwargs
):
    """to_maf 

    Args:
        vcf (_type_): _description_
        tokeep (_type_, optional): _description_. Defaults to TOKEEP_SMALL.
    """
    # dropping
    if drop_multi:
        #  drops 2% of the variants
        vcf = vcf[~vcf["multiallelic"]]

    loc = (
        # drops 30% of the variants
        (vcf["af"].astype(float) >= min_freq)
        & (vcf["ad"].str.split(",").str[1].astype(int) >= min_depth)
        # drops 90% of the variants
        & ~(
            vcf["map_qual"]
            | vcf["slippage"]
            | vcf["strand_bias"]
            | vcf["weak_evidence"]
            | vcf["clustered_events"]
            | vcf["base_qual"]
            # | vcf["fragments"]
        )
    )
    if only_coding:
        # drops 99.5% of the variants
        loc = loc & (vcf["is_coding"] | (vcf["variant_info"] == "SPLICE_SITE"))

    # we will drop 99.993% of the variants and 90% of the columns
    vcf = vcf[loc]

    # redefine somatic
    # drops 80% of the variants
    if only_somatic:
        vcf = vcf[
            (
                ~(vcf["germline"] | vcf["pon"])
                & (vcf["popaf"].astype(float) > max_log_pop_af)
            )
            | vcf["hotspot"]
            | vcf["ccle_deleterious"]
        ]
    print(len(vcf))
    # creating count columns
    vcf[["ref_count", "alt_count"]] = np.array(vcf["ad"].str.split(",").to_list())

    # subsettingnz
    vcf = vcf[list(tokeep.keys())]

    # setting the right type
    for k, v in tokeep.items():
        vcf[k] = vcf[k].astype(v)
        if v == "str":
            vcf[k] = vcf[k].replace(",", "%3B")
    vcf.to_csv(sample_id + "-maf-coding_somatic-subset.csv.gz", **kwargs)

