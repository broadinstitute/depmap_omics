"""Mutation postprocessing module."""

import re
from depmapomics import constants
from mgenepy.utils import helper as h
import os
import pandas as pd
from collections import Counter
from itertools import repeat
import multiprocessing
import subprocess
import pandera as pa
from tqdm import tqdm
import numpy as np


def annotateLikelyImmortalized(
    maf,
    sample_col=constants.SAMPLEID,
    genome_change_col="dna_change",
    chrom_col="chrom",
    pos_col="pos",
    hotspotcol="cosmic_hotspot",
    max_recurrence=constants.IMMORTALIZED_THR,
):
    """Annotate the maf file with the likely immortalized mutations

    Based on occurence accross samples

    Args:
        maf (pandas.DataFrame): the maf file with columns: sample_col, genome_change_col, TCGAlocs
        sample_col (str): the column name of the sample id
        genome_change_col (str, optional): the column name of the genome change. Defaults to "Genome_Change".
        TCGAlocs (list, optional): the column names of the counts that would make the mutation non immortalization induced. Defaults to ['TCGAhsCnt', 'COSMIChsCnt'].
        max_recurrence (float, optional): the maximum recurrence rate to call immortalize. Defaults to 0.05.
        min_tcga_true_cancer (int, optional): the minimum number of TCGA true cancer samples to not call immortalize. Defaults to 5.

    Returns:
        pandas.DataFrame: the maf file with the added column: immortalized
    """
    maf["is_likely_immortalization"] = False
    maf["combined_mut"] = (
        maf[chrom_col] + "_" + maf[pos_col].astype(str) + "_" + maf[genome_change_col]
    )
    leng = len(set(maf[sample_col]))
    maf[
        (maf[hotspotcol] != True)
        & (
            maf["combined_mut"].isin(
                [
                    k
                    for k, v in Counter(maf["combined_mut"].tolist()).items()
                    if v > max_recurrence * leng
                ]
            )
        )
    ]["LikelyImmortalized"] = True
    maf = maf.drop(columns=["combined_mut"])
    return maf


def addAnnotation(maf, NCBI_Build="37", Strand="+"):
    """
    adds NCBI_Build and Strand annotation on the whole maf file

    Args:
        maf (pandas.DataFrame): the maf file with columns: sample_col, genome_change_col, TCGAlocs
        NCBI_Build (str, optional): the NCBI build. Defaults to "37".
        Strand (str, optional): the strand. Defaults to "+".

    Returns:
        pandas.DataFrame: the maf file with the added columns: NCBI_Build, Strand
    """
    maf["NCBI_Build"] = NCBI_Build
    maf["Strand"] = Strand
    return maf


def makeMatrices(
    maf,
    homin=0.95,
    id_col=constants.SAMPLEID,
    hugo_col=constants.HUGO_COL,
    lof_col=constants.LIKELY_LOF_COL,
    hotspot_col=constants.HOTSPOT_COL,
):
    """generates genotyped hotspot, driver and damaging mutation matrices

    Returns:
        hotspot_mat (pd.DataFrame): genotyped hotspot mutation matrix. 0 == not damaging, 1 == heterozygous, 2 == homozygous
        lof_mat (pd.DataFrame): genotyped damaging mutation matrix. 0 == not damaging, 1 == heterozygous, 2 == homozygous
        driver_mat (pd.DataFrame): genotyped driver mutation matrix. 0 == not driver, 1 == heterozygous, 2 == homozygous
    """
    print("generating genotyped driver and damaging mutation matrix")
    gene_names = list(maf[hugo_col].unique())
    hotspot_mat = pd.DataFrame(columns=gene_names)
    lof_mat = pd.DataFrame(columns=gene_names)
    sample_ids = list(maf[id_col].unique())
    le = len(sample_ids)
    for j in range(le):
        h.showcount(j, le)
        sample = sample_ids[j]
        subset_maf = maf[maf[id_col] == sample]
        # hotspot
        hotspot = subset_maf[subset_maf[hotspot_col] == True]
        homhotspot = set(hotspot[hotspot["GT"] == "1|1"][hugo_col])
        for dup in h.dups(hotspot[hugo_col]):
            if hotspot[hotspot[hugo_col] == dup]["AF"].astype(float).sum() >= homin:
                homhotspot.add(dup)
        hethotspot = set(hotspot[hugo_col]) - homhotspot
        hotspot_mat.loc[sample, list(homhotspot)] = "2"
        hotspot_mat.loc[sample, list(hethotspot)] = "1"
        # damaging
        lof = subset_maf[(subset_maf[lof_col] == True)]
        homlof = set(lof[lof["GT"] == "1|1"][hugo_col])
        for dup in h.dups(lof[hugo_col]):
            if lof[lof[hugo_col] == dup]["AF"].astype(float).sum() >= homin:
                homlof.add(dup)
        hetlof = set(lof[hugo_col]) - homlof
        lof_mat.loc[sample, list(homlof)] = "2"
        lof_mat.loc[sample, list(hetlof)] = "1"
    hotspot_mat = hotspot_mat.dropna(axis="columns", how="all")
    lof_mat = lof_mat.dropna(axis="columns", how="all")
    hotspot_mat = hotspot_mat.fillna(0).astype(int)
    lof_mat = lof_mat.fillna(0).astype(int)

    return hotspot_mat, lof_mat


def managingDuplicates(samples, failed, datatype, tracker):
    """
    managingDuplicates manages the duplicates in the samples

    by only keeping the ones that are not old and did not fail QC

    Args:
        samples (list): the list of samples
        failed (list): the list of failed samples
        datatype (str): the data type to look at in the sample tracker
        tracker (pd.df): the sample tracker

    Returns:
        dict: the renaming dict
    """
    # selecting the right arxspan id (latest version)
    renaming = tracker.removeOlderVersions(
        names=samples,
        refsamples=tracker[tracker.datatype == datatype],
        priority="prioritized",
    )

    # reparing QC when we have a better duplicate
    ref = pd.DataFrame(tracker[tracker.datatype == datatype]["arxspan_id"])
    replace = 0
    for val in failed:
        if val in list(renaming.keys()):
            a = ref[ref.arxspan_id == ref.loc[val, "arxspan_id"]].index
            for v in a:
                if v not in failed:
                    renaming[v] = renaming.pop(val)
                    replace += 1
                    break
    print("could replace:")
    print(replace)
    return renaming


def aggregateMAFs(
    wm,
    sampleset="all",
    mafcol=constants.MAF_COL,
    keep_cols=constants.MUTCOL_DEPMAP,
    debug=False,
):
    """Aggregate MAF files from terra

    Args:
        refworkspace (str): the reference workspace
        sampleset (str, optional): the sample set to use. Defaults to 'all'.
        mutCol (str, optional): the MAF column name. Defaults to "somatic_maf".
        keep_cols (list, optional): which columns to keep in the aggregate MAF file. Defaults to constants.MUTCOL_DEPMAP

    Returns:
        aggregated_maf (df.DataFrame): aggregated MAF
    """
    sample_table = wm.get_samples()
    samples_in_set = wm.get_sample_sets().loc[sampleset]["samples"]
    sample_table = sample_table[sample_table.index.isin(samples_in_set)]
    print(mafcol)
    sample_table_valid = sample_table[~sample_table[mafcol].isna()]
    na_samples = set(sample_table.index) - set(sample_table_valid.index)
    print(str(len(na_samples)) + " samples don't have corresponding maf: ", na_samples)
    all_mafs = []
    counter = 0
    for name, row in tqdm(sample_table_valid.iterrows(), total=len(sample_table_valid)):
        # prints out progress bar
        maf = pd.read_csv(row[mafcol])
        maf[constants.SAMPLEID] = name
        # >1 because of the hess_signature typo in input mafs
        # can be 0 once the type is fixed upstream
        # TODO: replace hess_signature later
        if len(set(keep_cols.keys()) - set(maf.columns)) > 1:
            print(name + " is missing columns:")
            print(set(keep_cols.keys()) - set(maf.columns))
        all_mafs.append(maf)
        if debug:
            counter += 1
            if counter > 6:
                break
    all_mafs = pd.concat(all_mafs)
    return all_mafs


def aggregateSV(
    wm,
    sampleset="all",
    sv_colname=constants.SV_COLNAME,
    sv_header=constants.SV_HEADER,
    save_output=constants.WORKING_DIR,
    save_filename=constants.SV_FILENAME,
):
    """aggregate SVs pulled from a terra workspace

    Args:
        wm (dm.WorkspaceManager): workspace manager for the workspace where the data is stored
        sampleset (str, optional): name of the sample set to aggregate
        sv_colname (str, optional): name of the column in terra workspace where sample-level SVs are stored
        sv_header (list, optional): columns that are expected in the SV file
        save_output (str, optional): whether to save our data.
        save_filename (str, optional): name of the saved file

    Returns:
        all_svs (pd.DataFrame): aggregated SVs in the bedpe format

    """
    print("aggregating SVs")
    sample_table = wm.get_samples()
    samples_in_set = wm.get_sample_sets().loc[sampleset]["samples"]
    sample_table = sample_table[sample_table.index.isin(samples_in_set)]
    sample_table_valid = sample_table[~sample_table[sv_colname].isna()]
    na_samples = set(sample_table.index) - set(sample_table_valid.index)
    print(str(len(na_samples)) + " samples don't have corresponding sv: ", na_samples)
    all_svs = []
    for name, row in sample_table_valid.iterrows():
        sv = pd.read_csv(row[sv_colname], sep="\t")
        sv[constants.SAMPLEID] = name
        all_svs.append(sv)
    all_svs = pd.concat(all_svs)
    # making sure all expected columns are present in the dataframe
    assert set(sv_header) - set(all_svs.columns) == set()
    print("saving aggregated SVs")
    all_svs.to_csv(save_output + save_filename, index=False)
    return all_svs


def sv_internal_af_filter(bedpe, cutoff=constants.SV_INTERNAL_AF_CUTOFF):
    """In order to filter out artifacts, calculate the allele frequencies of SVs
    within the sample set. Remove ones that are above the cutoff and not rescued

    Args:
        bedpe (pd.DataFrame): aggregated SV table in bedpe format
        cutoff (float): max allele frequency allowed

    Returns:
        filtered_bedpe (pd.DataFrame): aggregated SV table in bedpe format, with SVs that pass the AF filter only
    """
    print("filtering SVs based on internal AF")
    from collections import Counter

    internal_afs = bedpe.loc[
        :,
        ["CHROM_A", "START_A", "END_A", "ALT_A", "CHROM_B", "START_B", "END_B", "TYPE"],
    ].apply(lambda x: ":".join(map(str, x)), axis=1)
    total_samples = bedpe[constants.SAMPLEID].unique().shape[0]
    internal_afs_ratio_dict = {}
    for k, v in Counter(internal_afs.tolist()).items():
        internal_afs_ratio_dict[k] = v / total_samples
    bedpe.loc[:, "internal_afs"] = internal_afs.map(internal_afs_ratio_dict)

    filtered_bedpe = bedpe[
        (bedpe.internal_afs <= cutoff) | (bedpe.Rescue == True)
    ].drop(["internal_afs"], axis=1)

    return filtered_bedpe


def generate_sv_matrix(
    df,
    id_col=constants.SAMPLEID,
    type_colname="TYPE",
    genea_colname="SYMBOL_A",
    geneb_colname="SYMBOL_B",
    del_colname="DEL_SYMBOLS",
    dup_colname="DUP_SYMBOLS",
    save_output=constants.WORKING_DIR,
    save_filename=constants.SV_MAT_FILENAME,
):
    """generate sample x gene matrix indicating which genes are affected by which type(s) of SVs

    Args:
        df (pd.DataFrame): aggregated SV table in bedpe format
        id_col (str, optional): name of the column in df that contains sample IDs
        type_colname (str, optional): name of the column in df that contains SV type annotation
        genea_colname (str, optional): name of the column in df that contains the name of the gene at breakpoint A
        geneb_colname (str, optional): name of the column in df that contains the name of the gene at breakpoint B
        del_colname (str, optional): name of the column in df that contains genes spanned by DELs
        dup_colname (str, optional): name of the column in df that contains genes spanned by DUPs
        save_output (str, optional): whether to save our data. Defaults to "".
        save_filename (str, optional): name of the saved file

    Returns:
        sv_mat (pd.DataFrame): sample x gene matrix
    """

    df[del_colname] = df[del_colname].fillna(".")
    df[dup_colname] = df[dup_colname].fillna(".")

    # gather all unique gene symbols in the df
    symbols = (
        df[genea_colname].str.split(", ").tolist()
        + df[geneb_colname].str.split(", ").tolist()
        + df[del_colname].str.split(", ").tolist()
        + df[dup_colname].str.split(", ").tolist()
    )
    symbol_flat_list = set([x for xs in symbols for x in xs])
    symbol_flat_list.remove(".")

    sample_ids = list(df[id_col].unique())
    le = len(sample_ids)

    # init list of dicts, each dict in this list accounts for one row (one sample)
    ds = []
    # iterate over all sample ids
    for j in range(le):
        h.showcount(j, le)
        sample = sample_ids[j]
        # init a gene: SV type dict for each sample
        d = {k: [] for k in symbol_flat_list}
        # subset df to only one sample
        subset_bedpe = df[df[id_col] == sample]

        # for BNDs, get genes at breakpoint A and breakpoint B
        bnds = subset_bedpe[subset_bedpe[type_colname] == "BND"]
        bnd_genes = (
            bnds[genea_colname].str.split(", ").tolist()
            + bnds[geneb_colname].str.split(", ").tolist()
        )
        bnd_genes = set([x for xs in bnd_genes for x in xs])
        bnd_genes.remove(".")

        for g in bnd_genes:
            d[g].append("BND")

        # for INSs, get genes at breakpoint A
        inss = subset_bedpe[subset_bedpe[type_colname] == "INS"]
        ins_genes = inss[genea_colname].str.split(", ").tolist()
        ins_genes = set([x for xs in ins_genes for x in xs])
        ins_genes.remove(".")

        for g in ins_genes:
            d[g].append("INS")

        # for DELs, get genes that intersect with the whole segment
        dels = subset_bedpe[subset_bedpe[type_colname] == "DEL"]
        del_genes = dels[del_colname].str.split(", ").tolist()
        del_genes = set([x for xs in del_genes for x in xs])
        del_genes.remove(".")

        for g in del_genes:
            d[g].append("DEL")

        # for DUPs, get genes that intersect with the whole segment
        dups = subset_bedpe[subset_bedpe[type_colname] == "DUP"]
        dup_genes = dups[dup_colname].str.split(", ").tolist()
        dup_genes = set([x for xs in dup_genes for x in xs])
        dup_genes.remove(".")

        for g in dup_genes:
            d[g].append("DUP")

        ds.append(d)

    sv_mat = pd.DataFrame(ds, index=sample_ids).applymap(lambda x: ", ".join(x))
    sv_mat.index.name = id_col
    sv_mat.to_csv(save_output + save_filename)

    return sv_mat


def aggregateGermlineMatrix(
    wm,
    sampleset="all",
    binary_mut_colname_dict=constants.BINARY_MUT_COLNAME_DICT,
    save_output="",
):
    """aggregate binary guide mutations pulled from a terra workspace

    Args:
        wm (str): terra workspace where the data is stored
        sampleset (str): name of the sample set in the workspace
        binary_mut_colname_dict (dict): dictionary mapping library name to corresponding column name
                                        in terra workspace where the binarized mutation calls are stored
        save_output (str, optional): whether to save our data. Defaults to "".

    Returns:
    """
    print("aggregating binary mutation matrices")
    sample_table = wm.get_samples()
    samples_in_set = wm.get_sample_sets().loc[sampleset]["samples"]
    sample_table = sample_table[sample_table.index.isin(samples_in_set)]
    all_guide_matrices = dict()
    for lib, colname in binary_mut_colname_dict.items():
        sample_table_valid = sample_table[~sample_table[colname].isna()]
        na_samples = set(sample_table.index) - set(sample_table_valid.index)
        print(
            str(len(na_samples))
            + " samples don't have corresponding binarized mutation calls for "
            + lib
            + ": ",
            na_samples,
        )
        all_muts = []
        header = False
        for name, row in tqdm(
            sample_table_valid.iterrows(), total=len(sample_table_valid)
        ):
            sample_mut = pd.read_csv(row[colname], names=["chrom", "start", "end", "sgRNA", name], sep="\t")
            if header == False:
                all_muts.append(sample_mut)
                header = True
            else:
                all_muts.append(sample_mut[[name]])
        all_muts = pd.concat(all_muts, axis=1)
        print("saving aggregated binary mutation matrices")
        all_guide_matrices[lib] = all_muts
        all_muts.to_csv(save_output + lib + "_binary_guide_mutations", index=False)
    return all_guide_matrices


def postProcess(
    wm,
    sampleset="all",
    mafcol=constants.MAF_COL,
    save_output=constants.WORKING_DIR,
    sv_col=constants.SV_COLNAME,
    sv_filename=constants.SV_FILENAME,
    sv_mat_filename=constants.SV_MAT_FILENAME,
    sv_header=constants.SV_HEADER,
    run_sv=True,
    sv_af_cutoff=constants.SV_INTERNAL_AF_CUTOFF,
    debug=False,
):
    """Calls functions to aggregate MAF files, annotate likely immortalization status of mutations,

    and aggregate structural variants (SVs)

    Args:
        wm (dalmatian.WorkspaceManager): workspace manager of the reference workspace
        sampleset (str, optional): the sample set to use. Defaults to 'all'.
        mutCol (str, optional): the mutation column name. Defaults to "mut_AC".
        save_output (str, optional): the output file name to save results into. Defaults to "".
        doCleanup (bool, optional): whether to clean up the workspace. Defaults to False.
        rename_cols (dict, optional): the rename dict for the columns.
            Defaults to {"i_ExAC_AF": "ExAC_AF",
                        "Tumor_Sample_Barcode": constants.SAMPLEID,
                        "Tumor_Seq_Allele2": "Tumor_Allele"}.

    Returns:
        pandas.DataFrame: the maf file with the added columns: variant_annotation
    """
    h.createFoldersFor(save_output)
    print("loading from Terra")
    # if save_output:
    # terra.saveConfigs(refworkspace, save_output + 'config/')
    mutations = aggregateMAFs(
        wm,
        sampleset=sampleset,
        mafcol=mafcol,
        keep_cols=constants.MUTCOL_DEPMAP,
        debug=debug,
    )

    print("further filtering and standardizing maf")
    mutations_with_standard_cols = postprocess_main_steps(mutations)

    print("saving somatic mutations (all)")
    #  /home/ubuntu/depmap_omics/depmapomics/mutations.py:314:71 - error: Argument of type "None" cannot be assigned to parameter "index" of type "_bool" in function "to_csv"
    #      Type "None" cannot be assigned to type "_bool" (reportGeneralTypeIssues)
    mutations_with_standard_cols.to_csv(
        save_output + "somatic_mutations_all.csv", index=False
    )
    print("done")

    svs = None
    sv_mat = None
    if run_sv:
        svs = aggregateSV(
            wm,
            sampleset=sampleset,
            sv_colname=sv_col,
            save_output=save_output,
            save_filename=sv_filename,
            sv_header=sv_header,
        )
        svs = sv_internal_af_filter(svs, cutoff=sv_af_cutoff)
        sv_mat = generate_sv_matrix(
            svs, save_output=save_output, save_filename=sv_mat_filename
        )

    return mutations_with_standard_cols, svs, sv_mat


def GetVariantClassification(
    vep_seq_ontology: str, var_type: str, inframe: bool
) -> str:
    """Map VEP sequence ontology into MAF variant classifications,
    VEP consequences is ordered by http://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
    """

    if re.match(
        r"^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)",
        vep_seq_ontology,
    ):
        return "Splice_Site"

    if re.match(r"^(stop_gained)", vep_seq_ontology):
        return "Nonsense_Mutation"

    if (
        re.match(r"^(frameshift_variant)", vep_seq_ontology)
        or (re.match(r"^(protein_altering_variant)", vep_seq_ontology) and not inframe)
    ) and (var_type == "DEL"):
        return "Frame_Shift_Del"

    if (
        re.match(r"^(frameshift_variant)", vep_seq_ontology)
        or (re.match(r"^(protein_altering_variant)", vep_seq_ontology) and not inframe)
    ) and (var_type == "INS"):
        return "Frame_Shift_Ins"

    if re.match(r"^(stop_lost)", vep_seq_ontology):
        return "Nonstop_Mutation"

    if re.match(r"^(initiator_codon_variant|start_lost)", vep_seq_ontology):
        return "Translation_Start_Site"

    if re.match(
        r"^(inframe_insertion|conservative_inframe_insertion|disruptive_inframe_insertion)",
        vep_seq_ontology,
    ) or (
        re.match(r"^(protein_altering_variant)", vep_seq_ontology)
        and inframe
        and (var_type == "INS")
    ):
        return "In_Frame_Ins"

    if re.match(
        r"^(inframe_deletion|disruptive_inframe_deletion|conservative_inframe_deletion)",
        vep_seq_ontology,
    ) or (
        re.match(r"^(protein_altering_variant)", vep_seq_ontology)
        and inframe
        and (var_type == "DEL")
    ):
        return "In_Frame_Del"

    if re.match(
        r"^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)",
        vep_seq_ontology,
    ):
        return "Missense_Mutation"

    if re.match(
        r"^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)",
        vep_seq_ontology,
    ):
        return "Intron"

    if re.match(
        r"^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant|start_retained_variant)",
        vep_seq_ontology,
    ):
        return "Silent"

    if re.match(
        r"^(splice_region_variant|splice_polypyrimidine_tract_variant|splice_donor_5th_base_variant|splice_donor_region_variant)",
        vep_seq_ontology,
    ):
        return "Splice_Region"

    if re.match(
        r"^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant|coding_transcript_variant)",
        vep_seq_ontology,
    ):
        return "RNA"

    if re.match(
        r"^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)",
        vep_seq_ontology,
    ):
        return "5'UTR"

    if re.match(r"^3_prime_UTR_variant", vep_seq_ontology):
        return "3'UTR"

    if re.match(r"^upstream_gene_variant", vep_seq_ontology):
        return "5'Flank"

    if re.match(r"^downstream_gene_variant", vep_seq_ontology):
        return "3'Flank"

    if re.match(
        r"^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)",
        vep_seq_ontology,
    ):
        return "IGR"

    if vep_seq_ontology == "":
        return "NoAnnotation"

    return "TargetedRegion"


def GetMafEndPosition(start: int, ref: str, alt: str) -> tuple:
    """Get the end position from the VCF start position and ref alt alleles

    Learn the complex InDel from https://github.com/qinqian/vcf2maf/blob/main/vcf2maf.pl#L706


    Return
    -----------
    (start, vartype, inframe): tuple
    """
    assert len(ref) > 0
    assert len(alt) > 0
    if len(ref) == len(alt):
        var_type_dict = {1: "SNP", 2: "DNP", 3: "TNP"}
        inframe = False
        if len(alt) > 3:
            var_type = "ONP"
        else:
            var_type = var_type_dict[len(alt)]
        return start, start + len(alt) - 1, var_type, inframe
    elif len(ref) < len(alt):
        # Insertion
        var_type = "INS"
        inframe = abs(len(ref) - len(alt)) % 3 == 0
        if ref == "-":
            return start - 1, start, var_type, inframe
        else:
            return start, start + len(ref) - 1, var_type, inframe
    else:
        # Deletion
        inframe = abs(len(ref) - len(alt)) % 3 == 0
        var_type = "DEL"
        return start, start + len(ref) - 1, var_type, inframe


def standardize_maf(maf: pd.DataFrame):
    """Standardize DepMap csv file into MAF 2.4 format

    Parameter
    --------------
    maf: pd.DataFrame
         a data frame that loads all the variants

    """
    formatted_coords = maf.loc[:, ["pos", "ref", "alt"]].apply(
        lambda x: GetMafEndPosition(*x), axis=1, result_type="expand"
    )

    maf.loc[:, "Strand"] = "+"
    maf.loc[:, "Start_Position"] = formatted_coords[0]
    maf.loc[:, "End_Position"] = formatted_coords[1]
    maf.loc[:, "Variant_Type"] = formatted_coords[2]
    maf.loc[:, "InFrame"] = formatted_coords[3]
    maf.loc[:, "Variant_Classification"] = maf.loc[
        :, ["variant_info", "Variant_Type", "InFrame"]
    ].apply(lambda x: GetVariantClassification(*x), axis=1)

    print((maf["pos"] - maf["Start_Position"]).sum())
    assert (
        maf["pos"] - maf["Start_Position"]
    ).sum() == 0, "Standardizing MAF shifted start position"

    maf["Hugo_Symbol"] = maf["hugo_symbol"]
    maf["Chromosome"] = maf["chrom"]
    maf["Reference_Allele"] = maf["ref"]
    maf["Alternate_Allele"] = maf["alt"]

    try:
        maf["Tumor_Sample_Barcode"] = maf[constants.SAMPLEID]
    except KeyError:
        maf["Tumor_Sample_Barcode"] = maf["CDS_ID"]
    maf["Protein_Change"] = maf["protein_change"]
    maf.loc[:, "NCBI_Build"] = "GRCh38"
    maf.loc[:, "Center"] = "DepMap"
    maf.loc[:, "Tumor_Seq_Allele1"] = maf.loc[:, "Reference_Allele"]
    maf.loc[:, "Tumor_Seq_Allele2"] = maf.loc[:, "Alternate_Allele"]
    reordered_columns = [
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
    ]
    reordered_columns += list(set(maf.columns) - set(reordered_columns))
    maf = maf.loc[:, reordered_columns]
    return maf


def patchEGFR(
    maf,
    hugo_col="Hugo_Symbol",
    protein_col="Protein_Change",
    inframe_col="InFrame",
    oncohotspot_col="oncokb_hotspot",
):
    """mark EGFR in frame deletions as hotspots"""
    topatch = maf[
        (maf[hugo_col] == "EGFR")
        & (maf[protein_col].str.endswith("del"))
        & (maf[inframe_col])
    ].index.tolist()
    maf.loc[topatch, oncohotspot_col] = True
    return maf


def convertProteinChange(
    maf,
    protein_cols=["Protein_Change", "protein_change"],
    protein_dict=constants.PROTEIN_DICT,
):
    """reformat the protein change column to exclude ensembl protein ids,
    and rename protein code from 3-letter to 1-letter"""
    for protein_col in protein_cols:
        maf[protein_col] = maf[protein_col].str.split(":").str[1]
        maf[protein_col] = maf[protein_col].fillna("")
        maf[protein_col] = maf[protein_col].replace(protein_dict, regex=True)
        maf[protein_col] = maf[protein_col].replace(r"^\s*$", np.nan, regex=True)

    return maf


def addLikelyLoF(row, vep_col="vep_impact", oncoimpact_col="oncokb_effect"):
    """add likely LoF column: true if a variant is high vep impact or likely lof according to oncoKB"""
    if (
        row[vep_col] == "HIGH"
        or row[oncoimpact_col] == "Likely Loss-of-function"
        or row[oncoimpact_col] == "Loss-of-function"
    ):
        return True
    else:
        return False


def addRescueReason(maf, rescue_reason_colname="rescue_reason"):
    """add column to indicate why variants are rescued

    Args:
        maf (df): MAF-like dataframe containing variants
        rescue_reason_colname (str): name of the rescue reason column

    Returns:
        maf (df): MAF-like dataframe containing variants, with rescue reason"""
    # initialize column as empty lists
    maf[rescue_reason_colname] = [[] for _ in range(maf.shape[0])]

    # go over each possible rescue reason in vcf_to_depmap and append them to the list
    maf.loc[
        (maf["oncokb_effect"].isin(["Loss-of-function", "Gain-of-function"]))
        | (maf["oncokb_oncogenic"] == "Oncogenic")
        | (maf["oncokb_hotspot"] == "Y"),
        "rescue_reason",
    ].apply(lambda x: x.append("OncoKB"))
    maf.loc[(maf["cosmic_tier"] == 1), "rescue_reason"].apply(
        lambda x: x.append("Cosmic")
    )
    maf.loc[(maf["brca1_func_score"].astype(float) <= -1.328), "rescue_reason"].apply(
        lambda x: x.append("BRCA1_score")
    )
    maf.loc[(maf["oncogene_high_impact"] == True), "rescue_reason"].apply(
        lambda x: x.append("Oncogene_high_impact")
    )
    maf.loc[(maf["tumor_suppressor_high_impact"] == True), "rescue_reason"].apply(
        lambda x: x.append("TS_high_impact")
    )
    maf.loc[(maf["hess_driver"] == "Y"), "rescue_reason"].apply(
        lambda x: x.append("Hess")
    )
    maf.loc[
        (
            (maf["Hugo_Symbol"] == "TERT")
            & (maf["pos"] >= 1295054)
            & (maf["pos"] <= 1295365)
        ),
        "rescue_reason",
    ].apply(lambda x: x.append("TERT"))
    maf.loc[
        (
            (maf["Hugo_Symbol"] == "MET")
            & (maf["pos"] >= 116771825)
            & (maf["pos"] <= 116771840)
        ),
        "rescue_reason",
    ].apply(lambda x: x.append("MET"))

    # join list of strings into one big string
    maf[rescue_reason_colname] = maf[rescue_reason_colname].apply(
        lambda x: ", ".join(x)
    )

    return maf


def postprocess_main_steps(
    maf: pd.DataFrame,
    adjusted_gnomad_af_cutoff: float = 1e-3,
    max_recurrence: float = 0.1,
) -> pd.DataFrame:
    """DepMap postprocessing steps after vcf_to_depmap

    Parameter
    ------------
    maf: pd.DataFrame
        a data frame from variants aggregation results

    """
    # force nan to be zero
    maf.loc[:, "gnomade_af"] = maf.loc[:, "gnomade_af"].fillna(0)
    maf.loc[:, "gnomadg_af"] = maf.loc[:, "gnomadg_af"].fillna(0)

    # step 1: filter the leftmost synonymous mutation
    maf = maf.query("~variant_info.str.contains('^synony', regex=True)")

    # step 2: less stringent cutoff for gnomad
    maf = maf[
        (maf["gnomade_af"] < adjusted_gnomad_af_cutoff)
        & (maf["gnomadg_af"] < adjusted_gnomad_af_cutoff)
    ]

    # step 3: remove all silent mutation classes
    #         remove variants without gene symbols
    maf = standardize_maf(maf)
    # add rescue for silent mutations that include TERT
    # because TERT belongs to 5'Flank
    maf = maf.loc[
        (
            ~maf.Variant_Classification.isin(
                [
                    "Silent",
                    "RNA",
                    "Intron",
                    "5'UTR",
                    "3'Flank",
                    "Splice_Region",
                    "5'Flank",
                ]
            )
        )
        | (maf.rescue),
        :,
    ]
    # maf = maf.loc[~maf.Variant_Classification.isin(['Silent', 'RNA', 'Intron', "3'Flank", 'Splice_Region']), :]
    maf = maf.loc[~maf.Hugo_Symbol.isnull(), :]
    maf = maf.sort_values(by=["Chromosome", "Start_Position", "End_Position"])

    # step 4: re-annotate missing EGFR hotspots
    maf = patchEGFR(maf)

    # step 5: convert protein change from 3-letter to 1-letter
    maf = convertProteinChange(maf)

    # step 6: add likely LoF column based on vep impact and oncokb mutation effect
    maf["likely_lof"] = maf.apply(addLikelyLoF, axis=1)

    # step 7: add hotspot column based on multiplt criteria
    maf["hotspot"] = False
    maf.loc[
        (
            (
                (maf[constants.HESS_COL] == "Y")
                | (maf[constants.ONCOKB_HOTSPOT_COL] == "Y")
                | (maf[constants.COSMIC_TIER_COL] == 1)
            ),
            "hotspot",
        )
    ] = True
    # manually classify two TERT promoter mutations
    maf.loc[
        (maf["chrom"] == "chr5") & (maf["pos"] == 1295135) & (maf["alt"] == "A"),
        "hotspot",
    ] = True
    maf.loc[
        (maf["chrom"] == "chr5") & (maf["pos"] == 1295113) & (maf["alt"] == "A"),
        "hotspot",
    ] = True

    # manually classify MET intron13 PPT mutations
    maf.loc[
        (maf["chrom"] == "chr7")
        & (maf["pos"] >= 116771825)
        & (maf["pos"] <= 116771840),
        "hotspot",
    ] = True
    print("unique hotspot genes: ")
    print(len(maf[maf["hotspot"] == True]["Hugo_Symbol"].unique()))

    # add rescue reason column
    maf = addRescueReason(maf)

    # step 8: remove high af from DepMap cohort
    internal_afs = maf.loc[
        :,
        [
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
        ],
    ].apply(lambda x: ":".join(map(str, x)), axis=1)
    total_samples = maf.Tumor_Sample_Barcode.unique().shape[0]
    # assume there are very few duplicated variants per sample
    # actually we have total 4 duplicated variants, it is trivial
    internal_afs_ratio_dict = {}
    for k, v in Counter(internal_afs.tolist()).items():
        internal_afs_ratio_dict[k] = v / total_samples
    maf.loc[:, "internal_afs"] = internal_afs.map(internal_afs_ratio_dict)
    maf = maf.loc[(maf.internal_afs <= max_recurrence) | (maf.rescue), :]
    return maf
