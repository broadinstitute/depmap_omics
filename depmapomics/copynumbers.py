# cn.py

from depmapomics import constants

from depmap_omics_upload import tracker as track

from IPython.display import Image, display
import dalmatian as dm
import pandas as pd
import os
from mgenepy import mutations as mut
from mgenepy.utils import helper as h
import pybedtools  # type: ignore
import subprocess


def renameColumns(df):
    """
    rename some of the main columns names from RSEM, GATK.. to more readable column names

    Args:
    -----
      df: the df to rename

    Returns:
    ------
      df the renamed df
    """
    return df.rename(columns=constants.COLRENAMING)


def loadFromGATKAggregation(
    refworkspace,
    setEntity="sample_set",
    sortby=[constants.SAMPLEID, "Chromosome", "Start", "End"],
    doCleanup=True,
    todrop=[],
    showPlots=False,
    colname="combined_seg_file",
    plotColname="modeled_segments_plot_tumor",
    tempFolder=constants.WORKING_DIR,
    toremove=[
        "readgroup_ubams",
    ],
    sampleset="all",
    colRenaming=constants.COLRENAMING,
):
    """fetching the data from the Terra and loading it into a dataframe

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [constants.SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        doCleanup (bool, optional): if do cleanup of Terra workspace unused output and logs. Defaults to True.
        todrop (list, optional): [description]. Defaults to [].
        showPlots (bool, optional): whether to show plot output from the GATK pipeline. Defaults to False.
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        plotColname (str, optional): the column on terra where the plots exist in sample. Defaults to "modeled_segments_plot_tumor".
        tempFolder (str, optional): where to put temp files. Defaults to "working/".
        toremove (list, optional): columns in Terra samples to remove (and delete corresponding data). Defaults to ["readgroup_ubams", ].
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to constants.COLRENAMING.

    Returns:
        pd.dataframe: dataframe containing the segments concatenated in a bed like format
    """
    wm = dm.WorkspaceManager(refworkspace)

    if doCleanup:
        print("cleaning up")
        res = wm.get_samples()
        for val in toremove:
            if val in res.columns.tolist():
                wm.disable_hound().delete_entity_attributes(
                    "sample", res[val], delete_files=True
                )

    segments = pd.read_csv(
        wm.get_entities(setEntity).loc[sampleset, colname], sep="\t"
    ).rename(columns=colRenaming)

    # removing the duplicates
    segments = segments[~segments[constants.SAMPLEID].isin(todrop)].reset_index(
        drop=True
    )
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.SegmentMean = 2**segments.SegmentMean
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments.loc[
        segments[segments.Chromosome.isin(["X", "Y"])].index, "SegmentMean"
    ] = (segments[segments.Chromosome.isin(["X", "Y"])]["SegmentMean"] / 2)
    segments = segments.sort_values(by=sortby)

    print("loading " + str(len(set(segments[constants.SAMPLEID]))) + " rows")
    if showPlots:
        # plotting results of CN calls for this new sample set
        for i, (k, val) in enumerate(
            wm.get_samples().loc[wm.get_sample_sets().loc[sampleset].samples].iterrows()
        ):
            plot = val[plotColname]
            os.system("gsutil cp " + plot + " " + tempFolder)
            print(k)
            print(val["arxspan_id"], val["sex"])
            if i > 30:
                continue
            display(Image(os.path.join(tempFolder, plot.split("/")[-1])))
    return segments


def managingDuplicates(samples, failed, datatype, tracker, newname="arxspan_id"):
    """removes duplicates and solves failed data

    Args:
        failed (list): list of samples that failed QC
        datatype (str): 'wes' or 'wgs'
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        samples (list): list of samples to filter
        newname (str): name of the new sample set

    Returns:
        dict: a dict of duplicates to keep and their newname name
    """
    # selecting the right arxspan id (latest version)
    renaming = track.removeOlderVersions(
        names=samples,
        refsamples=tracker[tracker.datatype == datatype],
        priority="prioritized",
    )

    # reparing QC when we have a better duplicate
    ref = pd.DataFrame(tracker[tracker.datatype == datatype][newname])
    replace = 0
    for val in failed:
        if val in renaming:
            a = ref[ref[newname] == ref.loc[val][newname]].index
            for v in a:
                if v not in failed:
                    renaming[v] = renaming.pop(val)
                    replace += 1
                    break
    print("could replace:")
    print(replace)
    return renaming


def arm_call(
    df,
    cn_colname="SegmentAbsoluteCN",
    width_colname="seg_width",
    ploidy_colname="Ploidy",
):
    """determines arm-level CNA. Outputs can be 1 (gain), 0 (neutral), or -1 (loss)"""
    df = df[~df.arm.isna()]
    df.sort_values(cn_colname, inplace=True)
    cumsum = df[width_colname].cumsum()
    cutoff = df[width_colname].sum() / 2.0
    median = df[cn_colname][cumsum >= cutoff].iloc[0].round().astype(int)

    ploidy = df[ploidy_colname].iloc[0].round().astype(int)
    status = 0
    if median > ploidy:
        status = 1
    elif median < ploidy:
        status = -1
    return status


def get_which_arm(df, start_colname="Start", end_colname="End"):
    """given a segment, determine which chromosome arm it is on"""
    df["seg_cent"] = 0.5 * (df[start_colname] + df[end_colname])
    df["arm"] = None
    df.loc[df["seg_cent"] < df["cent_start"], "arm"] = "p"
    df.loc[df["seg_cent"] > df["cent_end"], "arm"] = "q"
    return df


def get_cna_and_aneuploidy(
    seg,
    sig_table,
    cent_filename=constants.HG38_CENTROMERE,
    id_col=constants.SAMPLEID,
    ploidy_col="Ploidy",
    save_output="",
):
    """Arm-level CNA matrix and add aneuploidy scores to signature table"""
    print("generating arm-level CNA and aneuploidy score")
    # parse centromere file
    cent = pd.read_csv(cent_filename, sep="\t", index_col=False)
    cent = (
        cent[
            ~(cent["#region_name"].str.startswith("HET"))
            & (~cent["chr"].isin(["X", "Y"]))
        ]
        .drop(columns=["#region_name"])
        .rename(
            columns={"chr": "Chromosome", "start": "cent_start", "stop": "cent_end"}
        )
    )
    cent["cent_mid"] = (
        (0.5 * (cent["cent_start"] + cent["cent_end"])).round().astype(int)
    )
    cent["Chromosome"] = cent["Chromosome"].astype(int)

    seg["Chromosome"] = seg["Chromosome"].astype(int)
    seg["seg_width"] = seg["End"] - seg["Start"]
    merged_seg = seg.merge(cent, on=["Chromosome"], how="left")
    sig_table = sig_table.reset_index().rename(columns={"index": id_col})
    sig_table[ploidy_col] = sig_table[ploidy_col].astype("float")
    merged_seg = merged_seg.merge(
        sig_table[[id_col, ploidy_col]], on=[id_col], how="left"
    )

    seg_with_arm = get_which_arm(merged_seg)
    seg_with_arm["chrom_arm"] = (
        seg_with_arm["Chromosome"].astype(str) + seg_with_arm["arm"]
    )

    cna_table = (
        seg_with_arm.groupby([id_col, "chrom_arm"]).apply(arm_call).unstack(level=1)
    )
    cna_table = cna_table[constants.CNA_ARMS]

    aneuploidy = cna_table.abs().sum(axis=1).to_dict()
    sig_table["Aneuploidy"] = sig_table[id_col].map(aneuploidy)
    sig_table = sig_table.set_index(id_col)

    print("Saving arm-level CNA matrix and signature table with aneuploidy score")
    sig_table.to_csv(save_output + "globalGenomicFeaturesWithAneuploidy_all.csv")
    cna_table.to_csv(save_output + "arm_cna_all.csv")
    print("done")

    return cna_table, sig_table


def pureCNpostprocess(
    refworkspace,
    sortby=[constants.SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    lohcolname=constants.PURECN_LOH_COLNAME,
    failedcolname=constants.PURECN_FAILED_COLNAME,
    sampleset=constants.PURECN_SAMPLESET,
    colRenaming=constants.PURECN_COLRENAMING,
    lohvals=constants.PURECN_LOHVALUES,
    save_output="",
    mappingdf=None,
    min_gof=constants.PURECN_MIN_GOF,
    max_ploidy=constants.PURECN_MAX_PLOIDY,
):
    """fetching PureCN data from Terra, generate one matrix for absolute copy number, one matrix for LOH,
    and one cell line signature table

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [constants.SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        todrop (list, optional): [description]. Defaults to [].
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to constants.COLRENAMING.
        segmentsthresh (int, optional): above this threshold of number of segments the WGS sample is considered failed. Defaults to 2000.

    Returns:
        pd.dataframe: dataframe containing the segments concatenated in a bed like format
    """

    print("loading PureCN merged LOH file")
    wm = dm.WorkspaceManager(refworkspace)
    segments = pd.read_csv(
        wm.get_entities("sample_set").loc[sampleset, lohcolname]
    ).rename(columns=colRenaming)
    samples = wm.get_samples()
    failed = set(samples[samples[failedcolname] == "TRUE"].index)

    # filter out lines with low GOF (would require manual curation)
    samples["PureCN_gof"] = (
        samples.PureCN_comment.str.extract(r"([0-9]+)", expand=True)
        .fillna(100)
        .astype(int)
    )
    samples["Non_aberrant"] = samples.PureCN_comment.str.contains(
        "NON-ABERRANT"
    ).astype(bool)
    samples = samples[samples.PureCN_ploidy != "NA"]
    to_curate_idx = set(
        samples[
            ((samples.PureCN_gof < min_gof) & ~samples.Non_aberrant)
            | (samples.PureCN_ploidy.astype(float) > max_ploidy)
            | (samples.PureCN_curated == "TRUE")
        ].index
    )
    failed.update(to_curate_idx)

    # removing the duplicates
    segments = segments[~segments[constants.SAMPLEID].isin(todrop)].reset_index(
        drop=True
    )
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments = segments.sort_values(by=sortby)

    mappingdf = mappingdf[mappingdf["Chromosome"].isin(set(segments["Chromosome"]))]

    print("loading " + str(len(set(segments[constants.SAMPLEID]))) + " rows")

    # Generate gene-level absolute cn matrix
    absolute_genecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments),
        mappingdf,
        style="closest",
        value_colname="SegmentAbsoluteCN",
        gene_names_col="ensembl_gene_id",
    )

    segments = segments[
        ~segments[constants.SAMPLEID].isin(set(failed) | set(todrop))
    ].reset_index(drop=True)
    absolute_genecn = absolute_genecn[
        ~absolute_genecn.index.isin(set(failed) | set(todrop))
    ]

    print("summary of PureCN absolute cn data:")
    print(
        absolute_genecn.values.min(),
        absolute_genecn.values.mean(),
        absolute_genecn.values.max(),
    )

    print("PureCN: saving seg and gene cn files")
    segments.to_csv(save_output + "purecn_segments_all.csv", index=False)
    absolute_genecn.to_csv(save_output + "purecn_genecn_all.csv")
    print("done")

    # Generate gene-level LOH status matrix
    segments_binarized = segments.copy()
    segments_binarized["LoHStatus"] = segments_binarized["LoHStatus"].replace(
        lohvals, 1
    )
    segments_binarized["LoHStatus"] = segments_binarized["LoHStatus"].fillna(0)

    loh_status = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments_binarized),
        mappingdf,
        style="closest",
        value_colname="LoHStatus",
        gene_names_col="ensembl_gene_id",
    )

    loh_status = loh_status[~loh_status.index.isin(set(failed) | set(todrop))]
    loh_status = (loh_status > 0).astype(int)

    print("PureCN: saving LOH matrix")
    loh_status.to_csv(save_output + "purecn_loh_all.csv")
    print("done")

    return segments, absolute_genecn, loh_status, failed


def generateSigTable(
    refworkspace,
    todrop=[],
    purecncols=constants.SIGTABLE_TERRACOLS,
    misccols=constants.MISC_SIG_TERRACOLS,
    purecn_failed_col=constants.PURECN_FAILED_COLNAME,
    binary_cols=constants.SIGTABLE_BINARYCOLS,
    colrenaming=constants.SIGTABLE_RENAMING,
    save_output="",
):
    print("generating global genomic feature table")
    wm = dm.WorkspaceManager(refworkspace).disable_hound()
    samples = wm.get_samples()
    samples = samples[samples.index != "nan"]

    # discard lines that have failed PureCN for the PureCN columns
    purecn_passed = samples[samples[purecn_failed_col] != "TRUE"]
    purecn_table = purecn_passed[purecncols]
    purecn_table = purecn_table[~purecn_table.index.isin(set(todrop))]

    # use 1/0 for binary values
    purecn_table = purecn_table.replace({"TRUE": 1, "FALSE": 0})
    purecn_table[binary_cols] = purecn_table[binary_cols].astype("Int64")

    # table containing cols that are not generated by pureCN,
    # which may include samples that have failed pureCN
    misc_table = samples[misccols]
    sig_table = misc_table.merge(
        purecn_table, how="outer", left_index=True, right_index=True
    )

    # rename columns
    sig_table.index.rename(constants.SAMPLEID, inplace=True)
    sig_table = sig_table.rename(columns=colrenaming)

    print("saving global genomic feature table")
    sig_table.to_csv(save_output + "globalGenomicFeatures_all.csv")
    print("done")

    return sig_table


def exonUnion(df):
    """given a dataframe from biomart, take the union of intervals"""
    sorted_data = sorted(list(zip(df.exon_chrom_start, df.exon_chrom_end)))
    b = []
    for begin, end in sorted_data:
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    return b


def maskGenes(
    save_output="",
    maskthresh=constants.GENEMASKTHRESH,
    segdup_bed=constants.SEGDUP_BED,
    repeat_bed=constants.RM_BED,
    bedtoolspath=constants.BEDTOOLSPATH,
    rescue_list=constants.ONCOKB_ONCOGENE_ENSG_LIST,
):
    """given a bed file consisting of highly repeated/duplicated regions, mask
    genes that overlap with those regions (a gene is masked if the portion of its gene
    body length that overlaps with those regions is higher than maskthresh)

    Returns a list of geness to be masked, minus rescue (intersection with oncokb oncogenes)
    """
    mybiomart = h.generateGeneNames(
        ensemble_server=constants.ENSEMBL_SERVER_V,
        attributes=["start_position", "end_position", "chromosome_name"],
    )
    mybiomart = mybiomart.rename(
        columns={
            "start_position": "start",
            "end_position": "end",
            "chromosome_name": "Chromosome",
        }
    )
    mybiomart["Chromosome"] = mybiomart["Chromosome"].astype(str)
    mybiomart = mybiomart.drop_duplicates("ensembl_gene_id", keep="first")

    # sort and format biomart
    mybiomart["Chromosome"] = mybiomart["Chromosome"].replace(
        {"X": "23", "Y": "24", "MT": "25"}
    )
    mybiomart = mybiomart[mybiomart.Chromosome.isin(set(map(str, range(1, 26))))]
    mybiomart["Chromosome"] = mybiomart["Chromosome"].astype(int)
    mybiomart = mybiomart.sort_values(by=["Chromosome", "start", "end"])
    mybiomart["Chromosome"] = mybiomart["Chromosome"].replace(
        {23: "X", 24: "Y", 25: "MT"}
    )
    mybiomart["Chromosome"] = "chr" + mybiomart["Chromosome"].astype(str)

    mybiomart[["Chromosome", "start", "end", "ensembl_gene_id"]].to_csv(
        save_output + "biomart_cngenes.bed", sep="\t", header=False, index=False
    )

    subprocess.call(
        [
            bedtoolspath
            + "bedtools intersect -a "
            + save_output
            + "biomart_cngenes.bed -b "
            + segdup_bed
            + " > "
            + save_output
            + "mask_overlap_segdup.bed"
        ],
        shell=True,
    )

    overlap_segdup = pd.read_csv(
        save_output + "mask_overlap_segdup.bed",
        sep="\t",
        names=["chrom", "start", "end", "ensembl_gene_id"],
    )

    # download and reformat exon info
    biomart_exon = h.generateGeneNames(
        attributes=[
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_gene_id",
        ],
        default_attr=[],
    )

    biomart_exon_union = (
        pd.DataFrame(
            biomart_exon.groupby(["ensembl_gene_id"]).apply(exonUnion),
            columns=["exons"],
        )
        .reset_index()
        .explode("exons")
        .reset_index(drop=True)
    )
    biomart_exon_union = pd.concat(
        [
            biomart_exon_union[["ensembl_gene_id"]],
            pd.DataFrame(biomart_exon_union["exons"].tolist()).add_prefix("col"),
        ],
        axis=1,
    )
    biomart_exon_union = biomart_exon_union.merge(
        biomart_exon[["ensembl_gene_id", "chromosome_name"]],
        on=["ensembl_gene_id"],
        how="left",
    )
    biomart_exon_union["chromosome_name"] = biomart_exon_union[
        "chromosome_name"
    ].replace({"X": "23", "Y": "24", "MT": "25"})
    biomart_exon_union = biomart_exon_union[
        biomart_exon_union.chromosome_name.isin(set(map(str, range(1, 26))))
    ]
    biomart_exon_union["chromosome_name"] = biomart_exon_union[
        "chromosome_name"
    ].astype(int)
    biomart_exon_union = biomart_exon_union.sort_values(
        by=["chromosome_name", "col0", "col1"]
    )
    biomart_exon_union["chromosome_name"] = biomart_exon_union[
        "chromosome_name"
    ].replace({23: "X", 24: "Y", 25: "MT"})
    biomart_exon_union["chromosome_name"] = "chr" + biomart_exon_union[
        "chromosome_name"
    ].astype(str)

    biomart_exon_union[["chromosome_name", "col0", "col1", "ensembl_gene_id"]].to_csv(
        save_output + "biomart_exons.bed", sep="\t", header=False, index=False
    )

    subprocess.call(
        [
            bedtoolspath
            + "bedtools intersect -a "
            + save_output
            + "biomart_exons.bed -b "
            + repeat_bed
            + " > "
            + save_output
            + "mask_overlap_rm.bed"
        ],
        shell=True,
    )

    overlap_rm = pd.read_csv(
        save_output + "mask_overlap_rm.bed",
        sep="\t",
        names=["chrom", "start", "end", "ensembl_gene_id"],
    )

    gene_dict = (
        mybiomart[["Chromosome", "start", "end", "ensembl_gene_id"]]
        .set_index("ensembl_gene_id")
        .T.to_dict("list")
    )

    with open(rescue_list, "r") as f:
        to_rescue = f.read().splitlines()

    to_rescue = set(to_rescue) & set(mybiomart["ensembl_gene_id"])

    print("rescuing " + str(len(to_rescue)) + " genes from oncokb's oncogene list")

    # segdup
    masked_genes_segdup = []
    for g in overlap_segdup.ensembl_gene_id.unique().tolist():
        _, start, end = gene_dict[g]
        gene_length = end - start
        overlap_length = 0
        overlap_segments = overlap_segdup[overlap_segdup.ensembl_gene_id == g]
        for i, v in overlap_segments.iterrows():
            overlap_length += v["end"] - v["start"]
        if overlap_length / gene_length > maskthresh:
            masked_genes_segdup.append(g)

    # repeat masker
    masked_genes_rm = []
    for g in overlap_rm.ensembl_gene_id.unique().tolist():
        all_overlaps = overlap_rm[overlap_rm.ensembl_gene_id == g]
        exons = biomart_exon_union[biomart_exon_union.ensembl_gene_id == g]
        exon_length = exons["col1"].sum() - exons["col0"].sum()
        overlap_length = all_overlaps["end"].sum() - all_overlaps["start"].sum()
        if overlap_length / exon_length > maskthresh:
            masked_genes_rm.append(g)

    gene_list = list(set(masked_genes_segdup + masked_genes_rm) - set(to_rescue))
    with open(save_output + "genes_to_mask_minus_rescue_ensg.txt", "w") as f:
        for line in gene_list:
            f.write(f"{line}\n")

    return gene_list


def read_ms_repeats(
    refworkspace, sampleset="all", colname="ms_repeats", save_output=""
):
    """read aggregated microsatellite repeat table into a dataframe

    Args:
        refworkspace (str): name of the terra workspace
        sampleset (str): name of the sample set in refworkspace
        colname (str): name of the column in the sampleset table where the aggregated file is stored

    Returns:
        ms_df (pd.DataFrame): dataframe containing MS repeats"""

    wm = dm.WorkspaceManager(refworkspace)
    ms_df = pd.read_csv(wm.get_sample_sets().loc[sampleset, colname])
    ms_df.to_csv(save_output + "ms_repeats_all.csv", index=False)

    return ms_df


def postProcess(
    refworkspace,
    setEntity="sample_set",
    sampleset="all",
    purecnsampleset=constants.PURECN_SAMPLESET,
    save_output="",
    doCleanup=True,
    sortby=[constants.SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    priority=[],
    segmentsthresh=constants.SEGMENTSTHR,
    ensemblserver=constants.ENSEMBL_SERVER_V,
    source_rename={},
    useCache=False,
    maxYchrom=constants.MAXYCHROM,
):
    """post process an aggregated CN segment file, the CCLE way

    take a CN segment file from the Aggregate_WGS master terra workflow and post process it
    in the CCLE way.

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        sortby (list, optional): columns to sort df. Defaults to [constants.SAMPLEID, 'Chromosome', "Start", "End"].
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        genechangethresh (float, optional): above this threshold of variance of gene CN, the sample is considered failed. Defaults to 0.025.
        segmentsthresh (int, optional): above this threshold of number of segments the WGS sample is considered failed. Defaults to 1500.
        ensemblserver (str, optional): ensembl server biomart version . Defaults to constants.ENSEMBL_SERVER_V.
        source_rename (dict, optional): dict to rename the source column if needed. Defaults to {}.
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.

    Returns:
        pd.dataframe: a dataframe with the post processed segment file as a bed file concatenated over all samples.
        pd.dataframe: a dataframe with the post processed gene CN file as a bed file of samples x genes.
        list: list of samples that failed to be processed.
    """
    h.createFoldersFor(save_output)
    print("loading CN from Terra")
    segments = loadFromGATKAggregation(
        refworkspace,
        setEntity=setEntity,
        sampleset=sampleset,
        sortby=sortby,
        todrop=todrop,
        doCleanup=doCleanup,
    )
    print("making gene level copy number")

    mybiomart = h.generateGeneNames(
        ensemble_server=ensemblserver,
        useCache=useCache,
        attributes=["start_position", "end_position", "chromosome_name"],
    )
    mybiomart = mybiomart.rename(
        columns={
            "start_position": "start",
            "end_position": "end",
            "chromosome_name": "Chromosome",
        }
    )
    mybiomart["Chromosome"] = mybiomart["Chromosome"].astype(str)
    mybiomart = mybiomart.sort_values(by=["Chromosome", "start", "end"])
    mybiomart = mybiomart[mybiomart["Chromosome"].isin(set(segments["Chromosome"]))]
    mybiomart = mybiomart.drop_duplicates("ensembl_gene_id", keep="first")
    # drop Ychrom if > maxYchrom
    ychrom = segments[segments.Chromosome.str.contains("Y")]
    countYdrop = [
        i
        for i in set(ychrom[constants.SAMPLEID])
        if len(ychrom[ychrom[constants.SAMPLEID] == i]) > maxYchrom
    ]
    segments = segments[
        ~(
            (segments[constants.SAMPLEID].isin(countYdrop))
            & (segments.Chromosome == "Y")
        )
    ]
    genecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments),
        mybiomart,
        value_colname="SegmentMean",
        gene_names_col="ensembl_gene_id",
    )
    # validation step
    print("summary of the gene cn data:")
    print(genecn.values.min(), genecn.values.mean(), genecn.values.max())
    failed = mut.checkAmountOfSegments(segments, thresh=segmentsthresh)

    print("failed our QC")
    print(failed)

    if source_rename:
        segments = segments.replace({"Source": source_rename})
    if save_output:
        h.listToFile(failed, save_output + "failed.txt")
    # subsetting
    segments = segments[
        ~segments[constants.SAMPLEID].isin((set(failed) | set(todrop)) - set(priority))
    ].reset_index(drop=True)
    genecn = genecn[~genecn.index.isin((set(failed) | set(todrop)) - set(priority))]

    # saving
    print("saving files")
    segments.to_csv(save_output + "segments_all.csv", index=False)
    genecn.to_csv(save_output + "genecn_all.csv")
    print("done")
    purecn_segments, purecn_genecn, loh_status, purecn_failed = pureCNpostprocess(
        refworkspace,
        sampleset=purecnsampleset,
        mappingdf=mybiomart,
        sortby=sortby,
        todrop=todrop,
        save_output=save_output,
    )
    feature_table = generateSigTable(
        refworkspace, todrop=purecn_failed, save_output=save_output
    )
    cna_table, feature_table = get_cna_and_aneuploidy(
        purecn_segments,
        feature_table,
        id_col=constants.SAMPLEID,
        save_output=save_output,
    )
    ms_df = read_ms_repeats(refworkspace, save_output=save_output)
    return (
        segments,
        genecn,
        failed,
        purecn_segments,
        purecn_genecn,
        loh_status,
        feature_table,
        cna_table,
        ms_df,
    )
