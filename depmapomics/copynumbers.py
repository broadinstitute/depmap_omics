# cn.py

from depmap_omics_upload import tracker as track
from depmapomics.config import *
from IPython.display import Image, display
import dalmatian as dm
import pandas as pd
import os
from genepy import mutations as mut
from genepy.utils import helper as h
from genepy import terra


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
    return df.rename(columns=COLRENAMING)


def loadFromGATKAggregation(
    refworkspace,
    setEntity="sample_set",
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    save_output="",
    doCleanup=True,
    todrop=[],
    showPlots=False,
    colname="combined_seg_file",
    plotColname="modeled_segments_plot_tumor",
    tempFolder=WORKING_DIR,
    toremove=["readgroup_ubams",],
    sampleset="all",
    colRenaming=COLRENAMING,
):
    """fetching the data from the Terra and loading it into a dataframe

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        doCleanup (bool, optional): if do cleanup of Terra workspace unused output and logs. Defaults to True.
        todrop (list, optional): [description]. Defaults to [].
        showPlots (bool, optional): whether to show plot output from the GATK pipeline. Defaults to False.
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        plotColname (str, optional): the column on terra where the plots exist in sample. Defaults to "modeled_segments_plot_tumor".
        tempFolder (str, optional): where to put temp files. Defaults to "working/".
        toremove (list, optional): columns in Terra samples to remove (and delete corresponding data). Defaults to ["readgroup_ubams", ].
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to COLRENAMING.

    Returns:
        pd.dataframe: dataframe containing the segments concatenated in a bed like format
    """
    wm = dm.WorkspaceManager(refworkspace)
    if save_output:
        terra.saveConfigs(refworkspace, os.path.join(save_output, "terra/"))

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
    segments = segments[~segments[SAMPLEID].isin(todrop)].reset_index(drop=True)
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.SegmentMean = 2 ** segments.SegmentMean
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments.loc[
        segments[segments.Chromosome.isin(["X", "Y"])].index, "SegmentMean"
    ] = (segments[segments.Chromosome.isin(["X", "Y"])]["SegmentMean"] / 2)
    segments = segments.sort_values(by=sortby)

    print("loading " + str(len(set(segments[SAMPLEID]))) + " rows")
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


def pureCNpostprocess(
    refworkspace,
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    lohcolname=PURECN_LOH_COLNAME,
    failedcolname=PURECN_FAILED_COLNAME,
    sampleset=PURECN_SAMPLESET,
    colRenaming=PURECN_COLRENAMING,
    lohvals=PURECN_LOHVALUES,
    save_output="",
    mappingdf=None,
):
    """fetching PureCN data from Terra, generate one matrix for absolute copy number, one matrix for LOH,
    and one cell line signature table

    Args:
        refworkspace (str): workspace path
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        save_output (str, optional): location to save output. Defaults to ''.
        todrop (list, optional): [description]. Defaults to [].
        colname (str, optional): the column in Terra where the file is saved in sampleset. Defaults to "combined_seg_file".
        sampleset (str, optional): sample set to load in terra. Defaults to "all".
        colRenaming (dict, optional): segment renaming dict. Defaults to COLRENAMING.
        genechangethresh (float, optional): above this threshold of variance of gene CN, the sample is considered failed. Defaults to 0.025.
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
    failed = samples[samples[failedcolname] == "TRUE"].index

    # removing the duplicates
    segments = segments[~segments[SAMPLEID].isin(todrop)].reset_index(drop=True)
    if "chr" in segments["Chromosome"][0]:
        segments["Chromosome"] = [i[3:] for i in segments["Chromosome"]]
    # tranforming the df
    segments.Start = segments.Start.astype(int)
    segments.End = segments.End.astype(int)
    segments = segments.sort_values(by=sortby)

    mappingdf = mappingdf[mappingdf["Chromosome"].isin(set(segments["Chromosome"]))]

    print("loading " + str(len(set(segments[SAMPLEID]))) + " rows")

    # Generate gene-level absolute cn matrix
    absolute_genecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments),
        mappingdf,
        style="closest",
        value_colname="MajorAlleleAbsoluteCN",
    )

    segments = segments[
        ~segments[SAMPLEID].isin(set(failed) | set(todrop))
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
    )

    loh_status = loh_status[~loh_status.index.isin(set(failed) | set(todrop))]
    loh_status = (loh_status > 0).astype(int)

    print("PureCN: saving LOH matrix")
    loh_status.to_csv(save_output + "purecn_loh_all.csv")
    print("done")

    return segments, absolute_genecn, loh_status


def generateSigTable(
    refworkspace,
    todrop=[],
    purecncols=SIGTABLE_TERRACOLS,
    misccols=MISC_SIG_TERRACOLS,
    purecn_failed_col=PURECN_FAILED_COLNAME,
    binary_cols=SIGTABLE_BINARYCOLS,
    colrenaming=SIGTABLE_RENAMING,
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
    sig_table = sig_table.rename(columns=colrenaming)

    print("saving global genomic feature table")
    sig_table.to_csv(save_output + "globalGenomicFeatures_all.csv")
    print("done")

    return sig_table


def postProcess(
    refworkspace,
    setEntity="sample_set",
    sampleset="all",
    purecnsampleset=PURECN_SAMPLESET,
    save_output="",
    doCleanup=True,
    sortby=[SAMPLEID, "Chromosome", "Start", "End"],
    todrop=[],
    priority=[],
    genechangethresh=GENECHANGETHR,
    segmentsthresh=SEGMENTSTHR,
    ensemblserver=ENSEMBL_SERVER_V,
    source_rename={},
    useCache=False,
    maxYchrom=MAXYCHROM,
):
    """post process an aggregated CN segment file, the CCLE way

    take a CN segment file from the Aggregate_WGS master terra workflow and post process it
    in the CCLE way.

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        sortby (list, optional): columns to sort df. Defaults to [SAMPLEID, 'Chromosome', "Start", "End"].
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        genechangethresh (float, optional): above this threshold of variance of gene CN, the sample is considered failed. Defaults to 0.025.
        segmentsthresh (int, optional): above this threshold of number of segments the WGS sample is considered failed. Defaults to 1500.
        ensemblserver (str, optional): ensembl server biomart version . Defaults to ENSEMBL_SERVER_V.
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
    mybiomart = mybiomart.drop_duplicates("hgnc_symbol", keep="first")
    mybiomart["gene_name"] = [
        i["hgnc_symbol"] + " (" + str(i["entrezgene_id"]).split(".")[0] + ")"
        for _, i in mybiomart.iterrows()
    ]
    mybiomart = mybiomart[
        ~mybiomart.entrezgene_id.isna()
    ]  # dropping all nan entrez id cols
    # drop Ychrom if > maxYchrom
    ychrom = segments[segments.Chromosome.str.contains("Y")]
    countYdrop = [
        i
        for i in set(ychrom[SAMPLEID])
        if len(ychrom[ychrom[SAMPLEID] == i]) > maxYchrom
    ]
    segments = segments[
        ~((segments[SAMPLEID].isin(countYdrop)) & (segments.Chromosome == "Y"))
    ]
    genecn = mut.toGeneMatrix(
        mut.manageGapsInSegments(segments), mybiomart, value_colname="SegmentMean"
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
        ~segments[SAMPLEID].isin((set(failed) | set(todrop)) - set(priority))
    ].reset_index(drop=True)
    genecn = genecn[~genecn.index.isin((set(failed) | set(todrop)) - set(priority))]

    # saving
    print("saving files")
    segments.to_csv(save_output + "segments_all.csv", index=False)
    genecn.to_csv(save_output + "genecn_all.csv")
    print("done")
    purecn_segments, purecn_genecn, loh_status = pureCNpostprocess(
        refworkspace,
        sampleset=purecnsampleset,
        mappingdf=mybiomart,
        sortby=sortby,
        todrop=todrop,
        save_output=save_output,
    )
    feature_table = generateSigTable(refworkspace, save_output=save_output)
    return (
        segments,
        genecn,
        failed,
        purecn_segments,
        purecn_genecn,
        loh_status,
        feature_table,
    )
