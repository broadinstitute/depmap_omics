from genepy import terra
from genepy.utils import helper as h
from genepy import mutations as mut
import os
import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
import numpy as np
from collections import Counter

from depmapomics.config import *

import dalmatian as dm
import pandas as pd
from genepy.epigenetics import chipseq as chip
from itertools import repeat
import multiprocessing
import subprocess


def download_maf_from_workspace(
    refwm,
    sample_set_ids=["all_ice", "all_agilent"],
    output_maf="/tmp/mutation_filtered_terra_merged.txt",
):
    """
    TODO: javad to document
    """
    sample_sets = refwm.get_sample_sets()
    dfs = []
    for sample_set_id in sample_sets.index.intersection(sample_set_ids):
        cpFiles(
            sample_sets.loc[sample_set_id, "filtered_CGA_MAF_aggregated"],
            "/tmp/tmp.txt",
            payer_project_id="broad-firecloud-ccle",
            verbose=False,
        )
        df = pd.read_csv("/tmp/tmp.txt", sep="\t", low_memory=False)
        dfs.append(df)
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(output_maf, index=False, sep="\t")
    return dfs_concat


def annotateLikelyImmortalized(
    maf,
    sample_col=SAMPLEID,
    genome_change_col="Genome_Change",
    TCGAlocs=["TCGAhsCnt", "COSMIChsCnt"],
    max_recurrence=0.05,
    min_tcga_true_cancer=5,
):
    """annotateLikelyImmortalized annotates the maf file with the likely immortalized mutations

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
    leng = len(set(maf[sample_col]))
    tocheck = []
    for k, v in Counter(maf[genome_change_col].tolist()).items():
        if v > max_recurrence * leng:
            tocheck.append(k)
    for val in list(set(tocheck) - set([np.nan])):
        if (
            np.nan_to_num(maf[maf[genome_change_col] == val][TCGAlocs], 0).max()
            < min_tcga_true_cancer
        ):
            maf.loc[
                maf[maf[genome_change_col] == val].index, "is_likely_immortalization"
            ] = True
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


def add_variant_annotation_column(maf):
    """
    adds variant annotation column to the maf file

    Args:
        maf (pandas.DataFrame): the maf file with columns: sample_col, genome_change_col, TCGAlocs

    Returns:
        pandas.DataFrame: the maf file with the added column: variant_annotation
    """
    rename = {}
    for k, v in MUTATION_GROUPS.items():
        for e in v:
            rename[e] = k
    maf["Variant_annotation"] = [
        rename[i] for i in maf["Variant_Classification"].tolist()
    ]
    return maf


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


def postProcess(
    refworkspace,
    sampleset="all",
    mutCol="mut_AC",
    save_output="",
    doCleanup=False,
    rename_cols={
        "i_ExAC_AF": "ExAC_AF",
        "Tumor_Sample_Barcode": SAMPLEID,
        "Tumor_Seq_Allele2": "Tumor_Allele",
    },
    sv_col=SV_COLNAME,
    sv_filename=SV_FILENAME,
    sv_renaming=SV_COLRENAME,
):
    """post process an aggregated MAF file the CCLE way

    (usually a MAF file from the Aggregate_MAF Terra worklflow)

    Args:
        refworkspace (str): the reference workspace
        sampleset (str, optional): the sample set to use. Defaults to 'all'.
        mutCol (str, optional): the mutation column name. Defaults to "mut_AC".
        save_output (str, optional): the output file name to save results into. Defaults to "".
        doCleanup (bool, optional): whether to clean up the workspace. Defaults to False.
        rename_cols (dict, optional): the rename dict for the columns.
            Defaults to {"i_ExAC_AF": "ExAC_AF", 
                        "Tumor_Sample_Barcode": SAMPLEID,
                        "Tumor_Seq_Allele2": "Tumor_Allele"}.

    Returns:
        pandas.DataFrame: the maf file with the added columns: variant_annotation
    """
    h.createFoldersFor(save_output)
    print("loading from Terra")
    # if save_output:
    # terra.saveConfigs(refworkspace, save_output + 'config/')
    refwm = dm.WorkspaceManager(refworkspace)
    mutations = pd.read_csv(
        refwm.get_sample_sets().loc[sampleset, "filtered_CGA_MAF_aggregated"], sep="\t"
    )
    mutations = mutations.rename(columns=rename_cols).drop(
        columns=["Center", "Tumor_Seq_Allele1"]
    )

    mutations[mutCol] = [
        str(i[0]) + ":" + str(i[1])
        for i in np.nan_to_num(
            mutations[["t_alt_count", "t_ref_count"]].values, 0
        ).astype(int)
    ]
    mutations = mut.filterCoverage(mutations, loc=[mutCol], sep=":", cov=2)
    mutations = mut.filterAllelicFraction(mutations, loc=[mutCol], sep=":", frac=0.1)
    mutations = addAnnotation(mutations, NCBI_Build="37", Strand="+")
    mutations = annotateLikelyImmortalized(
        mutations,
        TCGAlocs=["TCGAhsCnt", "COSMIChsCnt"],
        max_recurrence=0.05,
        min_tcga_true_cancer=5,
    )

    mutations.to_csv(save_output + "somatic_mutations_all.csv", index=None)
    print("done")

    sampleids = refwm.get_sample_sets().loc[sampleset, "samples"]
    svs = aggregateSV(
        refworkspace,
        sampleids,
        all_sv_colname=sv_col,
        save_output=save_output,
        save_filename=sv_filename,
    )
    return mutations, svs


def aggregateSV(
    refworkspace,
    sampleids,
    all_sv_colname=SV_COLNAME,
    sv_renaming=SV_COLRENAME,
    save_output="",
    save_filename="",
):
    """aggregate SVs pulled from a terra workspace

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleids (list[str]): list of sequencing IDs
        all_sv_colname (str): name of column in terra workspace that contains sv output files. Defaults to "filtered_annotated_sv"
        save_output (str, optional): whether to save our data. Defaults to "".

    Returns:
        
    """
    print("aggregating SVs")
    wm = dm.WorkspaceManager(refworkspace).disable_hound()
    sample_table = wm.get_samples()
    sample_table = sample_table[sample_table.index.isin(sampleids)]
    all_svs = []
    for name, row in sample_table.iterrows():
        sv = pd.read_csv(row[all_sv_colname], sep="\t")
        sv[SAMPLEID] = name
        all_svs.append(sv)
    all_svs = pd.concat(all_svs)
    all_svs = all_svs.rename(columns=sv_renaming)
    print("saving aggregated SVs")
    all_svs.to_csv(save_output + save_filename, sep="\t", index=False)
    return all_svs


def mapBed(file, vcfdir, bed_location):
    """map mutations in one vcf file to regions in the guide bed file"""

    guides_bed = pd.read_csv(
        bed_location,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "foldchange"],
    )

    bed = pd.read_csv(
        vcfdir + file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "foldchange"],
    )
    bed["foldchange"] = 1
    name = file.split("/")[-1].split(".")[0].split("_")[1]
    if len(bed) == 0:
        return (name, None)
    val = chip.putInBed(guides_bed, bed, mergetype="sum")

    return (name, val)


def generateGermlineMatrix(
    vcfs,
    vcfdir,
    savedir=WORKING_DIR + SAMPLESETNAME + "/",
    filename="binary_mutguides.tsv.gz",
    bed_locations=GUIDESBED,
    cores=16,
):
    """generate profile-level germline mutation matrix for achilles' ancestry correction. VCF files are generated
    using the CCLE pipeline on terra
    Args:
        vcfs (list): list of vcf file locations (gs links)
        taiga_dataset (str, optional): taiga folder location. Defaults to TAIGA_CN.
        vcfdir (str, optional): directory where vcf files are saved.
        savedir (str, optional): directory where output germline matrices are saved.
        bed_location (str, optional): location of the guides bed file.
        vcf_colname (str, optional): vcf column name on terra.
        cores (int, optional): number of cores in parallel processing.
    
    Returns:
        sorted_mat (pd.DataFrame): binary matrix where each row is a region in the guide, and each column corresponds to a profile
    """
    # check if bcftools is installed
    print("generating germline matrix")
    print(
        "bcftools is required in order to generate the matrix. Checking if bcftools is installed..."
    )
    try:
        subprocess.call(["bcftools"])
    except FileNotFoundError:
        raise Exception("bcftools not installed!")

    h.createFoldersFor(savedir)
    # load vcfs from workspace using dalmatian
    h.createFoldersFor(vcfdir)

    # save vcfs from workspace locally,
    # and run bcftools query to transform vcfs into format needed for subsequent steps
    # only including regions in the guide bed file

    cmds = []
    for lib, _ in bed_locations.items():
        h.createFoldersFor(vcfdir + lib + "/")
    for sam in vcfs:
        cmd = (
            "gsutil cp "
            + sam
            + " "
            + vcfdir
            + sam.split("/")[-1]
            + " && gsutil cp "
            + sam
            + ".tbi"
            + " "
            + vcfdir
            + sam.split("/")[-1]
            + ".tbi && "
        )
        for lib, fn in bed_locations.items():
            cmd += (
                "bcftools query\
                    --exclude \"FILTER!='PASS'&GT!='mis'&GT!~'\.'\"\
                    --regions-file "
                + fn
                + " \
                    --format '%CHROM\\t%POS\\t%END\\t%ALT{0}\n' "
                + vcfdir
                + sam.split("/")[-1]
                + " > "
                + vcfdir
                + lib
                + "/"
                + "loc_"
                + sam.split("/")[-1].split(".")[0]
                + ".bed && "
            )
        cmd += "rm " + vcfdir + sam.split("/")[-1] + "*"
        cmds.append(cmd)

    print("running bcftools index and query")
    h.parrun(cmds, cores=cores)

    pool = multiprocessing.Pool(cores)
    binary_matrices = dict()
    for lib, fn in bed_locations.items():
        print("mapping to library: ", lib)
        res = pool.starmap(
            mapBed,
            zip(os.listdir(vcfdir + lib + "/"), repeat(vcfdir + lib + "/"), repeat(fn)),
        )
        guides_bed = pd.read_csv(
            fn, sep="\t", header=None, names=["chrom", "start", "end", "foldchange"],
        )
        sorted_guides_bed = guides_bed.sort_values(
            by=["chrom", "start", "end"]
        ).reset_index(drop=True)
        print("done pooling")
        for name, val in res:
            if val is not None:
                sorted_guides_bed[name] = val
        print("saving binary matrix for library: ", lib)
        sorted_guides_bed.to_csv(savedir + "_" + lib + "_" + filename)
        binary_matrices[lib] = sorted_guides_bed

    return binary_matrices
