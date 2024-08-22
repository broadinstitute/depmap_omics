import pytest

from depmapomics import constants
import dalmatian as dm
import pandas as pd
import os.path
import seaborn as sns
from tqdm import tqdm

from mgenepy import terra


def addToMainFusion(input_filenames, main_filename, sample_id=constants.SAMPLEID):
    """
    Given a tsv fusion files from RSEM algorithm, merge it to a tsv set of fusion data

    Args:
    ----
        input_filenames: a set of filepath to input the files should be c|tsv from Terra fusion pipeline
        main_filename: a filepath to input the files should be c|tsv from Terra aggregation pipeline
    """
    maindata = pd.read_csv(main_filename, sep="\t")
    if "." in maindata[sample_id][0]:
        maindata[sample_id] = [
            i[0] for i in maindata[sample_id].str.split(".").tolist()
        ]
    samples = set(maindata[sample_id].tolist())
    with open(main_filename, "a") as f:
        for input_filename in input_filenames:
            df = pd.read_csv(input_filename, sep="\t")
            input_filename = input_filename.split("/")[-1].split(".")[0]
            if input_filename in samples:
                print(input_filename + " is Already in main fusions")
            df[sample_id] = pd.Series(
                [input_filename] * len(df.index.tolist()), index=df.index
            )
            cols = df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            df = df[cols]
            df.to_csv(f, header=False, sep="\t", index=False)


def filterFusions(
    fusions,
    sampleCol,
    maxfreq=constants.FUSION_MAXFREQ,
    minffpm=constants.FUSION_MINFFPM,
    maxffpm=constants.FUSION_MAXFFPM,
    countCol="CCLE_count",
    red_herring=constants.FUSION_RED_HERRING,
    **kwargs
):
    """
    Given a fusion file from star fusion, filters it (will also filter Mitochrondria and HLA genes)

    We want to apply filters to the fusion table to reduce the number of artifacts in the dataset. Specifically, we filter the following:
    * Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes
    * Remove red herring fusions (from STAR-Fusion annotations column)
    * Remove recurrent in CCLE (>= 25 samples)
    * Remove fusion with (SpliceType=" INCL_NON_REF_SPLICE" and LargeAnchorSupport="No" and FFPM < 0.1)
    * Remove fusions with FFPM < 0.05 (STAR-Fusion suggests using 0.1, but looking at the translocation data, this looks like it might be too aggressive)

    Args:
        fusions (pd.df): the fusion data. Should contain: LeftBreakpoint, RightBreakpoint, FusionName, annots, SpliceType, LargeAnchorSupport, FFPM columns
        maxfreq (int): the max allowed frequency of that fusion across our samples. default is 0.1
        samplecol (str): colname for the sample ids. Should be in the fusions dataframe.
        countCol (str): colname where are stored counts of that fusion name across our samples. Default is "CCLE_count"
        minffpm (int): minimum ffpm freq to filter on. Default is 0.05
        red_herring (list[str]): flags to filter on. default is constants.FUSION_RED_HERRING

    Returns:
        (pd.df): the filtered fusion dataframe
    """
    fusions = fusions.copy()
    # remove recurrent
    fusions = fusions[fusions[countCol] < len(set(fusions[sampleCol])) * maxfreq]
    # (1) Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes,
    fusions = fusions[
        ~(
            fusions["LeftBreakpoint"].str.contains("chrM")
            & fusions["RightBreakpoint"].str.contains("chrM")
        )
    ]
    fusions = fusions[~fusions["FusionName"].str.contains("^HLA\\-")]
    # (2) Remove red herring fusions
    fusions = fusions[
        ~fusions["annots"].str.contains("|".join(red_herring), case=False)
    ]
    # (4) Removed fusion with (SpliceType=" INCL_NON_REF_SPLICE" and
    # LargeAnchorSupport="No" and minFAF<0.02), or
    fusions = fusions[
        ~(
            (fusions["SpliceType"] == "INCL_NON_REF_SPLICE")
            & (fusions["LargeAnchorSupport"] == "NO_LDAS")
            & (fusions["FFPM"] < maxffpm)
        )
    ]
    # STAR-Fusion suggests using 0.1, but after looking at the
    # translocation data, this looks like it might be too aggressive
    fusions = fusions[fusions["FFPM"] > minffpm]
    return fusions


def renameFusionGene(a):
    """
    Given a list of fusion names from star-fusion, renames and returns them

    Args:
        a (list[str]): list of fusion names from star-fusion

    Returns:
        (list[str]): list of renamed fusion names
    """
    return [
        str(i.split("^")).replace(", ", " (").replace("'", "")[1:-1] + ")" for i in a
    ]


def standardizeGeneNames(fusions):
    """
    converts [GENE_NAME]^[ENSG] --> [GENE_NAME] ([ENSG])

    Example: "SMAD4^ENSG00000141646.14" --> "SMAD4 (ENSG00000141646.14)"

    Args:
        fusions (pd.df): fusion dataframe

    Returns:
        (pd.df): fusion dataframe with standardized gene names
    """
    fusions[["LeftGene", "RightGene"]] = fusions[["LeftGene", "RightGene"]].applymap(
        lambda x: "{} ({})".format(*x.split(r"^"))
    )
    return fusions


def postProcess(
    refworkspace,
    sampleCol='TCGA_sample_id',
    samplesetToLoad="all",
    colnames=constants.FUSION_COLNAME,
    todrop=[],
    doplot=True,
    countCol="CCLE_count",
    save_output="",
    rnFunc=None,
    renaming=None,
    **kwargs
):
    """post process an aggregated fusion files in the CCLE way

    (usually from the aggregate_Fusion terra workflow)

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether and where to save our data. Defaults to "".
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        sampleCol (str, optional): column name for the sample id in the dataset. Defaults to "CCLE_sample_id".
        colnames (str, optional): column names where the fusion file is, on the workspace. Defaults to constants.FUSION_COLNAME.
        doplot (bool, optional): whether to plot the data. Defaults to True.
        countCol (str, optional): column name for the count of the fusion. Defaults to "CCLE_count".
        save_output (str, optional): whether and where to save our data. Defaults to "".
        rnFunc (function, optional): function to rename the sample names
        (takes a list of sample names and returns a list of sample names). Defaults to None.
        renaming (dict(str:str), optional): dictionary to rename the sample names otherwise. Defaults to None.

    Returns:
        (pd.df): fusion dataframe
        (pd.df): filtered fusion dataframe
    """

    """
R code for aggregating fusions:
library(plyr)
library(dplyr)
library(readr)
library(stringr)

## reads selected fields from maf files and aggregates them 
args <- commandArgs(TRUE)
print(args)

output_file_name = args[1]
input_file_names = args[2]

# Load files and ensure that we are only reading in files that exist
fns = as.character(read.csv(input_file_names, header=FALSE)[,])
fe = sapply(fns, file.exists)
samples = fns[which(fe)]
nsamples = length(samples)

# Iterate across samples, read tsv, add the DepMap_ID column and keep the rest of them
for(i in 1:nsamples){
  f <- samples[i]
  
  # Extract the arxspan id from the file name
  # sample <- stringr::str_extract(string = f, pattern = 'ACH\\-[0-9]+')
  
  # In this workspace, samples are not indexed by ARXSPAND ID but by CCLE_name. Will need to update script in the future
  sample <- gsub('\\.fusions.annotated$', '', gsub('.*/', '', f))
  segs <- read_tsv(f, col_names = TRUE, col_types = cols(.default = "c")) %>%
    mutate(DepMap_ID=sample) %>% 
    dplyr::select(DepMap_ID, everything())
  
  # Write to the output file
  write.table(segs,  file=output_file_name,  sep='\t', row.names = FALSE, quote = FALSE, append = (i>1), col.names = !(i>1))
}

    """
    refwm = dm.WorkspaceManager(refworkspace)
    samples = refwm.get_samples()
    fusions_list = []
    print("loading fusions")
    # pytest.set_trace()
    for i,s in tqdm(samples.iterrows()):
        sample_df = pd.read_csv(s.fusion_tsv,sep='\t')
        sample_df[sampleCol] = s.name
        fusions_list.append(sample_df)

    fusions = pd.concat(fusions_list)
    # pytest.set_trace()
    #aggregated = refwm.get_sample_sets().loc[samplesetToLoad]["fusions_star"]
    """
    fusions = pd.read_csv(
        aggregated, names=[sampleCol] + colnames, skiprows=1, sep="\t"
    )
    """
    fusions = fusions.rename(columns={'#FusionName':'FusionName'})
    fusions[sampleCol] = fusions[sampleCol].str.split(".").str[0]
    print("postprocessing fusions")
    fusions.RightGene = renameFusionGene(fusions.RightGene)
    fusions.LeftGene = renameFusionGene(fusions.LeftGene)

    count = (
        fusions[["LeftBreakpoint", "RightBreakpoint"]]
        .value_counts()
        .to_frame(name=countCol)
    )
    fusions = pd.merge(fusions, count, on=["LeftBreakpoint", "RightBreakpoint"])

    # removing failed
    fusions = fusions[~fusions[sampleCol].isin(todrop)]

    print("saving")
    fusions.to_csv(os.path.join(save_output, "fusions_all.csv"), index=False)
    if rnFunc is not None or renaming is not None:
        print("renaming")
        renaming = rnFunc(set(fusions[sampleCol])) if rnFunc is not None else renaming
        fusions = (
            fusions[fusions[sampleCol].isin(renaming.keys())]
            .replace({sampleCol: renaming})
            .reset_index(drop=True)
        )
        fusions.to_csv(os.path.join(save_output, "fusions_latest.csv"), index=False)

    fusions_filtered = filterFusions(
        fusions, sampleCol=sampleCol, countCol=countCol, **kwargs
    )
    if doplot:
        print("Plotting filtered fusions:")
        sns.kdeplot(fusions_filtered[countCol])
    fusions_filtered.to_csv(
        os.path.join(save_output, "filteredfusions_latest.csv"), index=False
    )

    print("done")
    return fusions, fusions_filtered
