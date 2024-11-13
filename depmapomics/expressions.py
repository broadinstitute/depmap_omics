from depmapomics import constants
import os.path

import dalmatian as dm
import pandas as pd
import numpy as np
from scipy.stats import zscore

from mgenepy.utils import helper as h
from mgenepy import rna
from depmapomics.qc import rna as myQC
from tqdm import tqdm

def addSamplesRSEMToMain(input_filenames, main_filename):
    """
    given a tsv RNA files from RSEM algorithm, merge it to a tsv set of RNA data

    Args:
    ----
        input_filenames: a list of dict like file path in Terra gs://, outputs from the rsem pipeline
        main_filename: a dict like file paths in Terra gs://, outputs from rsem aggregate
    """
    genes_count = pd.read_csv(
        constants.WORKING_DIR
        + main_filename["rsem_genes_expected_count"].split("/")[-1],
        sep="\t",
        compression="gzip",
    )
    transcripts_tpm = pd.read_csv(
        constants.WORKING_DIR + main_filename["rsem_transcripts_tpm"].split("/")[-1],
        sep="\t",
        compression="gzip",
    )
    genes_tpm = pd.read_csv(
        constants.WORKING_DIR + main_filename["rsem_genes_tpm"].split("/")[-1],
        sep="\t",
        compression="gzip",
    )

    for input_filename in input_filenames:
        name = input_filename["rsem_genes"].split("/")[-1].split(".")[0].split("_")[-1]
        rsem_genes = pd.read_csv(
            constants.WORKING_DIR + input_filename["rsem_genes"].split("/")[-1],
            sep="\t",
        )
        rsem_transcripts = pd.read_csv(
            constants.WORKING_DIR + input_filename["rsem_isoforms"].split("/")[-1],
            sep="\t",
        )
        genes_count[name] = pd.Series(
            rsem_genes["expected_count"], index=rsem_genes.index
        )
        transcripts_tpm[name] = pd.Series(
            rsem_transcripts["TPM"], index=rsem_transcripts.index
        )
        genes_tpm[name] = pd.Series(rsem_genes["TPM"], index=rsem_genes.index)

    genes_count.to_csv(
        constants.WORKING_DIR
        + main_filename["rsem_genes_expected_count"].split("/")[-1],
        sep="\t",
        index=False,
        index_label=False,
        compression="gzip",
    )
    transcripts_tpm.to_csv(
        constants.WORKING_DIR + main_filename["rsem_transcripts_tpm"].split("/")[-1],
        sep="\t",
        index=False,
        index_label=False,
        compression="gzip",
    )
    genes_tpm.to_csv(
        constants.WORKING_DIR + main_filename["rsem_genes_tpm"].split("/")[-1],
        sep="\t",
        index=False,
        index_label=False,
        compression="gzip",
    )


def solveQC(tracker, failed, save="", newname="arxspan_id"):
    """create a renaming dict to rename the columns of the QC file

    based on which samples have failed QC and which are blacklisted int he sample tracker

    Args:
        tracker (dataframe[datatype, prioritized, arxspan_id, index, ($newname)]): the sample tracker containing necessary info to compute which duplicates to keep
        failed (list): list of samples that failed QC
        save (str): save the renaming dict to a file
        newname (str): name of the column in the tracker to rename to
    Returns:
        dict: a dict to rename the samples to
    """
    newfail = []
    rename = {}
    # finding other replicates to solve failed ones
    for val in failed:
        if val not in tracker.index:
            continue
        a = tracker.loc[val][newname]
        res = tracker[
            (tracker.datatype == "rna")
            & (tracker[newname] == a)
            & (tracker.blacklist != 1)
        ]
        if len(res) >= 1:
            for k in res.index:
                if k not in failed:
                    rename[val] = k
                else:
                    newfail.append(val)
        else:
            newfail.append(val)
    print("samples that failed:")
    print(newfail)
    if save:
        h.listToFile(newfail, save + "rna_newfailed.txt")
    return rename


def loadFromRSEMaggregate(
    refworkspace,
    todrop=[],
    filenames=constants.RSEMFILENAME,
    sampleset="all",
    renamingFunc=None,
    rsemfilelocs=None,
    folder="",
):
    """Load the rsem aggregated files from Terra

    Args:
        refworkspace (str): the workspace where to load the files from
        todrop (list[str], optional): list of samples to drop. Defaults to [].
        filenames (list[str], optional): the filenames to load. Defaults to constants.RSEMFILENAME.
        sampleset (str, optional): the sample set to load. Defaults to 'all'.
        renamingFunc (function, optional): the function to rename the samples
        (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
        rsemfilelocs (pd.DataFrame, optional): locations of RSEM output files if refworkspace is not provided (no interaction with terra)

    Returns:
        dict(str: pd.df): the loaded dataframes
        dict: the renaming dict used to rename the dfs columns

    """
    files = {}
    renaming = {}

    refwm = dm.WorkspaceManager(refworkspace)
    rsemfilelocs = refwm.get_sample_sets().loc[sampleset]
    for val in filenames:
        file = pd.read_csv(
            rsemfilelocs[val],
            compression="gzip",
            header=0,
            sep="\t",
            quotechar='"',
            error_bad_lines=False,
        )
        if renamingFunc is not None:
            # removing failed version
            renaming = renamingFunc(file.columns[2:], todrop, folder=folder)
        else:
            renaming.update({i: i for i in file.columns[2:] if i not in todrop})
        renaming.update({"transcript_id(s)": "transcript_id"})
        # we remove the failed samples where we did not found anything else to replace them with
        files[val] = file[
            file.columns[:2].tolist()
            + [i for i in file.columns[2:] if i in set(renaming.keys())]
        ].rename(columns=renaming)
    return files, renaming


def subsetGenes(
    files,
    gene_rename,
    filenames=constants.RSEM_TRANSCRIPTS,
    drop=[],
    index_id="transcript_id",
):
    """
    Subset the rsem transcripts file to keep only the genes of interest

    Args:
        files (dict(str: pd.dfs)): the rsem transcripts dfs to subset samples x genes
        gene_rename (dict): the gene renaming dict (here we expect a dict of ensembl transcript ids: gene names)
        filenames (list[str], optional): the dict dfs to look at. Defaults to constants.RSEM_TRANSCRIPTS.
        drop (list[str], optional): the genes to drop. Defaults to [].
        index (str, optional): the index to use. Defaults to 'transcript_id'.

    Returns:
        dict(str: pd.df): the subsetted dfs
    """
    print("subsetting " + index_id + " columns")
    rename_transcript = {}
    missing = []
    for val in filenames:
        if len(rename_transcript) == 0 and index_id == "transcript_id":
            for _, v in files[val].iterrows():
                if v["gene_id"].split(".")[0] in gene_rename:
                    rename_transcript[v["transcript_id"].split(".")[0]] = (
                        gene_rename[v["gene_id"].split(".")[0]].split(" (")[0]
                        + " ("
                        + v.transcript_id.split(".")[0]
                        + ")"
                    )
                else:
                    missing.append(v.gene_id.split(".")[0])
            print("missing: " + str(len(missing)) + " genes")
        file = files[val].drop(columns=drop).set_index(index_id)
        if "effective_length" not in val:
            file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
        r = [i.split(".")[0] for i in file.index]
        dup = h.dups(r)
        par_y_drop = []
        for d in dup:
            idx = file[
                (file.index.str.startswith(d)) & (~file.index.str.endswith("_PAR_Y"))
            ].index[0]
            idx_par_y = file[
                (file.index.str.startswith(d)) & (file.index.str.endswith("_PAR_Y"))
            ].index[0]
            file.loc[idx] = file.loc[[idx, idx_par_y], :].sum(numeric_only=True)
            par_y_drop.append(idx_par_y)
        file = file.drop(index=par_y_drop)
        r = [i.split(".")[0] for i in file.index]
        dup = h.dups(r)
        if len(dup) > 0:
            print(dup)
            raise ValueError("duplicate " + index_id)
        file.index = r
        file = file.rename(
            index=rename_transcript if len(rename_transcript) != 0 else gene_rename
        ).T
        files[val] = file
    return files


def extractProtCod(
    files,
    mybiomart,
    protcod_rename,
    filenames=constants.RSEMFILENAME_GENE,
    dropNonMatching=False,
    rep=("genes", "proteincoding_genes"),
):
    """extracts protein coding genes from a merged RSEM gene dataframe and a biomart dataframe

    Args:
        files (dict(str: pd.dfs)): the rsem transcripts dfs to subset samples x genes
        mybiomart (pd.df): the biomart dataframe should contain the following columns:
        'ensembl_gene_id', 'entrezgene_id', 'gene_biotype'
        protcod_rename (dict(str, str)): the protein coding gene renaming dict
        (here we expect a dict of ensembl transcript ids: gene names)
        filenames (list[str], optional): the dict dfs to look at. Defaults to constants.RSEMFILENAME_GENE.
        rep (tuple, optional): how to rename the protein gene subseted df copies in the dict. Defaults to ('genes', 'proteincoding_genes').

    Raises:
        ValueError: if the biomart dataframe does not contain the required columns

    Returns:
        dict(str: pd.df): the subsetted dfs
    """
    for val in filenames:
        name = val.replace(rep[0], rep[1])
        files[name] = files[val].drop(columns="transcript_id").set_index("gene_id")
        files[name] = files[name][(files[name].sum(1) != 0) & (files[name].var(1) != 0)]
        r = [i.split(".")[0] for i in files[name].index]
        dup = h.dups(r)
        # STAR, RSEM and GENCODE were updated in 22Q2, indroducing dup gene ids ending with "_PAR_Y"
        # here we sum these entries with their corresponding _PAR_Y entries
        par_y_drop = []
        for d in dup:
            idx = files[name][
                (files[name].index.str.startswith(d))
                & (~files[name].index.str.endswith("_PAR_Y"))
            ].index[0]
            idx_par_y = files[name][
                (files[name].index.str.startswith(d))
                & (files[name].index.str.endswith("_PAR_Y"))
            ].index[0]
            files[name].loc[idx] = (
                files[name].loc[[idx, idx_par_y], :].sum(numeric_only=True)
            )
            par_y_drop.append(idx_par_y)
        files[name] = files[name].drop(index=par_y_drop)

        r = [i.split(".")[0] for i in files[name].index]
        dup = h.dups(r)
        if len(dup) > 0:
            print(dup)
            raise ValueError("duplicate genes")
        files[name].index = r
        files[name] = files[name][
            files[name].index.isin(set(mybiomart.ensembl_gene_id))
        ].rename(index=protcod_rename)
        # removing genes that did not match.. pretty unfortunate
        if dropNonMatching:
            files[name] = files[name].loc[[i for i in files[name].index if " (" in i]]
        # check: do we have any duplicates?
        # if we do, managing duplicates
        if len(set(h.dups(files[name].index.tolist()))) > 0:
            print("we have duplicate gene names!!")
            for dup in h.dups(files[name].index):
                a = files[name].loc[dup].sum()
                files[name].drop(index=dup)
                files[name].loc[dup] = a
        files[name] = files[name].T

    return files


async def ssGSEA(tpm_genes, geneset_file=constants.SSGSEAFILEPATH, recompute=True):
    """the way we run ssGSEA on the CCLE dataset

    Args:
        tpm_genes (pd.df): the tpm genes dataframe
        geneset_file (str, optional): the path to the geneset file. Defaults to constants.SSGSEAFILEPATH.

    Returns:
        pd.df: the ssGSEA results
    """

    tpm_genes = tpm_genes.copy()
    tpm_genes.columns = [i.split(" (")[0] for i in tpm_genes.columns]

    # summing the different exons/duplicates
    for i in h.dups(tpm_genes.columns):
        val = tpm_genes[i].sum(1)
        tpm_genes = tpm_genes.drop(columns=i)
        tpm_genes[i] = val

    # total size of data
    print("total size of data")
    print(len(set([val for val in tpm_genes.columns if "." not in val])))
    tpm_genes = pd.DataFrame(
        data=zscore(np.log2(tpm_genes + 1), nan_policy="omit"),
        columns=tpm_genes.columns,
        index=tpm_genes.index,
    )

    # MAYBE NOT NEEDED
    #### merging splicing variants into the same gene
    # counts_genes_merged, _, _= h.mergeSplicingVariants(counts_genes.T, defined='.')

    enrichments = (
        await rna.gsva(
            tpm_genes.T, geneset_file=geneset_file, method="ssgsea", recompute=recompute
        )
    ).T
    enrichments.index = [i.replace(".", "-") for i in enrichments.index]
    return enrichments


def saveFiles(files, folder=constants.WORKING_DIR, rep=("rsem", "expression")):
    """
    saves the files in the dict to the folder

    Args:
        files (dict(str: pd.df)): the dfs to save
        folder (str, optional): the folder to save the files. Defaults to constants.WORKING_DIR.
        rep (tuple, optional): how to rename (parts of) the files. Defaults to ('rsem', 'expression').
    """
    print("storing files in {}".format(folder))
    for k, val in files.items():
        val.to_csv(os.path.join(folder, k.replace(rep[0], rep[1]) + ".csv"))
        if "tpm" in k:
            val.apply(lambda x: np.log2(x + 1)).to_csv(
                os.path.join(folder, k.replace(rep[0], rep[1]) + "_logp1.csv")
            )


def load_rnaseqc(terra_path):
    """read rnaseqc2 gene counts table from a file path"""
    rnaseqc_count_df = pd.read_csv(terra_path, sep='\t', skiprows=2)
    rnaseqc_count_df = rnaseqc_count_df.set_index(rnaseqc_count_df.apply(lambda x: f"{x[1]} ({x[0].split('.')[0]})", axis=1))
    rnaseqc_count_df = rnaseqc_count_df.drop(["Name", "Description"], axis=1)
    return rnaseqc_count_df

def parse_rnaseqc_counts(refworkspace, samplesetname, colname):
    """parse rnaseqc gene counts to one a list of dataframes"""
    refwm = dm.WorkspaceManager(refworkspace)
    terra_rnaseq_df = refwm.get_samples()
    samplesinset = [i["entityName"]for i in refwm.get_entities("sample_set").loc[samplesetname].samples]
    terra_rnaseq_df = terra_rnaseq_df.loc[samplesinset, :]
    rnaseqc_count_dfs = []
    for _, row in tqdm(terra_rnaseq_df.iterrows(), total=len(terra_rnaseq_df)):
        rnaseqc_count_dfs.append(load_rnaseqc(row[colname]))
    return rnaseqc_count_dfs

async def postProcess(
    refworkspace,
    samplesetname,
    save_output="",
    doCleanup=False,
    colstoclean=[],
    ensemblserver=constants.ENSEMBL_SERVER_V,
    todrop=[],
    samplesetToLoad="all",
    priority=[],
    geneLevelCols=constants.RSEMFILENAME_GENE,
    trancriptLevelCols=constants.RSEMFILENAME_TRANSCRIPTS,
    ssGSEAcol="genes_tpm",
    renamingFunc=None,
    samplesinset=[],
    rsemfilelocs=None,
    rnaqclocs={},
    useCache=False,
    dropNonMatching=False,
    recompute_ssgsea=True,
    compute_enrichment=True,
    dry_run=False,
):
    """postprocess a set of aggregated Expression table from RSEM in the CCLE way

    (usually using the aggregate_RSEM terra worklow)

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        colstoclean (list, optional): the columns to clean in the terra workspace. Defaults to [].
        ensemblserver (str, optional): ensembl server biomart version . Defaults to constants.ENSEMBL_SERVER_V.
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        geneLevelCols (list, optional): the columns that contain the gene level
        expression data in the workspace. Defaults to constants.RSEMFILENAME_GENE.
        trancriptLevelCols (list, optional): the columns that contain the transcript
        level expression data in the workspacce. Defaults to constants.RSEMFILENAME_TRANSCRIPTS.
        ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
        rsemfilelocs (pd.DataFrame, optional): locations of RSEM output files if refworkspace is not provided (bypass interaction with terra)
        samplesinset (list[str], optional): list of samples in the sampleset if refworkspace is not provided (bypass interaction with terra)
        rnaqclocs (dict(str:list[str]), optional): dict(sample_id:list[QC_filepaths]) of rna qc file locations if refworkspace is not provided (bypass interaction with terra)
        renamingFunc (function, optional): the function to use to rename the sample columns
        (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
        compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
        dropNonMatching (bool, optional): whether to drop the non matching genes
        between entrez and ensembl. Defaults to False.
        recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
    """
    if not samplesetToLoad:
        samplesetToLoad = samplesetname
    refwm = dm.WorkspaceManager(refworkspace)
    # if save_output:
    #     terra.saveWorkspace(refworkspace, save_output + "terra/")

    samplesinset = [
        i["entityName"]
        for i in refwm.get_entities("sample_set").loc[samplesetname].samples
    ]
    if doCleanup:
        print("cleaninp up data")
        res = refwm.get_samples()
        for val in colstoclean:
            if val in res.columns.tolist():
                refwm.disable_hound().delete_entity_attributes(
                    "sample", res[val], delete_files=True
                )
            else:
                print(val + " not in the workspace's data")

    print("load QC and generate QC report")
    _, lowqual, failed = myQC.plot_rnaseqc_results(
        refworkspace,
        samplesinset,
        save=bool(save_output),
        output_path=save_output + "rna_qcs/",
    )

    all_qc_df = myQC.export_qc(refworkspace, selected_samples=[]).transpose()
    if not os.path.exists(save_output + "rna_qcs/"):
        os.mkdir(save_output + "rna_qcs/") # type: ignore
    all_qc_df.to_csv(save_output + "rna_qcs/all_qc.csv")

    failed = failed.index.tolist()
    print("those samples completely failed qc: ", failed)
    print("rescuing from whitelist list: ", set(failed) & set(priority))
    failed = list(set(failed) - set(priority))
    failed.extend(todrop)

    print("generating gene names")

    mybiomart = h.generateGeneNames(ensemble_server=ensemblserver, useCache=useCache)
    # creating renaming index, keeping top name first
    gene_rename = {}
    for _, i in mybiomart.iterrows():
        if i.ensembl_gene_id not in gene_rename:
            gene_rename.update(
                {i.ensembl_gene_id: i.hgnc_symbol + " (" + i.ensembl_gene_id + ")"}
            )
    protcod_rename = {}
    for _, i in mybiomart[
        (~mybiomart.entrezgene_id.isna()) & (mybiomart.gene_biotype == "protein_coding")
    ].iterrows():
        if i.ensembl_gene_id not in protcod_rename:
            protcod_rename.update(
                {
                    i.ensembl_gene_id: i.hgnc_symbol
                    + " ("
                    + str(int(i.entrezgene_id))
                    + ")"
                }
            )

    print("loading files")
    files, renaming = loadFromRSEMaggregate(
        refworkspace,
        todrop=failed,
        filenames=trancriptLevelCols + geneLevelCols,
        sampleset=samplesetToLoad,
        renamingFunc=renamingFunc,
        rsemfilelocs=rsemfilelocs,
        folder=save_output,
    )
    if save_output:
        h.dictToFile(renaming, save_output + "rna_sample_renaming.json")
        lowqual.to_csv(save_output + "rna_lowqual_samples.csv")
        h.listToFile(list(set(failed)), save_output + "rna_failed_samples.txt")
    print("renaming files")
    # gene level
    if len(geneLevelCols) > 0:
        # import pdb; pdb.set_trace()
        files = extractProtCod(
            files,
            mybiomart[mybiomart.gene_biotype == "protein_coding"],
            protcod_rename,
            dropNonMatching=dropNonMatching,
            filenames=geneLevelCols,
        )
        # assert {v.columns[-1] for k,v in files.items()} == {'ACH-000052'}
        files = subsetGenes(
            files,
            gene_rename,
            filenames=geneLevelCols,
            index_id="gene_id",
            drop="transcript_id",
        )
        # assert {v.columns[-1] for k,v in files.items()} == {'ACH-000052'}
    if len(trancriptLevelCols) > 0:
        files = subsetGenes(
            files,
            gene_rename,
            filenames=trancriptLevelCols,
            drop="gene_id",
            index_id="transcript_id",
        )

    if compute_enrichment:
        print("doing ssGSEA")
        enrichments = await ssGSEA(files[ssGSEAcol], recompute=recompute_ssgsea)
        print("saving files")
        enrichments.to_csv(save_output + "gene_sets_all.csv")
    saveFiles(files, save_output)
    print("done")

    return (
        files,
        failed,
        samplesinset,
        renaming,
        lowqual,
        enrichments if compute_enrichment else None,
    )

async def postProcessStranded(
    refworkspace,
    samplesetname,
    failed,
    save_output="",
    ensemblserver=constants.ENSEMBL_SERVER_V,
    samplesetToLoad="all",
    geneLevelCols=constants.RSEMFILENAME_GENE,
    trancriptLevelCols=constants.RSEMFILENAME_TRANSCRIPTS,
    renamingFunc=None,
    samplesinset=[],
    rsemfilelocs=None,
    useCache=False,
    dropNonMatching=False,
):
    """postprocess a set of aggregated Expression table from RSEM in the CCLE way

    (usually using the aggregate_RSEM terra worklow)

    Args:
        refworkspace (str): terra workspace where the ref data is stored
        sampleset (str, optional): sampleset where the red data is stored. Defaults to 'all'.
        save_output (str, optional): whether to save our data. Defaults to "".
        doCleanup (bool, optional): whether to clean the Terra workspaces from their unused output and lo. Defaults to True.
        colstoclean (list, optional): the columns to clean in the terra workspace. Defaults to [].
        ensemblserver (str, optional): ensembl server biomart version . Defaults to constants.ENSEMBL_SERVER_V.
        todrop (list, optional): if some samples have to be dropped whatever happens. Defaults to [].
        priority (list, optional): if some samples have to not be dropped when failing QC . Defaults to [].
        useCache (bool, optional): whether to cache the ensembl server data. Defaults to False.
        samplesetToLoad (str, optional): the sampleset to load in the terra workspace. Defaults to "all".
        geneLevelCols (list, optional): the columns that contain the gene level
        expression data in the workspace. Defaults to constants.RSEMFILENAME_GENE.
        trancriptLevelCols (list, optional): the columns that contain the transcript
        level expression data in the workspacce. Defaults to constants.RSEMFILENAME_TRANSCRIPTS.
        ssGSEAcol (str, optional): the rna file on which to compute ssGSEA. Defaults to "genes_tpm".
        rsemfilelocs (pd.DataFrame, optional): locations of RSEM output files if refworkspace is not provided (bypass interaction with terra)
        samplesinset (list[str], optional): list of samples in the sampleset if refworkspace is not provided (bypass interaction with terra)
        rnaqclocs (dict(str:list[str]), optional): dict(sample_id:list[QC_filepaths]) of rna qc file locations if refworkspace is not provided (bypass interaction with terra)
        renamingFunc (function, optional): the function to use to rename the sample columns
        (takes colnames and todrop as input, outputs a renaming dict). Defaults to None.
        compute_enrichment (bool, optional): do SSgSEA or not. Defaults to True.
        dropNonMatching (bool, optional): whether to drop the non matching genes
        between entrez and ensembl. Defaults to False.
        recompute_ssgsea (bool, optional): whether to recompute ssGSEA or not. Defaults to True.
    """
    print("postprocessing stranded samples")
    if not samplesetToLoad:
        samplesetToLoad = samplesetname
    refwm = dm.WorkspaceManager(refworkspace)
    # if save_output:
    #     terra.saveWorkspace(refworkspace, save_output + "terra/")

    samplesinset = [
        i["entityName"]
        for i in refwm.get_entities("sample_set").loc[samplesetname].samples
    ]

    print("generating gene names")

    mybiomart = h.generateGeneNames(ensemble_server=ensemblserver, useCache=useCache)
    # creating renaming index, keeping top name first
    gene_rename = {}
    for _, i in mybiomart.iterrows():
        if i.ensembl_gene_id not in gene_rename:
            gene_rename.update(
                {i.ensembl_gene_id: i.hgnc_symbol + " (" + i.ensembl_gene_id + ")"}
            )
    protcod_rename = {}
    for _, i in mybiomart[
        (~mybiomart.entrezgene_id.isna()) & (mybiomart.gene_biotype == "protein_coding")
    ].iterrows():
        if i.ensembl_gene_id not in protcod_rename:
            protcod_rename.update(
                {
                    i.ensembl_gene_id: i.hgnc_symbol
                    + " ("
                    + str(int(i.entrezgene_id))
                    + ")"
                }
            )

    print("loading files")
    files, _ = loadFromRSEMaggregate(
        refworkspace,
        todrop=failed,
        filenames=trancriptLevelCols + geneLevelCols,
        sampleset=samplesetToLoad,
        renamingFunc=renamingFunc,
        rsemfilelocs=rsemfilelocs,
        folder=save_output,
    )
    print("renaming files")
    # gene level
    if len(geneLevelCols) > 0:
        files = extractProtCod(
            files,
            mybiomart[mybiomart.gene_biotype == "protein_coding"],
            protcod_rename,
            dropNonMatching=dropNonMatching,
            filenames=geneLevelCols,
        )
        files = subsetGenes(
            files,
            gene_rename,
            filenames=geneLevelCols,
            index_id="gene_id",
            drop="transcript_id",
        )
    if len(trancriptLevelCols) > 0:
        files = subsetGenes(
            files,
            gene_rename,
            filenames=trancriptLevelCols,
            drop="gene_id",
            index_id="transcript_id",
        )
    saveFiles(files, save_output)
    print("done")

    return files
