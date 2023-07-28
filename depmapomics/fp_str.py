import io
import os
import pandas as pd
import gzip
import numpy as np
from depmapomics import constants
from depmap_omics_upload import tracker as track

from google.cloud import storage  # TODO: look for alternatives

STR_COL_RENAME_OMICS = {
    "D3S1358": "d3s1358",
    "TH01": "th01",
    "D21S11": "d21s11",
    "D18S51": "d18s51",
    "PentaE": "penta_e",
    "D5S818": "d5s818",
    "D13S317": "d13s317",
    "D7S820": "d7s820",
    "D16S539": "d16s539",
    "CSF1PO": "csf1po",
    "PentaD": "penta_d",
    "vWA": "vwa",
    "D8S1179": "d8s1179",
    "TPOX": "tpox",
    "FGA": "fga",
}


def read_vcf(path):
    """reads in a vcf in google cloud storage and converts it into a df"""
    storage_client = storage.Client()
    bucket = storage_client.bucket(path.split("/")[2])
    blob = bucket.blob("/".join(path.split("/")[3:]))
    if path.endswith(".gz"):
        data = io.BytesIO(blob.download_as_string())
        with gzip.open(data, "r") as f:
            lines = [l.decode("utf-8") for l in f if not l.startswith(b"#")]
    else:
        data = blob.download_as_string().decode("utf-8")
        f = data.split("\n")
        lines = [l + "\n" for l in f if not l.startswith("#")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "GT",
        ],
        sep="\t",
    )


def transformGB(row):
    if len(row) < 2:
        return np.nan
    else:
        split_gt = row.split(":")[1].split("|")
        return [int(i) for i in split_gt]


def gb2str(row):
    """transforms the GB field in hipstr's VCF to standard STR format"""
    if not isinstance(row["GB"], list):
        return "NA"
    else:
        strs = []
        for i in row["GB"]:
            integer = (row["REF"] * row["PERIOD"] + i) // row["PERIOD"]
            remainder = (row["REF"] * row["PERIOD"] + i) % row["PERIOD"]
            if remainder == 0:
                strs.append(str(int(integer)))
            else:
                strs.append(str(int(integer)) + "." + str(int(remainder)))
        return ", ".join(list(set(strs)))


def altAllele2str(row):
    """transforms the ALT field in gangstr's VCF to standard STR format"""
    if pd.isna(row["ALT"]) or row["ALT"] == ".":
        return "NA"
    else:
        strs = []
        for i in row["ALT"].split(","):
            l = len(i)
            integer = l // row["PERIOD"]
            remainder = l % row["PERIOD"]
            if remainder == 0:
                strs.append(str(int(integer)))
            else:
                strs.append(str(int(integer)) + "." + str(int(remainder)))
        return ", ".join(list(set(strs)))


def generateSTRRow(
    paths_df,
    method=constants.STR_METHOD,
    str_bed=constants.STR_BED,
    colname=constants.STR_COLNAME,
):
    """
    given a dataframe containing VCF locations for multiple samples,
    extract STR for all samples

    inputs:
        paths_df (pd.DataFrame): dataframe containing sample ids and their corresponding VCFs
        method (str): STR inferrence method. "hipstr" or "gangstr"
        str_bed (str): location of str bed file
        colname (str): name of column in paths_df that contains VCFs

    outputs:
        df containing STR for all samples
    """
    assert colname in set(paths_df.columns), (
        colname + " is not a column in input dataframe"
    )
    hg38_sites = pd.read_csv(
        str_bed,
        sep="\t",
        names=["CHROM", "START", "END", "PERIOD", "REF", "ID"],
    ).astype({"PERIOD": "int32"})
    str_rows = []
    for i, p in paths_df.iterrows():
        if pd.isna(p[colname]):
            print("no hipSTR vcf available for: " + i)
        else:
            df = read_vcf(p[colname])
            if method == "hipstr":
                df["GB"] = df.apply(lambda x: (transformGB(x["GT"])), axis=1)
                df = hg38_sites.merge(df[["ID", "GB"]], on="ID", how="left")
                df["STR"] = df.apply(lambda x: gb2str(x), axis=1)
            elif method == "gangstr":
                df = df.rename(columns={"POS": "START"})
                df = hg38_sites.merge(
                    df[["CHROM", "START", "ALT"]], on=["CHROM", "START"], how="left"
                )
                df["STR"] = df.apply(lambda x: altAllele2str(x), axis=1)
            # Maybe no need to convert to model IDs here yet?
            df["sample_id"] = i
            str_row = df.pivot(index="sample_id", columns="ID", values="STR")
            str_rows.append(str_row)
    return pd.concat(str_rows)


def computeTanabe(df1, idx1, df2, idx2, loci=constants.STR_LOCI_13):
    """compute tanabe similarity between two STR profiles"""
    match = 0
    total = 0
    for col in loci:
        # TODO: how to best handle NAs?
        a1 = set(df1.loc[idx1, col].split(", "))
        a2 = set(df2.loc[idx2, col].split(", "))
        if a1 != set(["NA"]) & a2 != set(["NA"]):
            match += len(set(a1) & set(a2))
            total += len(set(a1)) + len(set(a2))
    return 2 * match / total


def makeScoreMatrixDatabase(df_seqid, df_achid, loci=constants.STR_LOCI_13):
    """compute a match score matrix between STR profiles in two dataframes"""
    mytracker = track.SampleTracker()
    seq_table = mytracker.add_model_cols_to_seqtable(cols=[constants.MODEL_TABLE_INDEX])
    valid_achids = list(set(df_achid.index) - set([np.nan]))
    scoremat = pd.DataFrame(
        columns=valid_achids + [constants.MODEL_TABLE_INDEX], index=(df_seqid.index)
    )
    for i in df_seqid.index:
        scoremat.loc[i, constants.MODEL_TABLE_INDEX] = seq_table.loc[
            i, constants.MODEL_TABLE_INDEX
        ]
        for j in valid_achids:
            scoremat.loc[i, j] = computeTanabe(df_seqid, i, df_achid, j, loci=loci)
    return scoremat


def strCheck(
    paths_df,
    ref_df,
    method=constants.STR_METHOD,
    str_bed=constants.STR_BED,
    colname=constants.STR_COLNAME,
    loci=constants.STR_LOCI_13,
):
    """given a df containing STR VCF locations and a df containing reference STR profiles,
    generate a tanabe score matrix"""
    inferred_str = generateSTRRow(
        paths_df, method=method, str_bed=str_bed, colname=colname
    )
    score_table = makeScoreMatrixDatabase(inferred_str, ref_df, loci=loci)

    return score_table
