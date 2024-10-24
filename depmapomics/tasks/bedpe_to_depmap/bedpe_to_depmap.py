import pandas as pd
import gzip
import argparse
import itertools

COLS_TO_KEEP = [
    "CHROM_A",
    "START_A",
    "END_A",
    "ID",
    "STRAND_A",
    "TYPE",
    "FILTER",
    "REF_A",
    "ALT_A",
    "SVLEN_A",
    "MATEID_A",
    "SVINSLEN_A",
    "BND_DEPTH_A",
    "MATE_BND_DEPTH_A",
    "SYMBOL_A",
    "GENEID_A",
    "vep_SV_overlap_name_A",
    "vep_SV_overlap_AF_A",
    "CHROM_B",
    "START_B",
    "END_B",
    "STRAND_B",
    "REF_B",
    "ALT_B",
    "SVLEN_B",
    "MATEID_B",
    "SVINSLEN_B",
    "BND_DEPTH_B",
    "MATE_BND_DEPTH_B",
    "SYMBOL_B",
    "GENEID_B",
    "vep_SV_overlap_name_B",
    "vep_SV_overlap_AF_B",
    "DEL_SYMBOLS",
    "DUP_SYMBOLS",
    "PR",
    "SR",
    "Rescue",
]


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("bedpe_filename")
    parser.add_argument("gene_annotation_filename")
    parser.add_argument("annotated_overlap_del")
    parser.add_argument("annotated_overlap_dup")
    parser.add_argument("cosmic_fusion_pairs")
    parser.add_argument("sample_name")

    args = parser.parse_args()

    bedpe_filename = args.bedpe_filename
    sample_name = args.sample_name
    gene_annotation_filename = args.gene_annotation_filename
    annotated_overlap_del = args.annotated_overlap_del
    annotated_overlap_dup = args.annotated_overlap_dup
    cosmic_fusion_pairs = args.cosmic_fusion_pairs

    print("expanding INFO fields")
    bedpe_df = bedpe_to_df(bedpe_filename)

    print("reannotating gene symbols")
    bedpe_reannotated = reannotate_genes(
        bedpe_df, gene_annotation_filename, annotated_overlap_del, annotated_overlap_dup
    )
    bedpe_reannotated.to_csv(
        sample_name + ".svs.expanded.reannotated.bedpe", index=False, sep="\t"
    )

    print("filtering & rescuing")
    df_filtered = filter_svs(bedpe_reannotated, cosmic_fusion_pairs=cosmic_fusion_pairs)
    df_filtered.to_csv(
        sample_name + ".svs.expanded.reannotated.filtered.bedpe", index=False, sep="\t"
    )


def bedpe_to_df(
    path,
    additional_cols=[],
):
    """
    transforms a bedpe file into a dataframe file as best as it can
    largely borrowed from vcf_to_df in the same repo

    Args:
    -----
        path: str filepath to the vcf file
        additional_cols: list of flags that are boolean (present == True, absent == False)

    Returns:
    --------
      a dataframe for the bedpe
    """
    uniqueargs = [
        "IMPRECISE",
    ] + additional_cols

    VEP_CSQ_DESC = "Consequence annotations from Ensembl VEP."

    def read_comments(f):
        description = {}
        colnames = []
        rows = 0
        for l in f:
            l = l.decode("utf-8") if type(l) is not str else l
            if l.startswith("##"):
                rows += 1
                if "INFO" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    if res == "CSQ":
                        # print("parsing VEP CSQ")
                        for val in l.split("Description=")[1][:-3].split("|"):
                            val = "vep_" + val.split("Format: ")[-1]
                            description.update({val: VEP_CSQ_DESC})
                    elif res != "END":
                        desc = l.split("Description=")[1][:-2]
                        description.update({res: desc})
            elif l.startswith("#"):
                colnames = l[1:-1].split("\t")
                rows += 1
            else:
                break
        return description, colnames, rows

    if path.endswith(".gz"):
        with gzip.open(path, "r") as f:
            description, colnames, nrows_toskip = read_comments(f)
    else:
        with open(path, "r") as f:
            description, colnames, nrows_toskip = read_comments(f)

    colnames = [i for i in colnames]
    csvkwargs = {
        "sep": "\t",
        "index_col": False,
        "header": None,
        "names": colnames,
        "skiprows": nrows_toskip,
    }
    data = pd.read_csv(path, **csvkwargs)

    vep_fields_a = [k + "_A" for k, v in description.items() if VEP_CSQ_DESC in v]
    vep_fields_b = [k + "_B" for k, v in description.items() if VEP_CSQ_DESC in v]
    fields_a = {k + "_A": [] for k, _ in description.items()}
    fields_b = {k + "_B": [] for k, _ in description.items()}
    fields_combined = {}
    for suffix, side, fields, vep_fields in [
        ("_A", "INFO_A", fields_a, vep_fields_a),
        ("_B", "INFO_B", fields_b, vep_fields_b),
    ]:
        try:
            for j, info in enumerate(data[side].str.split(";").values.tolist()):
                res = {}
                # show speed
                if j % 10_000 == 0:
                    print(j, end="\r")
                for annot in info:
                    if annot == ".":
                        pass
                    elif annot in uniqueargs:
                        res.update({annot: True})
                    elif "=" in annot:
                        # taking care of the funcotator special fields
                        if "CSQ=" in annot:
                            annot = annot.replace("CSQ=", "")
                            res.update({name: [] for name in vep_fields})
                            for site in annot.split(","):
                                for i, sub_annot in enumerate(site.split("|")):
                                    res[vep_fields[i]].append(sub_annot)
                            for k in vep_fields:
                                if "".join(res[k]) != "":
                                    res[k] = ",".join(res[k])
                                else:
                                    res[k] = ""
                        elif "END" not in annot:
                            k, annot = annot.split("=")
                            res.update({k + suffix: annot})
                    else:
                        raise ValueError("unknown argument: " + annot)
                for k in list(fields.keys()):
                    fields[k].append(res.get(k, None))
        except ValueError:
            raise ValueError("unknown field")
        fields_combined.update(fields)

    data = pd.concat(
        [
            data.drop(columns=["INFO_A", "INFO_B"]),
            pd.DataFrame(data=fields_combined, index=data.index),
        ],
        axis=1,
    )

    samples = [i for i in colnames[21:]]
    sorting = data["FORMAT"][0].split(":")
    for sample in samples:
        res = data[sample].str.split(":").values.tolist()
        maxcols = max([len(v) for v in res])
        if maxcols - len(sorting) > 0:
            for i in range(maxcols - len(sorting)):
                sorting.append(sorting[-1] + "_" + str(i + 1))
        if len(samples) > 1:
            sorting = [sample + "_" + v for v in sorting]
        data = pd.concat(
            [
                data.drop(columns=sample),
                pd.DataFrame(data=res, columns=sorting, index=data.index),
            ],
            axis=1,
        )

    return data


def reannotate_genes(bedpe, annotation_path, del_annotation_path, dup_annotation_path):
    """since VEP can't reliably give the correct gene symbol annotation, redo it here

    Args:
        bedpe (df.DataFrame): SVs in bedpe format
        annotation_path (str): path to the gene reannotation file for breakpoints
        del_annotation_path (str): path to gene reannotation file for DELs only, considering the entire deleted intervals
        dup_annotation_path (str): path to gene reannotation file for DUPs only, considering the entire duplicated intervals

    Returns:
        merged (pd.DataFrame): SVs in bedpe format, with corrected gene annotation

    """
    # read in annotation files for all breakpoints, DELs, and DUPs separately
    gene_annotation = pd.read_csv(
        annotation_path,
        sep="\t",
        names=[
            "CHROM_A",
            "START_A",
            "END_A",
            "NAME_A",
            "GENE_A",
            "CHROM_B",
            "START_B",
            "END_B",
            "NAME_B",
            "GENE_B",
        ],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_A": "int64",
            "NAME_A": "string",
            "GENE_A": "string",
            "CHROM_B": "string",
            "START_B": "int64",
            "END_B": "int64",
            "NAME_B": "string",
            "GENE_B": "string",
        }
    )

    del_annotation = pd.read_csv(
        del_annotation_path,
        sep="\t",
        names=["CHROM_A", "START_A", "END_B", "NAME_A", "NAME_B", "DELGENES"],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_B": "int64",
            "NAME_A": "string",
            "NAME_B": "string",
            "DELGENES": "string",
        }
    )

    dup_annotation = pd.read_csv(
        dup_annotation_path,
        sep="\t",
        names=["CHROM_A", "START_A", "END_B", "NAME_A", "NAME_B", "DUPGENES"],
    ).astype(
        {
            "CHROM_A": "string",
            "START_A": "int64",
            "END_B": "int64",
            "NAME_A": "string",
            "NAME_B": "string",
            "DUPGENES": "string",
        }
    )

    # merge annotations to the SV table one by one
    merged = pd.merge(
        bedpe,
        gene_annotation[
            [
                "CHROM_A",
                "START_A",
                "END_A",
                "NAME_A",
                "GENE_A",
                "CHROM_B",
                "START_B",
                "END_B",
                "GENE_B",
            ]
        ],
        on=["CHROM_A", "START_A", "END_A", "NAME_A", "CHROM_B", "START_B", "END_B"],
        how="outer",
        indicator=True,
    )

    assert merged["_merge"].eq("both").all()

    merged = pd.merge(
        merged.drop(columns="_merge"),
        del_annotation[["CHROM_A", "START_A", "END_B", "NAME_A", "DELGENES"]],
        on=["CHROM_A", "START_A", "END_B", "NAME_A"],
        how="left",
    )

    merged = pd.merge(
        merged,
        dup_annotation[["CHROM_A", "START_A", "END_B", "NAME_A", "DUPGENES"]],
        on=["CHROM_A", "START_A", "END_B", "NAME_A"],
        how="left",
    )

    # in the case where there are multiple genes in one cell (HUGO1@ENSEMBL1; HUGO2@ENSEMBL2; ...)
    # split them in to separate comma-separated columns
    def split_multi(s):
        if pd.isna(s) or s == ".":
            return ".;."
        else:
            s = s.split(",")
            symbols = ", ".join([elem.split("@")[0] for elem in s])
            ids = ", ".join([elem.split("@")[1] for elem in s])
            return symbols + ";" + ids

    merged["GENE_A"] = merged["GENE_A"].map(split_multi)
    merged[["SYMBOL_A", "GENEID_A"]] = merged["GENE_A"].str.split(";", n=1, expand=True)
    merged["GENE_B"] = merged["GENE_B"].map(split_multi)
    merged[["SYMBOL_B", "GENEID_B"]] = merged["GENE_B"].str.split(";", n=1, expand=True)

    merged["DELGENES"] = merged["DELGENES"].map(split_multi)
    merged[["DEL_SYMBOLS", "DEL_GENEIDS"]] = merged["DELGENES"].str.split(
        ";", n=1, expand=True
    )
    merged["DUPGENES"] = merged["DUPGENES"].map(split_multi)
    merged[["DUP_SYMBOLS", "DUP_GENEIDS"]] = merged["DUPGENES"].str.split(
        ";", n=1, expand=True
    )

    # drop redundant columns
    merged = merged.drop(["GENE_A", "GENE_B", "DELGENES", "DUPGENES"], axis=1)

    return merged


def filter_svs(
    df,
    cosmic_fusion_pairs,
    sv_gnomad_cutoff=0.001,
    oncogene_list="/home/oncogene_oncokb.txt",
    ts_list="/home/tumor_suppressor_oncokb.txt",
    large_sv_size=1e9,
    cols_to_keep=COLS_TO_KEEP,
):
    """filter SVs in bedpe while rescuing important ones

    Args:
        df (pd.DataFrame): SVs in bedpe format
        cosmic_fusion_pairs (str): path to file containing cosmic fusion gene pairs
        sv_gnomad_cutoff (float): max gnomad allele frequency for an SV to be considered somatic
        oncogene_list (str): path to file containing list of oncogenes according to OncoKB
        ts_list (str): path to file containing list of tumor suppressor genes according to OncoKB
        large_sv_size (int): size of SVs beyond which is considered large and needs to be rescued
        cols_to_keep (list): list of columns to keep in the SV data

    Returns:
        df (pd.DataFrame): filtered SVs in bedpe format

    """

    df["SVLEN_A"] = df["SVLEN_A"].astype("Int64")

    # drop variants shorter than 50
    df = df.loc[df["SVLEN_A"].isna() | df["SVLEN_A"].abs().ge(50)]

    with open(oncogene_list) as f:
        oncogenes = [val[:-1] for val in f.readlines()]

    with open(ts_list) as f:
        tumor_suppressors = [val[:-1] for val in f.readlines()]

    oncogenes_and_ts = set(oncogenes + tumor_suppressors)

    cosmic = pd.read_csv(cosmic_fusion_pairs, dtype="string")
    cosmic_pairs = list(zip(cosmic["Gene_A"], cosmic["Gene_B"]))
    cosmic_pairs_sorted = set([tuple(sorted(elem)) for elem in cosmic_pairs])

    df["Rescue"] = False

    # rescue large SVs
    df.loc[df["SVLEN_A"].abs().ge(large_sv_size), "Rescue"] = True

    def onco_ts_overlap(s):
        l = s.split(", ")
        return len(set(l) & oncogenes_and_ts) > 0

    # rescue breakpoints that fall on oncogenes or tumor suppressors
    df["onco_ts_overlap_A"] = df["SYMBOL_A"].apply(onco_ts_overlap)
    df["onco_ts_overlap_B"] = df["SYMBOL_B"].apply(onco_ts_overlap)
    df.loc[df["onco_ts_overlap_A"] | df["onco_ts_overlap_A"], "Rescue"] = True

    # rescue gene pairs in cosmic
    def list_all_pairs(a, b):
        alist = a.split(", ")
        blist = b.split(", ")

        all_pairs = list(itertools.product(alist, blist))
        all_pairs = set([tuple(sorted(elem)) for elem in all_pairs])

        return len(all_pairs & cosmic_pairs_sorted) > 0

    df["pair_in_cosmic"] = df.apply(
        lambda row: list_all_pairs(row["SYMBOL_A"], row["SYMBOL_B"]), axis=1
    )
    df.loc[df["pair_in_cosmic"], "Rescue"] = True

    # gnomad AF parsing
    df["max_af"] = (
        df["vep_SV_overlap_AF_A"]
        .fillna("")
        .str.split("&")
        .apply(lambda x: max([float(e) if e != "" else 0 for e in x]))
    )

    # filter while keeping rescues
    df = df[df["Rescue"] | df["max_af"].lt(sv_gnomad_cutoff)]

    return df[cols_to_keep]


if __name__ == "__main__":
    main()
