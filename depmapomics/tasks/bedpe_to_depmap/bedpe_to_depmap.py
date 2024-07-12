import pandas as pd
import gzip
import argparse

COLS_TO_KEEP = ["CHROM_A",
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
                "vep_Consequence_A",
                "vep_IMPACT_A",
                "vep_SYMBOL_A",
                "vep_Gene_A",
                "vep_BIOTYPE_A",
                "vep_SV_overlap_name_A",
                "vep_SV_overlap_AF_A",
                "CHROM_B",
                "START_B",
                "END_B",
                "ID",
                "STRAND_B",
                "TYPE",
                "FILTER",
                "REF_B",
                "ALT_B",
                "SVLEN_B",
                "MATEID_B",
                "SVINSLEN_B",
                "BND_DEPTH_B",
                "MATE_BND_DEPTH_B",
                "vep_Consequence_B",
                "vep_IMPACT_B",
                "vep_SYMBOL_B",
                "vep_Gene_B",
                "vep_BIOTYPE_B",
                "vep_SV_overlap_name_B",
                "vep_SV_overlap_AF_B",
                "PR",
                "SR",
                "Rescue"]

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_filename")
    parser.add_argument("sample_name")

    args = parser.parse_args()

    vcf_filename = args.vcf_filename
    sample_name = args.sample_name

    print("expanding INFO fields")
    bedpe_df = bedpe_to_df(vcf_filename)
    bedpe_df.to_csv(sample_name + ".svs.expanded.bedpe", index=False)

    print("filtering & rescuing")
    df_filtered = filter_svs(bedpe_df)
    df_filtered = correct_bnd_gene(df_filtered)
    df_filtered.to_csv(sample_name + ".svs.expanded.filtered.bedpe", index=False)

def bedpe_to_df(
    path,
    additional_cols=[],
):
    """
    transforms a bedpe file into a dataframe file as best as it can

    Args:
    -----
        path: str filepath to the vcf file
        additional_filters: list[str] additional values added by the filtering tool looks for PASS, base_qual,
            clustered_events, fragment, germline, haplotype, map_qual, multiallelic,
            panel_of_normals, position, slippage, strand_bias, weak_evidence
        additional_cols: list[str] of additional colnames in the vcf already looks for 'DB',
            'SOMATIC', 'GERMLINE', "OVERLAP", "IN_PON", "STR", "ReverseComplementedAlleles"
        drop_null: bool if a column appears to be fully empty, will drop it
        force_keep: list[str] columns to force keep even if they are empty
        cols_to_drop: list[str] columns to drop even if they are not empty

    Returns:
    --------
      a dataframe fo the vcf
      a dict associating each column with its description (gathered from the vcf header)
      a list of the columns that have been dropped
    """
    uniqueargs = [
        "IMPRECISE",
    ] + additional_cols

    VEP_CSQ_DESC = "Consequence annotations from Ensembl VEP."

    dropped_cols = []

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
                        #print("parsing VEP CSQ")
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
    
    vep_fields_a = [k+"_A" for k, v in description.items() if VEP_CSQ_DESC in v]
    vep_fields_b = [k+"_B" for k, v in description.items() if VEP_CSQ_DESC in v]
    fields_a = {k+"_A": [] for k, _ in description.items()}
    fields_b = {k+"_B": [] for k, _ in description.items()}
    fields_combined = {}
    for suffix, side, fields, vep_fields in [("_A", "INFO_A", fields_a, vep_fields_a), ("_B", "INFO_B", fields_b, vep_fields_b)]:
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
                            res.update({k+suffix: annot})
                    else:
                        raise ValueError("unknown argument: " + annot)
                for k in list(fields.keys()):
                    fields[k].append(res.get(k, None))
        except ValueError:
            print(annot)
            raise ValueError("unknown field")
        fields_combined.update(fields)

    data = pd.concat(
        [data.drop(columns=["INFO_A", "INFO_B"]), pd.DataFrame(data=fields_combined, index=data.index)], axis=1
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
        print(sample)
        print(sorting)
        data = pd.concat(
            [
                data.drop(columns=sample),
                pd.DataFrame(data=res, columns=sorting, index=data.index),
            ],
            axis=1,
        )

    return data



def filter_svs(df, 
               sv_gnomad_cutoff = 0.001, 
               cosmic_fusion_pairs="gs://cds-cosmic/cosmic_fusion_gene_pairs_v100.csv",
               oncogene_list="/home/oncogene_oncokb.txt",
               ts_list="/home/tumor_suppressor_oncokb.txt",
               large_sv_size = 20000,
               cols_to_keep=COLS_TO_KEEP,
              ):
    
    # drop variants shorter than 50
    df = df[(df["SVLEN_A"].isna()) | 
       (df["SVLEN_A"].astype(float).astype('Int64') <= -50) | 
       (df["SVLEN_A"].astype(float).astype('Int64') >= 50)]
    
    with open(oncogene_list) as f:
        oncogenes = [val[:-1] for val in f.readlines()]
    
    with open(ts_list) as f:
        tumor_suppressors = [val[:-1] for val in f.readlines()]

    oncogenes_and_ts = set(oncogenes + tumor_suppressors)
    
    cosmic = pd.read_csv(cosmic_fusion_pairs)
    cosmic_pairs = list(zip(cosmic['Gene_A'], cosmic['Gene_B']))
    cosmic_pairs_sorted = set([tuple(sorted(elem)) for elem in cosmic_pairs])
    
    df["Rescue"] = False
    
    # rescue large SVs
    df.loc[(~df["SVLEN_A"].isna()) & 
       ((df["SVLEN_A"].astype(float).astype('Int64') <= -large_sv_size) | 
       (df["SVLEN_A"].astype(float).astype('Int64') >= large_sv_size)), "Rescue"] = True
    
    # rescue breakpoints that fall on oncogenes or tumor suppressors
    df.loc[(df["vep_SYMBOL_A"].isin(oncogenes_and_ts)) | (df["vep_SYMBOL_B"].isin(oncogenes_and_ts)), "Rescue"] = True
    
    # rescue gene pairs in cosmic
    df['vep_SYMBOL_A'] = df['vep_SYMBOL_A'].fillna("")
    df['vep_SYMBOL_B'] = df['vep_SYMBOL_B'].fillna("")
    df["pair"] = [tuple(sorted(elem)) for elem in list(zip(df['vep_SYMBOL_A'], df['vep_SYMBOL_B']))]
    df.loc[df["pair"].isin(cosmic_pairs_sorted), "Rescue"] = True
    
    # gnomad AF parsing
    df["max_af"] = df["vep_SV_overlap_AF_A"].fillna("").str.split("&").apply(lambda x: max([float(e) if e != "" else 0 for e in x]))
    
    # filter while keeping rescues
    df = df[(df["Rescue"] == True) |
           ((df["max_af"] < sv_gnomad_cutoff) &
           ((df["vep_BIOTYPE_A"].str.startswith("protein_coding")) |
           (df["vep_IMPACT_A"] == "HIGH") |
           ((df["vep_Consequence_A"].str.contains("splice")) & (df["vep_IMPACT_A"] == "MODERATE"))))]
    
    return df[cols_to_keep]

def correct_bnd_gene(bedpe):
    # if either breakend is intergenic OR only feature_truncation, remove gene label
    bedpe.loc[((bedpe.vep_Consequence_A.str.contains("intergenic_variant")) | (bedpe.vep_Consequence_A == "feature_truncation")) & (bedpe.TYPE == "BND"), ["vep_SYMBOL_A", "vep_Gene_A", "vep_BIOTYPE_A"]] = ""
    bedpe.loc[((bedpe.vep_Consequence_B.str.contains("intergenic_variant")) | (bedpe.vep_Consequence_B == "feature_truncation")) & (bedpe.TYPE == "BND"), ["vep_SYMBOL_B", "vep_Gene_B", "vep_BIOTYPE_B"]] = ""

    return bedpe
    

if __name__ == "__main__":
    main()