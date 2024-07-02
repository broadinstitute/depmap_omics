import pandas as pd
from mgenepy.utils import helper as h


def bedpe_to_df(
    path,
    additional_cols=[],
    **kwargs,
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
    
    tokeep = ['SVLEN', 'CIGAR', 'MATEID', 'EVENT', 'HOMLEN', 'HOMSEQ', 'SVINSLEN', 'SVINSSEQ', 
              'LEFT_SVINSSEQ', 'RIGHT_SVINSSEQ', 'BND_DEPTH', 'MATE_BND_DEPTH', 
              'vep_Consequence', 'vep_IMPACT', 'vep_SYMBOL', 'vep_Gene', 'vep_Feature_type', 'vep_Feature', 
              'vep_HGVSc', 'vep_HGVSp', 'vep_MANE_SELECT', 
              'vep_SIFT', 'vep_PolyPhen', 'vep_gnomAD_SV', 'vep_gnomAD_SV_AF', 'PR', 'SR']

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
        "skiprows": nrows_toskip + kwargs.get("skiprows", 0),
    }
    data = pd.read_csv(path, **{**kwargs, **csvkwargs})
    
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
                        if annot == "REF":
                            res.update({"REF_FLAG": True})
                        else:
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

    return data, description, dropped_cols


def filter_svs(df, 
               sv_gnomad_cutoff = 0.001, 
               cosmic_fusion_pairs="gs://cds-cosmic/cosmic_fusion_gene_pairs_v100.csv",
               oncogene_list="/home/xiaomeng/bin/depmap_omics/depmapomics/tasks/oncokb_dm/data/oncogene_oncokb.txt",
               ts_list="/home/xiaomeng/bin/depmap_omics/depmapomics/tasks/oncokb_dm/data/tumor_suppressor_oncokb.txt",
               large_sv_size = 20000,
              ):
    
    # drop variants shorter than 50
    df = df[(df["SVLEN_A"].isna()) | 
       (df["SVLEN_A"].astype(float).astype('Int64') <= -50) | 
       (df["SVLEN_A"].astype(float).astype('Int64') >= 50)]
    
    oncogenes = h.fileToList(oncogene_list)
    tumor_suppressors = h.fileToList(ts_list)
    oncogenes_and_ts = set(oncogenes + tumor_suppressors)
    
    cosmic = pd.read_csv(cosmic_fusion_pairs)
    cosmic_pairs = list(zip(cosmic['Gene_A'], cosmic['Gene_B']))
    cosmic_pairs_sorted = set([tuple(sorted(elem)) for elem in cosmic_pairs])
    
    df["rescue"] = False
    
    # rescue large SVs
    df.loc[(~df["SVLEN_A"].isna()) & 
       ((df["SVLEN_A"].astype(float).astype('Int64') <= -large_sv_size) | 
       (df["SVLEN_A"].astype(float).astype('Int64') >= large_sv_size)), "rescue"] = True
    
    # rescue breakpoints that fall on oncogenes or tumor suppressors
    df.loc[(df["vep_SYMBOL_A"].isin(oncogenes_and_ts)) | (df["vep_SYMBOL_B"].isin(oncogenes_and_ts)), "rescue"] = True
    
    # rescue gene pairs in cosmic
    df['vep_SYMBOL_A'] = df['vep_SYMBOL_A'].fillna("")
    df['vep_SYMBOL_B'] = df['vep_SYMBOL_B'].fillna("")
    df["pair"] = [tuple(sorted(elem)) for elem in list(zip(df['vep_SYMBOL_A'], df['vep_SYMBOL_B']))]
    df.loc[df["pair"].isin(cosmic_pairs_sorted), "rescue"] = True
    
    # gnomad AF parsing
    df["max_af"] = df["vep_gnomAD_SV_AF_A"].fillna("").str.split("&").apply(lambda x: max([float(e) if e != "" else 0 for e in x]))
    
    df = df[(df["rescue"] == True) |
           ((df["max_af"] < sv_gnomad_cutoff) &
           ((df["vep_BIOTYPE_A"].str.startswith("protein_coding")) |
           (df["vep_IMPACT_A"] == "HIGH")))]
    
    return df