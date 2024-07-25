import pandas as pd
import numpy as np
import seaborn as sns
from mgenepy.utils import helper as h
import gzip


def manageGapsInSegments(
    segtocp, Chromosome="Chromosome", End="End", Start="Start", cyto=None
):
    """
    extends the ends of segments in a segment file from GATK so as to remove all gaps ove the genome (works with multiple sample file)

    Args:
    ----
      segtocp: dataframe of segments from GATK CN pipeline
      Chromosome: str the value for the Chromosome columns
      End: str the value for the End columns
      Start: str the value for the Start columns
      cyto: dataframe with chrom;end; columns giving the size of each chromosome (else puts last segment to 1000000000)
    """
    prevchr = ""
    prevend = 0
    count = 0
    l = []
    segments = segtocp.copy()
    le = len(segments)
    for k, val in segments.iterrows():
        h.showcount(count, le)
        count += 1
        if val[Chromosome] != prevchr:  # we changed chromosome
            # we extend the previous segment (last of the prev chrom) to.. way enough
            if len(l) > 0:
                l[-1][2] = (
                    1000000000
                    if cyto is None
                    else cyto[cyto["chrom"] == prevchr]["end"].values[-1]
                )
            # we extend the first segment to 0
            l.append([val[Chromosome], 0, val[End]])
        else:
            if val[Start] > prevend + 1:  # we have a gap in the same chrom
                sizeofgap = val[Start] - prevend
                # we add to the previous one half of the gap
                l[-1][2] += (
                    int(sizeofgap / 2) if sizeofgap % 2 == 0 else int(sizeofgap / 2) + 1
                )
                # the rest to the other
                l.append([val[Chromosome], val[Start] - int(sizeofgap / 2), val[End]])
            elif val[Start] < prevend:  # this should never happen
                # import pdb; pdb.set_trace()
                raise ValueError("start comes after end")
            else:
                l.append([val[Chromosome], val[Start], val[End]])
        prevchr = val[Chromosome]
        prevend = val[End]
    # we extend the last one
    l[-1][2] = (
        1000000000 if cyto is None else cyto[cyto["chrom"] == prevchr]["end"].values[-1]
    )
    segments[[Chromosome, Start, End]] = l
    return segments.reset_index(drop=True)


def toGeneMatrix(
    segments,
    gene_mapping,
    style="weighted",
    missingchrom=["Y"],
    gene_names_col="gene_name",
    value_colname="Segment_Mean",
):
    """
    makes a geneXsample matrix from segment level copy number (works with multiple sample file)

    Args:
    ----
      style: str one of "weighted","mean","closest"
      segments: dataframe of segments containing: [Chromosome, Segment_Mean, Chromosome, start, end] columns
      gene_mapping: dataframe with symbol, ensembl_id columns for each gene
      missingchrom: list[str] chromosomes not to look into

    Returns:
    -------
      pd.dataframe: the matrix
    """
    samples = list(set(segments.DepMap_ID))
    data = np.zeros((len(samples), len(gene_mapping)))
    for i, sample in enumerate(samples):
        segs = segments[segments.DepMap_ID == sample][
            ["Chromosome", "Start", "End", value_colname]
        ].values
        hasmissing = set(missingchrom) - set(segs[:, 0])
        j = 0
        h.showcount(i, len(samples))
        for k, gene in enumerate(gene_mapping[["Chromosome", "start", "end"]].values):
            # print(i,j)
            if gene[0] in hasmissing:
                data[i, k] = np.nan
                continue
            try:
                while gene[0] != segs[j][0] or gene[1] >= segs[j][2]:
                    # print("went beyong",gene, segs[j])
                    j += 1
                # some genes are within other genes, we need to go back in the list of segment in that case
            except:
                raise ValueError("forgot to sort one of the DF?")
            while gene[1] < segs[j][1]:
                j -= 1
                # print("decrease gene",gene)
            # we are entirely within the segment
            c = 1
            if gene[2] <= segs[j][2]:
                data[i, k] = segs[j][3]
            else:
                # how much of the gene is covered by the segment
                coef = (segs[j][2] - gene[1]) / (gene[2] - gene[1])
                # print('coef',coef)
                val = segs[j][3] * coef if style == "weighted" else segs[j][3]
                end = segs[j][2]
                # until the end of a segments goes beyond the end of the gene (say if we have X segments within the gene)
                while end < gene[2]:
                    # pdb.set_trace()
                    j += 1
                    c += 1
                    nextend = segs[j][2] if segs[j][2] < gene[2] else gene[2]
                    # here, end (of prevsegment) is the next segment's start
                    ncoef = (nextend - end) / (gene[2] - gene[1])
                    # print('multi',gene, ncoef)
                    if style == "closest":
                        if ncoef > coef:
                            val = segs[j][3]
                        else:
                            # we switch it back (see line 894)
                            ncoef = coef
                    else:
                        val += segs[j][3] * ncoef if style == "weighted" else segs[j][3]
                    end = segs[j][2]
                    coef = ncoef
                data[i, k] = val if style == "weighted" else val / c
    return pd.DataFrame(data=data, index=samples, columns=gene_mapping[gene_names_col])


def checkAmountOfSegments(segmentcn, thresh=850, samplecol="DepMap_ID"):
    """
    if there is too many segments, something might be wrong (works with multiple sample file)

    will compute the number of segments for each samples from a df of segments from RSEM

    Args:
    ----
      segmentcn: segment dataframe
      thresh: max ok amount
    """
    failed = []
    celllines = set(segmentcn[samplecol].tolist())
    amounts = []
    for cellline in celllines:
        val = segmentcn[segmentcn[samplecol] == cellline].shape[0]
        amounts.append(val)
        if val > thresh:
            failed.append(cellline)
            print(cellline, val)
    sns.kdeplot(amounts)
    return failed


def vcf_to_df(
    path,
    additional_cols=[],
    additional_filters=[],
    parse_filter=False,
    drop_null=False,
    force_keep=[],
    cols_to_drop=[
        "clinvar_vcf_mc",
        "oreganno_build",
        "gt",
        "ad",
        "af",
        "dp",
        "f1r2",
        "f2r1",
        "fad",
        "sb",
        "pid",
    ],
    **kwargs,
):
    """
    transforms a vcf file into a dataframe file as best as it can

    Args:
    -----
        path: str filepath to the vcf file
        additional_filters: list[str] additional values added by the filtering tool looks for PASS, base_qual,
            clustered_events, fragment, germline, haplotype, map_qual, multiallelic,
            panel_of_normals, position, slippage, strand_bias, weak_evidence
        additional_cols: list[str] of additional colnames in the vcf already looks for 'DB',
            'SOMATIC', 'GERMLINE', "OVERLAP", "IN_PON", "STR", "ReverseComplementedAlleles"
        parse_filter: bool if true, will parse the filter field and add it to the dataframe
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
        "DB",
        "SOMATIC",
        "GERMLINE",
        "OVERLAP",
        "IN_PON",
        "STR",
        "ReverseComplementedAlleles",
    ] + additional_cols

    filters = [
        "PASS",
        "base_qual",
        "clustered_events",
        "fragment",
        "germline",
        "haplotype",
        "map_qual",
        "multiallelic",
        "panel_of_normals",
        "position",
        "slippage",
        "strand_bias",
        "weak_evidence",
    ] + additional_filters

    FUNCO_DESC = "Functional annotation from the Funcotator tool."
    VEP_CSQ_DESC = "Consequence annotations from Ensembl VEP."
    SNPEFF_ANN_DESC = "Functional annotations from SnpEff"
    LOF_DESC = "Predicted loss of function effects from VEP"

    dropped_cols = []

    def read_comments(f):
        description = {}
        colnames = []
        rows = 0
        for l in f:
            l = l.decode("utf-8") if type(l) is not str else l
            if l.startswith("##"):
                rows += 1
                if "FORMAT" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    desc = l.split("Description=")[1][:-2]
                    description.update({res: desc})
                if "INFO" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    if res == "FUNCOTATION":
                        #print("parsing funcotator special")
                        for val in l.split("Description=")[1][:-2].split("|"):
                            val = val.split("Funcotation fields are: ")[-1]
                            description.update({val: FUNCO_DESC})
                    elif res == "ANN":
                        #print("parsing funcitonal annotation from SnpEff")
                        l = l.replace(" / ", "_").replace(".", "_")
                        for val in l.split("Description=")[1][:-5].split(" | "):
                            val = "snpeff_" + val.split("Functional annotations: '")[-1]
                            description.update({val: SNPEFF_ANN_DESC})
                    elif res == "LOF":
                        #print("parsing predicted LOF status from SnpEff")
                        for val in l.split("Description=")[1][:-4].split(" | "):
                            val = "lof_" + val.split("Format: '")[-1]
                            description.update({val: LOF_DESC})
                    elif res == "CSQ":
                        #print("parsing VEP CSQ")
                        for val in l.split("Description=")[1][:-3].split("|"):
                            val = "vep_" + val.split("Format: ")[-1]
                            description.update({val: VEP_CSQ_DESC})
                    elif res == "REF":
                        description.update(
                            {"REF_FLAG": l.split("Description=")[1][:-2]}
                        )
                    else:
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
    #print(description)
    funco_fields = [k for k, v in description.items() if FUNCO_DESC in v]
    vep_fields = [k for k, v in description.items() if VEP_CSQ_DESC in v]
    snpeff_fields = [k for k, v in description.items() if SNPEFF_ANN_DESC in v]
    lof_fields = [k for k, v in description.items() if LOF_DESC in v]
    fields = {k: [] for k, _ in description.items()}
    try:
        for j, info in enumerate(data["INFO"].str.split(";").values.tolist()):
            res = {}
            # show speed
            if j % 10_000 == 0:
                print(j, end="\r")
            for annot in info:
                if annot in uniqueargs:
                    if annot == "REF":
                        res.update({"REF_FLAG": True})
                    else:
                        res.update({annot: True})
                elif "=" in annot:
                    # taking care of the funcotator special fields
                    if "FUNCOTATION=" in annot:
                        # for multi allelic site:
                        annot = annot.replace("FUNCOTATION=", "")[1:-1]
                        res.update({name: [] for name in funco_fields})
                        for site in annot.split("],["):
                            if "]#[" in site:
                                site = site.split("]#[")[0]
                            site = (
                                site.replace("_%7C_", " ")
                                .replace("_%20_", " ")
                                .replace("_%2C_", ",")
                                .replace("_%3D_", "=")
                                .split("|")
                            )
                            for i, sub_annot in enumerate(site):
                                res[funco_fields[i]].append(sub_annot)
                        for k in funco_fields:
                            res[k] = ",".join(res[k])
                    elif "ANN=" in annot:
                        annot = annot.replace("ANN=", "")
                        res.update({name: [] for name in snpeff_fields})
                        for site in annot.split(","):
                            for i, sub_annot in enumerate(site.split("|")):
                                res[snpeff_fields[i]].append(sub_annot)
                        for k in snpeff_fields:
                            if "".join(res[k]) != "":
                                res[k] = ",".join(res[k])
                            else:
                                res[k] = ""
                    elif "LOF=" in annot:
                        annot = annot.replace("LOF=(", "")
                        annot = annot.replace(")", "")
                        res.update({name: [] for name in lof_fields})
                        for site in annot.split(","):
                            for i, sub_annot in enumerate(site.split("|")):
                                res[lof_fields[i]].append(sub_annot)
                        for k in lof_fields:
                            if "".join(res[k]) != "":
                                res[k] = ",".join(res[k])
                            else:
                                res[k] = ""
                    elif "CSQ=" in annot:
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
                    else:
                        k, annot = annot.split("=")
                        res.update({k: annot})
                else:
                    raise ValueError("unknown argument: " + annot)
            for k in list(fields.keys()):
                fields[k].append(res.get(k, None))
    except ValueError:
        print(annot)
        raise ValueError("unknown field")

    data = pd.concat(
        [data.drop(columns="INFO"), pd.DataFrame(data=fields, index=data.index)], axis=1
    )
    if drop_null:
        to_drop = []
        for f in funco_fields:
            # drop columns that have the same value across all rows
            uniq = data[f].unique()
            if len(uniq) == 1 and f.lower() not in force_keep:
                to_drop.append(f)
                continue
            elif len(uniq) < 10:
                # checking multi allelic stuff
                multi = []
                for v in uniq:
                    multi += v.split(",")
                if len(set(multi)) == 1 and f.lower() not in force_keep:
                    to_drop.append(f)
        print("dropping uninformative columns:", to_drop)
        data = data.drop(columns=to_drop)
        dropped_cols += to_drop

    data.columns = [i.lower() for i in data.columns]
    samples = [i.lower() for i in colnames[9:]]
    if samples == ['normal', 'tumor']:
        pass
        data = data.drop('normal',axis = 1)
        colnames = [c for c in colnames if c != 'NORMAL']
        samples = [i.lower() for i in colnames[9:]]
    print("\nthe samples are:", samples)
    sorting = data["format"][0].split(":")
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

    # subsetting filters
    if parse_filter:
        data[filters] = False
        for f in filters:
            data.loc[data["filter"].str.contains(f), f] = True
        data = data.drop(columns="filter")
        dropped_cols.append("filter")

    # cleaning empty cols
    data = data.drop(columns="format")
    dropped_cols.append("format")

    todrop = []
    for val in cols_to_drop:
        if val in data.columns.tolist():
            todrop.append(val)
    data = data.drop(columns=todrop)

    if drop_null:
        empty = data.columns[data.isna().sum() == len(data)].tolist()
        empty = list(set(empty) - set(force_keep))
        print("dropping empty columns:", empty)
        data = data.drop(columns=empty)
        dropped_cols += empty

    # weird bug sometimes
    if "SB_1" in data.columns.tolist():
        loc = ~data.SB_1.isna()
        data.loc[loc, "PGT"] = data.loc[loc, "SB"]
        data.loc[loc, "SB"] = data.loc[loc, "SB_1_2_3"]
        data = data.drop(columns=["SB_1", "SB_1_2_3"])
        data = data.rename(columns={"SB_1_2": "PS", "SB_1": "PID"})
    elif 'SB' in data.columns.tolist():
        loc = data.SB.isna()
        data.loc[loc, "SB"] = data.loc[loc, "PGT"]
        data.loc[loc, "PGT"] = ""
    # sorting out issue with
    return data, description, dropped_cols
