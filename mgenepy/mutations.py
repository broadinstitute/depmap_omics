import pandas as pd
import numpy as np
import seaborn as sns
from mgenepy.utils import helper as h


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


def vcf_to_df(path, hasfilter=False, samples=["sample"], additional_cols=[]):
    """
    transforms a vcf file into a dataframe file as best as it can

    Args:
    -----
      path: str filepath to the vcf file
      hasfilter: bool whether or not the vcf has a filter column
      samples: list[str] colnames of the sample names.
      additional_cols: list[str] of additional colnames in the vcf already looks for 'DB', 'SOMATIC', 'GERMLINE', "OVERLAP", "IN_PON", "STR", "ReverseComplementedAlleles"

    Returns:
    --------
      a dataframe fo the vcf
      a dict associating each column with its description (gathered from the vcf header)
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

    def read_comments(f):
        fields = {}
        description = {}
        c = 0
        headerrow = 0
        for l in f:
            l = l.decode("utf-8") if type(l) is not str else l
            if l.startswith("##"):
                if "FORMAT" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    desc = l.split("Description=")[1][:-2]
                    description.update({res: desc})
                if "INFO" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    desc = l.split("Description=")[1][:-2]
                    description.update({res: desc})
                    fields.update({res: []})
                c += 1
            elif l.startswith("#CHROM"):
                headerrow = c + 1
            else:
                c += 1
                break
        return fields, description, headerrow

    if path.endswith(".gz"):
        with gzip.open(path, "r") as f:
            fields, description, headerrow = read_comments(f)
    else:
        with open(path, "r") as f:
            fields, description, headerrow = read_comments(f)
    names = ["chr", "pos", "id", "ref", "alt", "qual"]
    names += ["filter"] if hasfilter else ["strand"]
    names += ["data", "format"] + samples
    a = pd.read_csv(
        path, sep="\t", header=None, skiprows=headerrow, names=names, index_col=False
    )
    print(description)
    try:
        for j, val in enumerate(a.data.str.split(";").values.tolist()):
            res = dict(
                [
                    (v, True)
                    if v in uniqueargs
                    else (v, np.nan)
                    if "=" not in v
                    else tuple(v.split("="))
                    for v in val
                ]
            )
            for k in fields.keys():
                fields[k].append(res.get(k, None))
    except ValueError:
        print(val)
        raise ValueError("unknown field")
    a = pd.concat(
        [a.drop(columns="data"), pd.DataFrame(data=fields, index=a.index)], axis=1
    )
    for sample in samples:
        uniqformats = a.format.unique().tolist()
        formatcols = set()
        for i in uniqformats:
            formatcols.update(i.split(":"))
        sorting = list(formatcols)
        print("sorting column names: ", sorting)

        def make_format_list(row):
            f = row.format.split(":")
            v = row[sample].split(":")
            assert len(f) == len(v)
            l = []
            for s in sorting:
                if s in f:
                    l.append(v[f.index(s)])
                else:
                    l.append(pd.NA)
            return l

        res = a.apply(make_format_list, axis=1).tolist()

        if len(samples) > 1:
            sorting = [sample + "_" + v for v in sorting]
        a = pd.concat(
            [
                a.drop(columns=sample),
                pd.DataFrame(data=res, columns=sorting, index=a.index),
            ],
            axis=1,
        )
    return a.drop(columns="format"), description
