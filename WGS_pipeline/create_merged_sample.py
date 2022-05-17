import sys
import bgzip
import gzip
import numpy as np
from torch import zero_


print(sys.argv)

check_index = 0
multi_vcf = sys.argv[1]
sample_id = sys.argv[2]
rna_sample = sys.argv[3] if len(sys.argv) == 4 else "none"
min_af_tocall = float(sys.argv[4]) if len(sys.argv) == 5 else 0.05
min_read_good_qual = int(sys.argv[5]) if len(sys.argv) == 6 else 15
pseudo_count = float(sys.argv[6]) if len(sys.argv) == 7 else 1
likely_seq_error = 0
issues = 0


def solve_phasing(mut_s, geno, min_read_good_qual=15):
    """solve_phasing takes a list of mutations and a genotype and returns a list of mutations with the correct phasing

    Args:
        mut_s (_type_): _description_
        geno (_type_): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    conflicts = False
    pattern = set([])
    for gt in mut_s[0]:
        if "|" in gt:
            pattern = pattern | set(
                [tuple(j.split("|")) for j in gt.split("/") if "|" in j]
            )
    for val in pattern:
        if (val[1], val[0]) in pattern and val[1] != val[0]:
            print("Error in pattern:\n" + str(site[:5] + site[9:]))
            conflicts = True
            break
            # pattern.remove(val)
    rep = []
    if not conflicts:
        for val in pattern:
            if (val[0] in geno and val[1] in geno) and not (
                val[0] == val[1] and geno.count(val[0]) == 1
            ):
                if len(geno) == 3:
                    geno = "|".join(val)
                    continue
                if val[0] in rep or val[1] in rep:
                    print("Error in pattern:\n" + str(site[:5] + site[9:]))
                    conflicts = True
                    break
                if geno.index(val[0]) == 0:
                    geno = geno.replace(val[0] + "/", "")
                else:
                    geno = geno.replace("/" + val[0], "")
                if geno.index(val[1]) == 0:
                    geno = geno.replace(val[1] + "/", "")
                else:
                    geno = geno.replace("/" + val[1], "")
                rep.append(val[0])
                rep.append(val[1])
            geno = "|".join(val) + geno

    # managing conflicts:
    freqs = np.array([sam.split(",") for sam in mut_s[2]])
    freqs = freqs[mut_s[3].astype(int) > min_read_good_qual]
    if len(freqs) >= 2:
        for freq in freqs.T:  # .astype(float)
            for i, f in enumerate(freq[1:]):
                if freq[i] != "." and f != ".":
                    if abs(float(f) - float(freq[i])) > 0.4:
                        conflicts = True

    return geno, conflicts


def merging(mut_s, num, min_af_tocall=0.05, pseudo_count=0.95):
    """merging _summary_

    Args:
        mut_s (np.array): mutation information from vcf as a numpy array
        num (int): number of alleles 
        min_af_tocall (float, optional): minimum allele frequency to call a mutation. Defaults to 0.05.

    Returns:
        str: merged mutation information
        list: list of sites to drop
    """
    # ad: the sum of the depths across all samples for each alleles
    ad = np.zeros((num))
    ads = np.zeros((len((mut_s[1])), num))
    for i, depth in enumerate(mut_s[1]):
        if depth == ".":
            continue
        ads[i] = np.array([int(j) if j != "." else 0 for j in depth.split(",")])
        ad += ads[i]
    ad = ad.astype(int)

    # QCs if == 0 it is the ref allele, so we keep it
    loc = ad != 1
    if ~loc[0]:  # condition 103
        ad[0] = 0
        loc[0] = True

    to_drop = np.argwhere(~loc)[:, 0].tolist()
    # af: the frequency of each allele with a pseudocount of 0.95
    af = np.zeros((num - 1))
    for i in range(num - 1):
        af[i] = (ad[i + 1] + pseudo_count) / (sum(ad) + (pseudo_count * len(ad)))

    af = af[loc[1:]]

    ad = ad[loc]
    # depth
    depth = ad.sum()

    # F1R2: same as ad
    f1r2 = np.zeros((num), dtype=int)
    for sam in mut_s[4]:
        if sam == ".":
            continue
        sam = np.array([int(j) if j != "." else 0 for j in sam.split(",")])
        f1r2 += sam
    f1r2 = f1r2[loc]
    # F2R1: same as ad
    f2r1 = np.zeros((num), dtype=int)
    for sam in mut_s[5]:
        if sam == ".":
            continue
        sam = np.array([int(j) if j != "." else 0 for j in sam.split(",")])
        f2r1 += sam
    f2r1 = f2r1[loc]
    # FAD: same as ad
    fad = np.zeros((num), dtype=int)
    for sam in mut_s[6]:
        if sam == ".":
            continue
        sam = np.array([int(j) if j != "." else 0 for j in sam.split(",")])
        fad += sam
    fad = fad[loc]

    # creating string
    to_add = ["", "", "", "", "", "", "", ".,.,.,."]
    # gt
    if sum(af) < 1 - min_af_tocall:
        if not (sum(af) > 1 - (2 * min_af_tocall) and ad[0] == 0):
            to_add[0] += "0/"
        # else: condition 103
    for i, allele in enumerate(af):
        if allele >= 1 - min_af_tocall:
            # this allele on both strand
            to_add[0] += str(i + 1) + "/" + str(i + 1)
        elif allele >= min_af_tocall:
            to_add[0] += str(i + 1) + "/"
    if to_add[0] == "0/":
        to_add[0] += "0"
    if to_add[0][-1] == "/":
        if len(to_add[0]) > 2:
            to_add[0] = to_add[0][:-1]
        # in the case where we ran into condition 103
        else:
            to_add[0] += to_add[0][0]
    # rest
    to_add[1] = ",".join([str(i) for i in ad])
    to_add[2] = ",".join(["{:.3f}".format(i) for i in af])
    to_add[3] = str(depth)
    to_add[4] = ",".join([str(i) for i in f1r2])
    to_add[5] = ",".join([str(i) for i in f2r1])
    to_add[6] = ",".join([str(i) for i in fad])
    return to_add, to_drop


""" reads in a vcf file and adds the fake sample mutation information

Raises:
    ValueError: 
"""
with gzip.open(multi_vcf, "r+") as f:
    with open(sample_id + "_multi.vcf.gz", "wb") as raw_out:
        with bgzip.BGZipWriter(raw_out) as f_out:
            for i, line in enumerate(f):
                # print progress
                if i % 100_000 == 0:
                    print("Processed " + str(i) + " variants")
                original_string = line.decode("utf-8")
                # process header
                if original_string[:2] == "##":
                    # do nothing
                    to_print = original_string
                elif original_string[:2] == "#C":
                    # adding a descriptor for new info
                    to_print = '##INFO=<ID=CONFLICT, Number=0,Type=Flag, \
Description="0 or 1 for whether replicates are in conflicts">\n'
                    to_print += '##INFO=<ID=A_SPE_EXP, \
Number=0,Type=Flag,Description="Variant N is only/mostly expressed in the given \
RNAseq mutation file else "." or "nREF" for non ref">\n'
                    cols = original_string[1:-1].split("\t")
                    to_print += original_string + "\t" + sample_id + "\n"
                else:
                    site = original_string[:-1].split("\t")

                    mut_s = site[9:]

                    mut_s = np.array([mut.split(":") for mut in mut_s])

                    non_rna = mut_s[
                        [j for j, loc in enumerate(cols[9:]) if rna_sample not in loc],
                    ].T
                    rnaloc = [j for j, loc in enumerate(cols[9:]) if rna_sample in loc]
                    rna = mut_s[rnaloc].T
                    mut_s = mut_s.T

                    # speed up: don't do anything if only one sample has reads
                    loc = [i for i, af in enumerate(non_rna[1]) if af != "."]
                    if len(loc) == 1:
                        to_print = "\t".join(site + [site[9 + loc[0]]]) + "\n"
                        f_out.write(to_print.encode())
                        continue
                    elif len(loc) == 0:
                        print("RNA is creating a new mutation: ..skipping")
                        continue
                    non_rna_wval = non_rna[:, loc]
                    num = max(
                        [
                            len(non_rna_wval[1][i].split(","))
                            for i in range(len(mut_s[1]))
                        ]
                    )
                    # if num==2: single-allelic
                    # else: multi-allelic

                    # merging / summing
                    to_add, to_drop = merging(
                        non_rna_wval, num=num, min_af_tocall=min_af_tocall
                    )

                    if len(to_drop) > 0:
                        likely_seq_error += len(to_drop)

                    # re-phasing, using all available data:
                    # quick:
                    if len(set(non_rna_wval[0])) == 1:
                        geno = non_rna_wval[0][0]
                    else:
                        geno, conflicts = solve_phasing(non_rna_wval, to_add[0])
                        if conflicts:
                            site[6] += ";conflict"
                            issues += 1
                    to_add[0] = geno
                    # adding new info
                    site.append(":".join(to_add))

                    # allele_expression
                    if len(rnaloc) > 0:
                        if rna[2] != ".":
                            rna_geno = rna[0]
                            rna_geno = list(set(rna_geno.replace("|", "/").split("/")))
                            if len(rna_geno) == 1:
                                allele_specific = rna_geno[0]
                            else:
                                if "0" not in rna_geno:
                                    # case: 1/2, 1/3,
                                    allele_specific = "noREF"
                                else:
                                    # interesting case  0/2, 0/1/3
                                    allele_specific = "."
                            site[7] += ";A_SPE_EXP=" + allele_specific

                    to_print = "\t".join(site) + "\n"

                f_out.write(to_print.encode())
print("Processed " + str(i) + " variants")
print("likely_seq_error: " + str(likely_seq_error))
print("likely_seq_error/total: {:.3%}".format((likely_seq_error / i) * 100) + "%")
print("issues: " + str(issues))
print("issues/total: {:.3%}".format((issues / i) * 100) + "%")
with open(sample_id + "_seq_error.txt", "w") as f:
    f.write("{:.3%}".format((likely_seq_error / i) * 100) + "%")
with open(sample_id + "_issues.txt", "w") as f:
    f.write("{:.3%}".format((issues / i) * 100) + "%")

