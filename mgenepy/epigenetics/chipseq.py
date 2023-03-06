import numpy as np


# @jit(float32[:](float32[:, 3], float32[:, 4], int8, str))
def putInBed(conscensus, value, window=10, mergetype="mean"):
    """
    given a conscensus bed-like dataframe and another one, will merge the second one into the first

    Args:
    -----
      conscensus df[start,end,chrom] the conscensus (first one)
      value df[start, end, chrom,foldchange] the value one (second one)
      mergetype: flag: mean,first,last, on how to merge ttwo peaks that would fall on the same one
      window: int on max distance of second df peaks' to the first ones, to still merge them
        re.g. 0 is hard overap, 10000000 is nearest peak)

    Returns:
    -------
      np.array of values of the second dataframe to add to the first one
      (e.g. can do df1[newcol] = returned_array)
    """
    conscensus = conscensus.sort_values(by=["chrom", "start", "end"]).reset_index(
        drop=True
    )
    value = value.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
    locinvalue = 0
    loc = 0
    tot = 0
    num = []
    res = np.zeros(len(conscensus))
    not_end = True

    def add(res, num, not_end, loc):
        if len(num) > 0:
            if mergetype == "mean":
                res[loc] = np.mean(num)
            elif mergetype == "sum":
                res[loc] = np.sum(num)
            elif mergetype == "first":
                res[loc] = num[0]
            elif mergetype == "last":
                res[loc] = num[-1]
            else:
                raise ValueError("must be one of")
            num = []
        loc += 1
        if loc == len(conscensus):
            not_end = False
        return res, num, not_end, loc

    while not_end:
        # print(loc/len(conscensus),end="\r")
        a = conscensus.iloc[loc]
        b = value.iloc[locinvalue]
        if b.chrom < a.chrom:
            locinvalue += 1
            if locinvalue == len(value):
                not_end = False
        elif b.chrom > a.chrom:
            loc += 1
            if loc == len(conscensus):
                not_end = False
        elif b.start < a.start:
            if b.end + window > a.start:
                tot += 1
                num.append(b.foldchange)
                if b.end > a.end + window:
                    res, num, not_end, loc = add(res, num, not_end, loc)
                    continue
            locinvalue += 1
            if locinvalue == len(value):
                not_end = False
        elif b.start < a.end + window:
            tot += 1
            num.append(b.foldchange)
            if b.end > a.end + window:
                res, num, not_end, loc = add(res, num, not_end, loc)
                continue
            locinvalue += 1
            if locinvalue == len(value):
                not_end = False
        else:
            res, num, not_end, loc = add(res, num, not_end, loc)
    # print(str(tot) + " were merged into conscensus")
    return res
