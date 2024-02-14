import pandas as pd
from mgenepy.epigenetics import chipseq as chip
import argparse

def mapBed(file, name, guides_fn, lib_name):
    """map mutations in one vcf file to regions in the guide bed file"""

    bed = pd.read_csv(
        file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "foldchange"],
    )
    guide_df = pd.read_csv(
            guides_fn,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "foldchange"],
    )
    bed["foldchange"] = 1
    if len(bed) == 0:
        return (name, None)
    val = chip.putInBed(guide_df, bed, mergetype="sum")
    guide_df = guide_df.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True).rename(columns={"foldchange": "sgRNA"})
    guide_df[name] = val
    guide_df[name] = guide_df[name].astype(bool).astype(int)
    guide_df["end"] = guide_df["end"].astype(int)
    guide_df.to_csv(name + "_" + lib_name + "_mut_binary.csv", index=False)

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type=str)
    parser.add_argument("--bed_filenames", default=[], type=lambda x: x.split(","))
    parser.add_argument("--libraries", default=[], type=lambda x: x.split(","))
    parser.add_argument("--guides", default=[], type=lambda x: x.split(","))
    args = parser.parse_args()

    sample_id = args.sample_id
    bed_filenames = args.bed_filenames
    libraries = args.libraries
    guides = args.guides

    assert(len(bed_filenames) == len(libraries))
    assert(len(libraries) == len(guides))

    for i in range(len(libraries)):
        mapBed(bed_filenames[i], sample_id, guides[i], libraries[i])


if __name__ == "__main__":
    main()