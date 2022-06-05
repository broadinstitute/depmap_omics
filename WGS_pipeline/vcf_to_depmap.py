from depmapomics import vcf
from genepy import mutations
import sys
import pyarrow.parquet as pq
import pyarrow as pa


vcf_filename = sys.argv[1]
sample_name = sys.argv[2]
n_rows = int(sys.argv[3]) if len(sys.argv) > 3 else 500_000
use_multi = "True" == sys.argv[4] if len(sys.argv) > 4 else False
onco_kb = "True" == sys.argv[5] if len(sys.argv) > 5 else False
force_keep = sys.argv[6].split(",") if len(sys.argv) > 6 else []

print(
    "inputs: vcf_filename:",
    vcf_filename,
    ", sample_name:",
    sample_name,
    ", n_rows:",
    n_rows,
    ", use_multi:",
    use_multi,
    ", onco_kb:",
    onco_kb,
    ", force_keep:",
    force_keep,
)

tobreak = False
for i in range(10_000):
    # read in vcf as a df
    vcf_file, _, _ = mutations.vcf_to_df(
        vcf_filename,
        additional_cols=["PON"],
        parse_filter=True,
        force_keep=force_keep,
        drop_null=True,
        nrows=n_rows,
        skiprows=n_rows * i,
    )
    filen = len(vcf_file)
    if filen < n_rows:
        # we have reached the end:
        tobreak = True

    # improve
    vcf_file = vcf.improve(
        vcf_file,
        cols_to_drop=["clinvar_vcf_mc", "oreganno_build", "pid"],
        force_list=["oc_genehancer__feature_name"],
        with_onco_kb=onco_kb,
        split_multiallelic=use_multi,
        torename=vcf.TO_RENAME_OC if onco_kb else vcf.TO_RENAME,
    )

    # checking we have the same set of columns
    dropped = vcf_file.columns.tolist()
    if i == 0:
        prev_dropped = dropped
    if dropped != prev_dropped:
        raise ValueError(
            "we are removing different sets of columns",
            dropped,
            list(set(dropped) ^ set(prev_dropped)),
        )

    # save full
    # need pyarrows
    print("to parquet")
    pq.write_to_dataset(
        pa.Table.from_pandas(vcf_file), root_path=sample_name + "-maf-full.parquet"
    )

    # save maf
    print("saving maf")
    if i == 0:
        vcf.to_maf(vcf_file, sample_name, drop_multi=True)
    else:
        vcf.to_maf(vcf_file, sample_name, drop_multi=True, mode="a", header=False)
    del vcf_file

    if tobreak:
        break
print("finished, processed {} rows".format((n_rows * i) + filen))
