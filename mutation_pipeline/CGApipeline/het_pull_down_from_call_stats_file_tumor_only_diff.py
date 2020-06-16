import argparse
import csv
import math
import sys

MAX_ROWS = 10000
VALID_CHROMOSOMES = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                     "19", "20", "21", "22", "X", "Y"}


def het_pulldown(call_stats_file, db_variant_keys, het_writer, no_normal=False):
    header = None
    header_indexes = dict()
    for call_stats_line in call_stats_file:
        if call_stats_line[0] == "#":  # skip comments (assumes that the first character of the str is
            continue
        call_stats_line_tokens = call_stats_line.rstrip("\n").split("\t")
        if header is None:
            header = call_stats_line
            header_indexes = dict([(header_token, index)
                                   for index, header_token in enumerate(call_stats_line_tokens)])
            continue
        chrom = call_stats_line_tokens[header_indexes["contig"]]
        if chrom not in VALID_CHROMOSOMES:  # ignore chromosomes that are not needed
            continue
        pos = int(call_stats_line_tokens[header_indexes["position"]])
        ref_allele = call_stats_line_tokens[header_indexes["ref_allele"]]
        t_ref_allele_count = float(call_stats_line_tokens[header_indexes["t_ref_count"]])
        alt_allele = call_stats_line_tokens[header_indexes["alt_allele"]]
        t_alt_allele_count = float(call_stats_line_tokens[header_indexes["t_alt_count"]])
        if not no_normal:  # normal is present (you can determine germline hets)
            n_ref_allele_count = float(call_stats_line_tokens[header_indexes["n_ref_count"]])
            n_alt_allele_count = float(call_stats_line_tokens[header_indexes["n_alt_count"]])
            n_allele_count = n_ref_allele_count+n_alt_allele_count
            if n_allele_count > 0:
                n_alt_allelic_fraction = n_alt_allele_count/n_allele_count
                if math.fabs(n_alt_allelic_fraction-0.5) > 0.12:  # allelic fraction [0.38, 0.62]
                    continue
            else:  # no coverage in normal
                continue
        else:  # normal is missing
            # Remove all tumor homs as they're unlikely to be LOH events (unlikely that you'd get the same change in
            # base pair)
            t_allele_count = t_ref_allele_count+t_alt_allele_count
            if t_allele_count > 0 and t_ref_allele_count > 1 and t_alt_allele_count > 1:
                t_alt_allelic_fraction = t_alt_allele_count/t_allele_count
                t_ref_allelic_fraction = t_ref_allele_count/t_allele_count
                if t_allele_count < 10:
                    continue
            else:  # no coverage in tumor (defensive programming)
                continue


        t_ref_allele_count = int(t_ref_allele_count)
        t_alt_allele_count = int(t_alt_allele_count)
        variant_key = chrom + "_" + str(pos) + "_" + ref_allele + "_" + alt_allele
        if variant_key in db_variant_keys:
            row_2_write = [chrom, pos, ref_allele, alt_allele, t_ref_allele_count, t_alt_allele_count]
            het_writer.writerow(row_2_write)


def main():
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("--input-call-stats-filename", dest="input_call_stats_filename", action="store", required=True,
                        help="Input call-stats filename.")
    parser.add_argument("--db-variant-key-filename", dest="db_variant_key_filename", action="store", required=False,
                        help="Name of the db file that has a list of keys that are formed by the "
                             "concatenation of chrom, pos, ref and alt.",
                        default="/xchip/cga_home/mgupta/resources/dbsnp/dbsnp.variant_key.tsv")
    parser.add_argument("--no-normal", dest="no_normal", action="store", required=False,
                        help="No normal sample was provided.", default="False")
    parser.add_argument("--output-filename", dest="output_filename", action="store", required=False,
                        help="Name of the output file that is used as an input to Alleic CapSeg.", default="Tumor.cov")
    args, _ = parser.parse_known_args()

    with open(args.input_call_stats_filename, "r") as call_stats_file, \
            open(args.output_filename, "w") as output_tsv_file, \
            open(args.db_variant_key_filename, "r") as db_variant_key_file:

        het_writer = csv.writer(output_tsv_file, delimiter="\t")
        db_variant_keys = load_db_variant_keys(db_variant_key_file=db_variant_key_file)
        fieldnames = ["Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1", "i_t_ref_count",
                      "i_t_alt_count"]
        het_writer.writerow(fieldnames)
        het_pulldown(call_stats_file, db_variant_keys, het_writer, str2bool(v=args.no_normal))


def load_db_variant_keys(db_variant_key_file):
    sys.stdout.write("Begin reading variant keys.\n")
    db_variant_keys = set()
    for db_variant_key_line in db_variant_key_file:
        db_variant_keys.add(db_variant_key_line.rstrip("\n"))
    sys.stdout.write("Finished reading variant keys.\n")
    return db_variant_keys


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1", "y")

if __name__ == "__main__":
    status = main()
    sys.exit(status)
