from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
import os
import re
import sys
from GenericTsvReader import GenericTsvReader


def parseOptions():
    epilog = """This script is generic.
    It simply filters lines based on field names and regular expressions."""
    desc = "filter a tsv on fields."
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)

    parser.add_argument("-i","--input_tsv", type=str, help="Input aggregated tsv")
    parser.add_argument( "-f","--filter_cols", type=str, help="comma separated list of column names to filter (see values)")
    parser.add_argument( "-o","--output_tsv", type=str, help="Output file in MAF format")
    parser.add_argument( "-v","--filter_vals", type=str, help="comma separated list of REs to filter.  Keep only ones that match the RE.  No spaces between commas, please.")   

    args = parser.parse_args()
    return args


def main():
    args = parseOptions()

    # create compiled REs
    re_strings = args.filter_vals.split(",")
    print(str(re_strings))
    re_compiled = [re.compile(re_string) for re_string in re_strings]
    filter_cols = args.filter_cols.split(",")
    print(str(filter_cols))

    input_tsv_stub=os.path.basename(args.input_tsv)

    input_reader = GenericTsvReader(args.input_tsv)
    comments = input_reader.getComments()
    headers = input_reader.getFieldNames()

   
    # Initialize the file with comments and header.  Leave file pointer open to write (as DictWriter).
    fp = file(args.output_tsv, 'w')
    if comments is None or comments.endswith("\n") or comments.strip() == "":
        fp.write(comments + "# split from " + os.path.basename(args.input_tsv) + "\n")
    else:
        fp.write(comments + "\n# split from " + os.path.basename(args.input_tsv) + "\n")

    dict_writer = csv.DictWriter(fp, fieldnames=headers, delimiter="\t", lineterminator="\n")
    dict_writer.writeheader()

    # For each line in input tsv
    for line_dict in input_reader:

        is_skip_line = False

        # Check if the line should be kept.  If not, skip this iteration.
        for i,col in enumerate(filter_cols):
            val = line_dict[col]
            search_result = re_compiled[i].search(val)
            if search_result is None:
                is_skip_line = True
                break

        if is_skip_line:
            continue

        # Check the tumor and normal barcodes and map to a pair_id.

        # Write this line to the file pointer for this pair_id
        dict_writer.writerow(line_dict)

    fp.flush()
    fp.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())



