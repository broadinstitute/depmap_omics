import argparse


def parse_parameters():
    # The main argument parser
    parser = argparse.ArgumentParser(description="Script to ")
    parser.add_argument('--IN_MAF', '-m', dest='IN_MAF', help='input MAF file')
    parser.add_argument('--OUT_FILE', '-o', dest='OUT_FILE', help="File containing single word True (if SNP found) or False (if not)")
    args = parser.parse_args()
    if not args:
        parser.print_help()
    return args


def read_maf(IN_MAF):
    """

    """
    found_snv = False
    with open(IN_MAF, 'r') as reader:
        for line in reader:
            if line.startswith('Hugo_Symbol'):
                header = line.strip('\n').split('\t')
                var_type_idx = header.index('Variant_Type')
            elif line.startswith('#'):
                continue
            else:
                values = line.strip('\n').split('\t')
                type = values[var_type_idx]
                if type == 'SNP':
                    found_snv = True
                    return found_snv
    return found_snv


def write_file(found_snv, OUT_FILE):
    with open(OUT_FILE, 'w') as writer:
        writer.write(str(found_snv).lower())


def main():
    args = parse_parameters()
    found_snv = read_maf(args.IN_MAF)
    write_file(found_snv, args.OUT_FILE)


if __name__ == "__main__":
    main()
