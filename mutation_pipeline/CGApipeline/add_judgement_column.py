import sys
import argparse

def add_judgement_column(INPUT_MAF, OUTPUT_MAF, decision_column, pass_flag='1'):
	with open(INPUT_MAF, 'r') as reader, open(OUTPUT_MAF, 'w') as writer:
		for line in reader:
			if line.startswith('Hugo_Symbol'):
				header = line.strip('\n').split('\t')
				decision_col_idx = header.index(decision_column)
				writer.write('\t'.join(header + ['FILTER_JUDGEMENT']) + '\n')
			elif line.startswith('#'):
				writer.write(line)
			else:
				values = line.strip('\n').split('\t')
				decision = values[decision_col_idx]
				if decision == pass_flag or decision == '':
					writer.write('\t'.join(values + ['PASS']) + '\n')
				else:
					writer.write('\t'.join(values + ['REJECT']) + '\n')				


def parseOptions():
    description = '''
    Adding judgement column to MAF file after running through filter.
    Judgement column has two values PASS/REJECT 
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input',  metavar='input',  type=str, help='input MAF after filtering task')
    parser.add_argument('-o', '--output', metavar='output', type=str, help='output MAF with added judgement column')
    parser.add_argument('-c', '--column', metavar='column', type=str, help='column name that contains judgements for this filter')
    parser.add_argument('-p', '--pass_flag',   metavar='pass',   type=str, help='value for passing filter')
    args = parser.parse_args()
    return args


def main():	
	args = parseOptions()
	if args:
		add_judgement_column(args.input, args.output, args.column, args.pass_flag)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()