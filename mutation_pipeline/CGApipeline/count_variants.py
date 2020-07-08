from argparse import ArgumentParser

def count_variants(IN_FILE,OUT_FILE):	
	line_count = 0
	with open(IN_FILE, 'rb') as reader:		
		for line in reader:
			if (not line.startswith('#')) and (not line.startswith('contig') and (not line.startswith('Hugo_Symbol')) and (not line.startswith('chr'))):
				line_count += 1
	writer = open(OUT_FILE, 'w')
	writer.write(str(line_count)+'')


if __name__ == '__main__':
	parser = ArgumentParser(description='Count number of mutations in MAF file')
	parser.add_argument('IN_FILE',help='input MAF File')
	parser.add_argument('OUT_FILE',help='output file with count or mutations in input MAF file')
	args = parser.parse_args()
	if(args):
		print 'Picked up    IN_FILE  : ' + str(args.IN_FILE)
		print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
		count_variants(args.IN_FILE,args.OUT_FILE)
	else:
		parser.print_help()