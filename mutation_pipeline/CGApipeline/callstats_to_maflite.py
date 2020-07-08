import csv
from argparse import ArgumentParser

def convert_to_maflite(IN_FILE,OUT_FILE):
	with open(OUT_FILE, 'w') as writer:
		with open(IN_FILE, 'rb') as reader:
			for line in reader:
				if line.startswith('#'):
					writer.write(line)
				elif line.startswith('contig'):
					header = line.split('\t')
					idxPos = header.index('position')
					writer.write('\t'.join(['chr'] + header[1:idxPos] + ['start'] + ['end'] + header[idxPos+1:]) + '\n')
				else:
					values = line.strip().split('\t')
					writer.write('\t'.join(values[:idxPos+1] + [values[idxPos]] + values[idxPos+1:]) + '\n')


if __name__ == '__main__':
	parser = ArgumentParser(description='Convert MuTect1 Call Stats to MafLite for Oncotator')
	parser.add_argument('IN_FILE',help='input File (Mutect1 call stats)')
	parser.add_argument('OUT_FILE',help='output File (mafLite)')
	args = parser.parse_args()
	if(args):
		print 'Picked up    IN_FILE  : ' + str(args.IN_FILE)
		print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
		convert_to_maflite(args.IN_FILE,args.OUT_FILE)
	else:
		parser.print_help()