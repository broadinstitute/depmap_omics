# Merging call stats for Mutect1 and Mutect2 
#import csv
#import sys
from argparse import ArgumentParser

def filter_mutations(FILE_IN,FILE_OUT_REJECTED, FILE_OUT_PASSED, passed_flag):	
	rejected_writer = open(FILE_OUT_REJECTED, 'w')
	passed_writer = open(FILE_OUT_PASSED, 'w')
		
	with open(FILE_IN, 'rb') as in_file:			
		for line in in_file:
			if not(line.startswith('#')) and not(line.startswith('CONTIG')) and not(line.startswith('contig')):
				if passed_flag in line:
					passed_writer.write(line)
				else:
					rejected_writer.write(line)
			else:
				passed_writer.write(line)
				rejected_writer.write(line)
				
		in_file.close()			       
		rejected_writer.close()       
		passed_writer.close()       

if __name__ == '__main__':
	parser = ArgumentParser(description='Merge Mutect1 or Mutect2 Call Stats to one file')
	parser.add_argument('FILE_IN',help='input ')
	parser.add_argument('FILE_OUT_PASSED',help='output file with variants that passed filter')
	parser.add_argument('FILE_OUT_REJECTED',help="output file with variants that didn't pass filter")
	parser.add_argument('passed_flag',help="string to search for in file (KEEP, PASS)")
	
	args = parser.parse_args()
	if(args):
		print 'Picked up   FILE_IN  : '+str(args.FILE_IN)
		print 'Picked up   FILE_OUT_PASSED  : '+str(args.FILE_OUT_PASSED)
		print 'Picked up   FILE_OUT_REJECTED  : '+str(args.FILE_OUT_REJECTED)
		print 'Picked up   passed_flag  : '+str(args.passed_flag)

		filter_mutations(args.FILE_IN, args.FILE_OUT_REJECTED, args.FILE_OUT_PASSED, args.passed_flag)
	else:
		parser.print_help()