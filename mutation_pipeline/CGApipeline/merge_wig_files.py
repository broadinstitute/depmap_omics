# Merging Mutect1 coverage and power wiggle files
import csv
import sys
from argparse import ArgumentParser

def merge_files(FILES,OUT_FILE):
	# common header starters
	header_list = ['track']
	if len(FILES) == 0:
		sys.exit('Empty list of files provided for merging')	
	else:
		writer = open(OUT_FILE, 'w')
		# write header
		with open(FILES[0], 'r') as in_file:			
			for line in in_file:
				if any(substring in line for substring in header_list):
					writer.write(line)
			in_file.close()
	
		for IN_FILE in FILES:
			with open(IN_FILE, 'r') as in_file:				
				for line in in_file:
					if not any(substring in line for substring in header_list):
						writer.write(line)
				in_file.close()        
		writer.close()       

if __name__ == '__main__':
	parser = ArgumentParser(description='Merge Mutect1 coverage or power wiggle files')
	parser.add_argument('FILE_ARRAY',help='input array of wiggle files (coverage or power)')
	parser.add_argument('OUT_FILE',help='output File (merged wiggle file)')			
	args = parser.parse_args()
	if(args):
		print 'Picked up   FILE_ARRAY  : ' + str(args.FILE_ARRAY)
		print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
		FILES = list()
		with open(args.FILE_ARRAY, 'r') as file_list:
			for line in file_list:
				FILES.append(line.rstrip())
		merge_files(FILES,args.OUT_FILE)
	else:
		parser.print_help()