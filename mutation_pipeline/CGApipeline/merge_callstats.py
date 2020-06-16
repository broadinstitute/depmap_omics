# Merging call stats for Mutect1 and Mutect2 
import csv
import sys
from argparse import ArgumentParser

def merge_files(FILES,OUT_FILE):
	# common header starters
	header_list = ['##', '#CHROM', 'contig']
	if len(FILES) == 0:
		sys.exit('Empty list of files provided for merging')	
	else:
		writer = open(OUT_FILE, 'w')
		# write header
		with open(FILES[0], 'rb') as in_file:			
			for line in in_file:
				if any(substring in line for substring in header_list):
					writer.write(line)
			in_file.close()
	
		for IN_FILE in FILES:
			with open(IN_FILE, 'rb') as in_file:				
				for line in in_file:
					if not any(substring in line for substring in header_list):
						writer.write(line)
				in_file.close()        
		writer.close()       

if __name__ == '__main__':
	parser = ArgumentParser(description='Merge Mutect1 or Mutect2 Call Stats to one file')
	parser.add_argument('FILE_ARRAY',help='input array of CallStats Files (Mutect1 or Mutect2 call stats)')
	parser.add_argument('OUT_FILE',help='output File (callstat)')			
	args = parser.parse_args()
	if(args):
		print 'Picked up   FILE_ARRAY  : '+str(args.FILE_ARRAY)
		print 'Picked up   OUT_FILE  : '+str(args.OUT_FILE)
		FILES = list()
		with open(args.FILE_ARRAY, 'rb') as file_list:
			for line in file_list:
				FILES.append(line.rstrip())
		merge_files(FILES,args.OUT_FILE)
	else:
		parser.print_help()