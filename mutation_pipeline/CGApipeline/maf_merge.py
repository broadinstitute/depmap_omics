import sys
import os
import csv
from collections import defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter, FileType


_REQUIRED_MAF_HEADER   = ['Hugo_Symbol',
						  'Entrez_Gene_Id',
						  'Center',
						  'NCBI_Build',
						  'Chromosome',
						  'Start_position',
						  'End_position',
						  'Strand',
						  'Variant_Classification',
						  'Variant_Type',
						  'Reference_Allele',
						  'Tumor_Seq_Allele1',
						  'Tumor_Seq_Allele2',
						  'dbSNP_RS',
						  'dbSNP_Val_Status',
						  'Tumor_Sample_Barcode',
						  'Matched_Norm_Sample_Barcode',
						  'Match_Norm_Seq_Allele1',
						  'Match_Norm_Seq_Allele2',
						  'Tumor_Validation_Allele1',
						  'Tumor_Validation_Allele2',
						  'Match_Norm_Validation_Allele1',
						  'Match_Norm_Validation_Allele2',
						  'Verification_Status',
						  'Validation_Status',
						  'Mutation_Status',
						  'Sequencing_Phase',
						  'Sequence_Source',
						  'Validation_Method',
						  'Score',
						  'BAM_file',
						  'Sequencer']
_REQUIRED_NUM_MAF_COLS = len(_REQUIRED_MAF_HEADER)

PROGRAM_DESCRIPTION = 'This script generates a concatenated MAF.'


#===============================================================================
# Main
#===============================================================================
def main(args):

	parser = ArgumentParser(description=PROGRAM_DESCRIPTION, 
							formatter_class=RawDescriptionHelpFormatter)
	
	parser.add_argument('maf_In_1', help='File in MAF format.')
	parser.add_argument('maf_In_2', help='File in MAF format.') 
	parser.add_argument('maf_Out', help='File output')   
	parser.add_argument('-o', '--option', help='Intersection or union of the input maf files.', required = False, default = 'union')

	args = parser.parse_args()

	if(args):
		print 'Picked up   maf_In_1  : '+ str(args.maf_In_1)
		print 'Picked up   maf_In_2  : '+ str(args.maf_In_2)		
		print 'Picked up   maf_Out  : '+ str(args.maf_Out)		
		
		header_1 = _get_header(args.maf_In_1)				
		header_2 = _get_header(args.maf_In_2)	
		
		union_header = _get_header_union(header_1, header_2)				
		maf_data_1 = _maf_to_dictionary(args.maf_In_1, header_1, union_header)
		maf_data_2 = _maf_to_dictionary(args.maf_In_2, header_2, union_header)						
		_write_to_file([maf_data_1, maf_data_2], union_header, args.maf_Out)
	else:
		parser.print_help()


def _get_header(maf_file):
	''' Read the beginning of the MAF file to determine the headers. Skip starting
	lines that are empty or are comments (start with '#'). Return the array of
	header field names. 
	''' 	  	
	with open(maf_file, 'rb') as f:    
		for line in f:
			if line.startswith('Hugo'):
				header = line.strip('\n').split('\t')			
		if not all(field in header for field in _REQUIRED_MAF_HEADER):        
			raise Exception('MAF file (%s) does not contain required header fields.' % (maf_file)) 		                        
		f.close()
	return header


def _get_header_union (header_1, header_2):
	''' Return array of header fields that are in either maf except the Required MAF Headers
	'''
	header_1_no_required = [item for item in header_1 if item not in _REQUIRED_MAF_HEADER]
	header_2_no_required = [item for item in header_2 if item not in _REQUIRED_MAF_HEADER]	
	return list(set().union(header_1_no_required, header_2_no_required))


def _maf_to_dictionary (maf_file, header, union_header):
	''' Create dictionary from maf file 
	where key is the site: Chromosome_Start Position_End Position
	and value is the data for this site that is converted to dictionary
	(key - header field name, value - value for this field at the site) 		
	'''	
	data = list()		
	missing_fields = [field for field in union_header if field not in header]			
	with open(maf_file, 'rb') as f:    
		for line in f:			
			if (not line.startswith('#') and not line.startswith('Hugo')):			
				row =  line.strip('\n').split('\t')				
				maf_row = dict(zip(header, row))							
				for field in missing_fields:					
					maf_row[field] = ''																												
				data.append(maf_row)	
		f.close()
	return data


def _write_to_file(data, union_header, filename):
	'''
	 Write dictionary to MAF file
	'''    
	with open(filename, 'w') as f:
		f.write('\t'.join(_REQUIRED_MAF_HEADER + union_header) + '\n')
		for d in data:
			for row in d:
				values = list()
				for field in _REQUIRED_MAF_HEADER:
					values.append(row[field])				
				for field in union_header:
					values.append(row[field])								
				f.write('\t'.join(values) + '\n')
		f.close()


#===============================================================================
# Run Main
#===============================================================================
if __name__ == '__main__':
	main(sys.argv[1:])          