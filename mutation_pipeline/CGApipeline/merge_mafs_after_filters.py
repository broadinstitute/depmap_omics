import logging
import argparse
import collections
from time import gmtime, strftime


logging.basicConfig(filename='merge_maf_logging.log',level=logging.DEBUG)

def makehash():
    return collections.defaultdict(makehash)


class StoreInDict(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		d = getattr(namespace, self.dest)
		for opt in values:
			k,v = opt.split("=", 1)
			k = k.lstrip("-")
			if k in d:
				print  'Error duplicate key'
			else:
				d[k] = v
		setattr(namespace, self.dest, d)


def parse_mutation(line, header, filter_name):
	values = line.strip('\n').split('\t')
	if len(header) != len(values):
		error ="Error, line "+str(line_count)+" header/column length mismatch!"
		raise Exception(error)
	else:
		mutation = dict(zip(header, values))
		# find judgement column and add filter name to list of passed or failed filters
	return mutation


def read_mafs_into_dict(INPUT_DATA):
	"""		
		Reads in MAF files and merged all mutations into one dictionary
		where key = chrom, start_pos, end_pos, ref_allele, alt_allele
		value => all other fields in the MAFs
	"""
	comments = list()
	all_mutations = makehash()
	for filter_name, in_maf in INPUT_DATA.iteritems():
		with open(in_maf, 'r') as reader:
			print in_maf
			for line in reader:
				if line.startswith('Hugo_Symbol'):
					header = line.strip('\n').split('\t')
				elif line.startswith('#'):
					comments.append(line)
				else:
					mut = parse_mutation(line, header, filter_name)
					key = '-'.join([mut['Chromosome'], mut['Start_position'], mut['End_position'], 
									mut['Reference_Allele'], mut['Tumor_Seq_Allele2']])
					if key in all_mutations:
						mut_entry = all_mutations[key]
						for field, value in mut.iteritems():
							if not field in mut_entry:
								all_mutations[key][field] = value
					else:
						all_mutations[key] = mut
						all_mutations[key]['Passed_Filters'] = ''
						all_mutations[key]['Failed_Filters'] = ''

					if mut['FILTER_JUDGEMENT'] == 'PASS':
						all_mutations[key]['Passed_Filters'] += filter_name + ','
					else:
						all_mutations[key]['Failed_Filters'] += filter_name + ','	
	return all_mutations


def get_header(maf):
	""" 
		Returns list of column names from MAF file
		Searches for line that starts with Hugo_Symbol in MAF file
	"""
	with open(maf, 'r') as reader:
		for line in reader:
			if line.startswith('Hugo_Symbol'):
				return line.strip('\n').split('\t')
	error = "Error, did not find header in maf %s (e.g. line starts with Hugo_Symbol), check the correct formatting of maf!", str(maf)
	raise Exception(error)


def get_union_header(INPUT_DATA):
	"""
		Iterate over all headers and merge all fields in one header
		Read in MAF files and get list that is the union of the column names
	"""
	union_header = list()	
	for filter_name, in_maf in INPUT_DATA.iteritems():
		header = get_header(in_maf)
		for h in header:
			if not h in union_header:
				union_header.append(h)
	logging.info(" Unified header contains %d number of columns", len(union_header))
	return union_header


def sort_mutation_keys(keys):
	"""
		Sorts mutation keys first by chromosome and then by starting position
	"""
	sorted_keys = list()
	split_keys = [k.split('-') for k in keys]
	chromosomes = set([key[0] for key in split_keys])
	chrom_strs = list(filter(lambda x : not x.isdigit(),chromosomes))
	chrom_ints = [int(y) for y in list(filter(lambda x: x.isdigit(), chromosomes))]
	sorted_chromosomes = sorted(chrom_ints) + sorted(chrom_strs)
	for chrom in sorted_chromosomes:
		chrom_list = [key for key in split_keys if key[0] == str(chrom)]
		chrom_list.sort(key=lambda x: int(x[1]))
		for c in chrom_list:	
			sorted_keys.append('-'.join(c))
	return sorted_keys


def write_to_maf(all_mutations, header, pair_name):
	"""
		Writing mutations in order to an aggregated MAF
	"""
	union_maf = pair_name + '.merged_union.maf'
	intersection_maf = pair_name + '.passed_all_filters.maf'
	with open(union_maf, 'w') as union_writer, open(intersection_maf, 'w') as inter_writer:
		union_header = header + ['Passed_Filters', 'Failed_Filters']
		union_writer.write('\t'.join(union_header) + '\n')
		inter_header = header + ['Passed_Filters']
		inter_writer.write('\t'.join(inter_header) + '\n')
		# sort mutations by their keys
		keylist = sort_mutation_keys(all_mutations.keys())
		for key in keylist:
			mutation = all_mutations[key]			
			if mutation['Failed_Filters'] == '':
				inter_writer.write('\t'.join([mutation[k] if k in mutation else '' for k in inter_header]) + '\n')			
			union_writer.write('\t'.join([mutation[k] if k in mutation else '' for k in union_header]) + '\n')
			
	
def main():	
	logging.info(" Merging maf files %s", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	parser = argparse.ArgumentParser(description="Read MAF files, merge them based on the method provided.  NOTE that the merge/join is based on key='chrom,start,stop,ref,tum_allele_1,tum_allele_2'.  NOTE that 1) the set of columns if the set-union of the two and 2) if two columns of the same name disagree value in first input gets output", prefix_chars=' ')	
	parser.add_argument("IN_PARAMS", nargs="*", action=StoreInDict, default=dict(), help="input filter names and maf files")			
	args = parser.parse_args()
	INPUT_DATA = dict()	
	pair_name = 'None'

	if(args):
		logging.info(args.IN_PARAMS)
		for key in args.IN_PARAMS:
			if key.startswith('pair_name'):
				pair_name = args.IN_PARAMS[key]
			else:				
				INPUT_DATA[key] = args.IN_PARAMS[key]

		logging.info(" Reading data from %s" %str(INPUT_DATA))		
		
		union_header = get_union_header(INPUT_DATA)
		all_mutations = read_mafs_into_dict(INPUT_DATA)
				
		write_to_maf(all_mutations, union_header, pair_name)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()