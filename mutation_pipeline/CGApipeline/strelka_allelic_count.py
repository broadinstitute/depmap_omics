'''
Script for adding allelic depth and allelic fraction to the Strelka VCF file
Also for changing TUMOR and NORMAL to tumor and normal sample names
'''

from argparse import ArgumentParser

def _get_allelic_depth(fields, values):
	d = dict(zip(fields.split(':'), values.split(':')))
	# tier1RefCounts = First comma-delimited value from FORMAT/TAR
	ref_count = d['TAR'].split(',')[0]
	# tier1AltCounts = First comma-delimited value from FORMAT/TIR
	alt_count = d['TIR'].split(',')[0]
	return ref_count, alt_count 


def _compute_allelic_fraction(ref, alt):
	try:
		# Somatic allele frequency is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
		return float(alt) / (float(ref) + float(alt))
	except:
		print ('REF count (%s) or ALT count (%s) could not be converted into floats or both equal to 0.' %(ref, alt))


def _process_vcf_file(inVcf,outVcf, TUMOR_NAME, NORMAL_NAME):
	''' Function modifies Strelka indels VCF file 
		by adding alelic depth (AD) and alelic fraction (AF) to it
	'''
	tumor_idx = 10
	normal_idx = 9
	format_idx = 8
	add_ad_af_description = True
	with open (outVcf, 'w') as writer, open(inVcf, 'rb') as reader:
		for line in reader:
			if line.startswith('##'):
				if line.startswith('##FORMAT') and add_ad_af_description:
					writer.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">' + '\n')
					writer.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">' + '\n')
					add_ad_af_description = False
				writer.write(line)
			elif line.startswith('#CHROM'):
				pieces = line.strip('\n').split('\t')
				#print pieces
				tmr_idx = pieces.index('TUMOR')
				nrml_idx = pieces.index('NORMAL')
				pieces[tmr_idx] = TUMOR_NAME
				pieces[nrml_idx] = NORMAL_NAME
				writer.write('\t'.join(pieces) + '\n')
			else:
				pieces = line.strip('\n').split('\t')
				FORMAT = pieces[format_idx]
				NORMAL = pieces[normal_idx]
				TUMOR = pieces[tumor_idx]
				NORMAL_AD_REF, NORMAL_AD_ALT = _get_allelic_depth(FORMAT, NORMAL)
				NORMAL_AF = _compute_allelic_fraction(NORMAL_AD_REF, NORMAL_AD_ALT)
				TUMOR_AD_REF, TUMOR_AD_ALT = _get_allelic_depth(FORMAT, TUMOR)
				TUMOR_AF = _compute_allelic_fraction(TUMOR_AD_REF, TUMOR_AD_ALT)
				pieces[format_idx] = FORMAT +':AD:AF'
				pieces[normal_idx] = NORMAL + ':' + NORMAL_AD_REF + ',' + NORMAL_AD_ALT + ':' + str(NORMAL_AF)
				pieces[tumor_idx] = TUMOR + ':' + TUMOR_AD_REF + ',' + TUMOR_AD_ALT + ':' + str(TUMOR_AF)
				writer.write('\t'.join(pieces) + '\n')	


def main():
# The main argument parser
	parser = ArgumentParser(description="Script for modifying Strelka Indels VCF by adding alelic depth and fraction to it")
	parser.add_argument('IN_VCF', help='input Strelka VCF containing only indels')
	parser.add_argument('OUT_VCF',help="modified output Strelka VCF ")
	parser.add_argument('TUMOR_NAME',help="tumor sample name, string that will replace TUMOR in the output Strelka VCF ")
	parser.add_argument('NORMAL_NAME',help="normal sample name, string that will replace NORMAL in the output Strelka VCF ")
	args = parser.parse_args()
	args=parser.parse_args()
	if(args):
		_process_vcf_file(args.IN_VCF, args.OUT_VCF, args.TUMOR_NAME, args.NORMAL_NAME)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()
