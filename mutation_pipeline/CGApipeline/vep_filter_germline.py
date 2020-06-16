import argparse

def filter_vcf(IN_FILE, OUT_FILE):
    with open(IN_FILE, 'r') as reader, open(OUT_FILE, 'w') as writer:
        for line in reader:
            if line.startswith('##'):
                writer.write(line)
            elif line.startswith('#CHROM'):
                writer.write(line)
                header = line.strip('\n').split('\t')
                info_idx = header.index('INFO')
            else:
                values = line.strip('\n').split('\t')
                info = values[info_idx].split(';')
                found_germline_variant = False
                for i in info:
                    if i.startswith('gnomADg_AF='):
                        af = [float(x) for x in i.split('gnomADg_AF=')[1].split(',') if len(x)>1]
                        if len(af) > 0:
                            found_germline_variant = True
                            if max(af) < 0.001:
                                writer.write(line) 
                if found_germline_variant == False:
                    writer.write(line)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter out mutations that have AF in GNOMAD higher than 0.1%')
	parser.add_argument('IN_FILE',help='input VCF File (output from VEP) ')
	parser.add_argument('OUT_FILE',help='output VCF File (filtered Gnomad Germline variants)')
	args = parser.parse_args()
	if(args):
		print 'Picked up    IN_FILE  : ' + str(args.IN_FILE)
		print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
		filter_vcf(args.IN_FILE,args.OUT_FILE)
	else:
		parser.print_help()                