# Computer AF for Mutect2 call stats

import re
import traceback
import sys
from argparse import ArgumentParser

def processVCF(IN_VCF,OUT_VCF,TUMOR_STR,NORMAL_STR, CASE_NAME, CONTROL_NAME):
	reader=open(IN_VCF,'r')
	writer=open(OUT_VCF,'w')
	tumor_piece_idx=None
	for line in reader:
		sline=line.strip()
		if(sline.startswith('##')):			
			if 'normal_sample' in sline:
				pieces = sline.split('=')
				writer.write(''.join(pieces[:-1] + ['='] + [CONTROL_NAME]) + '\n')
			elif 'tumor_sample' in sline:
				pieces = sline.split('=')
				writer.write(''.join(pieces[:-1] + ['='] + [CASE_NAME]) + '\n')
			#if begins with '#' then pass-thru
			else:
				writer.write(sline+'\n')
		elif(sline.startswith('#CHROM')):
			pieces=sline.split('\t')
			tmr_idx=pieces.index(TUMOR_STR)
			nrml_idx=pieces.index(NORMAL_STR)
			pieces[tmr_idx] = CASE_NAME
			pieces[nrml_idx] = CONTROL_NAME
			writer.write('\t'.join(pieces)+'\n')
		elif(sline.startswith('#')):
			writer.write(sline+'\n')
		else:
			pieces=sline.split('\t')
			format_piece=pieces[8]
			format_pieces=format_piece.split(':')
			AD_idx=format_pieces.index('AD')
			AF_idx=format_pieces.index('AF')
			tmr_piece=pieces[tmr_idx]
			nrml_piece=pieces[nrml_idx]
			##FORMAT=<ID=AD,Number=R,Type=Integer,Description='Allelic depths for the ref and alt alleles in the order listed>
			tmr_pieces=tmr_piece.split(':')
			nrml_pieces=nrml_piece.split(':')
			tmr_ads=tmr_pieces[AD_idx]
			nrml_ads=nrml_pieces[AD_idx]
			tmr_ads_vals=tmr_ads.split(',')
			nrml_ads_vals=nrml_ads.split(',')
			tmr_af=0
			try: 
				tmr_af=float(tmr_ads_vals[1])/(float(tmr_ads_vals[1])+float(tmr_ads_vals[0]))
			except ZeroDivisionError:				
				pass
			nrml_af=0
			try: 
				nrml_af=float(nrml_ads_vals[1])/(float(nrml_ads_vals[1])+float(nrml_ads_vals[0]))
			except ZeroDivisionError:				
				pass
			tmr_pieces[AF_idx]='{:.4f}'.format(tmr_af)
			nrml_pieces[AF_idx]='{:.4f}'.format(nrml_af)
			pieces[tmr_idx]=':'.join(tmr_pieces)
			pieces[nrml_idx]=':'.join(nrml_pieces)
			writer.write('\t'.join(pieces)+'\n')
	writer.close()
	reader.close()

if __name__ == '__main__':
	parser = ArgumentParser(description='computer AF from AD, only AF for indicated TUMOR_STR sample')
	parser.add_argument('IN_VCF',help='input VCF')
	parser.add_argument('OUT_VCF',help='output MAF')
	parser.add_argument('TUMOR_NAME',help='name of tumor BAM (for acquiring proper sample-column)')
	parser.add_argument('NORMAL_NAME',help='name of the normal BAM (for acquiring proper sample-column')
	parser.add_argument('CASE_SAMPLE_NAME',help='name of tumor sample (for replacing TUMOR BAM NAME with it)')
	parser.add_argument('CONTROL_SAMPLE_NAME',help='name of the normal sample (for replacing NORMAL BAM NAME with it)')
	args = parser.parse_args()
	if(args):
		print 'Picked up    IN_VCF  : '+str(args.IN_VCF)
		print 'Picked up   OUT_VCF  : '+str(args.OUT_VCF)
		print 'Picked up TUMOR_NAME  : '+str(args.TUMOR_NAME)
		print 'Picked up NORMAL_NAME : '+str(args.NORMAL_NAME)
		print 'Picked up CASE_SAMPLE_NAME : ' + str(args.CASE_SAMPLE_NAME)
		print 'Picked up CONTROL_SAMPLE_NAME : ' + str(args.CONTROL_SAMPLE_NAME) 
		processVCF(args.IN_VCF,args.OUT_VCF,args.TUMOR_NAME,args.NORMAL_NAME, args.CASE_SAMPLE_NAME, args.CONTROL_SAMPLE_NAME)
	else:
		parser.print_help()