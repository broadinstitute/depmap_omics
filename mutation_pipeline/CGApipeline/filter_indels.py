#!/usr/bin/env python

from argparse import ArgumentParser


#routine wrapper for code readability
def isDot(s):
	if(s=="."):
		return True
	else:
		return False

#routine wrapper for core purpose and code readability
def isNonIndel(refStr,altStr):
	#if ref not dot and alt not dot, and both are length one, then is non-indel
	if(not(isDot(refStr)) and not(isDot(altStr))):
		#both not dot
		#if(len(refStr)==1 and len(altStr)==1):
		if len(refStr)==len(altStr):
			#is a SNV/non indel
			return True
		else:
			return False
	else:
		return False


#core function of program
def partition_vcf(inVcf,indelVcf,otherVcf):
	VCF=open(inVcf,'r')
	INDELS=open(indelVcf,'w')
	otherVcf=open(otherVcf,'w')
	for line in VCF:
		if(line.startswith("#")):
			#write to both
			INDELS.write(line)
			otherVcf.write(line)
		else:
			#get ref and alt strs and write based on result of isNonIndel
			pieces=line.split('\t')
			refStr=pieces[3]
			altStr=pieces[4]
			indelFlag=not(isNonIndel(refStr,altStr))
			if(indelFlag):
				#write to indel file
				INDELS.write(line)
			else:
				#write to non-indel file
				otherVcf.write(line)
	INDELS.close()
	otherVcf.close()





def main():
# The main argument parser
	parser = ArgumentParser(description="Script for partitioning a VCF into indel/non-indel variants")
	parser.add_argument('VCF', help='input VCF')
	parser.add_argument('INDELS',help="path for output VCF containing only indels")
	parser.add_argument('NON_INDELS',help="path for output VCF containing other variants")
	args = parser.parse_args()
	args=parser.parse_args()
	if(args):
		inVcf=args.VCF
		indelVcf=args.INDELS
		otherVcf=args.NON_INDELS
		partition_vcf(inVcf,indelVcf,otherVcf)
	else:
		parser.print_help()

if __name__ == "__main__":
	main()


