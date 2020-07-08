#Take in a MAF File, and then the following Bams with flags

import sys
import os
import pdb
import pandas as pd
import numpy as np

os.system(". /broad/tools/scripts/useuse; use .zlib-1.2.6")

import csv
import pysam

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

parser = ArgumentParser(description='Process input bam files')
parser.add_argument('--mafsnp', type=str, help='Maf file (snp)')
parser.add_argument('--mafindel', type=str, help='Maf file (indel)', default="None")
parser.add_argument('--wextumor',type=str, help='WEX tumor bam file path', default="None")
parser.add_argument('--wexnormal',type=str, help='WEX normal bam file path', default="None")
parser.add_argument('--wgstumor',type=str, help='WGS tumor bam file path', default="None")
parser.add_argument('--wgsnormal',type=str, help='WGS normal bam file path', default="None")
parser.add_argument('--rnatumor',type=str, help='RNA tumor bam file path', default="None")
parser.add_argument('--targetedtumor',type=str, help='targeted tumor bam file path', default="None")
parser.add_argument('--targetednormal',type=str, help='targeted normal bam file path', default="None")
parser.add_argument('--lowpasstumor',type=str, help='low pass tumor bam file path', default="None")
parser.add_argument('--lowpassnormal',type=str, help='low pass normal bam file path', default="None")
parser.add_argument('--othertumor',type=str, help="other validation tumor bam file path", default="None")
parser.add_argument('--othernormal',type=str, help="other validation normal bam file path", default="None")
parser.add_argument('--out', type=str, help='Output file', default='output_pileup_file' )
parser.add_argument('--rnatype', type=str, help='RNA alignment type', default='hg19' )
parser.add_argument('--minBaseQ', type=str, help='minimum base quality', default='5' )

args = parser.parse_args()
print(args)

##########################################
#DEFINE FUNCTIONS:
##########################################

#Convert maf file to a list of positions:
def maf_to_position(filename): #convert maf to a position file with Chromosome, Start_position, End_position.
    columns_care = ['Chromosome','Start_position','End_position','Start_Position','Variant_Type']
    count=0
    pos = []
    with open(filename, 'rb') as csvfile:
        fullfile = csv.reader(csvfile, delimiter='\t')
        for line in fullfile:
            #print(line)
            #tolerate empty gene names
            if len(line[0])==0 or line[0][0]!='#':
                if count==0:
                    header=line
                    ind = [int(i) for i, x in enumerate(header) if x in columns_care]
                    if len(ind)!=4:
                        raise Exception('Maf file must contain the following headers: Chromosome, Start_position, End_position, Variant_Type') #FIX THIS TO THROW ERROR!
                        exit()
                    count = 1
                else:
                    if line[ind[0]] in ('M','MT'):
                        print('Removing all mitochondria sites')
                    else:
                        pos.append([line[i] for i in ind])
    return pos

#Gernate bam soft links from an input bam file, tumor normal status, and the type of sequencing method: This is used becase pysam requires that the index format is .bam.bai whereas picard by default uses .bai
def create_soft_links(bamfile, sequencing, type_bam, bam_num):
    os.system("mkdir -p ./softlinked/")

    if not os.path.exists(bamfile):
        print(bamfile)
        print("Bam file does not exist! Please check your input files.")
        exit()
    if os.path.exists(bamfile + ".bai"):
        print("Bam file already indexed appropriately for pysam")
        # not sure why there was a ln -3s in this line. Seems to work better with ln -s.... (Chip S. 20 Aug 2015)
        #os.system("ln 3-s " + bamfile + " ./softlinked/" + sequencing + "." + type_bam + "." + bam_num + ".bam")
        os.system("ln -s " + bamfile + " ./softlinked/" + sequencing + "." + type_bam + "." + bam_num + ".bam")
        os.system("ln -s " + bamfile + ".bai" + " ./softlinked/" + sequencing +"." + type_bam + "." + bam_num + ".bam.bai")
    elif os.path.exists(bamfile):
        os.system("ln -s " + bamfile + " ./softlinked/" + sequencing +  "." + type_bam + "." + bam_num + ".bam")
        os.system("ln -s " + bamfile[:-4]  + ".bai ./softlinked/" + sequencing + "." + type_bam + "." + bam_num + ".bam.bai")

#Generate list of bam files from those that had symbolic link:
def dir_to_bam_list(dir_list):
    files = os.listdir(dir_list)
    fls = []
    for f in files:
        if f.endswith(".bam"):
            fls.append(os.path.join(dir_list,f))
    return fls
 

#######################################
#Functions to pull down the mutations:
def position_to_ins_del(read_start, position, cigar):
    #''' [(0, 87), (2, 1), (0, 62), (4, 1)] '''
    result = 0
    cigar_final = ""
    for tpl in cigar:
        cigar_final += str(tpl[0])*tpl[1]
    return cigar_final[(position-read_start)-1]

#Actual pileup:
def validate_mutation(bam_file,chrom,position,minBaseQ,Variant_Type):
    indel_counts = 0
    ins_counts = 0
    del_counts = 0
    for pileupcolumn in bam_file.pileup(chrom, position-2, position):
    # ,truncate=True):
        pos1=position
        # CS INS counts occur at previous position....
        if (Variant_Type == 'INS'):
            pos1=pos1-1
        if pileupcolumn.pos == pos1-1:
            readlist=[]
            for pileupread in pileupcolumn.pileups:
                if (readlist.count(pileupread.alignment.qname)>0):
                    continue
                if pileupread.alignment.is_duplicate or pileupread.alignment.mapq <=5 or pileupread.alignment.is_qcfail or pileupread.alignment.is_secondary or pileupread.alignment.is_unmapped:
                    #print "read failed quality metrics and therefore not counting in indels or snps"
                    continue
                if  pileupread.query_position is not None:
                    if pileupread.alignment.query_qualities[pileupread.query_position]<minBaseQ:
                        continue
                readlist.append(pileupread.alignment.qname)
                if pileupread.indel>0:
                        #print('im here')
                        indel_counts += 1
                        cnt = position_to_ins_del(pileupread.alignment.pos,position,pileupread.alignment.cigar)
                        if str(cnt) == "1":
                            ins_counts += 1
        # count DELs, and SNPS here
        if pileupcolumn.pos == position-1:
            bases = []
            #indel_counts = 0
            #ins_counts = 0
            #del_counts = 0
            #CS keep track of read names: skip the read seen for a second time to minimize double counting
            readlist=[]
            for pileupread in pileupcolumn.pileups:
                if (readlist.count(pileupread.alignment.qname)>0):
                    continue
                if pileupread.alignment.is_duplicate or pileupread.alignment.mapq <=5 or pileupread.alignment.is_qcfail or pileupread.alignment.is_secondary or pileupread.alignment.is_unmapped:
                    #print "read failed quality metrics and therefore not counting in indels or snps"
                    continue
                if  pileupread.query_position is not None:
                    if pileupread.alignment.query_qualities[pileupread.query_position]<minBaseQ:
                        continue
                readlist.append(pileupread.alignment.qname)
                if pileupread.is_del == 1 or pileupread.indel:
                        #print('im here')
                        indel_counts += 1

                        cnt = position_to_ins_del(pileupread.alignment.pos,position,pileupread.alignment.cigar)
                        #if str(cnt) == "1":
                        #    ins_counts += 1
                        if str(cnt) == "2":
                            del_counts += 1
                        continue


                bases.append(pileupread.alignment.seq[pileupread.query_position])
                if pileupread.indel:
                    indel_counts += 1
            return [pileupcolumn.n,bases.count("A"),bases.count("T"),bases.count("C"),bases.count("G"), ins_counts,del_counts]

#writes out the line in the preprocessing file based on the pileup for each position.
def validate_indel_mut_bam_position(sample,type_bam,bam,bam_number,position_list,output_file,rtype,minBaseQ):
    # pileup for indels and point mutations:
#    print(sample + type_bam
    for pos in position_list:
         #print pos
        samfile = pysam.Samfile(bam,"rb")
        tt = list(range(1,22))
        tt.append('X')
        tt.append('Y')
        pileup=None
        if (len(set(pos[0]).intersection(list('GL')))==0):
            if rtype== 'hg19-chr' and sample=='rna':
                chrom = 'chr' + pos[0]
                pileup = validate_mutation(samfile, chrom, int(pos[2]),minBaseQ,pos[3])
            else:
                pileup = validate_mutation(samfile, pos[0], int(pos[2]),minBaseQ,pos[3])

        #if (sample=='rna') & (len(set(pos[0]).intersection(str(tt)))>0) & (len(set(pos[0]).intersection(list('GL')))==0) & (rtype=="hg19-chr") :
        #    chrom = 'chr' + pos[0]
        #    pileup = validate_mutation(samfile, chrom, int(pos[2]),minBaseQ,pos[3])
        #elif (len(set(pos[0]).intersection(list('GL')))==0):
        #    pileup = validate_mutation(samfile, pos[0], int(pos[2]),minBaseQ,pos[3])
        if pileup==None:
            pileup = ['0']*7
        samfile.close()
        output_file.write(sample + "\t" + type_bam + "\t" + bam_number + "\t" + str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + "\t".join([str(x) for x in pileup]) + "\n")


##########################################################################
###################################
#SCRIPT to EXECUTE on INPUT FILES:
####################################

#Convert input maf files to position lists:
print('starting mutation validator preprocess')
print(args.mafsnp)
print(args.mafindel)
print(args.minBaseQ)
#pdb.set_trace()

if os.path.isfile(args.mafsnp) and not os.path.isfile(args.mafindel):
    print("snp only")
    positions = maf_to_position(args.mafsnp)
if not os.path.isfile(args.mafsnp) and os.path.isfile(args.mafindel):
	print("indel only")
	positions = maf_to_position(args.mafindel)
if os.path.isfile(args.mafsnp) and os.path.isfile(args.mafindel):
	print("snp and indel")
	print(args.mafsnp)
	print(args.mafindel)
	posmaf =  maf_to_position(args.mafsnp)
	posindel =  maf_to_position(args.mafindel)
	positions = posmaf + posindel

#Generate symbolic links to the bam files (required due to pysam .bam.bai formatting requirement). for normal samples, will take in either a single bam file or a list of bam files with .txt extension
if os.path.isfile(args.wextumor):
    #would only expect a single file for the tumors (only reason for combining normals is to model moise)
    print('wex tumor is a single bam file')
    create_soft_links(args.wextumor,'wex','tumor','0')

if os.path.isfile(args.wexnormal):
    if args.wexnormal[-3:]=='bam':
        print('wex normal is a single bam file')
        create_soft_links(args.wexnormal,'wex','normal','0')
    if args.wexnormal[-3:]=='txt':
        print('wex normal is a list of bam files')
        a = pd.read_csv(args.wexnormal, header=None, sep='\t')
        number_bams = a.shape[0]
        width = a.shape[1]
        for i in np.arange(number_bams):
            bam_name = a.iloc[i,width-1]
            create_soft_links(bam_name,'wex','normal',i.astype(str))

if os.path.isfile(args.wgstumor):
    create_soft_links(args.wgstumor,'wgs','tumor','0')

if os.path.isfile(args.wgsnormal):
    if args.wgsnormal[-3:]=='bam':
        print('wgs normal is a single bam file')
        create_soft_links(args.wgsnormal,'wgs','normal','0')
    if args.wgsnormal[-3:]=='txt':
        print('wgs normal is a list of bam files')
        a = pd.read_csv(args.wgsnormal, header=None, sep='\t')
        number_bams = a.shape[0]
        width = a.shape[1]
        for i in np.arange(number_bams):
            bam_name = a.iloc[i,width-1]
            print(bam_name)
            create_soft_links(bam_name,'wgs','normal',i.astype(str))

if os.path.isfile(args.rnatumor):
    create_soft_links(args.rnatumor,'rna','tumor','0')

if os.path.isfile(args.targetedtumor):
    create_soft_links(args.targetedtumor,'targeted','tumor','0')

if os.path.isfile(args.targetednormal):
    if args.targetednormal[-3:]=='bam':
        print('targeted normal is a single bam file')
        create_soft_links(args.targetednormal,'targeted','normal','0')
    if args.targetednormal[-3:]=='txt':
        print('targeted normal is a list of bam files')
        a = pd.read_csv(args.targetednormal, header=None, sep='\t')
        number_bams = a.shape[0]
        width = a.shape[1]
        for i in np.arange(number_bams):
            bam_name = a.iloc[i,width-1] #bam file must be in the last column of the tsv.
            create_soft_links(bam_name,'targeted','normal',i.astype(str))

if os.path.isfile(args.lowpasstumor):
    create_soft_links(args.lowpasstumor, 'lowpass','tumor','0')

if os.path.isfile(args.lowpassnormal):
    if args.lowpassnormal[-3:]=='bam':
        print('lowpass normal is a single bam file')
        create_soft_links(args.lowpassnormal,'lowpass','normal','0')
    if args.lowpassnormal[-3:]=='txt':
        print('lowpass normal is a list of bam files')
        a = pd.read_csv(args.lowpassnormal, header=None)
        number_bams = a.shape[0]
        width = a.shape[1]
        for i in np.arange(number_bams):
            bam_name = a.iloc[i,width-1] #bam file must be in the last column of the tsv.
            create_soft_links(bam_name,'lowpass','normal',i.astype(str))

if os.path.isfile(args.othertumor):
    create_soft_links(args.othertumor,'other','tumor','0')

if os.path.isfile(args.othernormal):
    if args.othernormal[-3:]=='bam':
        print('other normal is a single bam file')
        create_soft_links(args.othernormal,'other','normal','0')
    if args.othernormal[-3:]=='txt':
        print('other normal is a list of bam files')
        a = pd.read_csv(args.othernormal, header=None)
        number_bams = a.shape[0]
        width = a.shape[1]
        for i in np.arange(number_bams):
            bam_name = a.iloc[i,width-1] #bam file must be in the last column of the tsv.
            create_soft_links(bam_name,'other','normal',i.astype(str))

#list to pull down info on the softlinked bam files (preprocessed earlier)
files = dir_to_bam_list("./softlinked/")
#if len(files) =<2:
#	print('There are no validation bams in this set. Therefore running this script seems useless')
#	exit()

#Create a preprocessing file to write out pileups, then determine pileup info for each bam
output_file = args.out + ".pileup_preprocessing.txt"
results = open(output_file,"w")
header = ["Sample","Type","Bam_Number", "Chromosome", "Start_position", "End_position", "ref", "A","T","C","G","ins","del"]
results.write("\t".join([str(x) for x in header]) + "\n")

if len(positions)==0:
    print('Warning: You are passing in an empty maf. Only a header will be present in the pileup. You might want to check your input files.')
else:
    for fl in files:
        print(fl)
        bam_number = os.path.basename(fl).strip().split(".")[2]
        detailed_type = os.path.basename(fl).strip().split(".")[0] + "." + os.path.basename(fl).strip().split(".")[1]
        validate_indel_mut_bam_position(os.path.basename(fl).strip().split(".")[0],detailed_type,fl,bam_number,positions,results, args.rnatype,int(args.minBaseQ))
        #validate_indel_mut_bam_position(os.path.basename(fl).strip().split(".")[0],os.path.basename(fl).strip().split(".")[1],fl,positions,results, args.rnatype)

results.close()

raw_pileup = pd.read_csv(output_file, sep='\t')
summed_pileup = raw_pileup.groupby(["Sample", "Type","Chromosome","Start_position","End_position",]).sum()

summed_output = args.out + ".pileup_preprocessing.summed.txt"
summed_pileup.to_csv(summed_output, sep='\t')


os.system("rm -rf ./softlinked/")

#########################################
#EXAMPLE TO IMPLEMENT: THCA

#python validation_wrapper.py --mafsnp /xchip/cga/gdac-prod/cga/jobResults/filter_tsv_by_field/THCA-TCGA-BJ-A191-TP-NB-SM-1WK7Q-SM-1WK7S/4771644/THCA-TCGA-BJ-A191-TP-NB-SM-1WK7Q-SM-1WK7S.cross_center.SNP.maf --mafindel /xchip/cga/gdac-prod/genepattern/jobResults/6652540/THCA-TCGA-BJ-A191-Tumor-SM-1WK7Q.maf.annotated	--wextumor /seq/picard_aggregation/C499/TCGA-BJ-A191-01A-11D-A13W-08/v3/TCGA-BJ-A191-01A-11D-A13W-08.bam --wexnormal /seq/picard_aggregation/C499/TCGA-BJ-A191-10A-01D-A13W-08/v3/TCGA-BJ-A191-10A-01D-A13W-08.bam --wgstumor /seq/picard_aggregation/G32522/TCGA-BJ-A191-01A-11D-A13W-08/v1/TCGA-BJ-A191-01A-11D-A13W-08.bam --wgsnormal /seq/picard_aggregation/G32522/TCGA-BJ-A191-10A-01D-A13W-08/v1/TCGA-BJ-A191-10A-01D-A13W-08.bam --rnatumor /xchip/cga/gdac-prod/genepattern/jobResults/7105183/THCA-TCGA-BJ-A191-Tumor-SM-1WK7Q.recalibrated.bam


#python validation_wrapper_firehose_library_hack.py --mafsnp /local/cga-fh/cga/M2_validation_LUAD/Pair/LUAD-TCGA-05-4249-TP-NB-SM-1OIMZ-SM-1OINA/jobs/oncotator/m2/job.60104880/LUAD-TCGA-05-4249-TP-NB-SM-1OIMZ-SM-1OINA.m2_pass.maf --mafindel='None' --othernormal='None' --othertumor='None' --out='LUAD-TCGA-05-4249-TP-NB-SM-1OIMZ-SM-1OINA' --rnatumor='None' --targetednormal='/seq/picard_aggregation/C1524/TCGA-05-4249-10A-01D-1105-08/current/TCGA-05-4249-10A-01D-1105-08.bam' --targetedtumor='/seq/picard_aggregation/C1524/TCGA-05-4249-01A-01D-1105-08/current/TCGA-05-4249-01A-01D-1105-08.bam' --wexnormal='/seq/picard_aggregation/C347/TCGA-05-4249-10A-01D-1105-08/v6/TCGA-05-4249-10A-01D-1105-08.bam' --wextumor='/seq/picard_aggregation/C347/TCGA-05-4249-01A-01D-1105-08/v6/TCGA-05-4249-01A-01D-1105-08.bam' --wgsnormal='None' --wgstumor='None'

##################################################
