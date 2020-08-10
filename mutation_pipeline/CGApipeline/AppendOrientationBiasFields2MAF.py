'''
Created on Feb 15, 2016

@author: stewart
'''

import sys
import argparse
import csv
import os
import re
from math import log10
import subprocess

if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"

def parseOptions():
    description = '''
    Add annotation value to every line in maf file.
    write new maf with <fielname>.maf extension
    '''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--pair_id', metavar='pair_id', type=str, help='id : pair, sample, individual.... ')
    parser.add_argument('-m', '--maf_file', metavar='maf_file', type=str, help='maf_file.')
    parser.add_argument('-b', '--bam_file', metavar='bam_file', type=str, help='bam file.')
    parser.add_argument('-o', '--output', metavar='output', type=str, help='output area', default='.')
    #parser.add_argument('-l', '--lib_dir', metavar='lib_dir', type=str, help='lib_dir.')
    #parser.add_argument('-r', '--reference_fasta', metavar='reference_fasta', type=str, help='reference_fasta.')
    #parser.add_argument('-j', '--java', metavar='java', type=str, help='java version.', default='java')
    #parser.add_argument('-l', '--interval_file', metavar='interval_file', type=str, help='interval file derived from calls stats')
    parser.add_argument('-f', '--info_file', metavar='info_file', type=str, help='info file')    

    args = parser.parse_args()

    return args

def removeCommentLines(inputFP, commentPrepend="#"):
    ''' Removes lines starting with the comment prepend.
    Note that any space before the comment prepend character will throw this method off '''
    resetLocation = inputFP.tell()
    nextChar = inputFP.read(1)
    commentLines = ""

    # Get rid of blank lines
    while nextChar in ['\n', '\r']:
        resetLocation = inputFP.tell()
        nextChar = inputFP.read(1)

    while (nextChar == commentPrepend):
        commentLines = commentLines + (commentPrepend + inputFP.readline())
        resetLocation = inputFP.tell()
        nextChar = inputFP.read(1)

    # Go back one character to make sure that we have moved the file pointer to the beginning of the first non-comment line.
    inputFP.seek(resetLocation, os.SEEK_SET)
    return inputFP
        
def writeIntervalsFile(inputFile, outputIntervalFileName):
    outFP = file(outputIntervalFileName, 'w')
    inputFP = file(inputFile, 'r')
    inputFP = removeCommentLines(inputFP)
    inputTSVReader = csv.DictReader(inputFP, delimiter='\t')
    mouse = False
    # Read in the call stats file and for each chr, start row, add it to the interval file (of length 1)
    #    contig  position
    for line in inputTSVReader:
        #print line.keys()
        if 'contig' in line.keys():
            chrom = line['contig']
            position = line['position']
        elif 'chr' in line.keys():
            chrom = line['chr']
            position = line['start']
        else:
            chrom = line['Chromosome']
            position = line['Start_position']
            # mm8 MAFs need a 'chr' in the interval_list
            if 'NCBI_Build' in line.keys():
                mouse = line['NCBI_Build'] == 'mm9'
                if (mouse):
                    chrom = 'chr' + chrom
        # band-aid broken MT reference dictionary
        # MAF
        if ('M' in chrom) and (not mouse):
            chrom = 'MT'
        if ('judgement' not in line.keys()) or (line['judgement'] == 'KEEP'):
            outFP.write(chrom + ":" + position + "\n")
    
        #else:
            #outRejectFP.write(chrom + ":" + position + "\n")
    outFP.close()
    
def OrientationInfo(filename):
    ''' Given a filename of oxoGMetric data, create a dictionary of chr:position to prune.'''
    result = dict()
    
    oxoGFP = file(filename, 'r')
    oxoTSVReader = csv.DictReader(oxoGFP, delimiter='\t')

    print("Chromosome" + "\t" + "Start_position" + "\t" + "filter" + "\t" + "value" + "\t" +  "tumor_lod" + "\t" + "i_t_NaltArt" + "\t" + "i_t_NaltTot" + "\n")

    for line in oxoTSVReader:
        ref = line['ref']

        A_F1R2 = int(line['F1_A']) + int(line['R2_A'])
        A_F2R1 = int(line['F2_A']) + int(line['R1_A'])
        C_F1R2 = int(line['F1_C']) + int(line['R2_C'])
        C_F2R1 = int(line['F2_C']) + int(line['R1_C'])
        G_F1R2 = int(line['F1_G']) + int(line['R2_G'])
        G_F2R1 = int(line['F2_G']) + int(line['R1_G'])
        T_F1R2 = int(line['F1_T']) + int(line['R2_T'])
        T_F2R1 = int(line['F2_T']) + int(line['R1_T'])

        contig = line['contig']
        contig = contig.replace('chr','').strip(' \t\n\r')
        
        # band-aid dictionary hash key to match convention in main
        if 'M' == contig:
            contig='MT'
        
        result[contig + ":" + line['position']] = A_F1R2,A_F2R1,C_F1R2,C_F2R1,G_F1R2,G_F2R1,T_F1R2,T_F2R1
        #print contig + ":" + line['position'] 
    return result


if __name__ == "__main__":

    args = parseOptions()
    pair_id = args.pair_id
    maf_file = args.maf_file
    bam_file = args.bam_file    
    output = args.output
    info_file = args.info_file
    
    #lib_dir= args.lib_dir
    #reference_fasta=args.reference_fasta
    #java = args.java
    #output_interval_file = args.interval_file
    

    try:
        os.makedirs(output)
    except OSError:
        if not os.path.isdir(output):
            raise
    output_maf_file = output+'/'+pair_id+'.OrientationBiasInfo.maf'

    # These are the new columns added by this script.  Please note that order is very important. 
    additionalColumns = ['i_t_ALT_F1R2','i_t_ALT_F2R1','i_t_REF_F1R2','i_t_REF_F2R1']


    # load the input files
    maf_fileFP = file(maf_file, 'r')


    # Read preceding comment lines until all comment lines have been read.  Save comments to list
    numComments = 0
    commentLines=[]
    line = maf_fileFP.readline()
    while line.startswith('#'):
        numComments = numComments + 1
        commentLines.append(line)
        line = maf_fileFP.readline()
    maf_fileFP.seek(0,0)
    
    # Read until all comment lines have been read.  Discard comment lines
    for i in range(0,numComments):
        maf_fileFP.readline()
    
    inputTSVReader = csv.DictReader(maf_fileFP, delimiter='\t')

    # The output fields are going to be the input 
    maf_fields=inputTSVReader.fieldnames

    if additionalColumns in maf_fields:
        # make soft link and exit ok
        cmd="ln -s " + maf_file + " "+output_maf_file
        subprocess.call(cmd)

    else:
        maf_fileFP.close()
        # Create the output maf file
        outputFileFP = file(output_maf_file, 'w')

        # a) make MAF interval list
        #output_interval_file = output+'/'+pair_id+'.intervals'
        #output_info_file = output+'/'+pair_id+'.orientation_info.txt'

        #if not os.path.exists(output_info_file):

        #    writeIntervalsFile(maf_file, output_interval_file)

            # b) make oxoG info file
            #cmd="/usr/java/jre1.7.0_80/bin/java" #  -Xmx2g -jar /opt/broad_modules/oxoGMetrics/GenomeAnalysisTK.jar --analysis_type "OxoGMetrics" -R ${fasta} -I ${BAM} -L ${aliquot}.oxoG.interval_list -o "${aliquot}.oxoG.metrics.txt"
        #    cmd=java+" -Xmx2g -jar "+lib_dir+"/GenomeAnalysisTK.jar --analysis_type OxoGMetrics -R "+reference_fasta+" -I "+bam_file+" -L "+output_interval_file+" -o "+output_info_file
        #    cmds=cmd.split(" ")
        #    subprocess.call(cmds)

        # c) append fields to MAF
        Info=OrientationInfo(info_file)

        maf_fileFP = file(maf_file, 'r')


        # Read preceding comment lines until all comment lines have been read.  Save comments to list
        line = maf_fileFP.readline()
        maf_fileFP.seek(0,0)

        # Read until all comment lines have been read.  Discard comment lines
        #for i in range(0,numComments):
        #    maf_fileFP.readline()


        for line in commentLines:
            maf_fileFP.readline()
            outputFileFP.write(line)

        inputTSVReader = csv.DictReader(maf_fileFP, delimiter='\t')
        maf_fields=inputTSVReader.fieldnames

        maf_fields.extend(additionalColumns)

        tsvWriter = csv.DictWriter(outputFileFP, maf_fields, delimiter='\t', lineterminator='\n')
        tsvWriter.writeheader()

       # When a mutation is pruned out, rewrite the line.
        for line in inputTSVReader:
            ref_allele = line['Reference_Allele']
            alt_allele = line['Tumor_Seq_Allele2']
        

            i_t_ALT_F1R2=0
            i_t_ALT_F2R1=0
            i_t_REF_F1R2=0
            i_t_REF_F2R1=0

            # band-aid broken MT reference dictionary
            if 'M' is line['Chromosome']:
                line['Chromosome']='MT'

            key = line['Chromosome'] + ":" + line['Start_position']

            A_F1R2,A_F2R1,C_F1R2,C_F2R1,G_F1R2,G_F2R1,T_F1R2,T_F2R1 = Info[key]

            if (alt_allele == 'A'):
                i_t_ALT_F1R2=A_F1R2
                i_t_ALT_F2R1=A_F2R1
            if (alt_allele == 'C'):
                i_t_ALT_F1R2=C_F1R2
                i_t_ALT_F2R1=C_F2R1
            if (alt_allele == 'G'):
                i_t_ALT_F1R2=G_F1R2
                i_t_ALT_F2R1=G_F2R1
            if (alt_allele == 'T'):
                i_t_ALT_F1R2=T_F1R2
                i_t_ALT_F2R1=T_F2R1
            if (ref_allele == 'A'):
                i_t_REF_F1R2=A_F1R2
                i_t_REF_F2R1=A_F2R1
            if (ref_allele == 'C'):
                i_t_REF_F1R2=C_F1R2
                i_t_REF_F2R1=C_F2R1
            if (ref_allele == 'G'):
                i_t_REF_F1R2=G_F1R2
                i_t_REF_F2R1=G_F2R1
            if (ref_allele == 'T'):
                i_t_REF_F1R2=T_F1R2
                i_t_REF_F2R1=T_F2R1

            line1=line;
            line1['i_t_ALT_F1R2'] = str(i_t_ALT_F1R2)
            line1['i_t_ALT_F2R1'] = str(i_t_ALT_F2R1)
            line1['i_t_REF_F1R2'] = str(i_t_REF_F1R2)
            line1['i_t_REF_F2R1'] = str(i_t_REF_F2R1)


            #outputFileFP.write(line['Chromosome'] + "\t" + line['Start_position'] + "\t" + str(filter) + "\n")
            print((line1['Chromosome'] + "\t" + line1['Start_position'] +"\t" + ref_allele+alt_allele+  "\t" + line1['i_t_ALT_F1R2'] + "\t" + line1['i_t_ALT_F2R1']))

            tsvWriter.writerow(line1)

        outputFileFP.close()
        maf_fileFP.close()
