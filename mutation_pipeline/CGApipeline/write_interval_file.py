import sys
import argparse
import csv
import os
import re
from math import log10
import subprocess
from argparse import ArgumentParser

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
    #outFP.write(('chr:start-stop') + '\n')
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
        #if ('judgement' not in line.keys()) or (line['judgement'] == 'KEEP'):
        outFP.write(chrom + ":" + position + '-' + position + "\n")
    
        #else:
            #outRejectFP.write(chrom + ":" + position + "\n")
    outFP.close()


if __name__ == '__main__':
    parser = ArgumentParser(description='Read in the call stats file and for each chr, start row, add it to the interval file (of length 1)')
    parser.add_argument('IN_FILE',help='input File (call stats)')
    parser.add_argument('OUT_FILE',help='output File (intervals)')
    args = parser.parse_args()
    if(args):
        print 'Picked up    IN_FILE  : ' + str(args.IN_FILE)
        print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
        writeIntervalsFile(args.IN_FILE,args.OUT_FILE)
    else:
        parser.print_help()    