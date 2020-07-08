'''
Created on Apr 9, 2013

@author: stewart
'''

import sys
import argparse
import csv
import os
import re
from math import log10

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
    parser.add_argument('-f', '--fieldname', metavar='fieldname', type=str, help='fieldname.',default='x')
    parser.add_argument('-v', '--value', metavar='value', type=str, help='value.', default='CCG')
    parser.add_argument('-o', '--output', metavar='output', type=str, help='output area', default='.')

    args = parser.parse_args()

    return args

class maf():
    """
    load maf into container object
    store header lines identified with # in first column
    dict with key Chromosome:Start_position, date is maf line
    """
    def __init__(self, fileName):

        self.N = 0
        self.headline=[]
        self.fields=[]
        self.eventDict = {}
        self.sample = ''

        kAF = -1
        kALT = -1
        kREF = -1

        try:
            fh = open(fileName,'rt')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            raise
        
            
        for lines in fh:
            line=lines.strip('\n')
            if line[0] == '#':
                self.headline.append(line)
                continue

            mafL = line.split('\t')
            if mafL[0] == 'Hugo_Symbol':
                self.fields=mafL
                kALT = mafL.index('t_alt_count')
                kREF = mafL.index('t_ref_count')
                continue

            self.N=self.N+1
            
            chrom         = mafL[4]
            pstart        = mafL[5]
            pend          = mafL[6]
            obsAllele     = mafL[11]
            mutationType  = mafL[9]
            patientName   = mafL[15].replace('-Tumor','')
            #af            = float(mafL[kAF])
            altcount      = int(mafL[kALT])
            refcount      = int(mafL[kREF])
            af            = float(altcount)/(float(altcount)+float(refcount))

            key  = chrom + ":" + pstart + "-" + pend + ">" + obsAllele + "." + mutationType
            self.eventDict[key] = [line,patientName,af,altcount,refcount]

            if (self.N<2):
                self.sample = patientName

        fh.close()
        print "Loaded MAF:\t" + fileName
        print "sample:\t" + self.sample
        print "events:\t" + str(self.N)

if __name__ == "__main__":

    args = parseOptions()
    pair_id = args.pair_id
    maf_file= args.maf_file
    fieldname = args.fieldname
    value=args.value
    output = args.output

    if len(output)>1:
        if not os.path.exists(output):
            os.makedirs(output)

    A=maf(maf_file)

    out_maf_file = output  + '/' + pair_id + '.' + fieldname + '.maf.annotated'

    out_FP = open(out_maf_file, 'wt')

    for line in A.headline:
        out_FP.write(line+"\n")

    fields=A.fields
    fields.append('i_'+fieldname)
    out_FP.write("\t".join(fields) +"\n")

    for key in A.eventDict:
            line = A.eventDict[key][0]
            line = line + "\t" + value
            out_FP.write(line +"\n")


    out_FP.close()
    
    print('done')
    
