'''
Created on 28 November 2015

@author: stewart

extract Q value from pre_adapter_detail_metrics file
to produce an output file <OUT>/<ID>.annotate_orientationBiasQ with the phred Q value
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
    Parse pre_adapter_detail_metrics  file.
    link file to local area
    write Q value for context to output file
    '''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--sample_id', metavar='sample_id', type=str, help='sample_id.',default='')
    parser.add_argument('-m', '--metrics_file', metavar='metrics_file', type=str, help='metrics_file.')
    parser.add_argument('-a', '--altBase', metavar='altBase', type=str, help='alt Base.',default='T')
    parser.add_argument('-c', '--context', metavar='context', type=str, help='context.', default='.CG')
    parser.add_argument('-o', '--output', metavar='output', type=str, help='output area', default='.')

    args = parser.parse_args()

    return args


class metricsQ():

    def __init__(self, metrics_file,context,altBase):

        self.N = 0
     
        mFP = file(metrics_file, 'r')

        numComments = 0   
        head=mFP.readline()
        while not head.startswith('SAMPLE_ALIAS'):
            head = mFP.readline()
            numComments =  numComments +1 
        
        mFP.seek(0,0)   
        
        #for i in range(0,numComments):
        for _ in range(numComments):
            mFP.readline()
        
        mReader = csv.DictReader(mFP, delimiter='\t') #,fieldnames=oxoFN)
        self.NTOT=0
        self.NALTPRO=0
        self.NALTCON=0
        self.QAVE=0
        for line in mReader:
            CTXT=line['CONTEXT']
            ALT_BASE =line['ALT_BASE']
            if re.match(context, CTXT) and re.match(altBase,ALT_BASE):
                self.N = self.N + 1
                self.NTOT=self.NTOT+int(line['PRO_REF_BASES'])+int(line['PRO_ALT_BASES'])+int(line['CON_REF_BASES'])+int(line['CON_ALT_BASES'])
                self.NALTPRO=self.NALTPRO+int(line['PRO_ALT_BASES'])
                self.NALTCON=self.NALTCON+int(line['CON_ALT_BASES'])
                self.QAVE=self.QAVE+float(line['QSCORE'])
              
                
        if self.N>0:
            self.QAVE=self.QAVE/float(self.N)
            er=float(max(self.NALTPRO-self.NALTCON,1.0001))/float(self.NTOT)
            self.Q=-10.0*log10(abs(er)+1e-10)
    
       
        mFP.close()


if __name__ == "__main__":

    args = parseOptions()
    sample_id = args.sample_id
    metrics_file = args.metrics_file
    context = args.context
    altBase = args.altBase
    output = args.output
   

    if len(sample_id)<1:
        sample_id = os.path.split(metrics_file)[1].replace('.pre_adapter_detail_metrics','')

    if not os.path.exists(output):
        os.makedirs(output)

    if len(metrics_file)<1:
        sys.exit("need -m metrics file")

    outFile = output +'/'+sample_id+'.orientation_BiasQ.txt'
    outFP = file(outFile,'wt')

    m=metricsQ(metrics_file,context,altBase)

    outFP.write('%.2f'  % m.Q)
    print('Q:\n\t %.2f'  % m.Q)

    outFP.close()

    print('\ndone')
