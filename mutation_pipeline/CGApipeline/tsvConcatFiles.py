'''
Created on Oct 14, 2012

@author: lichtens
'''
import subprocess
import argparse
import sys
from tsvcat import tsvcat
import codecs

def parseOptions():
    description = '''Given a list of files as arguments and an output filename, append the input files into a bigger tsv and write it out to outputFilename'''
    epilog= '''
    This script passes a list of files through to tsvcat.py.  
    
    This script is necessary to get a large number of files processed in Firehose
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filenames', metavar='filenames', type=str, nargs='+',
                   help='filenames to concatenate.  Each as a separate parameter.')

    parser.add_argument('--outputFilename', metavar='outputFilename',  type=str, help ='tsv file name for output.')
    
    args = parser.parse_args()
    
    return args 

if __name__ == '__main__':
    args = parseOptions()
    inputFilenames = args.filenames
    outputFilename = args.outputFilename
    
    if outputFilename is None:
        raise StandardError('outfilename parameter is required (--outputFilename)')
    if outputFilename in inputFilenames:
        raise StandardError('outfilename parameter cannot also be one of the input files.')
    
    print('Appending ' + str(inputFilenames) + ' to ' + outputFilename)
    # fpOut = file(outputFilename, 'w')
    fpOut = codecs.open(outputFilename, 'w', encoding="iso-8859-1")
    if len(inputFilenames) > 0:
        tsvcat(inputFilenames, fpOut)
    pass
