import sys
import argparse

def parseOptions():
    description = '''
    Assign context, ref and alt allele depending on the
    name of the Orientation Bias Filter 
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--stub', metavar='stub', type=str, help='filter name (oxog, ffpe).')
    parser.add_argument('-i', '--sample_id', metavar='sample_id', type=str, help='sample_id.')
    args = parser.parse_args()

    return args

def write_to_file(sample_id, stub, description, text):    
    filename = '.'.join([sample_id, description, stub, 'txt'])
    with open(filename, 'w') as writer:
        writer.write(text)

if __name__ == "__main__":

    args = parseOptions()
    stub = args.stub
    sample_id = args.sample_id
    
    if stub == 'oxog':
        write_to_file(sample_id, stub, 'context', '.GG')
        write_to_file(sample_id, stub, 'ref_allele', 'G')
        write_to_file(sample_id, stub, 'artifact_allele', 'T')
        write_to_file(sample_id, stub, 'ref_allele_compliment', 'C')
        write_to_file(sample_id, stub, 'artifact_allele_compliment', 'A')
    elif stub == 'ffpe':
        write_to_file(sample_id, stub, 'context', '.CG')
        write_to_file(sample_id, stub, 'ref_allele', 'C')
        write_to_file(sample_id, stub, 'artifact_allele', 'T')
        write_to_file(sample_id, stub, 'ref_allele_compliment', 'G')
        write_to_file(sample_id, stub, 'artifact_allele_compliment', 'A')
    else:
        print "This stub is not recognized as a valid name for a filter. \n Acceptable values: [oxog, ffpe]"
        sys.exit(1)

