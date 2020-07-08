import sys
import os
import csv
from collections import defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter, FileType


_REQUIRED_MAF_HEADER   = ['Hugo_Symbol',
                          'Entrez_Gene_Id',
                          'Center',
                          'NCBI_Build',
                          'Chromosome',
                          'Start_position',
                          'End_position',
                          'Strand',
                          'Variant_Classification',
                          'Variant_Type',
                          'Reference_Allele',
                          'Tumor_Seq_Allele1',
                          'Tumor_Seq_Allele2',
                          'dbSNP_RS',
                          'dbSNP_Val_Status',
                          'Tumor_Sample_Barcode',
                          'Matched_Norm_Sample_Barcode',
                          'Match_Norm_Seq_Allele1',
                          'Match_Norm_Seq_Allele2',
                          'Tumor_Validation_Allele1',
                          'Tumor_Validation_Allele2',
                          'Match_Norm_Validation_Allele1',
                          'Match_Norm_Validation_Allele2',
                          'Verification_Status',
                          'Validation_Status',
                          'Mutation_Status',
                          'Sequencing_Phase',
                          'Sequence_Source',
                          'Validation_Method',
                          'Score',
                          'BAM_file',
                          'Sequencer']
_REQUIRED_NUM_MAF_COLS = len(_REQUIRED_MAF_HEADER)

PROGRAM_DESCRIPTION = 'This script generates a concatenated MAF.'


#===============================================================================
# Main
#===============================================================================
def main(args):

    parser = ArgumentParser(description=PROGRAM_DESCRIPTION, 
                            formatter_class=RawDescriptionHelpFormatter)
    
    parser.add_argument('strelka_maf', help='File containing Strelka passed indels in MAF format.')
    parser.add_argument('mutect2_maf', help='File containing Mutect2 passed indels in MAF format.') 
    parser.add_argument('maf_Out', help='File output')   
    
    args = parser.parse_args()

    if(args):
        print 'Picked up   strelka_maf  : '+ str(args.strelka_maf)
        print 'Picked up   mutect2_maf  : '+ str(args.mutect2_maf)      
        print 'Picked up   maf_Out  : '+ str(args.maf_Out)      
                
        strelka_header = _get_header(args.strelka_maf)
        strelka_data = _maf_to_dictionary(args.strelka_maf, strelka_header)
        mutect2_header = _get_header(args.mutect2_maf)
        mutect2_data = _maf_to_dictionary(args.mutect2_maf, mutect2_header)             
        merged_data = _merge_mafs(strelka_data, mutect2_data)
        _write_to_file(merged_data, strelka_header, args.maf_Out)
    else:
        parser.print_help()


def _get_header(maf_file):
    ''' Read the beginning of the MAF file to determine the headers. Skip starting
    lines that are empty or are comments (start with '#'). Return the array of
    header field names. 
    '''         
    with open(maf_file, 'rb') as f:    
        for line in f:
            if line.startswith('Hugo'):
                header = line.strip('\n').split('\t')                                   
        if not all(field in header for field in _REQUIRED_MAF_HEADER):        
            raise Exception('MAF file (%s) does not contain required header fields.' % (maf_file))                             
        f.close()
    return header


def _maf_to_dictionary (maf_file, header):
    ''' Create dictionary from maf file 
    where key is the site: Chromosome_Start Position_End Position
    and value is the data for this site that is converted to dictionary
    (key - header field name, value - value for this field at the site)         
    '''     
    data = defaultdict(dict)            
    with open(maf_file, 'rb') as f:    
        for line in f:          
            if (not line.startswith('#') and not line.startswith('Hugo')):              
                row =  line.strip('\n').split('\t')             
                maf_row = dict(zip(header, row))                                                                                    
                key = maf_row['Chromosome'] + '_' + maf_row['Start_position']+ '_' + maf_row['End_position']                                                    
                data[key] = maf_row 
        f.close()
    return data


def _getOverlap(a, b):
    '''
     Returns overlap of two intervals
    '''
    if a[0] == a[1] == b[0] == b[1]:
        return 1
    else:
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    


def _common_site(key, mutect2_data):
    '''
        Find all sites called by Mutect2 that overlap/the same that Strelka's site 
    '''
    common_sites = list()
    strelka_site = key.split('_')
    strelka_start, strelka_end = [int(x) for x in strelka_site[1:]]
    for k in mutect2_data.keys():
        mutect2_site = k.split('_')
        if strelka_site[0] == mutect2_site[0]:
            mutect2_start, mutect2_end = [int(x) for x in mutect2_site[1:]]
            if _getOverlap([strelka_start, strelka_end], [mutect2_start, mutect2_end]) != 0:
                common_sites.append(k)
    return common_sites


def _add_mutect2_annotation(common_sites, mutect2_data):
    '''
        Annotate MuTect2 site
    '''
    site_annotations = list()
    for mutect2_site in common_sites:
        annotation = 'Mutect2_' + mutect2_site + '_' + mutect2_data[mutect2_site]['Reference_Allele'] + '_' + mutect2_data[mutect2_site]['Tumor_Seq_Allele2']
        site_annotations.append(annotation)
    return ' '.join(site_annotations)


def _merge_mafs(strelka_data, mutect2_data):
    '''
     Take only sites called by Strelka but annotate if Mutect2 is also called the same site.
    '''  
    merged_data = defaultdict(dict)
    for k, v in strelka_data.items():
        # find if Mutect2 called the same site
        common_sites = _common_site(k, mutect2_data)
        if len(common_sites )!= 0:          
            v['Other_Variant_Caller'] = _add_mutect2_annotation(common_sites, mutect2_data)
        else:
            v['Other_Variant_Caller'] = ''
        merged_data[k] = v
    return merged_data


def _write_to_file(d, header, filename):
    '''
     Write dictionary to MAF file
    '''     
    other_fields = [field for field in header if field not in _REQUIRED_MAF_HEADER]
    other_fields.append('Other_Variant_Caller')
    with open(filename, 'w') as f:
        f.write('\t'.join(_REQUIRED_MAF_HEADER + other_fields) + '\n')
        for k,v in d.items():
            values = list()
            for field in _REQUIRED_MAF_HEADER:
                values.append(v[field])
            for field in other_fields:
                values.append(v[field])
            f.write('\t'.join(values) + '\n')
        f.close()


#===============================================================================
# Run Main
#===============================================================================
if __name__ == '__main__':
    main(sys.argv[1:])          