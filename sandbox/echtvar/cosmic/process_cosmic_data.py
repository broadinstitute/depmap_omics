from datetime import datetime
import sys
import gzip
import sh
import progressbar
import argparse


MAX_ALLELE_LEN = 1000

def process_cosmic_db(args):
    start_time = datetime.now()
    args = parse_args()
    cosmic_vcf = args.cosmic_vcf
    cosmic_export = args.cosmic_export

    # get genome and database version from VCF file
    genome_ver, db_ver = get_versions(cosmic_vcf)
    print(f"Genome version detected: {genome_ver}")
    print(f"COSMIC database version: {db_ver}")
    
    if args.exclude_large_vars:
        print("Variants >1000bp will not be included in the output VCF file.")
    
    # cosmic processed file
    cosmic_output_file = f"{cosmic_vcf[:-7]}_{genome_ver}_processed.vcf.gz"

    # process vcf file first
    unique_variants = {}
    process_vcf(cosmic_vcf, unique_variants, args.exclude_large_vars)

    # process mut export file
    sites = {}
    process_mut_export(cosmic_export, sites)

    # merge unique variants and sites count data and write to output file
    get_sites_per_variant(unique_variants, sites, cosmic_output_file, genome_ver, db_ver)
    end_time = datetime.now()
    print(f"Total processing time: {end_time-start_time}")

def process_vcf(vcf_file: str, unique_variants: dict, truncate_large_var: bool):
    gz_decomp = get_gzip_app()
    print('Calculating file size...')
    row_count = int(sh.bash("-c", f"{gz_decomp} {vcf_file} | wc -l"))
    print('Reading cosmic VCF file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for vcf_entry in read_vcf(vcf_file):
            temp_arr = vcf_entry.split('\t')
            chrom = temp_arr[0]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            variant = chrom + ':' + temp_arr[1] + ':' + temp_arr[3] + ':' + temp_arr[4]
            cosmic_id = temp_arr[2]
            # if truncation is activated filter variants with 
            # ref or alt fields > max allele length
            if truncate_large_var:
                if len(temp_arr[3]) <= MAX_ALLELE_LEN and len(temp_arr[4]) <= MAX_ALLELE_LEN:
                    unique_variants.setdefault(variant, []).append(cosmic_id)
            else:
                unique_variants.setdefault(variant, []).append(cosmic_id)
            counter += 1
            pbar.update(counter)
    print('Reading complete.')
    
def process_mut_export(filename: str, sites: dict):
    # run a shell command to create intermediate file
    sites_file = 'sites.tmp'
    print('Preprocessing cosmic mutant export file...')
    # check operating system
    gz_decomp = get_gzip_app()
    # version 95 from COSMIC has strange character encoding issue
    # using LC_ALL=C before the cut command resolves the issue with
    # illegal byte character error
    # sh.bash("-c", f"{gz_decomp} {filename} | LC_ALL=C cut -f 7,8,17 >{sites_file}")
    sh.bash("-c", f"{gz_decomp} {filename} | LC_ALL=C cut -f 7,8,17,29 >{sites_file}")
    print('Calculating file size...')
    row_count = int(sh.bash("-c", f"cat {sites_file} | wc -l"))
    
    # parse the intermediate file
    print ('Parsing intermediate file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for sites_item in read_sites_file(sites_file):
            temp_arr = sites_item.strip('\n').split('\t')
            cosmic_id = temp_arr[2]
            site = temp_arr[1]
            tumor_id = temp_arr[0]
            somatic_status = temp_arr[3]
            #if cosmic_id in sites:
            #    if site in sites[cosmic_id]:
            #        sites[cosmic_id][site].add(tumor_id)
            #    else:
            #        sites[cosmic_id][site] = {tumor_id}
            #else:
            #    temp_dict = {}
            #    temp_dict[site] = {tumor_id}
            #    sites[cosmic_id] = temp_dict
            if cosmic_id in sites:
                if site in sites[cosmic_id]:
                    sites[cosmic_id][site][0].add(tumor_id)
                    sites[cosmic_id][site][1].add(somatic_status)
                else:
                    sites[cosmic_id][site] = [{tumor_id}, {somatic_status}]
            else:
                temp_dict = {}
                temp_dict[site] = [{tumor_id},{somatic_status}]
                sites[cosmic_id] = temp_dict

            counter += 1
            pbar.update(counter)
    print('Parsing complete.')
    print('Removing temporary file...')
    sh.rm(sites_file)
    
def get_sites_per_variant(unique_variants: dict, sites: dict, output_file: str, genome_ver: str, db_ver: str):
    print ('Processing unique variants and associated sites...')
    with gzip.open(output_file, 'wb') as outf:
        # get vcf header and write it to VCF file
        print("Writing VCF header...")
        vcf_header = generate_vcf_header(genome_ver, db_ver)
        outf.write(vcf_header.encode('UTF-8'))
        # output buffer before writing to file
        print("Writing variants...")
        outbuffer = []
        buffer_cap = 10000
        row_count = 0
        # parse the unique variant list
        counter = 0
        var_count = len(unique_variants)
        with progressbar.ProgressBar(max_value=var_count) as pbar:
            missing_cosmic_ids = set()
            for variant, ids in unique_variants.items():
                # hold column element of each row in a list
                var_id_arr = variant.split(':')
                var_line = var_id_arr[:2]
                var_line.append(".")
                var_line.extend(var_id_arr[2:])
                var_line.extend(["32", "PASS"])
                info_field = []
                info_field.append(f"COSMIC_ID={'|'.join(list(set(ids)))}")
                # parse each of the cosmic ids to extract the site counts
                temp_site_info = {}
                for cosmic_id in ids:
                    if cosmic_id in sites:
                        #for site, tumor_ids in sites[cosmic_id].items():
                        #    if site in temp_site_info:
                        #        temp_site_info[site] = temp_site_info[site] | tumor_ids
                        #    else:
                        #        temp_site_info[site] = tumor_ids
                        for site, (tumor_ids, somatic_status) in sites[cosmic_id].items():
                            if site in temp_site_info:
                                temp_site_info[site][0] = temp_site_info[site][0] | tumor_ids
                                temp_site_info[site][1] = temp_site_info[site][1] | somatic_status
                            else:
                                temp_site_info[site] = [tumor_ids, somatic_status]
                    else:
                        missing_cosmic_ids.add(cosmic_id)
                # parse the temp_site_info
                site_list = []
                for site, (tumor_ids, somatic_status) in temp_site_info.items():
                    site_list.append(f"{site}:{str(len(tumor_ids))}:{'-'.join(list(somatic_status))}")
                if len(site_list) == 0:
                    site_list.append(".")
                info_field.append(f"TUMOR_TYPE={'|'.join(site_list)}")
                var_line.append(';'.join(info_field))
                # check if outbuffer is full
                if row_count > buffer_cap:
                    # write to file
                    outf.write(('\n'.join(outbuffer) + '\n').encode('UTF-8'))
                    # clear the buffer
                    del outbuffer[:]
                    # reset row_count
                    row_count = 0
                # add current line to buffer
                outbuffer.append('\t'.join(var_line))
                row_count += 1
                counter += 1
                pbar.update(counter)
                # print 'Progress: {0:.2f}%      \r'.format((counter/var_count) * 100),
        # write the last set of lines to the file
        outf.write(('\n'.join(outbuffer) + '\n').encode('UTF-8'))
        
        missing_cosmic_id_count = len(missing_cosmic_ids)
        with open("missing_cosmic_ids.txt", 'w') as f:
            if missing_cosmic_id_count > 0:
                print(f"{missing_cosmic_id_count} COSMIC IDs were not found in mutant export:")
                for cosmic_id in missing_cosmic_ids:
                    print(cosmic_id, file=f)

        # clear the buffer
        del outbuffer[:]
        print(f'Cosmic files processed. Cosmic mutation and site counts written to {output_file}.')
        
def get_gzip_app():
    gz_app = ''
    if sys.platform.startswith("linux"):
        gz_app = "zcat"
    elif sys.platform.startswith("darwin"):
        gz_app = "gunzip -c"
    else:
        pass
    return gz_app

def get_versions(vcf_file: str):
    genome_ver = ""
    db_ver = ""
    genome_ver_found = False
    db_ver_found = False

    with gzip.open(vcf_file, 'rb') as vcf:
        for bline in vcf:
            line = bline.decode('UTF-8')
            if line.startswith('#'):
                if 'source=COSMIC' in line:
                    db_ver = line.strip().split('=')[1]
                    db_ver_found = True
                elif 'reference=' in line:
                    genome_ver = line.strip().split('=')[1]
                    genome_ver_found = True
            else:
                return genome_ver, db_ver
            if genome_ver_found and db_ver_found:
                return genome_ver, db_ver

def read_vcf(vcf_file: str):
    with gzip.open(vcf_file, 'rb') as vcf:
        for line in vcf:
            linestr = line.decode('UTF-8')
            if not linestr.startswith('#'):
                yield linestr

def read_sites_file(filename: str):
    with open(filename, 'r') as input_file:
        for line in input_file:
            if not 'Primary' in line:
                yield line

def generate_vcf_header(genome_ver: str, db_ver: str):
    vcf_header = []
    vcf_header.append("##fileformat=VCFv4.2")
    vcf_header.append(f"##source={db_ver}|https://cancer.sanger.ac.uk/cosmic")
    vcf_header.append(f"##reference={genome_ver}")
    for c in list(range(1, 23)) + ["X", "Y", "MT"]:
        vcf_header.append(f"##contig=<ID=chr{c}>")
    
    vcf_header.append('##INFO=<ID=COSMIC_ID,Number=1,Type=String,Description="COSMIC unique variant identifier">')
    vcf_header.append('##INFO=<ID=TUMOR_TYPE,Number=1,Type=String,Description="Count of number of occurances of a variant in a set of tumor types">')
    vcf_header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    return '\n'.join(vcf_header) + '\n'

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
                        "cosmic_vcf",
                        type=str,
                        help="Gziped Cosmic VCF file. This can be the coding or non coding COSMIC VCF file (e.g., CosmicCoding.vcf.gz). Required input.")
    parser.add_argument(
                        "cosmic_export",
                        type=str,
                        help="Gziped COSMIC mutant export file (e.g., CosmicMutantExport.tsv.gz). Required input.")
    parser.add_argument(
                        "-l",
                        "--exclude-large-vars",
                        type=str2bool,
                        help="Exclude large variants (>1000bp long). Default is False.",
                        default=False,
                        const=True,
                        nargs="?")
    return parser.parse_args()

def str2bool(arg):
    '''
    This module handles boolean argument input
    Source: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    '''
    if isinstance(arg, bool):
        return arg
    elif arg.lower() in ['yes', 'y', '1', 'true', 't']:
        return True
    elif arg.lower() in ['no', 'n', '0', 'false', 'f']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':
    process_cosmic_db(sys.argv[1:])
