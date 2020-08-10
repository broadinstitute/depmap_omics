from argparse import ArgumentParser
from itertools import izip_longest
import h5py

def grouper(n, iterable, fillvalue=None):   
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def write_targets(targets, OUT_FILE):
    with open(OUT_FILE, 'w') as writer:
        header = ['contig', 'start', 'stop', 'name']
        writer.write('\t'.join(header) + '\n')
        for (chrom, start, stop, name) in targets:
            #name = 'chr' + chrom + ':' + start + '-' + stop
            writer.write('\t'.join([chrom, start, stop, name]) + '\n')

def run_version_7(pon_hdf5):
    original_data = pon_hdf5['original_data']
    raw_targets = original_data['intervals']
    values = raw_targets['transposed_index_start_end'].value
    contigs = [str(int(x) + 1) for x in values[0]]
    start_pos = [str(int(x)) for x in values[1]]
    stop_pos = [str(int(x)) for x in values[2]]
    targets_no_name = zip(contigs, start_pos, stop_pos)
    return [(chrom, start, stop, 'chr' + chrom + ':' + start + '-' + stop) for (chrom, start, stop) in targets_no_name]

        
def run_version_6(pon_hdf5):
    raw_targets = pon_hdf5['raw_targets']
    values = raw_targets['block0_values'].value
    index = raw_targets['index']
    counter = 0
    targets = list()
    print 'len index ', len(index)
    print 'len values ', len(values)
    for i in xrange(len(index)):
        name = index[i]        
        #counter += i
        chrom = values[counter]
        counter += 1
        start = values[counter]
        counter += 1
        stop = values[counter]
        counter += 1
        targets.append((chrom, start, stop, name))
    return targets

def get_pon_version(pon_hdf5):
    version_obj = pon_hdf5['version']
    try:
        return int(version_obj["value"].value[0])
    except:
        try:
            return int(version_obj["values"].value[0])
        except:
            print 'Version of PoN is not known. Accepted version 6 or 7'
            return None

if __name__ == '__main__':
    parser = ArgumentParser(description='Print targets used in the Copy Number PoN')
    parser.add_argument('IN_PON',help='input PON for retrieving targets')
    parser.add_argument('OUT_FILE',help='output File (padded_bed.bed)')
    args = parser.parse_args()
    if(args):
        print 'Picked up    IN_PON  : ' + str(args.IN_PON)
        print 'Picked up   OUT_FILE  : ' + str(args.OUT_FILE)
        try:
            pon_hdf5 = h5py.File(args.IN_PON, "r")
        except:
            print 'Incorrect format of PoN. Accepted version HDF5'
        version = get_pon_version(pon_hdf5)
        if version == 6 or version == 5:
            targets = run_version_6(pon_hdf5)
        elif version == 7:
            targets = run_version_7(pon_hdf5)
        else:
            print 'Version of PoN is not known'
        if not targets == None:
            write_targets(targets, args.OUT_FILE)
    else:
        parser.print_help()