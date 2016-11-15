#!/usr/bin/python
'''
###########################################################################
# Purpose:
#   Generates the master gap file which includes systematic gaps for specified samples 
#   Coverage files generated have the following format: 
#   [chr] [start] [end] [gene_name] [offset] [cov_depth]
#   A systematic gap can be treated as a base which has the folloiwng:
#   mean cov(across samples) - std deviation cov(across samples) < minimum threshold
#
# Usage:
#   See main() for arguments
#
# History:
#   1.0 2016-11-05  
###########################################################################
'''
import argparse
import os, sys
import gzip
import pickle
from collections import defaultdict
import numpy
numpy.set_printoptions(threshold='nan')

__author__          = 'Chris Ieng'
__division__        = 'VCGS'
__organisation__    = 'VCGS'
__date__            = '2016-11-07'
__version__         = '1.0'
__license__         = 'tba'
__description__     = "Generates the master gap file which includes systematic gaps for specified samples"

def check_path(path):
   if not os.path.exists(path):
        print >> sys.stderr, "Couldn't locate path %s, exiting" % path 
        sys.exit(1)

def get_id(cov):
    return cov.split('.cov.gz')[0]

def load_all_numpy_data(folder, chr_order):
    temp_dict = defaultdict(dict)
    for chrom in chr_order:
        print >> sys.stderr, "Loading chr: %s" % chrom
        data = numpy.load(os.path.join(folder, chrom+".npz"))
        for interval in data:
            start, end = interval.split('-')
            start = int(start)
            end = int(end)
            temp_dict[chrom][(start,end)] = data[interval]
    return temp_dict

def load_gzdata(folder, filename):
    print >> sys.stderr, "Loading pickle: %s" % filename
    with gzip.open(os.path.join(folder, filename), "rb") as gzfile:
        data = pickle.load(gzfile)
    gzfile.closed
    return data

def load_all_pickles(folder):
    print >> sys.stderr, "Loading from folder: %s" % folder
    chr_order = load_gzdata(folder, chr_order_filename)
    interval_order = load_gzdata(folder, interval_order_filename)
    batch_dict = load_gzdata(folder, batch_dict_filename)
    return batch_dict, chr_order, interval_order

def mkdir_p(folder):
    if not os.path.exists(folder):
        print >> sys.stderr, "Saving data into %s" % folder
        os.mkdir(folder)

def save_gzdata(folder, filename, data):
    mkdir_p(folder)
    print >> sys.stderr, "Saving pickle data: %s" % filename
    with gzip.open(os.path.join(folder, filename), 'wb') as file_write:
        pickle.dump(data, file_write)
    file_write.closed

def save_all_numpy_data(folder):
    mkdir_p(folder)
    print >> sys.stderr, "Saving numpy data"
    for chrom in t_dict.keys():
        save_dict = {}
        save_dict = {str(start)+"-"+str(end) : v for (start,end),v in t_dict[chrom].items()}
        filename = chrom+".npz"
        print >> sys.stderr, "\tSaving chr: %s data" % chrom
        numpy.savez_compressed(os.path.join(folder, filename), **save_dict)

def save_all(folder):
    mkdir_p(folder)
    print >> sys.stderr, "Saving to folder: %s" % folder
    save_gzdata(folder, chr_order_filename, chr_order)
    save_gzdata(folder, interval_order_filename, interval_order)
    save_all_numpy_data(folder)

def init_bstructs(batch_dict, cov_files):
    print >> sys.stderr, "Initializing batch related structures"
    fidx = len(batch_dict.keys()) 

    with open(cov_files, "rb") as all_coverage_files:
        for batchcov in all_coverage_files:
            check_path(batchcov.rstrip())
            directory, filename = os.path.split(batchcov)
            batch_id = get_id(filename)
            if batch_id in batch_dict:
                print >> sys.stderr, "[Batch: %s has been analysed, not adding]" % batch_id
                continue
            batch_dict[batch_id]['idx'] = fidx
            batch_dict[batch_id]['file_location'] = batchcov.rstrip()
            batch_dict[batch_id]['processed'] = False
            fidx += 1
    all_coverage_files.closed
    return batch_dict

def init_tstructs(target_bed, num_cov_files):
    print >> sys.stderr, "Initializing target bed structures"
    t_dict = defaultdict(dict)
    chr_order = []
    interval_order = defaultdict(list)

    with open(target_bed, "rb") as tb:
        for line in tb:
            toks = line.rstrip('\n').split('\t')
            chrom, start, end, name = toks[0], int(toks[1]), int(toks[2]), toks[3]
            if chrom not in chr_order:
                chr_order.append(chrom)
                print >> sys.stderr, "\tProcessed chr: %s from %s" % (chrom, target_bed)
            interval_order[chrom].append((start, end, name))
            t_dict[chrom][(start,end)] = numpy.zeros((end-start, num_cov_files))
    tb.closed
    return t_dict, chr_order, interval_order

def coerce_cov_type(tokens):
    try:
        chrom  = tokens[0]
        start  = int(tokens[1])
        end    = int(tokens[2])
        name   = tokens[3]
        offset = int(tokens[4])
        cov    = int(tokens[5])
    except:
        raise Exception("Type coersion error - "+'\t'.join(tokens))
    return chrom, start, end, name, offset, cov

def extend_tdict_columns(num_cov_files):
    print >> sys.stderr, "Extending t_dict struct by: %d" % num_cov_files
    for chrom in t_dict.keys():
        for interval in t_dict[chrom].keys():
            data = t_dict[chrom][interval]
            num_rows = data.shape[0]
            zero_array = numpy.zeros((num_rows, num_cov_files))
            t_dict[chrom][interval] = numpy.c_[data, zero_array] 

def process_coverage_files(batch_dict, reprocess=False):
    batches_to_process = []
    extend_count = 0
    for batch in batch_dict:
        if not batch_dict[batch]['processed']:
            batches_to_process.append(batch)
            extend_count+=1
    if reprocess:
        extend_tdict_columns(extend_count)

    for batch in batches_to_process:
        batch_file_location = batch_dict[batch]['file_location']
        idx = batch_dict[batch]['idx']
        with gzip.open(batch_file_location, 'rb') as open_batch_file:
            print >> sys.stderr, "Processing file: %s" % batch
            interval_diff = 0 
            chosen_chr = ""
            prev_chr = ""
            for cov_row in open_batch_file:
                if prev_chr != chosen_chr:  
                    chosen_chr = prev_chr
                    print >> sys.stderr, "\tChromosome: %s being processed" % chosen_chr 
                if interval_diff == 0:
                    chrom, start, end, name, offset, cov = coerce_cov_type(cov_row.rstrip('\n').split('\t'))
                    try:
                        chr_interval = t_dict[chrom][(start,end)]
                    except KeyError, e:
                        chr_interval = {} 
                        print >> sys.stderr, "====================== WARNING ======================"
                        print >> sys.stderr, '[Interval not found! Coverage files possibly analysed with a different target_bed!]' 
                        print >> sys.stderr, '[Cov file: %s]' % batch_file_location
                        print >> sys.stderr, '[Chromosome: %s]' % chrom
                        print >> sys.stderr, '[Interval: %s]' % str(e)
                        print >> sys.stderr, "====================================================="
                    interval_diff = len(chr_interval)
                    prev_chr = chrom
                else:
                    cov = cov_row[cov_row.rfind('\t')+1:].rstrip('\n')
                    offset += 1
                interval_diff -= 1 
                chr_interval[offset-1, idx] = cov
            # Set batch as processed
            print >> sys.stderr, "Setting batch: %s as processed" % batch
            batch_dict[batch]['processed'] = True
        open_batch_file.closed

    print >> sys.stderr, "Processed all .cov.gz files!" 

def process_gaps(chrom, fileptr):
    for start, end, name in interval_order[chrom]:
        diff = end-start
        data = t_dict[chrom][(start,end)]
        means = numpy.mean(data, axis=1)
        stdevs = numpy.std(data, axis=1)

        count_matrix = numpy.zeros((data.shape[0]))
        count_matrix[numpy.where(numpy.subtract(means,stdevs)<args.threshold)] = 1

        # gap_dict[index] = OFFSET_FROM_INDEX
        # This dictionary contains the intervals gaps (if any)
        # If an entry exists it contains the index of the gap and
        # number of offset indexes also contain gaps
        gap_dict = defaultdict(int)
        num_rows = len(count_matrix)
        find_start_idx = True

        for idx in range(num_rows):
            if count_matrix[idx]:
                if find_start_idx:
                    start_idx = idx
                    gap_dict[start_idx] = 0
                next_idx = idx +1
                if (next_idx < num_rows):
                    if count_matrix[next_idx]:
                        gap_dict[start_idx] += 1
                        find_start_idx = False
                        continue
                    else:
                        find_start_idx = True

        # Now print the bed file format
        # chr\tstart\tend\t\gene_name\t\mean
        if gap_dict:
            for idx in sorted(gap_dict.keys()):
                gap_start = idx
                gap_offset = gap_dict[idx]
                if gap_offset:
                    gap_end = gap_start + gap_offset
                    gap_mean = numpy.mean(means[gap_start:gap_end])
                else:
                    gap_mean = numpy.mean(means[gap_start])
                # Adjusting for 1-based BED notation
                interval_start = start + gap_start + 1
                interval_end = interval_start + gap_dict[idx]
                gap_line = chrom+"\t"+str(interval_start)+"\t"+str(interval_end)+"\t"+name+"\t"+str(gap_mean)+"\n"
                if fileptr:
                    fileptr.write(gap_line)
                else:
                    print gap_line 

def parse_args():
    parser = argparse.ArgumentParser(description='Generate the Master gap file')
    parser.add_argument('--target_bed', required=True, help='Target bed file')
    parser.add_argument('--covfile', required=True, help='File with full path to coverage files')
    parser.add_argument('--threshold', required=False, default=20, help='Minimum threshold for check')
    parser.add_argument('-s', '--save_data_folder', required=False, help='Specify folder to save data')
    parser.add_argument('-l', '--load_data_folder', required=False, help='Specify folder to load saved data')
    parser.add_argument('-o', '--output', required=False, help='Save Master Gap output to file')
    parser.add_argument('--skip', required=False, action='store_true', help='skip iinit datastruct steps')
    args = parser.parse_args()

    # Sanity Prints
    print >> sys.stderr, "====================== MASTER GAPS ======================"
    print >> sys.stderr, "Target Bed:", args.target_bed
    print >> sys.stderr, "Coverage File:", args.covfile
    print >> sys.stderr, "Threshold:", args.threshold
    if args.save_data_folder:  print >> sys.stderr, "Post-processing data saving to:", args.save_data_folder
    if args.load_data_folder:  print >> sys.stderr, "Loading existing data from:", args.load_data_folder
    if args.skip:  print >> sys.stderr, "Skip initial datastructure processing"

    return args

if __name__ == '__main__':
    args = parse_args()
    save_dir = ""
    load_dir = "" 

    chr_order_filename = "chr_order.gz"
    interval_order_filename = "interval_order.gz"
    batch_dict_filename = "batch_dict.gz"

    check_path(args.target_bed)
    check_path(args.covfile)

    if args.output:
        output_dir, output_filename = os.path.split(args.output)
        if not output_dir:
            output_dir = os.getcwd()
        print >> sys.stderr, "Master Gap Output directory:", output_dir 
        print >> sys.stderr, "Master Gap Output filename:", output_filename
        check_path(output_dir)
    else:
        print >> sys.stderr, "Master Gap Output: stdout"

    if args.save_data_folder: 
        save_dir = os.path.join(os.getcwd(), args.save_data_folder)

    if args.load_data_folder: 
        load_dir = os.path.join(os.getcwd(), args.load_data_folder)
        check_path(load_dir)
        batch_dict, chr_order, interval_order = load_all_pickles(load_dir)
        t_dict = load_all_numpy_data(load_dir, chr_order)
        if not args.skip:
            batch_dict = init_bstructs(batch_dict, args.covfile)
            process_coverage_files(batch_dict, True) 
    else:
        batch_dict = defaultdict(dict)
        batch_dict = init_bstructs(batch_dict, args.covfile)
        t_dict, chr_order, interval_order = init_tstructs(args.target_bed, len(batch_dict))
        process_coverage_files(batch_dict)

    if save_dir:
        save_gzdata(save_dir, batch_dict_filename, batch_dict)
        save_all(save_dir)
    
    # Main processing
    if args.output:
        with open(os.path.join(output_dir, output_filename), "wb") as master_gap_file:
            for chrom in chr_order:
                process_gaps(chrom, master_gap_file)
        master_gap_file.closed
    else:
        for chrom in chr_order:
            process_gaps(chrom, False) 


