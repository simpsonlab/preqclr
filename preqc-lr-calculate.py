#!/usr/bin/env python
"""Generate data to be used for the preqc-lr report
"""
from Bio import SeqIO
import subprocess, math, argparse
import matplotlib as MPL
MPL.use('Agg')
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pylab as plt
import sys, os, csv
import json, gzip
from operator import itemgetter
from scipy import stats
import time as time

def main():
    # --------------------------------------------------------
    # PART 0: Parse the input
    # --------------------------------------------------------
    parser = argparse.ArgumentParser(description='Display Pre-QC Long Read report')
    parser.add_argument('-i', '--input', action="store", required=True, dest="fa_filename", help="Fasta/q files containing reads")
    parser.add_argument('-o', '--output', action="store", dest="prefix", default="preqc-lr-output", help="Prefix for output pdf")
    parser.add_argument('-g', '--genome_size', action="store", required=True, dest="genome_size", help="Genome size. Default = 4641652 bps.")
    parser.add_argument('-t', '--type', action="store", required=True, dest="data_type", choices=['pb', 'ont'], help="Either pacbio or ONT.")
    parser.add_argument('-c', '--cov', action="store", required=True, dest="cov_reads_to_ref", help="Ref-SimReads coverage results in csv format with 2 columns: coverage, count")
    parser.add_argument('-n', '--sample_name', action="store", required=True, dest="sample_name", help="Sample name. You can use the name of species for example. No spaces allowed.")
    args = parser.parse_args()

    # --------------------------------------------------------
    # PART 1: Get all the information needed from fasta
    # --------------------------------------------------------
    print detect_filetype(args.fa_filename)
    
    # --------------------------------------------------------
    # PART 2: Initiate json obect, and record current info
    # --------------------------------------------------------
    data = {}
    data['sample_name'] = args.sample_name
    data['genome_size'] = args.genome_size
    data['data_type'] = args.data_type
    data['cov_reads_to_ref'] = args.cov_reads_to_ref
    data['fa_filename'] = args.fa_filename
    data['output_prefix'] = args.prefix
    
    # --------------------------------------------------------
    # PART 3: Create all folders
    # --------------------------------------------------------
    work_dir = './' + args.prefix
    overlaps_dir = './' + args.prefix + '/overlaps'
    downsample_dir = './' + args.prefix + '/downsample'
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(overlaps_dir):
        os.makedirs(overlaps_dir)
    if not os.path.exists(downsample_dir):
        os.makedirs(downsample_dir)

    # --------------------------------------------------------
    # PART 4: Create the preqclr report data
    # --------------------------------------------------------
    calculate_report(args.prefix, args.fa_filename, args.genome_size, args.data_type, args.cov_reads_to_ref, data)

class fasta_file:
    def __init__(self, fa_filename, read_seqs, read_lengths, genome_size, num_reads, mean_read_length, data_type):
        self.fa_filename = fa_filename
        self.read_seqs = read_seqs
        self.read_lengths = read_lengths
        self.genome_size = genome_size
        self.num_reads = num_reads
        self.mean_read_length = mean_read_length
        self.data_type = data_type

    def get_fa_filename(self):
        return self.fa_filename

    def get_read_seqs(self):
        return self.read_seqs

    def get_read_lengths(self):
        return self.read_lengths

    def get_genome_size(self):
        return self.genome_size

    def get_num_reads(self):
        return self.num_reads

    def get_mean_read_length(self):
        return self.mean_read_length

    def get_data_type(self):
        return self.data_type


def calculate_report(output_prefix, fa_filename, genome_size, data_type, cov_reads_to_ref, data):
    # --------------------------------------------------------
    # PART 0: Detect the file type, parse file, save to dict
    # --------------------------------------------------------
    file_type = detect_filetype(fa_filename)
    fa_sequences = dict(list())
    if ".gz" in file_type:
        with gzip.open(fa_filename, "rt") as handle:
            if "fasta.gz" in file_type:
                fa_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))		
            elif "fastq.gz" in file_type:
    		fa_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
    else:
        with open(fa_filename, "rt") as handle:
            if "fasta" in file_type:
                fa_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            elif "fastq" in file_type:
                fa_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))

    # --------------------------------------------------------
    # PART 1: Create all necessary info for fasta object
    # --------------------------------------------------------
    # read_seqs will contain all of the read sequences w/ read_ids
    # read_lengths will contain all the read lengths w/out read_ids
    # we will also need mean_read_length later for downsampling and for calculating max coverage for later on
    read_seqs = dict()
    read_lengths = list()
    total_num_bases = 0
    for read_id in fa_sequences:
        read_sequence, read_length = str(fa_sequences[read_id].seq), len(str(fa_sequences[read_id].seq))
        read_seqs[read_id] = read_sequence
        read_lengths.append(read_length)
        total_num_bases+=read_length
    mean_read_length = float(total_num_bases)/float(len(read_seqs))
    num_reads = len(read_seqs)
    # initialize a fasta object
    fasta = fasta_file(fa_filename, read_seqs, read_lengths, genome_size, num_reads, mean_read_length, data_type)

    # --------------------------------------------------------
    # PART 2: Let the calculations begin...
    # --------------------------------------------------------
    calculate_read_length(fasta, output_prefix, data)
    calculate_average_overlaps_vs_coverage(fasta, output_prefix, data)
    calculate_num_overlaps_per_read(fasta, output_prefix, data)
    calculate_estimated_coverage(fasta, output_prefix, cov_reads_to_ref, data)
    calculate_GC_content_per_read(fasta, output_prefix, data)
    calculate_expected_minimum_fractional_overlaps_vs_read_length(fasta, output_prefix, data)
#    print data['num_matching_residues']
    # --------------------------------------------------------
    # PART 3: After all calculations are done, store in json
    # --------------------------------------------------------
    preqclr_data = './' + output_prefix + '.preqclr'
    with open(preqclr_data, 'w') as outfile:  
    	json.dump(data, outfile, indent=4)

def detect_filetype(fa_filename):
    # gets the filetype
    path = fa_filename
    for ext in ['fastq.gz', 'fasta.gz', 'fastq', 'fasta']:
	if path.endswith(ext):
		return ext
    print "Must be either fasta, fastq, fasta.gz, fastq.gz"
    sys.exit(1)

def calculate_read_length(fasta, output_prefix, data):
    print "\n\n\n\n"
    print "Calculating read length distribution"
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Gets the list of pre-calculated read lengths
    # --------------------------------------------------------
    read_lengths = fasta.get_read_lengths()

    # --------------------------------------------------------
    # PART 1: Add to the data set
    # --------------------------------------------------------
    data['read_lengths'] = read_lengths

def calculate_average_overlaps_vs_coverage(fasta, output_prefix, data):
    print "\n\n\n\n"
    print "Calculating average overlaps vs coverage"
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    genome_size = fasta.get_genome_size()
    fa_filename = fasta.get_fa_filename()
    data_type = fasta.get_data_type()
    mean_read_length = fasta.get_mean_read_length()
    total_num_reads = fasta.get_num_reads()

    # --------------------------------------------------------
    # PART 1: Calculate the maximum coverage
    # --------------------------------------------------------
    n = int(total_num_reads)
    l = float(mean_read_length)
    g = float(genome_size)
    max_coverage = round(math.ceil((n * l) / g))
    
    # --------------------------------------------------------
    # PART 2: Calculate the coverage levels we will assess
    # --------------------------------------------------------
    # max_coverage is the maximum coverage
    # coverage_levels is a list that will hold the different coverage levels
    if max_coverage < 10: 
        print "Coverage too small for plot"
    else:
        coverage_levels = [2, 5]

        i = 10
        while i < max_coverage:
            coverage_levels.append(i)
            i = i * 2
    coverage_levels.append(max_coverage)

    # dictionary to store avg num of overlaps for each coverage level
    overlaps_per_coverage_level = dict()
    for c in coverage_levels:
    # --------------------------------------------------------
    # PART 3: Downsample files with seqtk
    # --------------------------------------------------------
    # num_reads will hold the number of reads needed to achieve coverage
    # this is calculated by the equation: num_reads = (coverage*genome_size)/read
    # downsample_file will be a fasta/q file containing the random subset of num_reads
        g = float(genome_size)
        l = float(mean_read_length)
        num_reads = round(math.ceil(float((c * g)/l)))

	downsample_file = create_downsample_file(fa_filename, output_prefix, c, num_reads)

    # --------------------------------------------------------
    # PART 4: Get overlap files with minimap2
    # --------------------------------------------------------
    # overlaps_filename will be the PAF file generated by minimap2
        overlaps_filename = create_overlaps_file(downsample_file, output_prefix, c, data_type)

    # --------------------------------------------------------
    # PART 5: Get the number of overlaps for each read
    # --------------------------------------------------------
    # reads_overlap_count will hold the number of overlaps per read
    # overlaps is the opened PAF file
    # each overlap is indicated as a new line in the PAF file
        reads_overlap_count = dict()
        overlaps = open(overlaps_filename, 'r')
        for line in overlaps:
            query_read_id = line.split('\t')[0]
            target_read_id = line.split('\t')[5]

            # do not count overlaps with self
            if target_read_id != query_read_id:    
                if query_read_id in reads_overlap_count:
                    reads_overlap_count[query_read_id] += 1
                else:
                    reads_overlap_count[query_read_id] = 1
        overlaps.close()

    # --------------------------------------------------------
    # PART 6: Get the average number of overlaps at cov level
    # --------------------------------------------------------
        total_overlaps = 0
        for query_read_id in reads_overlap_count:
            total_overlaps += reads_overlap_count[query_read_id]
        average_num_overlaps = total_overlaps/num_reads

        # record the information in dictionary
        overlaps_per_coverage_level[c] = average_num_overlaps

    # --------------------------------------------------------
    # PART 7: Add to the data set
    # --------------------------------------------------------
    data['coverage_avg_num_overlaps'] = overlaps_per_coverage_level

def calculate_num_overlaps_per_read(fasta, output_prefix, data):
    print "\n\n\n\n"
    print "Calculating number of overlaps per reads"
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    fa_filename = fasta.get_fa_filename()
    genome_size = fasta.get_genome_size()
    data_type = fasta.get_data_type()
    reads = fasta.get_read_seqs()
    mean_read_length = fasta.get_mean_read_length()
    num_reads = fasta.get_num_reads()

    # get overlaps file created with all reads
    max_coverage = round(math.ceil((int(num_reads) * float(mean_read_length)) / float(genome_size)))
    overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)

    # --------------------------------------------------------
    # PART 1: Get the num of overlaps per read
    # --------------------------------------------------------
    # overlaps is the opened PAF file
    # num_overlaps_per_read will, in the end, hold the total number of overlaps
    overlaps = open(overlaps_filename, 'r')    
    num_overlaps_per_read = dict()
    for line in overlaps:
        query_read_id = line.split('\t')[0]
        target_read_id = line.split('\t')[5]

        # do not count overlaps with self
        if target_read_id != query_read_id:
            if query_read_id in num_overlaps_per_read:
                num_overlaps_per_read[query_read_id] += 1
            else:
                num_overlaps_per_read[query_read_id] = 1
    overlaps.close()
        
    # --------------------------------------------------------
    # PART 2: Divide the num overlaps by the read length
    # --------------------------------------------------------
    for read_id in num_overlaps_per_read:
        read_length = float(len(reads[read_id]))
        num_overlaps = float(num_overlaps_per_read[read_id])
        num_overlaps_per_read[read_id] = num_overlaps / read_length 
    
    # --------------------------------------------------------
    # PART 3: Add to the data set
    # --------------------------------------------------------
    data['num_overlaps_per_read'] = num_overlaps_per_read

def calculate_estimated_coverage(fasta, output_prefix, cov_reads_to_ref, data):
    print "\n\n\n\n"
    print "Calculating average coverage per read"
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    fa_filename = fasta.get_fa_filename()
    genome_size = fasta.get_genome_size()
    data_type = fasta.get_data_type()
    reads = fasta.get_read_seqs()
    mean_read_length = fasta.get_mean_read_length()
    total_num_reads = fasta.get_num_reads()

    # --------------------------------------------------------
    # PART 1: get the PAF file calculated w all reads
    # --------------------------------------------------------
    # we get this by calculating the max_coverage
    n = int(total_num_reads)
    l = float(mean_read_length)
    g = float(genome_size)
    max_coverage = round(math.ceil((n * l) / g))
    overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)
    print overlaps_filename

    # --------------------------------------------------------
    # PART 2: get the overlap lengths for each overlap
    # --------------------------------------------------------
    sum_overlap_lengths = dict()
    overlap_accuracy = dict()			# key = read id, values = (sum_overlap_length, sum_matches )
    overlaps = open(overlaps_filename, "r")
    for line in overlaps:
        query_read_id = line.split('\t')[0]
        target_read_id = line.split('\t')[5]
        query_start_pos = int(line.split('\t')[2])
        query_end_pos = int(line.split('\t')[3])
        target_start_pos = int(line.split('\t')[7])
        target_end_pos = int(line.split('\t')[8])
        query_length = int(line.split('\t')[1])
        target_length = int(line.split('\t')[6])
        strand = line.split('\t')[4]
	num_matches = int(line.split('\t')[9])
        if not (query_read_id == target_read_id):
                query_prefix_len = query_start_pos
                query_suffix_len = query_length - query_end_pos

                target_prefix_len = target_start_pos
                target_suffix_len = target_length - target_end_pos

                # calculate length of overlap
                overlap_length = query_end_pos - query_start_pos

                left_clip = 0
                right_clip = 0
                if not (query_start_pos == 0) and not (target_start_pos == 0) :
                        if strand == "+":
                                left_clip = min(query_prefix_len, target_prefix_len)
                        else:
                                left_clip = min(query_prefix_len, target_suffix_len)
                if not (query_end_pos == 0) and not (target_end_pos == 0) :
                        if strand == "+":
                                right_clip = min(query_suffix_len, target_suffix_len)
                        else:
                                right_clip = min(query_suffix_len, target_prefix_len)
                overlap_length += left_clip + right_clip

    # --------------------------------------------------------
    # PART 3: add to the current recorded sum overlap length
    # --------------------------------------------------------
                if query_read_id in sum_overlap_lengths:
                        sum_overlap_lengths[query_read_id]+=overlap_length
			current_overlap_length = overlap_accuracy[query_read_id][0]
			new_overlap_length = current_overlap_length + overlap_length
			current_num_matches = overlap_accuracy[query_read_id][1]
			new_num_matches = current_num_matches + num_matches
			overlap_accuracy[query_read_id] = (new_overlap_length, new_num_matches)
                else:
                        sum_overlap_lengths[query_read_id]=overlap_length
			overlap_accuracy[query_read_id] = (overlap_length, num_matches)
                if target_read_id in sum_overlap_lengths:
                        sum_overlap_lengths[target_read_id]+=overlap_length
			current_overlap_length = overlap_accuracy[target_read_id][0]
                        new_overlap_length = current_overlap_length + overlap_length
                        current_num_matches = overlap_accuracy[target_read_id][1]
                        new_num_matches = current_num_matches + num_matches
                        overlap_accuracy[target_read_id] = (new_overlap_length, new_num_matches)
                else:
                        sum_overlap_lengths[target_read_id]=overlap_length
			overlap_accuracy[target_read_id] = (float(overlap_length), float(num_matches))
    overlaps.close()

    # --------------------------------------------------------
    # PART 4: for each coverage level count num reads w/ cov
    # --------------------------------------------------------
    # This part is also recording which reads fall under each coverage level,
    # this part is also generating data for scatter plot with read length vs estimated coverage per read.
    # per_cov_read_count will hold the read counts for each coverage level
    # data_seqs will hold the list of read ids for each coverage level
    # read_length_and_estimated_cov will hold the data points for the scatter plot for read length vs estimated coverage plot
    per_cov_read_count = dict()         	# key = avg coverage, value = number of reads with this avg coverage
    covs = list()
    accuracies = list()
    read_length_and_estimated_cov = list()    	# list of tuples (read_length, estimated_coverage)
    for read_id in sum_overlap_lengths:
    	seq = reads[read_id]
        length = len(seq)
        num_bases = float(len(seq))
        total_overlaps_length = float(sum_overlap_lengths[read_id])
        cov = int(round(float(total_overlaps_length) / float(num_bases)))
	accuracy = (overlap_accuracy[read_id][1])/float(overlap_accuracy[read_id][0])
	covs.append((cov, length))
	accuracies.append((accuracy, length))
    	if cov in per_cov_read_count:
        	per_cov_read_count[cov]+=1
	else:
        	per_cov_read_count[cov]=1
        read_length_and_estimated_cov.append((num_bases, cov))
    print "Finished coverage level count num reads with cov"    

    print "Length of read_lengths_estimated_cov: " +  str(len(read_length_and_estimated_cov))
    # --------------------------------------------------------
    # PART 5: Add the reads-to-ref data set by Hamza
    # --------------------------------------------------------
    reader = csv.reader(open(cov_reads_to_ref, 'r'))
    read_ref = {}
    next(reader) # skip first line
    for row in reader:
   	k, v = row
  	read_ref[int(k)] = int(v)

    # --------------------------------------------------------
    # PART 6: Filter outliers
    # --------------------------------------------------------
    max_cov = max(covs, key=itemgetter(0))[0]
    min_cov = min(covs, key=itemgetter(0))[0]
    print "Max cov: " + str(max_cov)
    print "Min cov: " + str(min_cov)
    q75, q25 = np.percentile([x[0] for x in covs], [75 ,25])
    IQR = float(q75) - float(q25)
    upperbound = q75 + IQR * 1.5
    lowerbound = q25 - IQR * 1.5
    print "Q75: " + str(q75)
    print "Q25: " + str(q25)
    print "Upperbound: " + str(upperbound)
    print "Lowerbound: " + str(lowerbound)
    temp = [i for i in covs if i[0] < upperbound]
    filtered_covs = [i for i in temp if i[0] > lowerbound]
    max_cov = max(filtered_covs, key=itemgetter(0))[0]
    min_cov = min(filtered_covs, key=itemgetter(0))[0]
    print "Max cov: " + str(max_cov)
    print "Min cov: " + str(min_cov)

    # --------------------------------------------------------
    # PART 7: Estimate genome size based off of coverage
    # --------------------------------------------------------
    mean_cov = float(np.mean([x[0] for x in filtered_covs]))
    median_cov = float(np.median([x[0] for x in filtered_covs]))
    mode_cov = stats.mode([x[0] for x in filtered_covs])
    mean_read_length = float(np.mean([x[1] for x in filtered_covs]))
    median_read_length = float(np.median([x[1] for x in filtered_covs]))
    mode_read_length = stats.mode([x[1] for x in filtered_covs])
    print "Cov (mean, median, mode): (" + str(mean_cov) + "," + str(median_cov) + "," + str(mode_cov) + ")"
    print "Read length (mean, median, mode): (" + str(mean_read_length) + "," + str(median_read_length) + "," + str(mode_read_length) + ")"   

    num_reads_used = float(len(filtered_covs))
    n = float(len(filtered_covs))
    l = float(mean_read_length)
    c = float(np.median([x[0] for x in filtered_covs]))
    estimated_genome_size = ( n * l ) / c

    # --------------------------------------------------------
    # PART 8: Estimate number of islands
    # --------------------------------------------------------
    T = 100.0			# amount of overlap in base pairs needed to detect overlap
    theta = T / l  		# expected minimum fractional overlap required between two clones
    print "T: " + str(T)
    print "Mean read length:" + str(l) 
    print "Theta: " + str(theta)
    sigma = 1 - theta
    estimated_num_islands = ( ( g * c ) / l ) * math.exp(-(1 - theta) * c)
    estimated_num_islands_1 = n/math.exp(sigma * c) - n/math.exp(2*sigma*c)
    print estimated_num_islands
    print estimated_num_islands_1
    # --------------------------------------------------------
    # PART 9: Add to the data set
    # --------------------------------------------------------
    data['estimated_coverage'] = (read_ref, [x[0] for x in filtered_covs])
    data['read_lengths_estimated_cov'] = read_length_and_estimated_cov	
    data['estimated_genome_size'] = estimated_genome_size
    data['overlap_accuracies'] = accuracies
    data['estimated_num_islands'] = estimated_num_islands

def calculate_GC_content_per_read(fasta, output_prefix, data):
    # ========================================================
    # Purpose: Calculate the mean GC content (%) per read,
    # this will be used in plot with mean GC content vs freq. 
    # ========================================================
    print "\n\n\n\n"
    print "Calculating GC-content per read"
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    reads = fasta.get_read_seqs()

    # --------------------------------------------------------
    # PART 1: Get the GC content per read, and store
    # --------------------------------------------------------
    read_counts_per_GC_content = dict()
    for read_id in reads:
        seq = reads[read_id]
        length = len(seq)
        count_C = seq.count('C')
        count_G = seq.count('G')
    	count_GC = count_C + count_G	
    	GC_content = math.ceil((float(count_GC)/float(length))*100)
    	if GC_content in read_counts_per_GC_content:
		read_counts_per_GC_content[GC_content]+=1
    	else:
		read_counts_per_GC_content[GC_content]=1

    data["read_counts_per_GC_content"] = read_counts_per_GC_content			

'''def calculate_homopolymer_distribution(fasta, output_prefix, data):
    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    reads = fasta.get_read_seqs()

    for read_id in reads:
	seq = reads[read_id]
	length = len(seq)
	for i in range(20):
		key = 'A' + str(i)
		j = 0
 		hp = ''
		while j < i
			hp = hp + 'A'
'''		
	

def calculate_expected_minimum_fractional_overlaps_vs_read_length(fasta, output_prefix, data):
    # get the distribution of matching residues
    # get the PAF file with the max number of 
    print "\n\n\n\n"
    print "Calculating num of matches per overlap distribution "
    print "___________________________________________"

    # --------------------------------------------------------
    # PART 0: Get all the information needed from fasta
    # --------------------------------------------------------
    fa_filename = fasta.get_fa_filename()
    genome_size = fasta.get_genome_size()
    data_type = fasta.get_data_type()
    reads = fasta.get_read_seqs()
    mean_read_length = fasta.get_mean_read_length()
    total_num_reads = fasta.get_num_reads()

    # --------------------------------------------------------
    # PART 1: get the PAF file calculated w all reads
    # --------------------------------------------------------
    # we get this by calculating the max_coverage
    n = int(total_num_reads)
    l = float(mean_read_length)
    g = float(genome_size)
    max_coverage = round(math.ceil((n * l) / g))
    overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)
    print overlaps_filename
    overlaps = open(overlaps_filename, "r")
    num = list()
    for line in overlaps:
        num_matching_residues = int(line.split('\t')[9])
	num.append(num_matching_residues)

    data['num_matching_residues'] = num
    
def create_downsample_file(fa_filename, output_prefix, coverage, num_reads):
        downsample_filename = "./" + output_prefix + "/downsample/" + output_prefix + "_" + str(coverage) + "X.fasta"
        downsample_command = "seqtk sample " + fa_filename + " " + str(num_reads) + " > " + downsample_file
        print downsample_command

	# create downsample file if it doesn't already exist
    	if not (os.path.exists(downsample_filename) and os.path.getsize(downsample_filename) > 0):
        	subprocess.call(downsample_command, shell=True)

	return downsample_filename

def create_overlaps_file(fa_filename, output_prefix, coverage, data_type):
    # use minimap2 to calculate long read overlaps
    overlaps_filename = "./" + output_prefix + "/overlaps/" + output_prefix + "_" + str(coverage) + "X_overlaps.paf"

    # create overlaps file if it doesn't already exist
    if not (os.path.exists(overlaps_filename) and os.path.getsize(overlaps_filename) > 0):
        if data_type == "ont":
            	minimap2_command = "minimap2 -x ava-ont " + fa_filename + " " + fa_filename + " > " + overlaps_filename
		print minimap2_command
        	subprocess.call(minimap2_command, shell=True)
        elif data_type == "pb": 
            	minimap2_command = "minimap2 -x ava-pb " + fa_filename + " " + fa_filename + " > " + overlaps_filename
            	subprocess.call(minimap2_command, shell=True)
        else:
            	print "Error: Long read sequencing data not recognized. Either ONT or PB data only."
            	sys.exit()

    return overlaps_filename

if __name__ == "__main__":
    main()


