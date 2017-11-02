#!/usr/bin/env python
"""Generate data to be used for the preqc-lr report
"""
try:
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
	import time
except ImportError:
	print('Missing package(s)')

paf_given=False
store_csv=False

def main():
	start = time.clock()
	# --------------------------------------------------------
    	# PART 0: Parse the input
    	# --------------------------------------------------------
    	parser = argparse.ArgumentParser(prog='preqc-lr',description='Calculate Pre-QC Long Read report')
    	parser.add_argument('-i', '--input', action="store", required=True, 
			dest="fa_filename", help="Fasta, fastq, fasta.gz, or fastq.gz files containing reads.")
    	parser.add_argument('-t', '--type', action="store", required=True, 
			dest="data_type", choices=['pb', 'ont'], help="Either pacbio (pb) or oxford nanopore technology data (ont).")
    	parser.add_argument('-n', '--sample_name', action="store", required=True, 
			dest="sample_name", help="Sample name; you can use the name of species for example. This will be used as output prefix.")
    	parser.add_argument('-p', '--paf', action="store", required=False, 
			dest="paf", help="Minimap2 pairwise alignment file (PAF). This is produced using 'minimap2 -x ava-ont sample.fastq sample.fastq'.")
    	parser.add_argument('--csv', action="store_true", required=False,
			dest="store_csv", default=False, help="Use flag to save comma separated values (CSV) files for each plot.")
    	parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.9')
    	args = parser.parse_args()

    	# --------------------------------------------------------
    	# PART 1: Check input integrity
    	# --------------------------------------------------------
    	# check sequences input
    	if not os.path.exists(args.fa_filename) or not os.path.getsize(args.fa_filename) > 0 or not os.access(args.fa_filename, os.R_OK):
		print "Fasta/fastq file does not exist, is empty or is not readable."
		raise InputFileError('Fasta/fastq')
		
    	# check output directory, report if already exists, see if user would like to change output directory
   	output_prefix = args.sample_name
    	working_dir = './' + output_prefix
    	while os.path.exists(working_dir):
    		use_working_dir = ''
		while not use_working_dir == "n" and not use_working_dir == "y":
			print "Output directory \'" + working_dir + "\' already exists. Do you still want to use this directory? [y/n]"
			use_working_dir = raw_input()
			if not use_working_dir == "n" and not use_working_dir == "y":
				print "Only 'y' or 'n' as response."
        	if use_working_dir == "n":
			print "Please input new output_prefix:"
			output_prefix = raw_input()
		else:
	       		working_dir = './' +  output_prefix
			print "Output directory: " + working_dir
			break

    	# check paf if given
    	global paf_given
    	if args.paf:
        	if not os.path.exists(args.paf) or not os.path.getsize(args.paf) > 0 or not os.access(args.paf, os.R_OK):
                	raise InputFileError('PAF')
		else:	
			paf_given = True	

    	global store_csv
    	if args.store_csv:
		store_csv=True

    	# --------------------------------------------------------
    	# PART 2: Initiate json object, and record input info
    	# --------------------------------------------------------
    	data = {}
    	data['sample_name'] = args.sample_name
    	data['data_type'] = args.data_type
    	data['fa_filename'] = args.fa_filename
    	data['csv_storage'] = args.store_csv
    	if paf_given:
    		data['paf'] = args.paf

    	# --------------------------------------------------------
    	# PART 3: Create work dir, and csv dir if requested
    	# --------------------------------------------------------
    	csv_dir = './' + output_prefix + '/csv'
    	work_dir = working_dir
    	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
    	if args.store_csv and not os.path.exists(csv_dir):
		os.makedirs(csv_dir)

    	# --------------------------------------------------------
    	# PART 4: Create the preqclr report data
    	# --------------------------------------------------------
    	if paf_given:
    		calculate_report(output_prefix, args.fa_filename, args.data_type, data, args.paf) 
    	else:
		calculate_report(output_prefix, args.fa_filename, args.data_type, data)
	end = time.clock()
	total_time = end - start
	print "Total time: " + str(total_time) + " seconds"

class fasta_file:
    def __init__(self, fa_filename, read_seqs, num_reads, mean_read_length, data_type, paf = ''):
        self.fa_filename = fa_filename
        self.read_seqs = read_seqs
        self.num_reads = num_reads
        self.mean_read_length = mean_read_length
        self.data_type = data_type
	self.paf = paf

    def get_fa_filename(self):
        return self.fa_filename

    def get_read_seqs(self):
        return self.read_seqs

    def get_num_reads(self):
        return self.num_reads

    def get_mean_read_length(self):
        return self.mean_read_length

    def get_data_type(self):
        return self.data_type

    def get_paf(self):
        if os.path.exists(self.paf) and os.path.getsize(self.paf) > 0 and os.access(self.paf, os.R_OK):
		return self.paf
	else:
                raise InputFileError('PAF')

class Error(Exception):
	"""Base class for other exceptions"""
	pass

class InputFileError(Error):
	"""Raised if PAF does not exist, or exists but it is empty or not readable"""
	def __init__(self, filetype):
		self.filetype = filetype
		self.msg =  self.filetype + " does not exist, or exists but it is empty or not readable"

def reduce_represent(read_seq):
	#reduced_rep = read_seq.replace('A','00')
        #reduced_rep = reduced_rep.replace('C','01')
        #reduced_rep = reduced_rep.replace('G', '10')
        #reduced_rep = reduced_rep.replace('T', '11')
        #int_rep = int( reduced_rep , 2)
        count_C = read_seq.count('C')
        count_G = read_seq.count('G')
        count_A = read_seq.count('A')
        count_T = read_seq.count('T')
        read_length = int(count_A) + int(count_C) + int(count_G) + int(count_T)
        count_rep = "A" + str(count_A) + "C" + str(count_C) + "G" + str(count_G) + "T" + str(count_T)
	return (count_rep, read_length)

def calculate_report(output_prefix, fa_filename, data_type, data, paf=''):
    	# --------------------------------------------------------
    	# PART 0: Detect the file type, parse file, save to dict
    	# --------------------------------------------------------
	print '\n \n'
	print "Parsing sequences file"
	print "-------------------------------------"
	start = time.clock()
    	file_type = detect_filetype(fa_filename)
    	fa_sequences = dict(list())
    	if ".gz" in file_type:
        	with gzip.open(fa_filename, "rt") as handle:
			print "unzipped fine!"
            		if "fasta.gz" in file_type:
				for record in SeqIO.parse(handle, "fasta"):
					read_id, read_seq = record.id, record.seq
					fa_sequences[read_id] = reduce_represent(read_seq)
            		elif "fastq.gz" in file_type:
				for record in SeqIO.parse(handle, "fastq"):
                        		read_id, read_seq = record.id, str(record.seq)
                        		fa_sequences[read_id] = reduce_represent(read_seq)
    	else:
        	with open(fa_filename, "rt") as handle:
            		if "fasta" in file_type:
                		for record in SeqIO.parse(handle, "fasta"):
                        		read_id, read_seq = record.id, str(record.seq)
                        		fa_sequences[read_id] = reduce_represent(read_seq)
            		elif "fastq" in file_type:
                		for record in SeqIO.parse(handle, "fastq"):
                        		read_id, read_seq = record.id, record.seq
                        		fa_sequences[read_id] = reduce_represent(read_seq)
	end = time.clock()
	total_time = end - start
	print "Time elapsed: " + str(total_time) + " seconds"

    	# --------------------------------------------------------
    	# PART 1: Create all necessary info for fasta object
    	# --------------------------------------------------------
    	total_num_bases = 0
    	for read_id in fa_sequences:
        	seq_length = fa_sequences[read_id][1]
        	total_num_bases+=seq_length
    	mean_read_length = float(total_num_bases)/float(len(fa_sequences))
    	num_reads = len(fa_sequences)

    	global paf_given
    	if not paf_given:
		paf = create_overlaps_file(fa_filename, output_prefix, data_type)

    	fasta = fasta_file(fa_filename, fa_sequences, num_reads, mean_read_length, data_type, paf)	

    	# add the fasta file information to data
    	data['num_reads'] = num_reads
    	data['mean_read_length_no_filter'] = mean_read_length
    	data['total_num_bases'] = total_num_bases
    	data['paf'] = paf
  
    	# --------------------------------------------------------
    	# PART 2: Let the calculations begin...
   	# --------------------------------------------------------
	start = time.clock()
    	calculate_read_length(fasta, output_prefix, data)
	end = time.clock()
	print "Time elapsed: " + str(end - start) + " seconds"

	start = time.clock()
    	calculate_num_overlaps_per_read(fasta, output_prefix, data)
	end = time.clock()
	print "Time elapsed: " + str(end - start) + " seconds"

	start = time.clock()
    	calculate_estimated_coverage(fasta, output_prefix, data)
        end = time.clock()
        print "Time elapsed: " + str(end - start) + " seconds"

        start = time.clock()
    	calculate_GC_content_per_read(fasta, output_prefix, data)
        end = time.clock()
        print "Time elapsed: " + str(end - start) + " seconds"

        start = time.clock()
    	calculate_expected_minimum_fractional_overlaps_vs_read_length(fasta, output_prefix, data)
        end = time.clock()
        print "Time elapsed: " + str(end - start) + " seconds"

        start = time.clock()
	calculate_total_num_bases_vs_min_read_length(fasta, output_prefix, data)
        end = time.clock()
        print "Time elapsed: " + str(end - start) + " seconds"

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
    	reads = fasta.get_read_seqs()
	read_lengths = list()
        for read_id in reads:
		read_lengths.append(int(reads[read_id][1]))

        # --------------------------------------------------------
        # PART 1: Write to csv if requested
        # --------------------------------------------------------
	if store_csv:
		csv_filename = 'read_lengths.csv'
		write_to_csv( csv_filename, output_prefix, read_lengths)

    	# --------------------------------------------------------
    	# PART 1: Add to the data set
    	# --------------------------------------------------------
    	data['per_read_read_length'] = read_lengths

def calculate_num_overlaps_per_read(fasta, output_prefix, data):
    	print "\n\n\n\n"
    	print "Calculating number of overlaps per reads"
    	print "___________________________________________"

    	# --------------------------------------------------------
    	# PART 0: Get all the information needed from fasta
    	# --------------------------------------------------------
    	fa_filename = fasta.get_fa_filename()
    	data_type = fasta.get_data_type()
    	reads = fasta.get_read_seqs()
    	mean_read_length = fasta.get_mean_read_length()
    	num_reads = fasta.get_num_reads()

    	# get overlaps file created with all reads
    	overlaps_filename = fasta.get_paf()

    	# --------------------------------------------------------
    	# PART 1: Get the num of overlaps per read
    	# --------------------------------------------------------
    	# overlaps is the opened PAF file
    	# num_overlaps_per_read will, in the end, hold the total number of overlaps
    	num_overlaps_per_read = dict()
    	with open(overlaps_filename, 'r') as overlaps:
    		for line in overlaps:
        		query_read_id = line.split('\t')[0]
        		target_read_id = line.split('\t')[5]

        		# do not count overlaps with self
        		if target_read_id != query_read_id:
            			if query_read_id in num_overlaps_per_read:
                			num_overlaps_per_read[query_read_id] += 1
            			else:
                			num_overlaps_per_read[query_read_id] = 1
 
    	# --------------------------------------------------------
    	# PART 2: Divide the num overlaps by the read length
    	# --------------------------------------------------------
    	for read_id in num_overlaps_per_read:
        	read_length = float(reads[read_id][1])
        	num_overlaps = float(num_overlaps_per_read[read_id])
        	num_overlaps_per_read[read_id] = num_overlaps / read_length 

        # --------------------------------------------------------
        # PART 3: Write to csv if requested
        # --------------------------------------------------------
        if store_csv:
                csv_filename = 'num_overlaps_per_read.csv'
                write_to_csv( csv_filename, output_prefix, num_overlaps_per_read)
 
    	# --------------------------------------------------------
    	# PART 4: Add to the data set
    	# --------------------------------------------------------
    	data['per_read_overlap_count'] = num_overlaps_per_read

def calculate_estimated_coverage(fasta, output_prefix, data):
    	print "\n\n\n\n"
    	print "Calculating average coverage per read"
   	print "___________________________________________"

    	# --------------------------------------------------------
    	# PART 0: Get all the information needed from fasta
    	# --------------------------------------------------------
    	fa_filename = fasta.get_fa_filename()
    	data_type = fasta.get_data_type()
    	reads = fasta.get_read_seqs()
    	mean_read_length = fasta.get_mean_read_length()
    	total_num_reads = fasta.get_num_reads()
    	overlaps_filename = fasta.get_paf()

    	# --------------------------------------------------------
    	# PART 1: Get the sum of overlap lengths for each read
    	# --------------------------------------------------------
    	sum_overlap_lengths = dict()
    	overlap_accuracy = dict()			# key = read id, values = (sum_overlap_length, sum_matches )
    	with open(overlaps_filename, "r") as overlaps:
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
        	num_bases = float(reads[read_id][1])
        	total_overlaps_length = float(sum_overlap_lengths[read_id])
        	cov = int(round(float(total_overlaps_length) / float(num_bases)))
		accuracy = (overlap_accuracy[read_id][1])/float(overlap_accuracy[read_id][0])
		covs.append((cov, num_bases))
		accuracies.append((accuracy, num_bases))
    		if cov in per_cov_read_count:
        		per_cov_read_count[cov]+=1
		else:
        		per_cov_read_count[cov]=1
        	read_length_and_estimated_cov.append((num_bases, cov))

    	# --------------------------------------------------------
    	# PART 5: Filter outliers
    	# --------------------------------------------------------
    	print "\n\n"
    	print "Pre-Filter"
    	print "-------------------------------------"
    	max_cov = max(covs, key=itemgetter(0))[0]
    	min_cov = min(covs, key=itemgetter(0))[0]
    	print "Max cov: " + str(max_cov)
    	print "Min cov: " + str(min_cov)
    	mean_cov = float(np.mean([x[0] for x in covs]))
    	median_cov = float(np.median([x[0] for x in covs]))
    	mode_cov = stats.mode([x[0] for x in covs])
    	mean_read_length = float(np.mean([x[1] for x in covs]))
    	median_read_length = float(np.median([x[1] for x in covs]))
    	mode_read_length = stats.mode([x[1] for x in covs])
    	print "Cov (mean, median, mode): (" + str(mean_cov) + "," + str(median_cov) + "," + str(mode_cov) + ")"
    	print "Read length (mean, median, mode): (" + str(mean_read_length) + "," + str(median_read_length) + "," + str(mode_read_length) + ")"
    	num_reads = float(len(covs))
    	pre_filter = (num_reads, max_cov, min_cov, median_cov, mean_read_length) 

    	print "\n\n"
    	print "Filters based on"
    	print "-------------------------------------"
    	q75, q25 = np.percentile([x[0] for x in covs], [75 ,25])
    	IQR = float(q75) - float(q25)
    	upperbound = q75 + IQR * 1.5
    	lowerbound = q25 - IQR * 1.5
    	print "Q75: " + str(q75)
    	print "Q25: " + str(q25)
    	print "Upperbound limit: " + str(upperbound)
    	print "Lowerbound limit: " + str(lowerbound)

    	# filtering
    	temp = [i for i in covs if i[0] < upperbound]
    	filtered_covs = [i for i in temp if i[0] > lowerbound]

    	print "\n\n"
    	print "Post-Filter"
    	print "-------------------------------------"
    	max_cov = max(filtered_covs, key=itemgetter(0))[0]
    	min_cov = min(filtered_covs, key=itemgetter(0))[0]
    	print "Max cov: " + str(max_cov)
    	print "Min cov: " + str(min_cov)
    	mean_cov = float(np.mean([x[0] for x in filtered_covs]))
    	median_cov = float(np.median([x[0] for x in filtered_covs]))
    	mode_cov = stats.mode([x[0] for x in filtered_covs])
    	mean_read_length = float(np.mean([x[1] for x in filtered_covs]))
    	median_read_length = float(np.median([x[1] for x in filtered_covs]))
    	mode_read_length = stats.mode([x[1] for x in filtered_covs])
    	print "Cov (mean, median, mode): (" + str(mean_cov) + "," + str(median_cov) + "," + str(mode_cov) + ")"
    	print "Read length (mean, median, mode): (" + str(mean_read_length) + "," + str(median_read_length) + "," + str(mode_read_length) + ")"
    	num_reads = float(len(filtered_covs))
        q75, q25 = np.percentile([x[0] for x in filtered_covs], [75 ,25])
    	post_filter = (num_reads, max_cov, min_cov, median_cov, mean_read_length, upperbound, lowerbound)

    	# --------------------------------------------------------
    	# PART 6: Estimate genome size based off of filtered cov
    	# --------------------------------------------------------
    	n = float(len(filtered_covs))
    	l = float(mean_read_length)
    	c = float(np.median([x[0] for x in filtered_covs]))
    	estimated_genome_size = ( n * l ) / c
    
   	# --------------------------------------------------------
    	# PART 7: Estimate number of islands
    	# --------------------------------------------------------
    	T = 100.0			# amount of overlap in base pairs needed to detect overlap
    	theta = T / l  		# expected minimum fractional overlap required between two clones
    	print "T: " + str(T)
    	print "Mean read length:" + str(l) 
    	print "Theta: " + str(theta)
    	sigma = 1 - theta
        g = estimated_genome_size
    	estimated_num_islands = ( ( g * c ) / l ) * math.exp(-(1 - theta) * c)
    	estimated_num_islands_1 = n/math.exp(sigma * c) - n/math.exp(2*sigma*c)
    	estimated_num_islands_2 = n / math.exp(sigma * c)
    	print estimated_num_islands
    	print estimated_num_islands_1
    	print estimated_num_islands_2

        # --------------------------------------------------------
        # PART 8: Write to csv if requested
        # --------------------------------------------------------
        if store_csv:
                csv_filename = 'per_read_estimated_coverage.csv'
                write_to_csv( csv_filename, output_prefix, [x[0] for x in covs])
		
		csv_filename = 'read_lengths_estimated_cov.csv'
		write_to_csv( csv_filename, output_prefix, read_length_and_estimated_cov)

    	# --------------------------------------------------------
    	# PART 9: Add to the data set
    	# --------------------------------------------------------
    	data['estimated_coverage_stats_pre_filter'] = pre_filter
    	data['estimated_coverage_stats_post_filter'] = post_filter
    	data['per_read_estimated_coverage'] = ([x[0] for x in covs], upperbound, n, q25, q75)
    	data['read_lengths_estimated_cov'] = read_length_and_estimated_cov	
    	data['estimated_genome_size'] = estimated_genome_size
    	data['overlap_accuracies'] = accuracies
    	data['estimated_num_islands'] = estimated_num_islands

def extract_from_count_rep(count_rep, nt):
	next_nt = ""
	if nt == "A":
		next_nt = "C"
	elif nt == "C":
		next_nt = "G"
	elif nt == "G":
		next_nt = "T"
	else:
		return count_rep.split("T")[1]
	start = count_rep.index(nt) + 1
	end = count_rep.index(next_nt, start)
	return count_rep[start:end]

def calculate_GC_content_per_read(fasta, output_prefix, data):
    	print "\n\n\n\n"
    	print "Calculating GC-content per read"
    	print "___________________________________________"

    	# --------------------------------------------------------
    	# PART 0: Get all the information needed from fasta
    	# --------------------------------------------------------
    	reads = fasta.get_read_seqs()
        read_counts_per_GC_content = dict()
	for read_id in reads:
		read_length = reads[read_id][1]
		count_rep = reads[read_id][0]	# this is a string containing info about nucleotide counts see reduce_represent()
		count_C = extract_from_count_rep(count_rep, 'C')
		count_G = extract_from_count_rep(count_rep, 'G')
		count_GC = int(count_G) + int(count_C)
		GC_content = (float(count_GC)/float(read_length))*100.0
    		if round(GC_content) in read_counts_per_GC_content:
			read_counts_per_GC_content[round(GC_content)]+=1
    		else:
			read_counts_per_GC_content[round(GC_content)]=1

        # --------------------------------------------------------
        # PART 1: Write to csv if requested
        # --------------------------------------------------------
        if store_csv:
                csv_filename = 'read_counts_per_GC_content.csv'
                write_to_csv( csv_filename, output_prefix, read_counts_per_GC_content)

        # --------------------------------------------------------
        # PART 2: Add to the data set 
        # --------------------------------------------------------
    	data["read_counts_per_GC_content"] = read_counts_per_GC_content			

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
    	data_type = fasta.get_data_type()
    	reads = fasta.get_read_seqs()
    	mean_read_length = fasta.get_mean_read_length()
    	total_num_reads = fasta.get_num_reads()

    	# --------------------------------------------------------
    	# PART 1: get the PAF file calculated w all reads
    	# --------------------------------------------------------
    	# we get this by calculating the max_coverage
    	overlaps_filename = fasta.get_paf()
    	print overlaps_filename
    	overlaps = open(overlaps_filename, "r")
    	num = list()
    	for line in overlaps:
       		num_matching_residues = int(line.split('\t')[9])
		num.append(num_matching_residues)

    	data['num_matching_residues'] = num

def calculate_total_num_bases_vs_min_read_length(fasta, output_prefix, data):
    	print "\n\n\n\n"
    	print "Calculating total number of bases as a function of min read length"
    	print "___________________________________________"

    	# --------------------------------------------------------
    	# PART 0: Get all the information needed from fasta
    	# --------------------------------------------------------
    	fa_filename = fasta.get_fa_filename()
    	data_type = fasta.get_data_type()
    	reads = fasta.get_read_seqs()
    	mean_read_length = fasta.get_mean_read_length()
    	total_num_reads = fasta.get_num_reads()

	# get read lengths
        read_lengths = dict()
        for read_id in reads:
		#read_lengths.append(int(reads[read_id][1]))
		l = int(reads[read_id][1])
		if l in read_lengths:
			read_lengths[l] += 1
		else:
			read_lengths[l] = 1

	# get 100th, 90th, 80th, 70th, 60th, ..., 0th percentile
	min_read_lengths = list()

        unique_lengths = read_lengths.keys()
	percentile = 90
	while percentile > 0:
		cutoff = np.percentile(unique_lengths, percentile)
		min_read_lengths.append(int(cutoff))
		percentile-=5

	total_num_bases_vs_min_read_length = dict()
        total_num_bases = 0
	while not len(min_read_lengths) == 0:
		# starting from the maximum read length cut off
		cutoff = min_read_lengths.pop(0)
		# starting from the largest read length value then going down ...
		unique_lengths.sort()
		l = unique_lengths.pop(len(unique_lengths)-1)
		while l > cutoff:
			num_reads_w_length = read_lengths[l]
			total_num_bases+=l*num_reads_w_length
			l = unique_lengths.pop(len(unique_lengths)-1)
		total_num_bases_vs_min_read_length[cutoff] = total_num_bases

	print total_num_bases_vs_min_read_length	

        # --------------------------------------------------------
        # PART 1: Write to csv if requested
        # --------------------------------------------------------
        if store_csv:
                csv_filename = 'total_num_bases_vs_min_read_length.csv'
                write_to_csv( csv_filename, output_prefix, total_num_bases_vs_min_read_length)


	data['total_num_bases_vs_min_read_length'] = total_num_bases_vs_min_read_length

def write_to_csv(csv_filename, output_prefix, data):
	csv_filepath = './' + output_prefix + '/csv/' + csv_filename 
	if isinstance(data, list):
		with open(csv_filepath, "w") as outfile:
        		for entry in data:
            			outfile.write(str(entry))
            			outfile.write("\n")
	elif isinstance(data, dict):
		with open(csv_filepath, 'wb') as outfile:
    			writer = csv.writer(outfile)
    			for key, value in data.items():
       				writer.writerow([key, value])
 
def create_overlaps_file(fa_filename, output_prefix, data_type, coverage=''):
    # use minimap2 to calculate long read overlaps
    if not coverage == "":
    	overlaps_filename = "./" + output_prefix + "/" + output_prefix + "_" + str(coverage) + "X_overlaps.paf"
    else:
	overlaps_filename = "./" + output_prefix + "/" + output_prefix + "_overlaps.paf"

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
