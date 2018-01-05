#!/usr/bin/env python
# ========================================================
# preqc-lr calculate:
# Generates all the information used to create plots with
# preqc-lr report script.
# ========================================================
try:
	from Bio import SeqIO
	import subprocess, math, argparse
	import numpy as np
	import sys, os, csv
	import json, gzip
	from operator import itemgetter
	import time
except ImportError:
	print('Missing package(s)')
	quit()

gfa_given=False
paf_given=False
store_csv=False
verbose=False
log=list()
est_genome=0

def main():
	start = time.clock()	
	# --------------------------------------------------------
	# PART 0: Parse the input
	# --------------------------------------------------------
	parser = argparse.ArgumentParser(prog='preqc-lr',description='Calculate Pre-QC Long Read report')
	parser.add_argument('-r', '--reads', action="store", required=True, 
	dest="fa_filename", help="Fasta, fastq, fasta.gz, or fastq.gz files containing reads.")
	parser.add_argument('-t', '--type', action="store", required=True, 
	dest="data_type", choices=['pb', 'ont'], help="Either pacbio (pb) or oxford nanopore technology data (ont).")
	parser.add_argument('-n', '--sample_name', action="store", required=True, 
	dest="sample_name", help="Sample name; you can use the name of species for example. This will be used as output prefix.")
	parser.add_argument('-p', '--paf', action="store", required=False, 
	dest="paf", help="Minimap2 pairwise alignment file (PAF). This is produced using 'minimap2 -x ava-ont sample.fastq sample.fastq'.")
	parser.add_argument('--csv', action="store_true", required=False,
	dest="store_csv", default=False, help="Use flag to save comma separated values (CSV) files for each plot.")
	parser.add_argument('-g', '--gfa', action="store", required=False, dest="gfa", help="Graph Fragment Assembly (GFA) file created by miniasm.")
	parser.add_argument('--verbose', action="store_true", required=False, dest="verbose", help="Use to output preqc-lr progress.")
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2')
	args = parser.parse_args()

	# initialize list that will log all output
	global verbose
	global log
	verbose = args.verbose

	custom_print( "========================================================")
	custom_print( "RUNNING PREQC-LR CALCULATE" ) 
	custom_print( "========================================================")
	custom_print( "========================================================")
	custom_print( "CHECKING INPUT" )
	custom_print( "========================================================")

	# --------------------------------------------------------
	# PART 1: Check input integrity
	# --------------------------------------------------------
	# check sequences input
	if not os.path.exists(args.fa_filename) or not os.path.getsize(args.fa_filename) > 0 or not os.access(args.fa_filename, os.R_OK):
		print( "ERROR: Fasta/fastq file does not exist, is empty or is not readable." )
		sys.exit(1)

	# set output directory
	output_prefix = args.sample_name
	working_dir = './' + output_prefix

	# check gfa if given
	global gfa_given
	if args.gfa:
		if not os.path.exists(args.gfa) or not os.path.getsize(args.gfa) > 0 or not os.access(args.gfa, os.R_OK):
			print( "ERROR: GFA file does not exist, is empty or is not readable." )
			sys.exit(1)
		else:
			gfa_given = True

	# check paf if given
	global paf_given
	if args.paf:
		if not os.path.exists(args.paf) or not os.path.getsize(args.paf) > 0 or not os.access(args.paf, os.R_OK):
			print( "ERROR: PAF file does not exist, is empty or is not readable." )
			sys.exit(1)
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
	custom_print( "[ Input ]" )
	custom_print( "[+] Sample name: " + args.sample_name )
	custom_print( "[+] Data type: " + args.data_type )
	custom_print( "[+] FASTA/Q: " + args.fa_filename )
	custom_print( "[+] Storing csv: " + str(args.store_csv) )
	if paf_given:
		custom_print( "[+] PAF: " + str(args.paf) )
	else:
 		custom_print( "[+] PAF: None given, will run minimap2." )
	if gfa_given:
		custom_print( "[+] GFA: " + str(args.gfa) )
	else:
		custom_print( "[+] GFA: None given, will not run NGX calculations." )

	# --------------------------------------------------------
	# PART 3: Create work dir, and csv dir if requested
	# --------------------------------------------------------
	custom_print( "[ Output ]" )
	csv_dir = './' + output_prefix + '/csv'
	work_dir = working_dir
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	if args.store_csv and not os.path.exists(csv_dir):
		os.makedirs(csv_dir)
	custom_print( "[+] Output directory: " + work_dir )
	if args.store_csv:
		custom_print( "[+] Storing csv in: " + csv_dir )

	# --------------------------------------------------------
	# PART 4: Create the preqclr report data
	# --------------------------------------------------------
	if paf_given and gfa_given:
		calculate_report(output_prefix, args.fa_filename, args.data_type, data, args.gfa, args.paf) 
	elif paf_given and not gfa_given:
		calculate_report(output_prefix, args.fa_filename, args.data_type, data, '', args.paf)
	elif not paf_given and gfa_given:
		calculate_report(output_prefix, args.fa_filename, args.data_type, data, args.gfa, '')
	else:
		calculate_report(output_prefix, args.fa_filename, args.data_type, data)
	end = time.clock()
	total_time = end - start
	custom_print( "[ Done ]" )
	custom_print( "[+] Total time: " + str(total_time) + " seconds" )

	# --------------------------------------------------------
	# Final: Store log in file if user didn't specify verbose
	# --------------------------------------------------------
	outfile = output_prefix + "_preqclr.log"
	with open(outfile, 'wb') as f:
		for o in log:
			f.write(o + "\n")

def calculate_report(output_prefix, fa_filename, data_type, data, gfa='', paf=''):
	# ========================================================
	# Preprocess and calls calculation functions:
	# --------------------------------------------------------
	# Performs all the preprocessing of input for each calc,
	# then calls all the calculation functions.
	# Input:	FASTA/Q file, data type, paf if given, and 
	#		   initiated json structure called data.
	# Output:   json file holding all information needed to
	#		   create plots.
	# ========================================================

	custom_print( "========================================================" )
	custom_print( "PRE-PROCESS INPUT" )
	custom_print( "========================================================" )
	# --------------------------------------------------------
	# PART 0: Detect the reads file type, parse, save to dict
	# --------------------------------------------------------
	start = time.clock()
	file_type = detect_filetype(fa_filename)
	fa_sequences = parse_fa(fa_filename)
	end = time.clock()
	total_time = end - start
	custom_print( "[+] Time elapsed: " + str(total_time) + " seconds" )

	# --------------------------------------------------------
	# PART 1: Create all necessary info for fasta object
	# --------------------------------------------------------
	# count total number of bases, and mean read length
	total_num_bases = 0
	for read_id in fa_sequences:
		seq_length = fa_sequences[read_id][1]
		total_num_bases+=seq_length
	mean_read_length = float(total_num_bases)/float(len(fa_sequences))
	num_reads = len(fa_sequences)

	# --------------------------------------------------------
	# PART 2: Create PAF file if not given
	# --------------------------------------------------------
	global paf_given
	if not paf_given:
		start = time.clock()
		paf = create_overlaps_file(fa_filename, output_prefix, data_type)
		end = time.clock()
		total_time = end - start
		custom_print( "[+] Time elapsed: " + str(total_time) + " seconds" ) 

	# --------------------------------------------------------
	# PART 3: Parse PAF file, extract only necessary info
	# --------------------------------------------------------
	start = time.clock()
	paf_records = parse_paf(paf)
	end = time.clock()
	total_time = end - start
	custom_print( "[+] Time elapsed: " + str(total_time) + " seconds" )

	# store all info about fasta file and PAF file into object
	fasta = fasta_file(fa_filename, fa_sequences, num_reads, mean_read_length, data_type, paf_records)  

	# add the fasta file information to data
	data['num_reads'] = num_reads
	data['mean_read_length_no_filter'] = mean_read_length
	data['total_num_bases'] = total_num_bases
	data['paf'] = paf

	custom_print( "========================================================" )
	custom_print( "CALCULATE INFO" )
	custom_print( "========================================================" )
	# --------------------------------------------------------
	# PART 4: Let the calculations begin...
	# --------------------------------------------------------

	start = time.clock()
	calculate_read_length(fasta, output_prefix, data)
	end = time.clock()
	custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )

	start = time.clock()
	calculate_num_overlaps_per_read(fasta, output_prefix, data)
	end = time.clock()
	custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )

	start = time.clock()
	calculate_est_cov(fasta, output_prefix, data)
	end = time.clock()
	custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )

	global gfa_given
	if gfa_given:
		start = time.clock()
		calculate_ngx(gfa, output_prefix, data)
		end = time.clock()
		custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )

	start = time.clock()
	calculate_GC_content_per_read(fasta, output_prefix, data)
	end = time.clock()
	custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )
	
	start = time.clock()
	calculate_total_num_bases_vs_min_read_length(fasta, output_prefix, data)
	end = time.clock()
	custom_print( "[+] Time elapsed: " + str(end - start) + " seconds" )

	custom_print( "========================================================" )
	custom_print( "STORE INFO" )
	custom_print( "========================================================" )

	# --------------------------------------------------------
	# FINAL: After all calculations are done, store in json
	# --------------------------------------------------------
	preqclr_data = './' + output_prefix + '.preqclr'
	with open(preqclr_data, 'w') as outfile:  
		json.dump(data, outfile, indent=4)
	custom_print( "[+] Preqc-lr output stored in: " + preqclr_data )

def calculate_read_length(fasta, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating read length distribution ]" )
	# --------------------------------------------------------
	# Input:	Dictionary of reads
	# Output:   List of read lengths
	# ========================================================

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
	# PART 2: Add to the data set
	# --------------------------------------------------------
	read_lengths.sort()
	data['per_read_read_length'] = read_lengths

def calculate_num_overlaps_per_read(fasta, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating number of overlaps per read ]" )
	# --------------------------------------------------------
	# Input: FASTA/Q filename
	# Output: Either ['fastq.gz', 'fasta.gz', 'fastq', 'fasta']
	# ========================================================

	# --------------------------------------------------------
	# PART 0: Get all the information needed from fasta
	# --------------------------------------------------------
	fa_filename = fasta.get_fa_filename()
	data_type = fasta.get_data_type()
	reads = fasta.get_read_seqs()
	mean_read_length = fasta.get_mean_read_length()
	num_reads = fasta.get_num_reads()
	paf_records = fasta.get_paf_records()

	# --------------------------------------------------------
	# PART 1: Number of overlaps already calculated ...
	# --------------------------------------------------------
	# --------------------------------------------------------
	# PART 2: Divide the num overlaps by the read length
	# --------------------------------------------------------
	per_read_overlap_count = dict() 
	for read_id in paf_records:
			read_length = float(paf_records[read_id].get_read_length())
			num_overlaps = float(paf_records[read_id].get_total_num_overlaps())
			per_read_overlap_count[read_id] = num_overlaps / read_length 

	# --------------------------------------------------------
	# PART 3: Write to csv if requested
	# --------------------------------------------------------
	if store_csv:
			csv_filename = 'per_read_overlap_count.csv'
			write_to_csv( csv_filename, output_prefix, per_read_overlap_count)

	# --------------------------------------------------------
	# PART 4: Add to the data set
	# --------------------------------------------------------
	data['per_read_overlap_count'] = per_read_overlap_count

def calculate_est_cov(fasta, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating est cov per read, and est genome size ]" )
	# --------------------------------------------------------
	# Uses for each read, length and length of all overlaps,
	# to calculate the est cov. It then performs
	# filtering to calculate genome size estimates.
	# Input: parse_paf() output dictionary paf_records
	# Output: Dictionary with overlap info:
	#		 key = read_id
	#		 value = read(read_id, length, total OLs length, 
	#				 number of OLs)
	# ========================================================

	# --------------------------------------------------------
	# PART 0: Get all the information needed from fasta
	# --------------------------------------------------------
	fa_filename = fasta.get_fa_filename()
	data_type = fasta.get_data_type()
	reads = fasta.get_read_seqs()
	mean_read_length = fasta.get_mean_read_length()
	total_num_reads = fasta.get_num_reads()
	paf_records = fasta.get_paf_records()

	# --------------------------------------------------------
	# PART 1: Calculate cov for each read
	# --------------------------------------------------------
	covs = list()
	for read_id in paf_records:
		total_len_overlaps = float(paf_records[read_id].get_total_len_overlaps())   
		read_len = float(paf_records[read_id].get_read_length())
		read_cov =  round ( total_len_overlaps / read_len ) 
		covs.append((read_cov, read_len))

	# --------------------------------------------------------
	# PART 2: Filter outliers
	# --------------------------------------------------------
	custom_print( "[+] Pre-filter stats: " )
	max_cov = max(covs, key=itemgetter(0))[0]
	min_cov = min(covs, key=itemgetter(0))[0]
	custom_print( "[-] 	Max cov: " + str(max_cov) )
	custom_print( "[-] 	Min cov: " + str(min_cov) )
	mean_cov = float(np.mean([x[0] for x in covs]))
	median_cov = float(np.median([x[0] for x in covs]))
	mean_read_length = float(np.mean([x[1] for x in covs]))
	median_read_length = float(np.median([x[1] for x in covs]))
	custom_print( "[-] 	Cov (mean, median): (" + str(mean_cov) + "," + str(median_cov)  + ")" )
	custom_print( "[-] 	Read length (mean, median): (" + str(mean_read_length) + "," + str(median_read_length)  + ")" )
	num_reads = float(len(covs))
	pre_filter = (num_reads, max_cov, min_cov, median_cov, mean_read_length) 

	custom_print( "[+] Filter parameters: " )
	q75, q25 = np.percentile([x[0] for x in covs], [75 ,25])
	IQR = float(q75) - float(q25)
	upperbound = q75 + IQR * 1.5
	lowerbound = q25 - IQR * 1.5
	custom_print( "[-] 	Coverage upperbound limit: " + str(upperbound) )
	custom_print( "[-] 	Coverage lowerbound limit: " + str(lowerbound) )

	# filtering
	temp = [i for i in covs if i[0] < upperbound]
	filtered_covs = [i for i in temp if i[0] > lowerbound]

	custom_print( "[+] Post-filter stats: " )
	max_cov = max(filtered_covs, key=itemgetter(0))[0]
	min_cov = min(filtered_covs, key=itemgetter(0))[0]
	custom_print( "[-] 	Max cov: " + str(max_cov) )
	custom_print( "[-] 	Min cov: " + str(min_cov) )
	mean_cov = float(np.mean([x[0] for x in filtered_covs]))
	median_cov = float(np.median([x[0] for x in filtered_covs]))
	mean_read_length = float(np.mean([x[1] for x in filtered_covs]))
	median_read_length = float(np.median([x[1] for x in filtered_covs]))
	custom_print( "[-] 	Cov (mean, median): (" + str(mean_cov) + "," + str(median_cov) + ")" )
	custom_print( "[-]	Read length (mean, median): (" + str(mean_read_length) + "," + str(median_read_length) + ")" )
	num_reads = float(len(filtered_covs))
	q75, q25 = np.percentile([x[0] for x in filtered_covs], [75 ,25])
	IQR = int( q75 - q25 )
	post_filter = (num_reads, max_cov, min_cov, median_cov, mean_read_length, upperbound, lowerbound)

	# --------------------------------------------------------
	# PART 3: Estimate genome size based off of filtered cov
	# --------------------------------------------------------
	n = float(len(filtered_covs))
	l = float(mean_read_length)
	c = float(np.median([x[0] for x in filtered_covs]))
	est_genome_size = ( n * l ) / c

	# --------------------------------------------------------
	# PART 4: Write to csv if requested
	# --------------------------------------------------------
	if store_csv:
		csv_filename = 'per_read_est_cov.csv'
		write_to_csv( csv_filename, output_prefix, [x[0] for x in covs])

	# --------------------------------------------------------
	# PART 5: Add to the data set
	# --------------------------------------------------------
	data['est_cov_stats_pre_filter'] = pre_filter
	data['est_cov_stats_post_filter'] = post_filter
	data['est_cov_post_filter_info'] = (lowerbound, upperbound, num_reads, IQR)
	data['per_read_est_cov_and_read_length'] = covs   
	data['est_genome_size'] = est_genome_size

	global est_genome
	est_genome = est_genome_size	

def calculate_GC_content_per_read(fasta, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating GC-content per read ]" )
	# --------------------------------------------------------
	# Uses for reduced representation of each read to sum 
	# the number of Gs and Cs, then divides by total read
	# length.
	# Input:	Dict of reads
	# Output:   Dictionary with GC content info:
	#		   	key = GC content level (0-100%)
	#		   	value = # of reads with GC content level
	# ========================================================

	# --------------------------------------------------------
	# PART 0: Get all the information needed from fasta
	# --------------------------------------------------------
	reads = fasta.get_read_seqs()
	read_counts_per_GC_content = dict()
	for read_id in reads:
		read_length = reads[read_id][1]
		count_rep = reads[read_id][0]   # this is a string containing info about nucleotide counts see reduce_represent()
		count_C = translate_reduce_represent(count_rep, 'C')
		count_G = translate_reduce_represent(count_rep, 'G')
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

def calculate_total_num_bases_vs_min_read_length(fasta, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating total number of bases as a function of min read length ]" )
	# --------------------------------------------------------
	# Uses for reduced representation of each read to sum 
	# the number of Gs and Cs, then divides by total read
	# length.
	# Input:	Dict of reads
	# Output:   Dictionary with GC content info:
	#		   key = GC content level (0-100%)
	#		   value = # of reads with GC content level
	# ========================================================

	# --------------------------------------------------------
	# PART 0: Get all the information needed from fasta
	# --------------------------------------------------------
	fa_filename = fasta.get_fa_filename()
	data_type = fasta.get_data_type()
	reads = fasta.get_read_seqs()
	mean_read_length = fasta.get_mean_read_length()
	total_num_reads = fasta.get_num_reads()

	# --------------------------------------------------------
	# PART 1: Bin the reads by read length
	# --------------------------------------------------------
	# read_lengths is a dictionary where:
	#	   key = unique read length
	#	   value = number of reads that have this read length
	read_lengths = dict()
	for read_id in reads:
		l = int(reads[read_id][1])
		if l in read_lengths:
			read_lengths[l] += 1
		else:
			read_lengths[l] = 1

	# --------------------------------------------------------
	# PART 2: Sort read lengths
	# --------------------------------------------------------
	min_read_lengths = list()
	unique_lengths = read_lengths.keys()
	unique_lengths.sort()											   # sort read lengths
	unique_lengths_desc = sorted(unique_lengths, key=int, reverse=True) # sort read lengths in descending order (from largest to smallest)
	total_num_bases_vs_min_read_length = list()
	longest_read_length = unique_lengths_desc.pop(0)
	num_reads_w_length = read_lengths[longest_read_length]
	total_num_bases = longest_read_length*num_reads_w_length
	read_length = longest_read_length

	# --------------------------------------------------------
	# PART 3: Traverse through unique read lengths from 
	# largest value, sum the length of reads that are greater
	# than the current read length. 
	# --------------------------------------------------------
	while not len(unique_lengths_desc) == 0:
		# starting from the largest read length value then going down ...
		total_num_bases_vs_min_read_length.append((read_length, total_num_bases))
		read_length = unique_lengths_desc.pop(0)
		num_reads_w_length = read_lengths[read_length]
		total_num_bases+=read_length*num_reads_w_length

	# --------------------------------------------------------
	# PART 4: Write to csv if requested
	# --------------------------------------------------------
	if store_csv:
			csv_filename = 'total_num_bases_vs_min_read_length.csv'
			write_to_csv( csv_filename, output_prefix, total_num_bases_vs_min_read_length)

	# --------------------------------------------------------
	# PART 5: Add to the data set
	# --------------------------------------------------------
	data['total_num_bases_vs_min_read_length'] = total_num_bases_vs_min_read_length

def calculate_ngx(gfa, output_prefix, data):
	# ========================================================
	custom_print( "[ Calculating NG(X) ]" )
	# --------------------------------------------------------
	# Uses GFA information to evaluate the assembly quality
	# Input:	Graphical fragment assembly
	# Output:	NG(X) values
	# ========================================================
	# --------------------------------------------------------
    # PART 0: Use the GFA file to extract info on contig lens
    # --------------------------------------------------------
	contig_lengths = list()
	sum_lengths = 0
	with open(gfa, "r") as file:
		for line in file:
			info = line.split("\t")
			type = info.pop(0).rstrip()
			if type == 'S':
				ln_flag = info.pop(2)
				length = int(ln_flag.split(':')[2].rstrip())
				contig_lengths.append(length)	

	# --------------------------------------------------------
	# PART 1: Sort the list of contig lengths in decr. order
	# --------------------------------------------------------
	contig_lengths_sorted = sorted(contig_lengths, reverse = True)

	# --------------------------------------------------------
	# PART 2: Calculate the xth percent of genome size est.
	# --------------------------------------------------------
	NGX_p_val = dict()
	x = 0
	while x <= 100:
		NGX_p_val[str(x)] = 0
		x+=1

	global est_genome
	est_genome_size = est_genome
	sum_contig_len = sum(contig_lengths)
	for x in NGX_p_val:
		NGX_p_val[str(x)] = float(est_genome_size)*(float(x)/100.0)	

	# --------------------------------------------------------
	# PART 3: Find NGX values
	# --------------------------------------------------------
	NGX = dict()
	curr_sum = 0
	start = 0
	end = 0

	for l in contig_lengths_sorted:
		# add unused longest contig length to the curr sum
		end+=l
		# for all percentile values that are less than the curr sum
		for x in NGX_p_val:
			p = NGX_p_val[x]
			if ( float(p) >= float(start) ) and ( float(p) <= float(end) ):
				NGX[str(x)] = l
		start+=l

	# if the genome size estimate is >>>>>> sum of contig lengths just assign all percentile values as last contig
	if len(NGX) == 0:
		for x in NGX_p_val:
			p = NGX_p_val[x]
			NGX[str(x)] = l

	# --------------------------------------------------------
	# PART 4: Write to csv if requested
	# --------------------------------------------------------	
	if store_csv:
		csv_filename = 'ngx_values.csv'
		write_to_csv( csv_filename, output_prefix, NGX)

	# --------------------------------------------------------
	# PART 5: Add to the dataset
	# --------------------------------------------------------
	data["NGX_values"] = NGX

class fasta_file:
	def __init__(self, fa_filename, read_seqs, num_reads, mean_read_length, data_type, paf_records):
		self.fa_filename = fa_filename
		self.read_seqs = read_seqs
		self.num_reads = num_reads
		self.mean_read_length = mean_read_length
		self.data_type = data_type
		self.paf_records = paf_records

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

	def get_paf_records(self):
		return self.paf_records

class read:
	def __init__(self, read_id, read_length, total_len_overlaps, total_num_overlaps):
		self.read_id = read_id
		self.read_length = read_length 
		self.total_len_overlaps = total_len_overlaps
		self.total_num_overlaps = total_num_overlaps

	def update_total_len_overlaps(self, ol_length):
		self.total_len_overlaps += int(ol_length)

	def update_total_num_overlaps(self):
		self.total_num_overlaps += 1

	def get_total_len_overlaps(self):
		return int(self.total_len_overlaps)

	def get_total_num_overlaps(self):
		return int(self.total_num_overlaps)

	def get_read_length(self):
		return int(self.read_length)

def reduce_represent(read_seq):
	# ========================================================
	# Represents reads in a reduced form: 
	# --------------------------------------------------------
	# We only need information on counts of nucleotides CG
	# Input:	FASTA/Q file
	# Output:   Dictionary of sequences:
	#		   	key = read_id
	#		   	value = (reduced rep. of seq, read length)
	# ========================================================
	count_C = read_seq.count('C')
	count_G = read_seq.count('G')
	read_length = int(len(read_seq))
	#count_rep = "A" + str(count_A) + "C" + str(count_C) + "G" + str(count_G) + "T" + str(count_T)
	count_rep = "C" + str(count_C) + "G" + str(count_G)
	return (count_rep, read_length)

def translate_reduce_represent(count_rep, nt):
	# ========================================================
	# Translates the read in the reduced form: 
	# --------------------------------------------------------
	# Retrieves the count for G or C nucleotides (nt) in
	# a given read.
	# Input: Reduced representation of read, and one of the nt
	#		{C, G} you want counts for.
	# Output: Number of occurences of indicated nt in read.
	# ========================================================
	next_nt = ""
	if nt == "C":
		next_nt = "G"
	elif nt == "G":
		return count_rep.split("G")[1]
	else:
		custom_print( "ERROR: translating reduce representation of read. Nucleotide not recognized." )
	start = count_rep.index(nt) + 1
	end = count_rep.index(next_nt, start)
	return int(count_rep[start:end])

def detect_filetype(fa_filename):
	# ========================================================
	# Detects filetype of sequences input
	# --------------------------------------------------------
	# Input: FASTA/Q filename
	# Output: Either ['fastq.gz', 'fasta.gz', 'fastq', 'fasta']
	# ========================================================
	path = fa_filename
	for ext in ['fastq.gz', 'fasta.gz', 'fastq', 'fasta']:
		if path.endswith(ext):
			return ext
	custom_print( "Must be either fasta, fastq, fasta.gz, fastq.gz" )
	sys.exit(1)

def parse_fa(fa_filename):
	# ========================================================
	custom_print( "[ Parse FASTA/Q file ]" )
	# --------------------------------------------------------
	# Detects file type then reduces seq info to only needed 
	# information
	# Input:	FASTA/Q file
	# Output:   Dictionary of sequences:
	#		   	key = read_id
	#		   	value = (reduced rep. of seq, read length)
	# ========================================================
	
	file_type = detect_filetype(fa_filename)
	fa_sequences = dict(list())
	if ".gz" in file_type:
		with gzip.open(fa_filename, "rt") as handle:
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

	return fa_sequences

def parse_paf(overlaps_filename):
	# ========================================================
	custom_print( "[ Parse PAF file ]" )
	# --------------------------------------------------------
	# Gets for each read, length, length of all overlaps (OLs), 
	# and number of OLs.
	# Input:	minimap2 output PAF file
	# Output:   Dictionary with overlap info:
	#		   	key = read_id
	#		   	value = read(read_id, length, total OLs length, 
	#				 number of OLs)
	# ========================================================
	paf_records = dict()
	with open(overlaps_filename, "r") as overlaps:
		for line in overlaps:
			# get information for each overlap
			ol = list()
			ol = line.split('\t')
			query_read_id = ol.pop(0)
			query_length = int(ol.pop(0))
			query_start_pos = int(ol.pop(0))
			query_end_pos = int(ol.pop(0))
			strand = ol.pop(0)
			target_read_id = ol.pop(0)
			target_length = int(ol.pop(0))
			target_start_pos = int(ol.pop(0))
			target_end_pos = int(ol.pop(0))
			num_matches = int(ol.pop(0))
			alignment_length = int(ol.pop(0))
			overlap_length = 0
			read_length_diff = abs(target_length - query_length)
			min_len = min(target_length, query_length)
			aln_len_min_len = float(alignment_length)/float(min_len)

			# calculate overlap length and account for soft clipping
			if not (query_read_id == target_read_id):
				query_prefix_len = query_start_pos
				query_suffix_len = query_length - query_end_pos

				target_prefix_len = target_start_pos
				target_suffix_len = target_length - target_end_pos

				# calculate length of overlap
				overlap_length = query_end_pos - query_start_pos

				# minimap2 might be performing softclipping
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

				# add sum overlaps, increment overlap counter
				if query_read_id in paf_records:
					paf_records[query_read_id].update_total_len_overlaps(int(overlap_length))
					paf_records[query_read_id].update_total_num_overlaps()
				else:
					new_read = read(query_read_id, int(query_length), int(overlap_length), int(1))
					paf_records[query_read_id] = new_read

				if target_read_id in paf_records:
					paf_records[target_read_id].update_total_len_overlaps(int(overlap_length))
					paf_records[target_read_id].update_total_num_overlaps()
				else:
					new_read = read(target_read_id, int(target_length), int(overlap_length), int(1))
					paf_records[target_read_id] = new_read
				percent_id = float(num_matches)/float(alignment_length)
	return paf_records

def write_to_csv(csv_filename, output_prefix, data):
	# ========================================================
	# Writes data sets to a csv file:
	# --------------------------------------------------------
	# Detects data structure type, then writes to csv.
	# Input:	data structure, csv filename
	# Output:   csv file with data in csv directory
	# ========================================================
	csv_filepath = './' + output_prefix + '/csv/' + csv_filename 
	if isinstance(data, list):
		with open(csv_filepath, "w") as outfile:
			csv_out=csv.writer(outfile)
			for entry in data:
				if isinstance(entry, tuple):
					csv_out.writerow(entry)
				else:
					outfile.write(str(entry))
					outfile.write("\n")
	elif isinstance(data, dict):
		with open(csv_filepath, 'wb') as outfile:
			writer = csv.writer(outfile)
			for key, value in data.items():
				writer.writerow([key, value])
 
def create_overlaps_file(fa_filename, output_prefix, data_type):
	# ========================================================
	custom_print( "[ Running minimap2 to achieve all-vs-all overlaps ]" )
	# --------------------------------------------------------
	# Input:	FASTA/Q filename and data type {ont, pb}
	# Output:   PAF file which holds all info about overlaps
	# ========================================================
	overlaps_filename = "./" + output_prefix + "/" + output_prefix + "_overlaps.paf"

	# check if minimap2 is installed and in PATH
	program="minimap2"
	minimap2_works=False
	fpath, fname = os.path.split(program)
	if fpath and os.path.isfile(fpath) and os.access(fpath, os.X_OK):
		minimap2_works=True
	else:
		# check each path in PATH environment variable to see if minimap2 is installed in any of these directories
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if ( os.path.isfile(exe_file) and os.access(exe_file, os.X_OK) ):
				minimap2_works=True
	if minimap2_works:
		custom_print("[+] Minimap2 was succesfully found in PATH and is executable.")
	else:
		print("ERROR: Minimap2 was unsuccessfully found in PATH or is NOT executable.")
		sys.exit(1)

	# create overlaps file if it doesn't already exist
	if not (os.path.exists(overlaps_filename) and os.path.getsize(overlaps_filename) > 0):
		if data_type == "ont":
			minimap2_command = "minimap2 -x ava-ont " + fa_filename + " " + fa_filename + " > " + overlaps_filename
			custom_print( "[x] " + minimap2_command )
			try:
				subprocess.call(minimap2_command, shell=True)
			except subprocess.CalledProcessError:
				print( "ERROR: minimap2 run was unsuccessful." )
				sys.exit(1)
		else: 
			minimap2_command = "minimap2 -x ava-pb " + fa_filename + " " + fa_filename + " > " + overlaps_filename
			custom_print( "[x] " + minimap2_command )
			try:
				subprocess.call(minimap2_command, shell=True)
			except subprocess.CalledProcessError:
				print( "ERROR: minimap2 run was unsuccessful." )
				sys.exit(1)
	else:
		custom_print( "[+] Non-empty PAF file already exists: " + overlaps_filename )

	return overlaps_filename

class Error(Exception):
	"""Base class for other exceptions"""
	pass

class InputFileError(Error):
	"""Raised if PAF does not exist, or exists but it is empty or not readable"""
	def __init__(self, filetype):
		self.filetype = filetype

def custom_print(s):
	global verbose
	global log
	if verbose:
		print s
	log.append(s)

if __name__ == "__main__":
	main()
