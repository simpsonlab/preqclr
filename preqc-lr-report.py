#!/usr/bin/env python
"""Generate a readable report from preqc output.
"""
from Bio import SeqIO
import subprocess, math, argparse
import matplotlib as MPL
MPL.use('Agg')
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pylab as plt
import sys, os, csv

def main():
        # process input
        parser = argparse.ArgumentParser(description='Display Pre-QC Long Read report')
        parser.add_argument('-i', '--input', action="store", required=True, dest="fa_filename", help="Fasta/q files containing reads")
        parser.add_argument('-o', '--output', action="store", dest="prefix", default="preqc-lr-output", help="Prefix for output pdf")
        parser.add_argument('-g', '--genome_size', action="store", required=True, dest="genome_size", help="Genome size. Default = 4641652 bps.")
        parser.add_argument('-t', '--type', action="store", required=True, dest="data_type", choices=['pb', 'ont'], help="Either pacbio or ONT.")
        parser.add_argument('-c', '--cov', action="store", required=True, dest="cov_reads_to_ref", help="Ref-SimReads coverage results in csv format with 2 columns: coverage, count")
        args = parser.parse_args()
        
        work_dir = './' + args.prefix
        print work_dir
        overlaps_dir = './' + args.prefix + '/overlaps'
	downsample_dir = './' + args.prefix + '/downsample'
	plot_dir = './' + args.prefix + '/plot_data'
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        if not os.path.exists(overlaps_dir):
            os.makedirs(overlaps_dir)
	if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        if not os.path.exists(downsample_dir):
            os.makedirs(downsample_dir)
        # parse the fasta file
        create_report(args.prefix, args.fa_filename, args.genome_size, args.data_type, args.cov_reads_to_ref)

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


def create_report(output_prefix, fa_filename, genome_size, data_type, cov_reads_to_ref):

        # configure the PDF file 
        MPL.rc('figure', figsize=(8,10.5)) # in inches
        MPL.rc('font', size=11)
        MPL.rc('xtick', labelsize=6)
        MPL.rc('ytick', labelsize=6)
        MPL.rc('legend', fontsize=6)
        MPL.rc('axes', titlesize=10)
        MPL.rc('axes', labelsize=8)
        MPL.rcParams['lines.linewidth'] = 1.5

        # TO DO : we need to adjust how many plots should be on the page
        '''
            insert code here
        '''

        output_pdf= output_prefix + ".pdf"
        pp = PdfPages(output_pdf)

        # add figure
        fig = plt.figure()
        fig.suptitle("Preqc Long Read Results : " + output_prefix )
        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

        # six subplots per pdf page
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)

        # get reads and some read info from fasta file
        fasta_sequences = SeqIO.parse(open(fa_filename),'fastq')
        read_seqs = dict()
        read_lengths = list()
        total_num_bases = 0
        for read in fasta_sequences:
            read_id, read_sequence, read_length = read.id, str(read.seq), len(str(read.seq))
            read_seqs[read_id] = read_sequence
            read_lengths.append(read_length)
            total_num_bases+=read_length
        mean_read_length = float(total_num_bases)/float(len(read_seqs))
        num_reads = len(read_seqs)

        fasta = fasta_file(fa_filename, read_seqs, read_lengths, genome_size, num_reads, mean_read_length, data_type)

        print "Extracting data from: " + fa_filename
        print "Total number of reads: " + str(num_reads)

        plot_read_length_distribution(ax1, fasta, output_prefix)
        plot_average_overlaps_vs_coverage_distribution(ax2, fasta, output_prefix)
        plot_num_overlaps_per_read_distribution(ax3, fasta, output_prefix)
        plot_estimated_coverage(ax4, fasta, output_prefix, cov_reads_to_ref)
        fig.savefig(pp, format='pdf')
        pp.close()

def plot_read_length_distribution(ax, fasta, output_prefix):
        print "\n\n\n\n"
        print "Plotting read length distribution"
        print "___________________________________________"

        read_lengths = fasta.get_read_lengths()
	
	read_length_data = "./" + output_prefix + "/plot_data/" + output_prefix + "_read_length_data.csv"
	with open(read_length_data, 'wb') as myfile:
    		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    		wr.writerow(read_lengths)

        ax.hist(read_lengths)
	print read_lengths
        ax.set_title('Read length distribution')
        ax.set_xlabel('Read lengths (bps)')
        ax.set_ylabel('Frequency')
	#ax.set_xticks(np.arange(0, max(read_lengths)+1000, 1000))
        ax.grid(True, linestyle='-', linewidth=0.3)
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	ax.get_xaxis().get_major_formatter().set_useOffset(False)

def plot_average_overlaps_vs_coverage_distribution(ax, fasta, output_prefix):
        print "\n\n\n\n"
        print "Plotting average overlaps vs coverage"
        print "___________________________________________"

        genome_size = fasta.get_genome_size()
        fa_filename = fasta.get_fa_filename()
        data_type = fasta.get_data_type()
        mean_read_length = fasta.get_mean_read_length()
        num_reads = fasta.get_num_reads()

        max_coverage =  round(math.ceil((int(num_reads) * float(mean_read_length)) / float(genome_size)))
        print "Max coverage : " + str(max_coverage)
        
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
            data = dict() 

            # downsampling, and getting overlaps
            for c in coverage_levels:
                    # calculate the number of reads needed to achieve coverage
                    num_reads = round(math.ceil(float((c * int(genome_size)))/float(mean_read_length)))

                    # use seqtk to downsample fastq files
                    downsample_file = "./" + output_prefix + "/downsample/" + output_prefix + "_" + str(c) + "X.fasta"
                    downsample_command = "seqtk sample " + fa_filename + " " + str(num_reads) + " > " + downsample_file
             	    print downsample_command
                    subprocess.call(downsample_command, shell=True)

                    # use minimap2 to calculate long read overlaps
                    overlaps_filename = create_overlaps_file(downsample_file, output_prefix, c, data_type)

                    # get the number of overlaps for each read
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

                    # get average # of overlaps per coverage
                    total_overlaps = 0
                    for query_read_id in reads_overlap_count:
                        total_overlaps += reads_overlap_count[query_read_id]
                    average_num_overlaps = total_overlaps/num_reads

                    data[c] = average_num_overlaps

            # plotting the avg number of overlaps vs coverage
            # sorted by ke, return a list of tuples
            lists = sorted(data.items())
            # unpack a list of pairs into two tuples
            x, y = zip(*lists)

	    # write plot data to csv file
	    overlaps_vs_coverage_data = "./" + output_prefix + "/plot_data/" + output_prefix + "_overlaps_vs_coverage_data.csv"
            with open(overlaps_vs_coverage_data, 'wb') as myfile:
		writer = csv.writer(myfile)
		for key, value in data.items():
			writer.writerow([key, value])
            myfile.close()

            ax.scatter(x, y)
            ax.set_title('Average number of overlaps vs Coverage')
            ax.set_xlabel('Coverage')
            ax.set_ylabel('Average number of overlaps')
	    ax.grid(True, linestyle='-', linewidth=0.3)

def plot_num_overlaps_per_read_distribution(ax, fasta, output_prefix):
        print "\n\n\n\n"
        print "Plotting number of overlaps per reads"
        print "___________________________________________"

        fa_filename = fasta.get_fa_filename()
        genome_size = fasta.get_genome_size()
        data_type = fasta.get_data_type()
        reads = fasta.get_read_seqs()
        mean_read_length = fasta.get_mean_read_length()
        num_reads = fasta.get_num_reads()

        reads_overlap_count = dict()

        # get overlaps file created with all reads
        max_coverage = round(math.ceil((int(num_reads) * float(mean_read_length)) / float(genome_size)))
        overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)

        # get the number of overlaps for each read
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
        
        # divide number of overlaps by the length of each sequence
        for read_id in reads_overlap_count:
	    read_length = float(len(reads[read_id]))
	    num_overlaps = float(reads_overlap_count[read_id])
	    # taking into account the read lengths
            reads_overlap_count[read_id] = num_overlaps / read_length 

	# write plot data to csv
	read_overlaps_count_data = "./" + output_prefix + "/plot_data/" + output_prefix + "_read_overlaps_count_data.csv"
        with open(read_overlaps_count_data, 'wb') as myfile:
        	writer = csv.writer(myfile)
        	for key, value in reads_overlap_count.items():
                	writer.writerow([key, value])
        myfile.close()

        # plotting the number of overlaps/read
        ax.hist(reads_overlap_count.values(), bins=120)
        ax.set_title('Number of overlaps/read/length distribution')
        ax.set_xlabel('Number of overlaps/read/length')
        ax.set_ylabel('Frequency')    
        ax.grid(True, linestyle='-', linewidth=0.3)

def plot_estimated_coverage(ax, fasta, output_prefix, cov_reads_to_ref):
	print "\n\n\n\n"
	print "Plotting average coverage per base per read"
	print "___________________________________________"


        # get all information on fasta file needed
        mean_read_length = fasta.get_mean_read_length()
        genome_size = fasta.get_genome_size()
        num_reads = fasta.get_num_reads()
        fa_filename = fasta.get_fa_filename()
        data_type = fasta.get_data_type()
	reads = fasta.get_read_seqs()

	# getting the overlaps file created with all reads
        max_coverage = round(math.ceil((int(num_reads) * float(mean_read_length)) / float(genome_size)))
        overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)

        overlaps = open(overlaps_filename, 'r')
        reads_overlaps = dict()
        for line in overlaps:
            query_read_id = line.split('\t')[0]
	    target_read_id = line.split('\t')[5]
            query_start_pos = line.split('\t')[2]
            query_end_pos = line.split('\t')[3]
	    target_start_pos = line.split('\t')[7]
	    target_end_pos = line.split('\t')[8]
	    target_overlap_length = int(target_end_pos) - int(target_start_pos)
            query_overlap_length = int(query_end_pos) - int(query_start_pos)
	    if not (query_read_id == target_read_id):
            	if query_read_id in reads_overlaps: 
            		reads_overlaps[query_read_id]+=query_overlap_length
            	else:
                	reads_overlaps[query_read_id]=query_overlap_length
		if target_read_id in reads_overlaps:
                        reads_overlaps[target_read_id]+=target_overlap_length
		else:
                        reads_overlaps[target_read_id]=target_overlap_length
        overlaps.close()

        # initialize list of coverage_per_base_per_read values
        # get the reads list from the fasta file object so we can get length of reads
        reads = fasta.get_read_seqs()
        data = dict() # key = avg coverage, value = number of reads with this avg coverage
	data_seqs = dict() # key = avg coverage, value = list of read ids that have this avg coverage
        for read_id in reads_overlaps:
            seq = reads[read_id]
            num_bases = len(seq)
	    cov_per_base_per_read = int(round(float(reads_overlaps[read_id]) / float(num_bases)))
	    if cov_per_base_per_read in data:
		data[cov_per_base_per_read]+=1
		data_seqs[cov_per_base_per_read].append(read_id)
	    else:
		data[cov_per_base_per_read]=1
		read_ids = list()
		read_ids.append(read_id)
		data_seqs[cov_per_base_per_read]=read_ids

	# get outlier sequences to check if they are repetitive regions
	# outlier definition = any data point outside the interval [Q1 - 1.5*IQR, Q3 + 1.5*IQR]
	q1 = np.percentile(data.values(), 25)
	q3 = np.percentile(data.values(), 75)
	iqr = q3 - q1
        upper_bound = q3 +  1*iqr
        lower_bound = q1 - 1*iqr
	print "Q1: " + str(q1) + ", Q3: " + str(q3) + ", IQR: " + str(iqr) + ", Upper Bound : " + str(upper_bound) + ", Lower Bound : " + str(lower_bound) 


	for key in data:
		print key
		if (key < lower_bound) or (key > upper_bound):
			print str(key) + " true"
			for x in data_seqs[key]:	
				print reads[x]
							

	# write plot data to csv file
        estimated_coverage_data = "./" + output_prefix + "/plot_data/" + output_prefix + "_estimated_coverage_data.csv"
        with open(estimated_coverage_data, 'wb') as myfile:
                writer = csv.writer(myfile)
                for key, value in data.items():
                        writer.writerow([key, value])
	myfile.close()


        # plotting the cov per base per read
	
	# superimposing hamza's data
	reader = csv.reader(open(cov_reads_to_ref, 'r'))
	counted_data = {}
	next(reader) # skip first line
	for row in reader:
   		k, v = row
  		counted_data[int(k)] = int(v)
	l="SimReads-SimReads (" + output_prefix + ")"
	bins = np.linspace(0, 130, 130)
	print data.keys()
	ax.hist(counted_data.keys(), weights=counted_data.values(), alpha=0.5, label="Ref-SimReads", bins=bins)
        ax.hist(data.keys(), weights=data.values(), alpha=0.5, label=l, bins=bins)
        ax.set_title('Estimated coverage distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xticks(np.arange(0, 130, 10.0))
        ax.set_xlabel('Avg. number of overlaps per base per read')
        ax.set_ylabel('Frequency')   
	ax.legend(loc='upper right')

def create_overlaps_file(fa_filename, output_prefix, coverage, data_type):
        # use minimap2 to calculate long read overlaps
        overlaps_filename = "./" + output_prefix + "/overlaps/" + output_prefix + "_" + str(coverage) + "X_overlaps.paf"

        # create overlaps file if it doesn't already exist
        if not (os.path.exists(overlaps_filename) and os.path.getsize(overlaps_filename) > 0):
            if data_type == "ont":
                minimap2_command = "minimap2 -x ava-ont " + fa_filename + " " + fa_filename + " > " + overlaps_filename
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

