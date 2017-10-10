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

        args = parser.parse_args()
        
        work_dir = './' + args.prefix
        print work_dir
        overlaps_dir = './' + args.prefix + '/overlaps'
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
            print "HELLO"
        if not os.path.exists(overlaps_dir):
            os.makedirs(overlaps_dir)

        # parse the fasta file
        create_report(args.prefix, args.fa_filename, args.genome_size, args.data_type)

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


def create_report(output_prefix, fa_filename, genome_size, data_type):

        # configure the PDF file 
        MPL.rc('figure', figsize=(8,10.5)) # in inches
        MPL.rc('font', size=11)
        MPL.rc('xtick', labelsize=6)
        MPL.rc('ytick', labelsize=6)
        MPL.rc('legend', fontsize=10)
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
        fig.suptitle("Preqc Long Read Results :")
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
        plot_average_coverage_per_base_per_read(ax4, fasta, output_prefix)
        fig.savefig(pp, format='pdf')
        pp.close()

def plot_read_length_distribution(ax, fasta, output_prefix):
        print "\n\n\n\n"
        print "Plotting read length distribution"
        print "___________________________________________"

        read_lengths = fasta.get_read_lengths()
	
	read_length_data = output_prefix + "_read_length_data.csv"
	with open(read_length_data, 'wb') as myfile:
    		wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    		wr.writerow(read_lengths)

        ax.hist(read_lengths)
        ax.set_title('Read length distribution')
        ax.set_xlabel('Read lengths (bps)')
        ax.set_ylabel('Frequency')


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
                    downsample_file = "./" + output_prefix + "/" + output_prefix + "_" + str(c) + "X.fasta"
                    downsample_command = "seqtk sample " + fa_filename + " " + str(num_reads) + " > " + downsample_file
          #          print downsample_command
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

	    overlaps_vs_coverage_data = output_prefix + "_overlaps_vs_coverage_data.csv"
            with open(overlaps_vs_coverage_data, 'wb') as myfile:
		writer = csv.writer(myfile)
		for key, value in data.items():
			writer.writerow([key, value])
            myfile.close()

            ax.scatter(x, y)
            ax.set_title('Average number of overlaps vs Coverage')
            ax.set_xlabel('Coverage')
            ax.set_ylabel('Average number of overlaps')


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
#            print reads_overlap_count[read_id]
	    read_length = float(len(reads[read_id]))
	    num_overlaps = float(reads_overlap_count[read_id])
	    # taking into account the read lengths
            reads_overlap_count[read_id] = num_overlaps / read_length 
       #     print reads_overlap_count[read_id]

	read_overlaps_count_data = output_prefix + "_read_overlaps_count_data.csv"
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

def plot_average_coverage_per_base_per_read(ax, fasta, output_prefix):
	print "\n\n\n\n"
	print "Plotting average coverage per base per read"
	print "___________________________________________"


        # get overlaps file for max_coverage total data size
        # get overlaps file created with all reads
        mean_read_length = fasta.get_mean_read_length()
        genome_size = fasta.get_genome_size()
        num_reads = fasta.get_num_reads()
        fa_filename = fasta.get_fa_filename()
        data_type = fasta.get_data_type()

	# getting the overlaps file created with all reads
        max_coverage = round(math.ceil((int(num_reads) * float(mean_read_length)) / float(genome_size)))
        overlaps_filename = create_overlaps_file(fa_filename, output_prefix, max_coverage, data_type)

        overlaps = open(overlaps_filename, 'r')
        reads_overlaps = dict()
        for line in overlaps:
            query_read_id = line.split('\t')[0]
            query_start_pos = line.split('\t')[2]
            query_end_pos = line.split('\t')[3]
            query_length = int(query_end_pos) - int(query_start_pos)
            if query_read_id in reads_overlaps: 
                reads_overlaps[query_read_id]+= query_length
            else:
                reads_overlaps[query_read_id] = 0
        overlaps.close()

        # initialize list of coverage_per_base_per_read values
        # get the reads list from the fasta file object so we can get length of reads
        reads = fasta.get_read_seqs()
        data = list()
        for read_id in reads_overlaps:
            seq = reads[read_id]
            num_bases = len(seq)
            cov_per_base_per_read = reads_overlaps[read_id] / num_bases
            data.append(tuple((read_id, cov_per_base_per_read)))

        estimated_coverage_data = output_prefix + "_estimated_coverage_data.csv"
        with open(estimated_coverage_data, 'wb') as myfile:
        	writer = csv.writer(myfile)
                writer.writerow(['read_id','avg_cov_per_base'])
		for row in data:
			writer.writerow(row)
        myfile.close()


        # plotting the cov per base per read
	values_only = list()
	values_only = [x[1] for x in data]
        ax.hist(values_only,100)
        ax.set_title('Estimated coverage distribution')
        ax.set_xlabel('Avg. number of overlaps per base per read')
        ax.set_ylabel('Frequency')   

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

