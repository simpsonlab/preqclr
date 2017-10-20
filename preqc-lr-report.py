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
import json

def main():
        # process input
        parser = argparse.ArgumentParser(description='Display Pre-QC Long Read report')
        parser.add_argument('-i', '--input', action="store", required=True, dest="preqc_file", help="PreQC file")
        parser.add_argument('-o', '--output', action="store", dest="prefix", default="preqc-lr-output", help="Prefix for output pdf")
        args = parser.parse_args()
        
        create_report(args.prefix, args.preqc_file)

def create_report(output_prefix, preqclr_file):

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

        output_pdf = output_prefix + ".pdf"
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
        with open(preqclr_file) as json_file:  
    		data = json.load(json_file)
		read_lengths = data['read_lengths']					# histogram
		per_coverage_avg_num_overlaps = data['coverage_avg_num_overlaps']	# scatter plot
		per_read_overlap_count = data['num_overlaps_per_read']			# histogram
		estimated_coverage_via_avg_overlaps = data['estimated_coverage'] 	# histogram
		read_lengths_estimated_cov = data['read_lengths_estimated_cov'] 	# scatter plot
		per_read_GC_content = data['read_counts_per_GC_content'] 		# histogram
	
	plot_read_length_distribution(ax1, read_lengths, output_prefix)
        plot_average_overlaps_vs_coverage_distribution(ax2, per_coverage_avg_num_overlaps, output_prefix)
        plot_num_overlaps_per_read_distribution(ax3, per_read_overlap_count, output_prefix)
        plot_estimated_coverage(ax4, estimated_coverage_via_avg_overlaps, output_prefix)
	plot_estimated_coverage_per_read(ax5, read_lengths_estimated_cov, output_prefix)
	plot_per_read_GC_content(ax6, per_read_GC_content, output_prefix)

	fig.savefig(pp, format='pdf')
        pp.close()

def plot_read_length_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting read length distribution"
        print "___________________________________________"

	read_lengths = data

        binwidth = 1000.0
        bins = np.arange(0, max(read_lengths) + binwidth, binwidth)
        ax.hist(read_lengths, bins=bins)
        ax.set_title('Read length distribution')
        ax.set_xlabel('Read lengths (bps)')
        ax.set_ylabel('Frequency')
        ax.grid(True, linestyle='-', linewidth=0.3)
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	ax.get_xaxis().get_major_formatter().set_useOffset(False)

def plot_average_overlaps_vs_coverage_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting average overlaps vs coverage"
        print "___________________________________________"

	# sorted by key, return a list of tuples
	lists = sorted(data.items())
	# unpack a list of pairs into two tuples
	x, y = zip(*lists)

        ax.scatter(x, y)
	ax.set_title('Average number of overlaps vs coverage')
	ax.set_xlabel('Coverage')
	ax.set_ylabel('Average number of overlaps')
	ax.grid(True, linestyle='-', linewidth=0.3)

def plot_num_overlaps_per_read_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting number of overlaps per reads"
        print "___________________________________________"

	reads_overlap_count = data

        # plotting the number of overlaps/read
        binwidth = 0.01
        bins = np.arange(0, max(reads_overlap_count.values()) + binwidth, binwidth)
        ax.hist(reads_overlap_count.values(), bins=bins)
        ax.set_title('Number of overlaps/base distribution')
        ax.set_xlabel('(Number of overlaps/read)/length')
        ax.set_ylabel('Frequency')    
        ax.grid(True, linestyle='-', linewidth=0.3)

def plot_estimated_coverage(ax, data, output_prefix):
	print "\n\n\n\n"
	print "Plotting average coverage per base per read"
	print "___________________________________________"

        # plotting the cov per base per read
	ref_reads = {}
	reads_reads = {}
	ref_reads_1 = {}
	for key in data[0].keys():
		ref_reads[int(key)]= data[0][key]
	for key in data[1].keys():
                reads_reads[int(key)]= data[1][key]
	print ref_reads
	binwidth = 5.0	
	bins = np.arange(0, max(ref_reads.keys()) + binwidth, binwidth)

	ax.hist(ref_reads.keys(), weights=ref_reads.values(), alpha=0.5, label="ref (cov)", bins=bins)
        ax.hist(reads_reads.keys(), weights=reads_reads.values(), alpha=0.5, label="reads", bins=bins)
	#ax.hist(ref_reads_1.keys(), weights=ref_reads_1.values(),  alpha=0.5, label="ref (cov1)", bins=bins)
        ax.set_title('Estimated coverage distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xticks(np.arange(0, max(ref_reads.keys()), 50.0))
        ax.set_xlabel('Avg. number of overlaps per base per read')
        ax.set_ylabel('Frequency')   
	ax.legend(loc='upper right')

def plot_estimated_coverage_per_read(ax, data, output_prefix):
	
	# get x and y values from list of tuples
	x,y = zip(*data)
	# scatter plot with read length on x axis, and estimted coverage for read on y axis
	ax.scatter(x, y, alpha=0.4, marker="+")
        ax.set_title('Estimated coverage vs read length')
        ax.grid(True, linestyle='-', linewidth=0.3)
        #ax.set_xticks(np.arange(0, max(x), 50.0))
        #ax.set_yticks(np.arange(0, max(y) + 500, 50)
        ax.set_xlabel('Read length')
        ax.set_ylabel('Coverage')

def plot_per_read_GC_content(ax, data, output_prefix):
	per_read_GC_content = {}
	print "Per read GC content"
	print data
	for GC_content_level in data:
		per_read_GC_content[int(float(GC_content_level))] = int(data[GC_content_level])
        binwidth = 1.0
	ax.hist(per_read_GC_content.keys(), weights=per_read_GC_content.values(), bins=np.arange(0, 100, binwidth), alpha=0.6)
        ax.set_title('Per read GC content')
        ax.set_xlabel('% GC content')
        ax.set_ylabel('Frequency')
        ax.grid(True, linestyle='-', linewidth=0.3)
		
if __name__ == "__main__":
    main()

