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

	fig1 = plt.figure()
	fig1.suptitle("Preqc Long Read Results : " + output_prefix )
        fig1.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

        # six subplots per pdf page
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)

        ax7 = fig1.add_subplot(321)
        ax8 = fig1.add_subplot(322)
        ax9 = fig1.add_subplot(323)
        ax10 = fig1.add_subplot(324)
        ax11 = fig1.add_subplot(325)
        ax12 = fig1.add_subplot(326)


        # get reads and some read info from fasta file
        with open(preqclr_file) as json_file:  
    		data = json.load(json_file)
		read_lengths = data['read_lengths']					# histogram
		per_coverage_avg_num_overlaps = data['coverage_avg_num_overlaps']	# scatter plot
		per_read_overlap_count = data['num_overlaps_per_read']			# histogram
		estimated_coverage_via_avg_overlaps = data['estimated_coverage'] 	# histogram
		read_lengths_estimated_cov = data['read_lengths_estimated_cov'] 	# scatter plot
		per_read_GC_content = data['read_counts_per_GC_content'] 		# histogram
		num_matching_residues = data['num_matching_residues']		# histogram
	plot_read_length_distribution(ax1, read_lengths, output_prefix)
        plot_average_overlaps_vs_coverage_distribution(ax2, per_coverage_avg_num_overlaps, output_prefix)
        plot_num_overlaps_per_read_distribution(ax3, per_read_overlap_count, output_prefix)
        plot_estimated_coverage(ax4, estimated_coverage_via_avg_overlaps, output_prefix)
	plot_estimated_coverage_per_read(ax5, read_lengths_estimated_cov, output_prefix)
	plot_per_read_GC_content(ax6, per_read_GC_content, output_prefix)
	plot_num_matching_residues(ax7, num_matching_residues, output_prefix)
	fig.savefig(pp, format='pdf')
        fig1.savefig(pp, format='pdf')

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
        ax.set_title('Per-base coverage distribution')
        ax.set_xlabel('Coverage')
        ax.set_ylabel('Frequency')    
        ax.grid(True, linestyle='-', linewidth=0.3)

def plot_estimated_coverage(ax, data, output_prefix):
	print "\n\n\n\n"
	print "Plotting average coverage per base per read"
	print "___________________________________________"

        # plotting the cov per base per read
	# the json file has reads-reads info stored as a list with all the coverage levels for each read,
	# whereas for the ref-reads info, it is already pre-counted....This is temporary.
        reads_reads = data[1]
	ref_reads = {}
	for key in data[0].keys():
		ref_reads[int(key)]= data[0][key]

	# adjusting bins
	q25, q75 = np.percentile(reads_reads, [25, 75])
	IQR = float(q75 - q25)
	n = float(len(reads_reads))
	h = int(float(2 * IQR) / float(n ** (1/3)))
	binwidth = int(float(max(reads_reads)-min(reads_reads))/float(h))
	num_bins = int(float(max(reads_reads)-min(reads_reads))/binwidth)
	bins = np.arange(0, max(reads_reads) + binwidth, binwidth)

	# plotting
	ax.hist(ref_reads.keys(), weights=ref_reads.values(), alpha=0.5, label="ref (cov)", bins=bins)
	ax.hist(reads_reads, alpha=0.5, label="reads", bins=bins)
        ax.set_title('Estimated coverage distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xticks(np.arange(0, max(reads_reads), int(5 * round(float( num_bins/2)/5))))
        ax.set_xlabel('Avg. number of overlaps per base per read')
        ax.set_ylabel('Frequency')   
	ax.legend(loc='upper right')

def plot_estimated_coverage_per_read(ax, data, output_prefix):
	# get x and y values from list of tuples
	x,y = zip(*data)
	# scatter plot with read length on x axis, and estimted coverage for read on y axis
	ax.scatter(x, y, alpha=0.4, marker='o', s=2, edgecolors='#DB5461', linewidth=0.2, facecolors='#DB5461')
        ax.set_title('Estimated coverage vs read length')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Read length')
        ax.set_ylabel('Coverage')

def plot_num_matching_residues(ax, data, output_prefix):
        binwidth = 100.0
        bins = np.arange(0, max(data) + binwidth, binwidth)
	ax.hist(data, bins=bins)
	ax.set_title('Number of matching residues per overlap distributions')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Number of matching residues')
        ax.set_ylabel('Frequency')

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

def getshape(d):
    if isinstance(d, dict):
        return {k:getshape(d[k]) for k in d}
    else:
        # Replace all non-dict values with None.
        return None

'''def filter_outliers(data, thresh=1.5):
    points = ''
    if isinstance(data, dict):
	points = data.copy()
	q75, q25 = np.percentile(points.values(), [75,25])
    elif isinstance(data, list):
	points = data[:]
	q75, q25 = np.percentile(points, [75,25])
    else:
	return points
    IQR = q75 - q25
    upperbound = float(q75 + IQR * thresh)
    lowerbound = float(q25 - IQR * thresh)

    if isinstance(points, dict):
	for key in points:
		if points[key] > upperbound:
			points.pop(key, None)
		if points[key] < lowerbound:
			points.pop(key, None)
    elif isinstance(points, list):
	for value in points:
		if value > upperbound:
			points.remove(value)
		elif value < lowerbound:
			points.remove(value)

    return points'''	
	
if __name__ == "__main__":
    main()

