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
        parser.add_argument('-i', '--input', action="store", required=True, dest="preqc_file", nargs='+', help="preqclr file(s)")
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

	estimated_genome_sizes = dict()
	read_lengths = dict()
	per_coverage_avg_num_overlaps = dict()
	per_read_overlap_count = dict()
	estimated_coverage_via_avg_overlaps = dict()
	read_lengths_estimated_cov = dict()
	per_read_GC_content = dict()
	num_matching_residues = dict()

	markers = ['s', 'o', '^', 'p', '+']
	colors = ['#E84855', '#3C91E6', '#FFFD82', '#FF9B71', '#1B998B'] 
	for sample_preqclr_file in preqclr_file:
		# get sample name
		color = colors.pop(0)
		marker = markers.pop(0)
        	with open(sample_preqclr_file) as json_file:
			data = json.load(json_file)
			sample = data['sample_name']
			estimated_genome_sizes[sample] = (color, data['estimated_genome_size'], marker)
			read_lengths[sample] = (color, data['read_lengths'], marker)
			per_coverage_avg_num_overlaps[sample] = (color, data['coverage_avg_num_overlaps'], marker)		# scatter plot
			per_read_overlap_count[sample] = (color, data['num_overlaps_per_read'], marker)				# histogram
			estimated_coverage_via_avg_overlaps[sample] = (color, data['estimated_coverage'], marker) 		# histogram
			read_lengths_estimated_cov[sample] = (color, data['read_lengths_estimated_cov'], marker) 		# scatter plot
			per_read_GC_content[sample] = (color, data['read_counts_per_GC_content'], marker) 			# histogram
			num_matching_residues[sample] = (color, data['num_matching_residues'], marker)				# histogram

	plot_estimated_genome_size(ax1, estimated_genome_sizes, output_prefix)
	plot_read_length_distribution(ax2, read_lengths, output_prefix)
        plot_average_overlaps_vs_coverage_distribution(ax3, per_coverage_avg_num_overlaps, output_prefix)
        plot_num_overlaps_per_read_distribution(ax4, per_read_overlap_count, output_prefix)
        plot_estimated_coverage(ax5, estimated_coverage_via_avg_overlaps, output_prefix)
	plot_estimated_coverage_per_read(ax6, read_lengths_estimated_cov, output_prefix)
	plot_per_read_GC_content(ax7, per_read_GC_content, output_prefix)
	plot_num_matching_residues(ax8, num_matching_residues, output_prefix)
	fig.savefig(pp, format='pdf')
        fig1.savefig(pp, format='pdf')

        pp.close()

class sample(object):
	name = ""
	color = ""
	genome_size_estimate = ""
	read_lengths = list()
	per_coverage_avg_num_overlaps = list()
	per_read_overlap_count = dict()
	estimated_coverage_via_avg_overlaps = dict()
	read_lengths_estimated_cov = list()
	per_read_GC_content = list()
	num_matching_residues = list()
 
	def __init__(self, read_lengths, per_coverage_avg_num_overlaps, per_read_overlap_count, estimated_coverage_via_avg_overlaps, read_lengths_estimated_cov, per_read_GC_content, num_matching_residues):
		self.read_lengths =  read_lengths
		self.per_coverage_avg_num_overlaps = per_coverage_avg_num_overlaps
		self.per_read_overlap_count = per_read_overlap_count
		self.estimated_coverage_via_avg_overlaps = estimated_coverage_via_avg_overlaps
		self.read_lengths_estimated_cov = read_lengths_estimated_cov
		self.per_read_GC_content = per_read_GC_content
		self.num_matching_residues = num_matching_residues
	
	def get_read_lengths(self):
		return self.read_lengths

	def get_per_coverage_avg_num_overlaps(self):
		return self.per_coverage_avg_num_overlaps

	def get_per_read_overlap_count(self):
		return self.per_read_overlap_count

	def get_estimated_coverage_via_avg_overlaps(self):
		return self.estimated_coverage_via_avg_overlaps

	def get_read_lengths_estimated_cov(self):
		return self.read_lengths_estimated_cov

	def get_per_read_GC_content(self):
		return self.per_read_GC_content

	def get_genome_size_estimate(self):
		return self.genome_size_estimate

def plot_read_length_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting read length distribution"
        print "___________________________________________"

	read_lengths = data

	# get the maximum read length across all the samples
	max_read_length = 0
	for sample in read_lengths:
		sample_name = sample
		sample_color = read_lengths[sample][0]
		sample_read_lengths = read_lengths[sample][1]
		sample_max_read_length = max(sample_read_lengths)
		if sample_max_read_length > max_read_length:
			max_read_length = sample_max_read_length
        binwidth = 1000.0
        bins = np.arange(0, max_read_length + binwidth, binwidth)

	# now start plotting
	for sample in read_lengths:
                sample_name = sample
                sample_color = read_lengths[sample][0]
                sample_read_lengths = read_lengths[sample][1]
                ax.hist(sample_read_lengths, color=sample_color, label=sample_name, alpha=0.5, bins=bins)

        ax.set_title('Read length distribution')
        ax.set_xlabel('Read lengths (bps)')
        ax.set_ylabel('Frequency')
        ax.grid(True, linestyle='-', linewidth=0.3)
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc='upper right')

def plot_average_overlaps_vs_coverage_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting average overlaps vs coverage"
        print "___________________________________________"

	for sample in data:
		sample_name = sample
		sample_color = data[sample][0]
		sample_marker = data[sample][2]
		sample_data = data[sample][1]
		# sorted by key, return a list of tuples
		lists = sorted(sample_data.items())
		# unpack a list of pairs into two tuples
		x, y = zip(*lists)
        	ax.scatter(x, y, color=sample_color, marker=sample_marker, alpha=0.6, label=sample_name)

	ax.set_title('Average number of overlaps vs coverage')
	ax.set_xlabel('Coverage')
	ax.set_ylabel('Average number of overlaps')
	ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')

def plot_num_overlaps_per_read_distribution(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting number of overlaps per reads"
        print "___________________________________________"

        max_num_overlaps = 0
        for sample in data:
                sample_data = data[sample][1]
                sample_max_num_overlaps = max(sample_data.values())
                if sample_max_num_overlaps > max_num_overlaps:
                        max_num_overlaps = sample_max_num_overlaps
        binwidth = 0.01
        bins = np.arange(0, float(max_num_overlaps) + binwidth, binwidth)

	# now start plotting
        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
	        ax.hist(sample_data.values(), bins=bins, alpha=0.5, label=sample_name, color=sample_color)

        # plotting the number of overlaps/read
        ax.set_title('Per-base coverage distribution')
        ax.set_xlabel('Coverage')
        ax.set_ylabel('Frequency')    
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')

def plot_estimated_coverage(ax, data, output_prefix):
	print "\n\n\n\n"
	print "Plotting average coverage per base per read"
	print "___________________________________________"

        max_cov = 0
	binwidth = 0
	num_bins = 0
        for sample in data:
                sample_data = data[sample][1]
		sample_data_reads = sample_data[1]
                sample_max_cov = max(sample_data_reads)
                if sample_max_cov > max_cov:
                        max_cov = sample_max_cov
		        q25, q75 = np.percentile(sample_data_reads, [25, 75])
        		IQR = float(q75 - q25)
        		n = float(len(sample_data_reads))
        		h = int(float(2 * IQR) / float(n ** (1/3)))
        		binwidth = int(float(max(sample_data_reads)-min(sample_data_reads))/float(h))
        		num_bins = int(float(max(sample_data_reads)-min(sample_data_reads))/binwidth)
        
	bins = np.arange(0, max_cov + binwidth, binwidth)

	for sample in data:
		sample_name = sample
		sample_color = data[sample][0]
                sample_data = data[sample][1]
                sample_data_reads = sample_data[1]
		sample_data_ref = sample_data[0]

		# add read-reference info
        	for key in sample_data[0].keys():
                	sample_data_ref[int(key)] = sample_data[0][key]
		#ax.hist(sample_data_ref.keys(), weights=sample_data_ref.values(), alpha=0.3, label="ref", bins=bins)
		ax.hist(sample_data_reads, alpha=0.5, color=sample_color, label=sample_name, bins=bins)

        ax.set_title('Coverage distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xticks(np.arange(0, max_cov, int(5 * round(float( num_bins/2)/5))))
        ax.set_xlabel('Avg. number of overlaps per base per read')
        ax.set_ylabel('Frequency')   
	ax.legend(loc='upper right')

def plot_estimated_coverage_per_read(ax, data, output_prefix):

        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		sample_marker = data[sample][2]
        	# get x and y values from list of tuples
        	x,y = zip(*sample_data)
	        ax.scatter(x, y, alpha=0.5, marker=sample_marker, s=4, edgecolors=sample_color, linewidth=0.5, facecolors='None')
	
	# scatter plot with read length on x axis, and estimted coverage for read on y axis
        ax.set_title('Coverage vs read length')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Read length')
        ax.set_ylabel('Coverage')
        ax.legend(loc='upper right')

def plot_num_matching_residues(ax, data, output_prefix):

	max_num_matching_residues = 0
	for sample in data:
                sample_data = data[sample][1]
		max_sample_data = max(sample_data)	
		if max_sample_data > max_num_matching_residues:
			max_num_matching_residues = max_sample_data

        binwidth = 100.0
        bins = np.arange(0, max_num_matching_residues + binwidth, binwidth)

        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		ax.hist(sample_data, bins=bins, alpha=0.5, label=sample_name, color=sample_color)	

	ax.set_title('Number of matching residues per overlap distributions')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Number of matching residues')
        ax.set_ylabel('Frequency')

def plot_per_read_GC_content(ax, data, output_prefix):
        per_read_GC_content = {}
        binwidth = 1.0
	for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		for GC_content_level in sample_data:
			per_read_GC_content[int(float(GC_content_level))] = int(sample_data[GC_content_level])
		ax.hist(per_read_GC_content.keys(), weights=per_read_GC_content.values(), color=sample_color, bins=np.arange(0, 100, binwidth), alpha=0.6)
        ax.set_title('Per read GC content')
        ax.set_xlabel('% GC content')
        ax.set_ylabel('Frequency')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')

def plot_estimated_genome_size(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting genome size estimates"
        print "___________________________________________"

	genome_sizes = list()
	sample_names = list()
	colors = list()
        # now start plotting
        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		genome_sizes.append(sample_data)
		sample_names.append(sample_name)
		colors.append(str(sample_color))
		print sample_data
		print sample_name

	print sample_names
        # plotting the number of overlaps/read
	y_pos = np.arange(len(genome_sizes))
	ax.bar(y_pos, genome_sizes, align='center', alpha=0.5, color=colors)
	ax.set_xticks(y_pos)
	ax.set_xticklabels(sample_names)
        ax.set_title('Estimated genome size')
        ax.set_xlabel('Samples')
        ax.set_ylabel('Genome size (bps)')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')
	

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

