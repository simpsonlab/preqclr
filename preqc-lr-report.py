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
import collections

plots_available = ['est_genome_size', 'read_length_dist', 'start_pos_per_read_dist', 'est_cov_dist', 'est_cov_vs_read_length', 'per_read_GC_content_dist']

def main():
        global plots_available

        # process input
        parser = argparse.ArgumentParser(description='Display Pre-QC Long Read report')
        parser.add_argument('-i', '--input', action="store", required=True, dest="preqc_file", nargs='+', help="preqclr file(s)")
        parser.add_argument('-o', '--output', action="store", dest="prefix", default="preqc-lr-output", help="Prefix for output pdf")
	parser.add_argument('--plot', action="store", required=False, dest="plots_requested", nargs='+', choices=plots_available, help="List of plots wanted by name.")
	parser.add_argument('--list_plots', action="store_true", dest="list_plots", default=False, help="Use to see the plots available")
        args = parser.parse_args()

	# list plots if requested
	if args.list_plots:
		i = 1
		print "List of available plots: "
		for p in plots_available:
			print '\t' + str(i) + ". " + p
			i+=1 
		print "Use --plot and names to choose which plots to create."
		sys.exit(0)

	# list plots to make
	plots = list()
	if args.plots_requested:
		plots = args.plots_requested
	else:
		# if user did not specify plots, make all plots
		plots = plots_available

	try:
		if len(args.preqc_file) > 6:
			raise ValueError
	except ValueError:
		print "Warning: Large amount of samples may not display as well"

        create_report(args.prefix, args.preqc_file, plots)

def create_report(output_prefix, preqclr_file, plots_requested):
	# calculate number of plots to create
	# total number of plots = a + b
	# a = number of plots in plots_requested not including est. cov vs read length
	# b = est. cov vs read length * number of samples
	# if est. cov vs read length requested, we need to create one for each sample

	num_samples = len(preqclr_file)
	# calculate a
	a = len(plots_requested)

	# if est cov vs read length requested
	if 'est_cov_vs_read_length' in plots_requested:
		b = num_samples - 1
	else:
		b = 0

	num_plots = a + b
	if num_samples > 7:
		print "Too many samples."
		sys.exit(1)

        # configure the PDF file 
        MPL.rc('figure', figsize=(8,10.5)) # in inches
        MPL.rc('font', size=11)
        MPL.rc('xtick', labelsize=6)
        MPL.rc('ytick', labelsize=6)
        MPL.rc('legend', fontsize=6)
        MPL.rc('axes', titlesize=10)
        MPL.rc('axes', labelsize=8)
        MPL.rcParams['lines.linewidth'] = 1.5

	# set PDF name
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
        subplots = list()
        subplots.extend((ax1, ax2, ax3, ax4, ax5, ax6))

	if num_plots > 6:
		fig1 = plt.figure()
		fig1.suptitle("Preqc Long Read Results : " + output_prefix )
        	fig1.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

        	ax7 = fig1.add_subplot(321)
        	ax8 = fig1.add_subplot(322)
        	ax9 = fig1.add_subplot(323)
        	ax10 = fig1.add_subplot(324)
        	ax11 = fig1.add_subplot(325)
        	ax12 = fig1.add_subplot(326)
		
		subplots.extend((ax7, ax8, ax9, ax10, ax11, ax12))

	estimated_genome_sizes = dict()
	per_read_read_length = dict()
	per_coverage_avg_num_overlaps = dict()
	per_read_overlap_count = dict()
	per_read_estimated_coverage = dict()
	read_lengths_estimated_cov = dict()
	per_read_GC_content = dict()
	num_matching_residues = dict()
	overlap_accuracies = dict()

	markers = ['s', 'o', '^', 'p', '+', '*', 'v']
	colors = ['#E84855', '#3C91E6', '#FFFD82', '#FF9B71', '#1B998B', '#68A691' ] 
	for sample_preqclr_file in preqclr_file:
		# get sample name
		color = colors.pop(0)
		marker = markers.pop(0)
        	with open(sample_preqclr_file) as json_file:
			data = json.load(json_file)
			sample = data['sample_name']
			estimated_genome_sizes[sample] = (color, data['estimated_genome_size'], marker)
			per_read_read_length[sample] = (color, data['per_read_read_length'], marker)
			#per_coverage_avg_num_overlaps[sample] = (color, data['coverage_avg_num_overlaps'], marker)		# scatter plot
			per_read_overlap_count[sample] = (color, data['per_read_overlap_count'], marker)				# histogram
			per_read_estimated_coverage[sample] = (color, data['per_read_estimated_coverage'], marker) 		# histogram
			read_lengths_estimated_cov[sample] = (color, data['read_lengths_estimated_cov'], marker) 		# scatter plot
			per_read_GC_content[sample] = (color, data['read_counts_per_GC_content'], marker) 			# histogram
			num_matching_residues[sample] = (color, data['num_matching_residues'], marker)				# histogram
			overlap_accuracies[sample] = (color, data['overlap_accuracies'], marker)				# scatter plot


	if 'est_genome_size' in plots_requested:
		ax = subplots.pop(0)
		plot_estimated_genome_size(ax, estimated_genome_sizes, output_prefix)
	if 'read_length_dist' in plots_requested:
                ax = subplots.pop(0)
		plot_read_length_distribution(ax, per_read_read_length, output_prefix)
	if 'start_pos_per_read_dist' in plots_requested:
		ax = subplots.pop(0)
	        plot_num_overlaps_per_read_distribution(ax, per_read_overlap_count, output_prefix) 
	if 'est_cov_dist' in plots_requested:
		ax = subplots.pop(0)
       	 	plot_estimated_coverage(ax, per_read_estimated_coverage, output_prefix)
        if 'per_read_GC_content_dist' in plots_requested:
                ax = subplots.pop(0)
                plot_per_read_GC_content(ax, per_read_GC_content, output_prefix)
	if 'est_cov_vs_read_length' in plots_requested:
        	for sample in read_lengths_estimated_cov:
                	ax = subplots.pop(0)
			plot_estimated_coverage_vs_read_length(ax, read_lengths_estimated_cov, sample, output_prefix)

	#	plot_mean_overlap_accuracies_per_read(ax, overlap_accuracies, output_prefix)
	#plot_num_matching_residues(ax8, num_matching_residues, output_prefix)
        #plot_average_overlaps_vs_coverage_distribution(ax3, per_coverage_avg_num_overlaps, output_prefix)
	fig.savefig(pp, format='pdf', dpi=1000)
	if num_plots > 6:
	        fig1.savefig(pp, format='pdf')

        pp.close()

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
			x_lim = np.percentile(sample_read_lengths, 99)
        binwidth = 1000.0
        bins = np.arange(0, max_read_length + binwidth, binwidth)

	# now start plotting
	for sample in read_lengths:
                sample_name = sample
                sample_color = read_lengths[sample][0]
                sample_data = read_lengths[sample][1]
		base = 100
                sample_data_rounded = [ int(base * round(float(x)/base)) for x in sample_data ]
                labels, values = zip(*sorted(collections.Counter(sorted(sample_data_rounded)).items()))
	        ax.plot(labels, [float(i) for i in values], color=sample_color, label=sample_name)

        ax.set_title('Read length distribution')
        ax.set_xlabel('Read lengths (bps)')
        ax.set_ylabel('Frequency')
	ax.set_xlim(0, x_lim)
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
        print "Plotting number of overlaps per read"
        print "___________________________________________"

        max_num_overlaps = 0
        for sample in data:
                sample_data = data[sample][1]
                sample_max_num_overlaps = max(sample_data.values())
                if sample_max_num_overlaps > max_num_overlaps:
                        max_num_overlaps = sample_max_num_overlaps
			x_lim = np.percentile(sample_data.values(), 99)
        binwidth = 0.001
        bins = np.arange(0, float(max_num_overlaps) + binwidth, binwidth)

	# now start plotting
        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		sample_data_rounded = [ round(x, 4) for x in sample_data.values() ]
                labels, values = zip(*sorted(collections.Counter(sorted(sample_data_rounded)).items()))
		print labels
	        ax.plot(labels, values, label=sample_name, color=sample_color)

        # plotting the number of overlaps/read
        ax.set_title('Start positions per read distribution')
        ax.set_xlabel('Number of overlaps')
        ax.set_ylabel('Frequency')    
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')
	ax.set_xlim(0, float(x_lim))

def plot_estimated_coverage(ax, data, output_prefix):
	print "\n\n\n\n"
	print "Plotting estimated coverage per read"
	print "___________________________________________"

        max_cov = 0
	binwidth = 0
	num_bins = 0
        for sample in data:
                sample_data = data[sample][1]
		sample_data_reads = sample_data[0]
                sample_data_upperbound = sample_data[1]		# precalculated est. cov. upperbound 
		sample_data_num_reads = sample_data[2]		# number of reads after filtering
		sample_data_q25 = sample_data[3]		# 1st quartile cov after filtering
		sample_data_q75 = sample_data[4]		# 3rd quartile cov after filtering
                if sample_data_upperbound > max_cov:
                        max_cov = sample_data_upperbound
        		IQR = float(sample_data_q75 - sample_data_q25)
        		n = float(sample_data_num_reads)
        		binwidth = int(float(2 * IQR) / float(n ** (1/3)))
			print binwidth
        		num_bins = int(float(max_cov-min(sample_data_reads))/binwidth)
        print max(sample_data_reads)
	bins = np.arange(0, max_cov + binwidth, num_bins)
	print num_bins
	for sample in data:
		sample_name = sample
		sample_color = data[sample][0]
                sample_data = data[sample][1]
                sample_data_reads = sample_data[0]
		sample_data_upperbound = sample_data[1]
		labels, values = zip(*collections.Counter(sample_data_reads).items())
		ax.plot(labels, values, color=sample_color, label=sample_name)

        ax.set_title('Estimated coverage distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	#ax.set_xticks(np.arange(0, max_cov + binwidth, num_bins*5.0))
        ax.set_xlabel('Estimated coverage')
        ax.set_ylabel('Frequency')   
	ax.set_xlim(0, max_cov)
	ax.legend(loc='upper right')

def plot_estimated_coverage_vs_read_length(ax, data, sample, output_prefix):
        print "\n\n\n\n"
        print "Plotting estimated coverage vs read length scatter plot"
        print "___________________________________________"

	# get the specific sample's data
        sample_name = sample
        sample_color = data[sample][0]
        sample_data = data[sample][1]
	sample_marker = data[sample][2]
        # get x and y values from list of tuples
        x,y = zip(*sample_data)
	ax.scatter(x, y, alpha=0.2, marker=sample_marker, s=4, edgecolors=sample_color, linewidth=0.5, facecolors=sample_color, rasterized=True)
	
	# scatter plot with read length on x axis, and estimted coverage for read on y axis
        ax.set_title('Estimated coverage vs read length')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Read length (bps)')
        ax.set_ylabel('Estimated coverage')
        ax.legend(loc='upper right')
	#ax.set_yscale('log')

def plot_num_matching_residues(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting number of matching residues"
        print "___________________________________________"

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
        print "\n\n\n\n"
        print "Plotting GC content"
        print "___________________________________________"

        per_read_GC_content = {}
        binwidth = 1.0
	for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
		print sample_data
		for GC_content_level in sample_data:
			per_read_GC_content[int(float(GC_content_level))] = int(sample_data[GC_content_level])
#		ax.plot(per_read_GC_content.keys(), weights=per_read_GC_content.values(), color=sample_color, bins=np.arange(0, 100, binwidth), alpha=0.
#6)
		ax.plot(per_read_GC_content.keys(), per_read_GC_content.values(), color=sample_color, label=sample_name)
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
		genome_size_in_megabases = round(float(sample_data)/float(1000000), 2)
		genome_sizes.append(genome_size_in_megabases)
		sample_names.append(sample_name)
		colors.append(str(sample_color))
		print sample_data
		print sample_name

	print sample_names
        # plotting the number of overlaps/read
	y_pos = np.arange(len(genome_sizes))
	ax.barh(y_pos, genome_sizes, align='center')
	# add the actual estimated genome size values to the bars
	rects = ax.patches

	if len(sample_names) == 1 :
		height = 0.3
	else:
		height = 1.0/float(len(sample_names))

	for rect, value, sample in zip(rects, genome_sizes, sample_names):
		rect.set_facecolor(data[sample][0])
		rect.set_height(height)
		#genome_size_in_megabases = round(float(value)/float(1000000), 2)
		t = ax.text(0.05, rect.get_y() + rect.get_height()/2.0, sample_name + ": " + str(value) + "Mbp", ha='left', va='center', fontsize='8')
		t.set_bbox(dict(facecolor='#FFFFFF', alpha=0.5, edgecolor='#FFFFFF'))

	ax.set_yticks([])
        ax.set_title('Estimated genome size')
        ax.set_xlabel('Genome size (Mbp)')
	ax.set_ylabel(' \n \n ')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.legend(loc='upper right')

def plot_mean_overlap_accuracies_per_read(ax, data, output_prefix):
        print "\n\n\n\n"
        print "Plotting mean overlap accuracies per read vs read length"
        print "___________________________________________"

        for sample in data:
                sample_name = sample
                sample_color = data[sample][0]
                sample_data = data[sample][1]
                sample_marker = data[sample][2]
                # get x and y values from list of tuples
                y,x = zip(*sample_data)
                ax.scatter(x, y, alpha=0.1, marker=sample_marker, s=4, edgecolors=sample_color, linewidth=0.5, facecolors=sample_color, rasterized=True)

        # scatter plot with read length on x axis, and estimted coverage for read on y axis
        ax.set_title('Mean accuracy vs read length')
        ax.grid(True, linestyle='-', linewidth=0.3)
        ax.set_xlabel('Read length')
        ax.set_ylabel('Mean accuracy')
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

