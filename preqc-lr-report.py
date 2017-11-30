#!/usr/bin/env python
# ========================================================
# preqc-lr report:
# Generates plots and saves to a preqc-lr report in PDF
# ========================================================
try:
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
	from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
	custom_print('Missing package(s)')	

plots_available = ['est_genome_size', 'read_length_dist', 'start_pos_per_read_dist', 'est_cov_dist', 'est_cov_vs_read_length', 'per_read_GC_content_dist', 'total_num_bases_vs_min_read_length']
save_png=False

def main():
	print "========================================================"
	print "RUNNING PREQC-LR REPORT"
	print "========================================================"
	global plots_available

	# process input
	parser = argparse.ArgumentParser(description='Display Pre-QC Long Read report')
	parser.add_argument('-i', '--input', action="store", required=True, dest="preqc_file", nargs='+', help="preqclr file(s)")
	parser.add_argument('-o', '--output', action="store", dest="output_prefix", default="preqc-lr-output", help="Prefix for output pdf")
	parser.add_argument('--plot', action="store", required=False, dest="plots_requested", nargs='+', choices=plots_available, help="List of plots wanted by name.")
	parser.add_argument('--list_plots', action="store_true", dest="list_plots", default=False, help="Use to see the plots available")
	parser.add_argument('--save_png', action="store_true", dest="save_png", default=False, help= "Use to save png for each plot.")
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

	# set global variable of save png
	global save_png
	if args.save_png:
		save_png=True
		png_dir = "./" + args.output_prefix + "/png/"
		if not os.path.exists(png_dir):
			os.makedirs(png_dir)

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
		print "Warning: Large amount of ss may not display as well"

	create_report(args.output_prefix, args.preqc_file, plots)

def create_report(output_prefix, preqclr_file, plots_requested):
	# calculate number of plots to create
	# total number of plots = a + b
	# a = number of plots in plots_requested not including est. cov vs read length
	# b = est. cov vs read length * number of ss
	# if est. cov vs read length requested, we need to create one for each s

	num_ss = len(preqclr_file)
	# calculate a
	a = len(plots_requested)

	# if est cov vs read length requested
	if 'est_cov_vs_read_length' in plots_requested:
		b = num_ss - 1
	else:
		b = 0

	num_plots = a + b

	if num_ss > 7:
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
	MPL.rcParams['lines.linewidth'] = 1.0

	# set PDF name
	output_pdf = output_prefix + ".pdf"
	pp = PdfPages(output_pdf)

	i = 0
	figs = list()
	subplots = list()

	while i < num_plots:
		# add figure
		fig = plt.figure()
		fig.suptitle("Preqc Long Read Results : " + output_prefix )
		fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)

		# six subplots per pdf page
		ax1 = fig.add_subplot(321)
		ax2 = fig.add_subplot(322)
		ax3 = fig.add_subplot(323)
		ax4 = fig.add_subplot(324)
		ax5 = fig.add_subplot(325)
		ax6 = fig.add_subplot(326)

		subplots.extend((ax1, ax2, ax3, ax4, ax5, ax6))
		figs.append(fig)

		i+=6

	# extract and save data from json file
	est_genome_sizes = dict()
	per_read_read_length = dict()
	per_read_overlap_count = dict()
	per_read_est_cov_and_read_length = dict()
	est_cov_post_filter_info = dict()
	per_read_GC_content = dict()
	total_num_bases_vs_min_read_length = dict()
	ngx_values = dict()
	nx_values = dict()

	# each s will be represented with a unique marker and color
	markers = ['s', 'o', '^', 'p', '+', '*', 'v']
	colors = ['#FC614C', '#2DBDD8', '#B4E23D', '#F7C525', '#5B507A', '#0A2463' ] 

	# start reading the preqclr file(s)
	for s_preqclr_file in preqclr_file:
		color = colors.pop(0)
		marker = markers.pop(0)
		with open(s_preqclr_file) as json_file:
			data = json.load(json_file)
			s = data['sample_name']
			print s
			est_genome_sizes[s] = (color, data['est_genome_size'], marker)													# bar graph
			per_read_read_length[s] = (color, data['per_read_read_length'], marker)											# histogram
			per_read_overlap_count[s] = (color, data['per_read_overlap_count'], marker)										# histogram
			per_read_est_cov_and_read_length[s] = (color, data['per_read_est_cov_and_read_length'], marker) 				# histogram
			est_cov_post_filter_info[s] = data['est_cov_post_filter_info']
			per_read_GC_content[s] = (color, data['read_counts_per_GC_content'], marker) 									# histogram
			total_num_bases_vs_min_read_length[s] = (color, data['total_num_bases_vs_min_read_length'], marker)				# line
			ngx_values[s] = (color, data['NGX_values'], marker)
			nx_values[s] = (color, data['NX_values'], marker)

	print per_read_read_length.keys()
	# parameters for saving each subplot as a png
	expand_x = 1.4
	expand_y = 1.28

	if 'est_genome_size' in plots_requested:
		ax = subplots.pop(0)
		ax_temp = plot_est_genome_size(ax, est_genome_sizes, output_prefix)
		if save_png:
			temp_fig = ax_temp.get_figure()
			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
			ax_png_file = "./" + output_prefix + "/png/plot_est_genome_size.png"
			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)	
	if 'read_length_dist' in plots_requested:
		ax = subplots.pop(0)
		ax_temp = plot_read_length_distribution(ax, per_read_read_length, output_prefix)
		if save_png:
			temp_fig = ax_temp.get_figure()
			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
			ax_png_file = "./" + output_prefix + "/png/plot_read_length_distribution.png"
			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)
#	if 'start_pos_per_read_dist' in plots_requested:
#		ax = subplots.pop(0)
#		ax_temp = plot_num_overlaps_per_read_distribution(ax, per_read_overlap_count, output_prefix) 
#		if save_png:
#			temp_fig = ax_temp.get_figure()
#			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
#			ax_png_file = "./" + output_prefix + "/png/plot_start_pos_per_read_dist.png"
#			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)
	if 'est_cov_dist' in plots_requested:
		ax = subplots.pop(0)
		ax_temp = plot_est_cov(ax, per_read_est_cov_and_read_length, est_cov_post_filter_info, output_prefix)
		if save_png:
			temp_fig = ax_temp.get_figure()
			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
			ax_png_file = "./" + output_prefix + "/png/plot_est_cov.png"
			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)
	if 'per_read_GC_content_dist' in plots_requested:
		ax = subplots.pop(0)
		ax_temp = plot_per_read_GC_content(ax, per_read_GC_content, output_prefix)
		if save_png:
			temp_fig = ax_temp.get_figure()
			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
			ax_png_file =  "./" + output_prefix + "/png/plot_per_read_GC_content.png"
			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)
	if 'est_cov_vs_read_length' in plots_requested:
		for s in per_read_est_cov_and_read_length:
			ax = subplots.pop(0)
			ax_temp = plot_per_read_est_cov_vs_read_length(ax, per_read_est_cov_and_read_length, s, output_prefix)
			if save_png:
				temp_fig = ax_temp.get_figure()
				extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
				ax_png_file = "./" + output_prefix + "/png/plot_est_cov_vs_read_length_" + s +".png" 
				temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)
	if 'total_num_bases_vs_min_read_length' in plots_requested:
		ax = subplots.pop(0)
		ax_temp = plot_total_num_bases_vs_min_read_length(ax, total_num_bases_vs_min_read_length, output_prefix)
		if save_png:
			temp_fig = ax_temp.get_figure()
			extent = ax_temp.get_window_extent().transformed(temp_fig.dpi_scale_trans.inverted())
			ax_png_file = "./" + output_prefix + "/png/plot_total_num_bases_vs_min_read_length.png"
			temp_fig.savefig(ax_png_file, bbox_inches=extent.expanded(expand_x, expand_y), dpi=700)

	ax = subplots.pop(0)
	ax_temp = plot_assembly_quality_metrics(ax, nx_values, output_prefix, "NX")

	ax = subplots.pop(0)
	ax_temp = plot_assembly_quality_metrics(ax, ngx_values, output_prefix, "NGX")

	for f in figs:
		f.savefig(pp, format='pdf', dpi=1000)

	pp.close()

def plot_read_length_distribution(ax, data, output_prefix):
	print "[ Plotting read length distribution ]"

	read_lengths = data

	# get the maximum read length across all the ss
	max_read_length = 0
	for s in read_lengths:
		s_name = s
		s_color = read_lengths[s][0]
		s_read_lengths = read_lengths[s][1]
		s_max_read_length = max(s_read_lengths)
		if s_max_read_length > max_read_length:
			max_read_length = s_max_read_length
			x_lim = np.percentile(s_read_lengths, 95)
		binwidth = 1000.0
		bins = np.arange(0, max_read_length + binwidth, binwidth)

	# now start plotting
	max_y = 0
	for s in read_lengths:
		s_name = s
		s_color = read_lengths[s][0]
		sd = read_lengths[s][1]
		base = 100
		sd_rounded = [ int(base * round(float(x)/base)) for x in sd ]
		labels, values = zip(*sorted(collections.Counter(sorted(sd_rounded)).items()))
		# normalize labels
		s = sum(values)
		nlabels = list()
		nvalues = list()
		for v in values:
			nv = float(v)/float(s)
			nvalues.append(nv)
		if max(nvalues) > max_y:
			max_y = max(nvalues)
		ax.plot(labels, [float(i) for i in nvalues], color=s_color, label=s_name)

	ax.set_title('Read length distribution')
	ax.set_xlabel('Read lengths (bp)')
	ax.set_ylabel('Proportion')
	ax.set_xlim(0, x_lim)
	ax.set_ylim(0, max_y)
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.legend(loc='upper right')
	return ax

def plot_average_overlaps_vs_cov_distribution(ax, data, output_prefix):
	print "[ Plotting average overlaps vs cov ]"

	for s in data:
		s_name = s
		s_color = data[s][0]
		s_marker = data[s][2]
		sd = data[s][1]
		# sorted by key, return a list of tuples
		lists = sorted(sd.items())
		# unpack a list of pairs into two tuples
		x, y = zip(*lists)
		# normalize the data
		s = sum(y)
		ny = list()
		for i in y:
			ni = float(i)/float(s)
			ny.append(ni)

		ax.scatter(x, ny, color=s_color, marker=s_marker, alpha=0.6, label=s_name)

	ax.set_title('Average number of overlaps vs cov')
	ax.set_xlabel('Coverage')
	ax.set_ylabel('Average number of overlaps')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.legend(loc='upper right')
	return ax

def plot_num_overlaps_per_read_distribution(ax, data, output_prefix):
	print "[ Plotting per read number of overlaps distribution ]"

	max_num_overlaps = 0
	for s in data:
		sd = data[s][1]
		s_max_num_overlaps = max(sd.values())
		if s_max_num_overlaps > max_num_overlaps:
			max_num_overlaps = s_max_num_overlaps
			x_lim = np.percentile(sd.values(), 90)
			binwidth = 0.001
			bins = np.arange(0, float(max_num_overlaps) + binwidth, binwidth)

	# now start plotting
		for s in data:
				s_name = s
				s_color = data[s][0]
				sd = data[s][1]
				sd_rounded = [ round(x, 4) for x in sd.values() ]
				labels, values = zip(*sorted(collections.Counter(sorted(sd_rounded)).items()))
				s = sum(values)
				nvalues = list()
				for v in values:
					nv = float(v)/float(s)
					nvalues.append(nv)
				ax.plot(labels, nvalues, label=s_name, color=s_color)

	# plotting the number of overlaps/read
	ax.set_title('Start positions per read distribution')
	ax.set_xlabel('Number of overlaps')
	ax.set_ylabel('Proportion')	
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.legend(loc='upper right')
	ax.set_xlim(0, float(x_lim))
	return ax

def plot_est_cov(ax, data, est_cov_post_filter_info, output_prefix):
	print "[ Plotting est cov per read ]"

	max_cov = 0
	binwidth = 0
	num_bins = 0
	for s in data:
		sd = data[s]														# list of tuples (read_cov, read_len)
		sd_est_cov_read_length = data[s][1]
		sd_est_cov = [i[0] for i in data[s][1]]							# retrieving only coverage information from list of tuples
		sd_est_cov_post_filter_info = est_cov_post_filter_info[s]			# tuple (lowerbound cov, upperbound cov, num reads, IQR covs)
		print sd_est_cov_post_filter_info
		sd_upperbound = sd_est_cov_post_filter_info[1]				# precalculated est. cov. upperbound 
		sd_num_reads = sd_est_cov_post_filter_info[2]					# number of reads after filtering
		sd_IQR = sd_est_cov_post_filter_info[3]						# IQR cov after filtering
		if sd_upperbound > max_cov:
			max_cov = sd_upperbound
			IQR = float(sd_IQR) 
			n = float(sd_num_reads)
			binwidth = int(float(2 * IQR) / float(n ** (1/3)))
			num_bins = int(float(max_cov)/binwidth)
	bins = np.arange(0, max_cov + binwidth, num_bins)
	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		sd_est_cov = [i[0] for i in sd]
		sd_est_cov_post_filter_info = est_cov_post_filter_info[s]
		sd_upperbound = sd_est_cov_post_filter_info[1]
		labels, values = zip(*collections.Counter(sd_est_cov).items())
		ax.plot(labels, values, color=s_color, label=s_name)

	ax.set_title('Est. cov. distribution')
	ax.grid(True, linestyle='-', linewidth=0.3)
	#ax.set_xticks(np.arange(0, max_cov + binwidth, num_bins*5.0))
	ax.set_xlabel('Est. cov.')
	ax.set_ylabel('Frequency')   
	ax.set_xlim(0, max_cov)
	ax.legend(loc='upper right')
	return ax

def plot_per_read_est_cov_vs_read_length(ax, data, s, output_prefix):
	print "[ Plotting per read est cov vs read length scatter plot ]"

	# get the specific s's data
	s_name = s
	s_color = data[s][0]
	sd = data[s][1]
	s_marker = data[s][2]

	# get x and y values from list of tuples
	y,x = zip(*sd)

	# setting limits
	mx = np.mean(x)
	my = np.mean(y)
	len_lim = mx*3
	cov_lim = my*3
	print len_lim
	print cov_lim

	# clors for heatmap
	# discrete color scheme

	# creating heatmap
	heatmap, xedges, yedges = np.histogram2d(x, y, bins=(20,20), range=[[0, len_lim], [0, cov_lim]])
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]	
	im = ax.imshow(heatmap.T, extent=extent, interpolation='nearest', origin='lower', aspect='auto')
#	ax.scatter(x, y, alpha=0.2, marker=s_marker, s=4, label=s_name,edgecolors=s_color, linewidth=0.5, facecolors=s_color, rasterized=True)
	
	# scatter plot with read length on x axis, and estimted cov for read on y axis
	ax.set_title('Est. cov vs read length (' + s + ')')
	ax.grid(True, linestyle='-', linewidth=0.3)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	ax.figure.colorbar(im, cax=cax)
	ax.set_xlabel('Read length (bps)')
	ax.set_ylabel('Estimated cov')
	ax.legend(loc='upper right')
	return ax

def plot_num_matching_residues(ax, data, output_prefix):
	print "[ Plotting number of matching residues ]"

	max_num_matching_residues = 0
	for s in data:
		sd = data[s][1]
		max_sd = max(sd)	
		if max_sd > max_num_matching_residues:
			max_num_matching_residues = max_sd

	binwidth = 100.0
	bins = np.arange(0, max_num_matching_residues + binwidth, binwidth)
	
	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		ax.hist(sd, bins=bins, alpha=0.5, label=s_name, color=s_color)	

	ax.set_title('Number of matching residues per overlap distributions')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xlabel('Number of matching residues')
	ax.set_ylabel('Frequency')
	return ax

def plot_per_read_GC_content(ax, data, output_prefix):
	print "[ Plotting GC content ]"
	binwidth = 1.0
	for s in data:
		per_read_GC_content = {}
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		s = sum(sd.values())
		for GC_content_level in sd:
			value = sd[GC_content_level]
			nvalue = float(value)/float(s)
			per_read_GC_content[int(float(GC_content_level))] = nvalue
		ax.plot(per_read_GC_content.keys(), per_read_GC_content.values(), color=s_color, label=s_name)

	ax.set_title('Per read GC content')
	ax.set_xlabel('% GC content')
	ax.set_ylabel('Proportion')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.legend(loc='upper right')
	return ax

def plot_est_genome_size(ax, data, output_prefix):
	print "[ Plotting genome size estimates ]"

	genome_sizes = list()
	s_names = list()
	colors = list()
	# now start plotting
	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		genome_size_in_megabases = round(float(sd)/float(1000000), 2)
		genome_sizes.append(genome_size_in_megabases)
		s_names.append(s_name)
		colors.append(str(s_color))

	# plotting the number of overlaps/read
	y_pos = np.arange(len(genome_sizes))
	ax.barh(y_pos, genome_sizes, align='center')
	
	# add the actual est genome size values to the bars
	rects = ax.patches
	if len(s_names) == 1 :
		height = 0.3
	else:
		height = 1.0/float(len(s_names))

	for rect, value, s in zip(rects, genome_sizes, s_names):
		rect.set_facecolor(data[s][0])
		rect.set_height(height)
		t = ax.text(0.05, rect.get_y() + rect.get_height()/2.0, s + ": " + str(value) + "Mbp", ha='left', va='center', fontsize='8')
		t.set_bbox(dict(facecolor='#FFFFFF', alpha=0.5, edgecolor='#FFFFFF'))

	ax.set_yticks([])
	ax.set_title('Est. genome size')
	ax.set_xlabel('Genome size (Mbp)')
	ax.set_ylabel(' \n \n ')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.legend(loc='upper right')
	return ax

def plot_mean_overlap_accuracies_per_read(ax, data, output_prefix):
	print "[ Plotting mean overlap accuracies per read vs read length ]"

	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		s_marker = data[s][2]
		# get x and y values from list of tuples
		y,x = zip(*sd)
		ax.scatter(x, y, alpha=0.1, marker=s_marker, s=4, edgecolors=s_color, linewidth=0.5, facecolors=s_color, rasterized=True)

	# scatter plot with read length on x axis, and estimted cov for read on y axis
	ax.set_title('Mean accuracy vs read length')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xlabel('Read length')
	ax.set_ylabel('Mean accuracy')
	ax.legend(loc='upper right')	
	return ax

def plot_total_num_bases_vs_min_read_length(ax, data, output_prefix):
	print "[ Plotting total number of bases as a function of minimum read length ]"

	max_min_read_length = 0
	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]		# dictionary: key = min read length cut off, value = total bases
		s_marker = data[s][2]
		# get x and y values from list of tuples
		x,y=zip(*sd)

		# let's normalize the y values
		ny = list()
		s = max(y)
		for i in y:
			ny.append(float(i)/float(s))	
		
		# record if current highest x-value
		mx = max(x)
		if mx > max_min_read_length:
			max_min_read_length = mx
		ax.plot(x, ny, label=s_name, color=s_color)

	# set x limit
	x_lim = max_min_read_length*0.90

	# scatter plot with read length on x axis, and estimted cov for read on y axis
	ax.set_title('Total number of bases')
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xlabel('Minimum read length (bp)')
	ax.set_ylabel('Proportion')
	ax.legend(loc='upper right')
	ax.set_xlim(0, x_lim)
	return ax

def plot_assembly_quality_metrics(ax, data, output_prefix, type):
	if type == "NGX":
		print "[ Plotting NGX values ]"
		ax.set_title('NG(X)')
	else:
		print "[ Plotting NX values ]"
		ax.set_title('N(X)')

	for s in data:
		s_name = s
		s_color = data[s][0]
		sd = data[s][1]
		s_marker = data[s][2]
		lists = list()
		for key in sd:
			lists.append((int(key), sd[key]))
		lists.sort(key=lambda x: x[0])
		x, y = zip(*lists)
		ax.plot(x, y, label=s_name, color=s_color)

	# scatter plot with read length on x axis, and estimted cov for read on y axis
	ax.grid(True, linestyle='-', linewidth=0.3)
	ax.set_xlabel('X (%)')
	ax.set_ylabel('Contig length (bp)')
	ax.legend(loc='upper right')
	return ax

if __name__ == "__main__":
	main()


