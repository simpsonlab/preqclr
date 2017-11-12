.. preqc-lr documentation master file, created by
   sphinx-quickstart on Tue Nov  7 14:20:14 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

preqc-lr
====================================
preqc-lr is a software tool that performs quality control checks to long read sequencing data. 

Install
====================================
To install the latest version: ::

   git clone https://github.com/simpsonlab/preqc-lr.git

What is preqc-lr?
====================================
With the emergence of new long read sequencing technology such as Pacbio Single Molecule, Real-Time (SMRT) Sequencing technology and Oxford Nanopore Technologies (ONT), there is a need for a method that assesses sequencing quality prior to analyses. With preqc-lr we enable users to visualize metrics of quality.

There are two components to preqc-lr:
1. calculate
2. report

Calculate
----------------
The first tool will calculate all the datasets needed to create plots using overlap information provided by minimap2. 

Report
----------------
The second tool reads the calculated output and generates a pdf with the following plots:
	1. Estimated genome size
	2. Read length distribution
	3. Estimated coverage distribution
	4. Per read GC content distribution
	5. Estimated coverage vs read length scatter plot 
	6. Total number of bases as a function of minimum read length

Quickstart
====================================
First run preqc-lr-calculate ::

    python preqc-lr-calculate.py \
        -i [FASTA/Q] \
        -t {pb, ont} \
        -n [sample_name] \
        --csv --verbose \
        --paf [paf file]

This will produce a json format file: ``[sample_name].preqclr``

Then run preqc-lr-report: ::

	python preqc-lr-report.py \
		-i [sample_name].preqclr \
		-o [output prefix] \
		--list_plots \
		--save_png

This will produce a pdf file: ``[output_prefix].pdf``


Command reference
====================================


.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
