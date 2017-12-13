.. _quickstart:

Quickstart
============

**Time:** 10 minutes 

preqc-lr generates a PDF report containing several plots such as estimated genome size and coverage. This report can be used to evaluate the quality of your sequencing data. Here, we provide a step-by-step tutorial to get you started!

**Requirements:**

* `preqc-lr <https://github.com/simpsonlab/preqc-lr>`_ 
* `minimap2 <https://github.com/lh3/minimap2>`_
* `miniasm <https://github.com/lh3/miniasm/>`_

Download example dataset
""""""""""""""""""""""""""

You can download the example dataset we will use here: ::

    wget http://s3.climb.ac.uk/nanopolish_tutorial/ecoli_100kb_region.tar.gz
    tar -xf ecoli_100kb_region.tar.gz
	cd ecoli_100kb_region/

**Details:**

This is a subset of reads that align to a 100kb region of *E. coli*.  The data was produced using Oxford Nanopore Technologies (ONT) MinION sequencer.

* Sample :    E. coli str. K-12 substr. MG1655
* Instrument : ONT MinION sequencing R9.4 chemistry
* Basecaller : Albacore v2.0.1
* Region length: Approx. 100kb
* Number of reads: 1581

Generate overlap information with minimap2
""""""""""""""""""""""""""""""""""""""""""""""""

We use minimap2 to find overlaps between our ONT long reads: ::

   minimap2 -x ava-ont reads.fasta reads.fasta > overlaps.paf 

If we take a peak at the first few lines of the Pairwise mApping Format (PAF) file, we see the following: ::

    00ffc7ec-ce0c-472d-a1f4-9a9ab65353fd	6379	78	6343	-	93d821e7-f696-4fab-b727-e18d4f179747	10079	3330	9672	2988	6467	0	tp:A:S	cm:i:509	s1:i:2946
    00ffc7ec-ce0c-472d-a1f4-9a9ab65353fd	6379	132	6320	+	b3f12a71-53ed-4245-9f72-67593fd00ec9	10347	2847	9059	2950	6352	0	tp:A:S	cm:i:539	s1:i:2912
    00ffc7ec-ce0c-472d-a1f4-9a9ab65353fd	6379	230	6343	+	044c27f2-7c69-4485-bea3-a5c75b7459f7	8566	70	6195	2643	6274	0	tp:A:S	cm:i:456	s1:i:2597
    00ffc7ec-ce0c-472d-a1f4-9a9ab65353fd	6379	127	6343	-	5bba56d6-dcd3-4e5a-8705-ffa9fe5049b2	8434	1222	7417	2442	6373	0	tp:A:S	cm:i:374	s1:i:2395
    00ffc7ec-ce0c-472d-a1f4-9a9ab65353fd	6379	102	6343	-	083d17cc-7c65-4404-85d2-baa6a36b4183	9218	1062	7319	2372	6407	0	tp:A:S	cm:i:376	s1:i:2329

You can find more information about the format of the PAF file `here <https://github.com/lh3/miniasm/blob/master/PAF.md>`_.

Generate assembly graph with miniasm
"""""""""""""""""""""""""""""""""""""""""""""""""

We use miniasm to get an assembly graph in the `Graphical Fragment Assembly <https://github.com/GFA-spec/GFA-spec/blob/master/GFA-spec.md>`_ format: ::

   miniasm -f reads.fasta overlaps.paf > layout.gfa

Perform calculations
""""""""""""""""""""""""

We now have the necessary files to run preqc-lr (``reads.fasta``, ``overlaps.paf``, and ``layout.gfa``). 
To generate the data needed for the reports we first run preqc-lr-calculate ::

    python /path/to/preqc-lr-calculate.py \
        --reads reads.fasta \
        --type ont \
        --sample_name ecoli_100kb.ONT \
        --paf overlaps.paf \
        --gfa layout.gfa

This will produce a JSON formatted file: ``ecoli_100kb.ONT.preqclr``.

Generate report
"""""""""""""""""""

Now we are read to run preqc-lr-report to generate a PDF file describing quality metrics of your sequencing data: ::

    python /path/to/preqc-lr-report.py \
        -i ecoli_100kb.ONT.preqclr \
        -o ecoli_100kb.ONT

This will produce a PDF file: ``ecoli.ONT.pdf``.

Example report
"""""""""""""""""""

The report produces plots as seen below.

**Plot 0:**

.. figure:: _static/plot_est_genome_size.png
  :scale: 80%
  :alt: plot_est_genome_size

**Plot 1:**

.. figure:: _static/plot_read_length_distribution.png
  :scale: 80%
  :alt: plot_read_length_distribution

**Plot 2:**

.. figure:: _static/plot_est_cov.png
  :scale: 80%
  :alt: plot_est_cov

**Plot 3:**

.. figure:: _static/plot_per_read_GC_content.png
  :scale: 80%
  :alt: plot_per_read_GC_content

**Plot 4:**

.. figure:: _static/plot_est_cov_vs_read_length.png
  :scale: 80%
  :alt: plot_est_cov_vs_read_length

**Plot 5:**

.. figure:: _static/plot_total_num_bases.png
  :scale: 80%
  :alt: plot_total_num_bases

**Plot 6:**

.. figure:: _static/plot_NGX.png
  :scale: 80%
  :alt: plot_NGX.png
