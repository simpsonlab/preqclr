.. _manual

Manual
==================

Calculate
------------


Overview
"""""""""""""""""""""""

Generates data needed to create plots in preqc-lr-report.

Input
"""""""""""""""""""""""

    * READS file: long read sequences in fasta or fastq format
    * PAF file: information on overlaps between reads in READS file
    * GFA file: graph assembly that contains contig information

Output
"""""""""""""""""""""""

    * JSON file containing data needed to generate plots in preqclr-report
    * Log file summarizing statistics calculated, input, and output

Usage example
"""""""""""""""""""""""

::

   ./preqclr [-h/--help] -r/--reads <fasta|fastq|fasta.gz|fastq.gz> \
           -n/--sample_name sample_name \
           -p/--paf <PAF> -g/--gfa <GFA> \
           --rlen_cutoff INT \
           --verbose -v/--version  

.. list-table:: 
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``-r``, ``--reads``
     - Y
     - NA
     - Fasta, fastq, fasta.gz, or fastq.gz files containing reads.

   * - ``-n``, ``--sample_name``
     - Y
     - NA
     - Sample name; you can use the name of species for example. This will be used as output prefix.

   * - ``-p``, ``--paf``
     - N
     - NA
     - Minimap2 Pairwise mApping Format (PAF) file. This is produced using ``minimap2 -x ava-ont sample.fastq sample.fasta``.

   * - ``-g``, ``--gfa``
     - N
     - NA
     - Miniasm Graph Fragment Assembly (GFA) file. This is produced using ``miniasm -f reads.fasta overlaps.paf > layout.gfa``. This is required only if user wants to generate an NGX plot. If not given, it will **NOT CALCULATE NGX STATISTICS**.

   * - ``--verbose``
     - N
     - False
     - Use to output preqc-lr progress to stdout.

Report
---------


Overview
"""""""""""""""""""""""

Generates a report with plots describing QC metrics for long read data sets.

Input
"""""""""""""""""""""""

    * JSON file(s) containing data for sample(s) needed to generate plots created in preqclr calculate 

Output
"""""""""""""""""""""""

    * PDF file report

Plots:

1. Estimated genome size
   This is a bar plot that shows the estimated genome size for one or more samples. As coverage was inferred from overlap information, we can use this to calculate genome size with Lander-Waterman statistics. 
2. Read length distribution
   This is the distribution of read lengths calculated from the READS file. preqclr imposes an x-limit such that 90% of all of the read lengths falls under this limit. This was done to avoid extremely long tails.
3. Estimated coverage distribution
   This shows the distribution of coverage for each read inferred from the overlap information file (PAF). 
4. Per read GC content distribution
   In this plot we show the distribution of GC content per read for a sample of 40% of reads. To calculate this for each read, we summed the number of C and G nucleotides then divided by the read length.
5. Total number of bases vs minimum read length
   We show the total number of bases with reads of a minimum length of x.
6. NGX
   This shows the contigiuity of the data. Miniasm produces contigs from your sequencing data. To interpret this let's look at x=50 and it's NG(50) value on the y-axis. The contig length on the y-axis describes the length at which 50% of the genome size estimate is capture in contigs with length greater than or equal to the NG(50) value.


Usage example
"""""""""""""""""""""""

::

   python preqclr-report.py [-h/--help] -i/--input <*.preqclr> \
        --save_png --list_plots -o/--output <output_prefix> --plot <list of user specified plots> \
        --verbose 

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``-i``, ``--input``
     - Y
     - NA
     - Output of preqclr calculate step. JSON formatted file with '.preqclr' extension.

   * - ``-o``, ``--output``
     - N
     - If only one preqclr file given, it will infer from prefix. Else if multiple, prefix will be "preqc-lr-output".
     - Prefix for output PDF.

   * - ``--plot``
     - N
     - NA
     - Users can specify which plots they want. To do so, use ``--list_plots`` and use the names of plots.

   * - ``--list_plots``
     - N
     - NA
     - Use this to see which plots are available. Note that NGX plots are also dependent on whether or not it was calculated in preqc-lr-calculate step and this depends on whether or not miniasm's GFA file was passed as input.

   * - ``--save_png``
     - N
     - False
     - Use to save each subplot as a png.

   * - ``--verbose``
     - N
     - False
     - Use to print progress to stdout.
