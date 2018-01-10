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

    * READ file: long read sequences
    * PAF file: information on overlaps between reads in READ file
    * GFA file: graph assembly

Output
"""""""""""""""""""""""

    * JSON file containing data needed to generate plots in preqc-lr-report
    * Log file summarizing statistics calculated, input, and output
    * Directory of csv files optionally given that contains all calculations

Usage example
"""""""""""""""""""""""

::

   python preqc-lr-calculate.py [-h/--help] -r/--reads <fasta|fastq|fasta.gz|fastq.gz> \
           -t/--type {pb, ont} -n/--sample_name sample_name \
           -p/--paf <PAF> -g/--gfa <GFA> \
           --csv --verbose -v/--version  

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

   * - ``-t``, ``--type``
     - Y
     - NA
     - Type of sequencing performed to achieve reads file. Either pacbio (pb) or oxford nanopore technology data (ont).

   * - ``-n``, ``--sample_name``
     - Y
     - NA
     - Sample name; you can use the name of species for example. This will be used as output prefix.

   * - ``-p``, ``--paf``
     - N
     - NA
     - Minimap2 pairwise alignment file (PAF). This is produced using ``minimap2 -x ava-ont sample.fastq sample.fastq``. This is **REQUIRED** for analysis. However, if users do not pass a PAF file as input, preqc-lr will attempt to run minimap2. This requires minimap2 to be installed and found via PATH variable. We recommend running minimap2 on your own and passing the PAF file.

   * - ``-g``, ``--gfa``
     - N
     - NA
     - Miniasm graph gragment assembly (GFA) file. This is produced using ``miniasm -f reads.fasta overlaps.paf``. This is required only if user wants to generate an NGX plot. If not given, it will **NOT CALCULATE NGX STATISTICS**.

   * - ``--verbose``
     - N
     - False
     - Use to output preqc-lr prorgress to stdout.

Report
---------


Overview
"""""""""""""""""""""""

Generates a report with plots describing QC metrics for long read data sets.

Input
"""""""""""""""""""""""

    * JSON file(s) containing data for a sample needed to generate plots created by preqc-lr-calculate 

Output
"""""""""""""""""""""""

    * PDF file report

Plots:

1. Estimated genome size
   This is a bar plot that shows the estimated genome size for one or more samples. As coverage was inferred from overlap information, we can use this to calculate genome size with Lander-waterman statistics. 
2. Read length distribution
   This is the distribution of read lengths. preqc-lr imposes an x-limit such that 90% of all of the read lengths falls under this limit. This was done to avoid extremely long tails.
3. Estimated coverage distribution
   This shows the distribution of coverage for each read inferred from the overlap information file (PAF). 
4. Per read GC content distribution
   In this plot we show the distribution of GC content per read. To calculate this for each read, we summed the number of C and G nucleotides then divided by the read length.

Usage example
"""""""""""""""""""""""

::

   python preqc-lr-report.py [-h/--help] -i/--input <*.preqclr> \
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
     - Output of preqc-lr-calculate. JSON formatted file with '.preqclr' extension.

   * - ``-o``, ``--output``
     - N
     - If only one preqc-lr file given, it will infer from prefix. Else if multiple, prefix will be "preqc-lr-output".
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
