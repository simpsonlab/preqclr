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
