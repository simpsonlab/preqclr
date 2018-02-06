# preqc-lr
preqc-lr is a software tool that reports on quality control metrics for long read sequencing data without the use of a reference genome.

## Dependencies

For the calculation step:
* C++ compiler with C++11 support

For the report generation step:
* Python2.7
* matplotlib
* BioPython
* setuptools (to download report script dependencies)

## Install
To install from source:

```bash
git clone https://github.com/simpsonlab/preqc-lr.git
cd preqc-lr
make

# download report script dependencies
# create virtual environment
virtualenv preqc-lr-venv
source preqc-lr-venv/bin/activate
python setup.py install
```

## What is preqc-lr?

With the emergence of new long read sequencing technology such as Pacbio Single Molecule, Real-Time (SMRT) Sequencing technology and Oxford Nanopore Technologies (ONT), there is a need for a method that assesses sequencing quality prior to analyses. With preqc-lr we enable users to visualize metrics of quality.

There are two components to preqc-lr:

    1. calculate
    2. report

## Quick usage

```bash
#STEP 1: calculate data for plots
./preqclr 	-r reads.fq \
			-p overlaps.paf \
			-t pb \
			-n ecoli_sample.pacbio \
			--verbose

#STEP 2: create a PDF report
python preqc-lr-report.py -i ecoli_sample.pacbio.preqclr --verbose 
```

## Learn

* Documentation [here](http://preqc-lr.readthedocs.io/en/latest/)
* Quickstart tutorial [here](http://preqc-lr.readthedocs.io/en/latest/quickstart.html)

## Acknowledgements

preqc-lr uses overlap and assembly information from minimap2 and miniasm, respectively. To parse the output PAF file, and efficiently read fasta files we used kseq and PAF parser in miniasm. Therefore, I would like to thank Heng Li for developing these tools.
