# preqc-lr
preqc-lr is a software tool that reports on quality control metrics for long read sequencing data without the use of a reference genome.

## Dependencies

* Python2.7
* matplotlib
* BioPython

## Install
To install from source:

```bash
git clone https://github.com/simpsonlab/preqc-lr.git
cd preqc-lr
python setup.py install
```

## What is preqc-lr?

With the emergence of new long read sequencing technology such as Pacbio Single Molecule, Real-Time (SMRT) Sequencing technology and Oxford Nanopore Technologies (ONT), there is a need for a method that assesses sequencing quality prior to analyses. With preqc-lr we enable users to visualize metrics of quality.

There are two components to preqc-lr:

    1. calculate
    2. report

## Learn

* Documentation [here](http://preqc-lr.readthedocs.io/en/latest/)
* Quickstart tutorial [here](http://preqc-lr.readthedocs.io/en/latest/quickstart.html)
