.. _installation.rst

Installation
================

Requirements
-------------

To generate files needed:
* minimap2 - required to create required input PAF file
* miniasm - required if NG(X) plots requested

For the calculation step:
* C++ compiler with C++11 support

For the report generation step:
* Python 2.7.11
* setuptools - required for installation
* BioPython
* matplotlib v2.0.0

Installing the latest code from github
----------------------------------------
::

    git clone https://github.com/simpsonlab/preqc-lr.git
    cd preqc-lr
    make

Installing dependencies
---------------------------------------

First we need to make sure we have everything to properly use pip or the setup.py script.

::

    # create virtual environment
    virtualenv preqclr-venv
    source preqclr-venv/bin/activate

    # check that you are using correct environment
    which pip

    # check that setuptools is installed
    pip list    

    # update setuptools if needed
    python -m pip install --upgrade pip setuptools

    # check that we are using the right version of python
    python -V

Okay, we are ready to install dependencies.

::   

    # download report script dependencies 
    python setup.py install

    # or we can use pip
    pip install preqc-lr

    pip list

    # update setuptools
    python -m pip install --upgrade pip setuptools wheel

    # download report script dependencies
    pip install preqc-lr

This will install the following packages with all their dependencies:

* matplotlib v2.0.0
* BioPython

