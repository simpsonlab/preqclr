.. _installation.rst

Installation
================

Requirements
-------------

To generate files needed:

- minimap2 (required to create required input PAF file)
- miniasm (required if NG(X) plots requested)

For the calculation step:

- C++ compiler with C++11 support

For the report generation step:

- Python 2.7.11
- setuptools - required for installation
- BioPython
- matplotlib v2.0.0

Installing the latest code from github
----------------------------------------
::

    git clone --recursive  https://github.com/simpsonlab/preqc-lr.git
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
    pip freeze    

    # update setuptools if needed
    python -m pip install --upgrade pip setuptools

    # check that we are using the right version of python (2.7.11+)
    python -V

Okay, we are ready to install dependencies.

::   

    # download report script dependencies 
    python setup.py install

    # OR we can use pip
    pip install preqc-lr

To check that we have installed all the packages and the right versions we run `pip freeze` and we should see the following:

::

    biopython==1.70
    cycler==0.10.0
    functools32==3.2.3.post2
    gevent==1.3a1
    greenlet==0.4.13
    matplotlib==2.0.0
    numpy==1.14.0
    preqc-lr==2.0
    pyparsing==2.2.0
    python-dateutil==2.6.1
    pytz==2017.3
    six==1.11.0
    subprocess32==3.5.0rc1 


