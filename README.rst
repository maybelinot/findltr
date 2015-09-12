======================
    findltr
======================

Command-line tool for de novo identification of LTR retrotransposons.

Requirements
============

+ Basic Local Alignment Search Tool ([BLAST](http://www.ncbi.nlm.nih.gov/books/NBK279690/))::
    
    $ sudo apt-get install ncbi-blast+


Installations
=============
Example install, using VirtualEnv::

    # install/use python virtual environment
    virtualenv ~/virtenv_scratch --no-site-packages

    # activate the virtual environment
    source ~/virtenv_scratch/bin/activate

    # upgrade pip in the new virtenv
    pip install -U pip setuptools

    # install this package in DEVELOPMENT mode
    python setup.py develop

    # or simply install
    # python setup.py install


Options
~~~~~~~

The list of available options:

General
-------

--input
    Input **\*.fasta** or **\*.fa** file

--output
    Output **\*.gff** file

Search parameters
-----------------

--algorithm
    Algorithm, [kmp|salcp|original] (default: original)

--pattern_len
   Length of pattern for searching algorithms (default: 20)

--min_distance
    Minimal distance between LTRs (default: 1000)

--max_distance
    Maximal distance between LTRs (default: 20000)

--min_ltr_len
    Minimal length of LTRs (default: 100)

--max_ltr_len
    Maximal length of LTRs (default: 1000)

Other
-----

--verbose
    Include more details (like modified git directories)
--debug
    Turn on debugging output, do not catch exceptions

See ``findltr --help`` for complete list of available options.
