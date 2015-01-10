VUT_project
===========

De novo identification of LTR retrotransposons

==========
programs/toolboxes

Python 3.4
Biopython 1.64 [https://pypi.python.org/pypi/biopython/1.64]
bcbio-gff      [https://github.com/chapmanb/bcbb/tree/master/gff]
blastn(add to path): 2.2.30+ (Package: blast 2.2.30, build Oct 27 2014 16:53:52) [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST]

==========
To run the program, type in the console(run from current folder):
python run.py -i <inputfile> -a <algorithm> <parametrs>

-i input file (full path without spaces or .\input.fasta if file is located in current directory)
-a one of 3 algorithms: kmp (KnuthMorrisPratt), salcp (Suffix array with long common prefix), original. By default original.

Parametrs:
--min_pattern_len Minimal length of pattern for pattern searching algorithms (by default 40),
--min_distance    Minimal distance between LTRs (by default 1000),
--max_distance    Maximal distance between LTRs (by default 20000),
--max_ltr_len     Maximal length of one LTR (by default 1000),
--min_ltr_len 	  Minimal length of one LTR (by default 100).

Example:

python run.py -a original -i .\example.fasta --min_pattern_len 50 --max_ltr_len 1500