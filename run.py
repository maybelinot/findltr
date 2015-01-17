__author__ = 'Eduard Trott'

from Bio import SeqIO
import shelve, os, sys
import algorithms, arguments

class GenomeClass:
    """
    Class to represent genome
    """

    def __init__(self, args):
        handle = open(args['input_file'], "rU")
        self.data = next(SeqIO.parse(handle, "fasta"))
        self.algorithm = args['algorithm']
        self.outputfile = args['output_file']
        self.min_pattern_len = args['min_pattern_len']
        self.min_distance = args['min_distance']
        self.max_distance = args['max_distance']
        self.max_ltr_len = args['max_ltr_len']
        self.min_ltr_len = args['min_ltr_len']

    def __str__(self):
        """
        Human readable representation of genome
        """

    def run(self):
        self.de_novo()

    def de_novo(self):
        """
        De novo identification of young intact LTR retroelements
        """

        seq = str(self.data.seq).upper()

        # Pattern searching
        db = shelve.open('LCP.db', writeback=True)
        # db = [[start_of_pattern, start_of_appropriate_pattern], ...] where pattern has a length equal min_pattern_len
        # and distance between patterns is in range (min_distance : max_distance)
        db['young_lcp_parts'] = getattr(algorithms, self.algorithm)(seq = seq,\
                                                                    min_pattern_len = self.min_pattern_len,\
                                                                    min_distance = self.min_distance,\
                                                                    max_distance = self.max_distance)
        db.close()

        # de_novo_second_step
        # grouping of patterns to LTRs
        os.system("grouping.py %d %d %d %d" % (self.max_ltr_len, self.min_ltr_len, self.min_pattern_len, self.min_distance))

        # creation of sequences records
        algorithms.gff_writing(seq)



if __name__ == "__main__":
    genome = GenomeClass(arguments.arguments(sys.argv[1:]))
    genome.run()