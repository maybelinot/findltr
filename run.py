__author__ = 'Eduard Trott'

from Bio import SeqIO

FILENAME = "Homo sapiens chromosome X genomic scaffold, GRCh38 Primary Assembly HSCHRX_CTG3.fasta"


class GenomeClass:
    """
    Class to represent genome
    """
    def __init__(self, filename):
        handle = open(filename, "rU")
        self.data = next(SeqIO.parse(handle, "fasta"))
        self.__x = 5
        self.x = 6

    def __str__(self):
        """
        Human readable representation of genome
        """

    def run(self):
        self.de_novo(self)

    @staticmethod
    def de_novo(self):
        """
        De novo identification of young intact LTR retroelements
        """
        # for one fasta sequence

        # de_novo_first_step(binary searching with LCP array)
        min_seq_len = 40
        min_distances = 1000
        max_distances = 20000
        print(len(self.data.seq))
        for idx in range(len(self.data.seq) - (min_seq_len + min_distances)):
            pass

        # de_novo_second_step

        # de_novo_last_step

genome = GenomeClass(FILENAME)
genome.run()
# print(genome.__x)
# self.data[0].seq

