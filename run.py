__author__ = 'Eduard Trott'

from Bio import SeqIO

FILENAME = "Homo sapiens chromosome X genomic scaffold, GRCh38 Primary Assembly HSCHRX_CTG3.fasta"


class GenomeClass:
    """
    Class to represent genome
    """
    def __init__(self, filename):
        handle = open(filename, "rU")
        self.data = list(SeqIO.parse(handle, "fasta"))

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
        def de_novo_first_step(name):
            """

            """
            print(name)
        de_novo_first_step('Vasya')

genome = GenomeClass(FILENAME)
genome.de_novo(genome)
# self.data[0].seq

