__author__ = 'Eduard Trott'

from Bio import SeqIO

FILENAME = "Homo sapiens chromosome X genomic scaffold, GRCh38 Primary Assembly HSCHRX_CTG3.fasta"

def get_height(text, suffix_array):
    lcp = [0]*len(suffix_array)
    rank = [0]*len(suffix_array)
    n = len(suffix_array)
    for idx in range(n):
        rank[suffix_array[idx]] = idx
    h = 0
    for idx in range(n):
        k = rank[idx]
        if k == 0:
            lcp[k] = -1
        else:
            j = suffix_array[k-1]
            while idx+h < n and j+h < n and text[idx+h] == text[j+h]:
                h += 1
            lcp[k] = h
        if h > 0:
            h -= 1
    return lcp


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
        seq = '$'+str(self.data.seq[0:10])
        satupes = sorted([(seq[i:],i) for i in range(0,len(seq))])
        suffix_array = list(map(lambda x: x[1], satupes))
        lcp = get_height(seq, suffix_array)
        print(seq)
        print(suffix_array)
        print(lcp)
        for idx in range(len(self.data.seq) - (min_seq_len + min_distances)):
            pass

        # de_novo_second_step

        # de_novo_last_step

genome = GenomeClass(FILENAME)
genome.run()




# print(genome.__x)
# self.data[0].seq

