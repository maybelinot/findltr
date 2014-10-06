__author__ = 'Eduard Trott'

from Bio import SeqIO
import time

FILENAME = "Homo sapiens chromosome X genomic scaffold, GRCh38 Primary Assembly HSCHRX_CTG3.fasta"

def get_lcp(text, suffix_array):
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


def count_lcp(text, pattern):
    min_len = min(len(text), len(pattern))
    if min_len < 1:
        return 0
    for i in range(min_len):
        if text[i] != pattern[i]:
            i -= 1
            break
    return i+1

print(count_lcp('', 'a'))


def find_pattern(t,p,lcp):
    L = 0
    R = len(t)
    M = len(t)//2

    pass


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

        # seq = str(self.data.seq[0:10]) + '$'
        seq = 'parampampam$'
        start_time = time.time()
        satupes = sorted([(seq[i:], i) for i in range(0, len(seq))])
        suffix_array = list(map(lambda x: x[1], satupes))
        lcp = get_lcp(seq, suffix_array)
        print(seq)
        print(suffix_array)
        print(lcp)
        print("--- %s seconds ---" % (time.time() - start_time))
        pattern = 'pam'
        L = 0
        # for idx in range(len(self.data.seq) - (min_seq_len + min_distances)):
        #     pass

        # de_novo_second_step

        # de_novo_last_step

genome = GenomeClass(FILENAME)
genome.run()




# print(genome.__x)
# self.data[0].seq

