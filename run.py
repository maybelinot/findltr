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


def find_pattern(text, pattern, lcp, sa):
    left = 0
    right = len(text)
    middle = (right + left) // 2
    # print('left middle right\n', left, middle, right)
    while True:
        if abs(right - left) <= 1:
            if pattern == text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                return sa[middle][1], len(pattern)
            return None
        # print(sa)
        pattern_lcp = max([count_lcp(pattern, text[sa[idx][1]:]) for idx in range(left, middle)])
        # print('lcp', lcp)
        # print('sa', sa)
        # print(pattern_lcp, lcp[middle])
        if pattern_lcp > lcp[middle]:
            right = middle
            middle = (left + right) // 2
            # print(left, middle, right)
        elif pattern_lcp < lcp[middle]:
            left = middle
            middle = (left + right) // 2
        elif len(pattern) <= len(text[sa[middle][1]:]):
            if pattern > text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                left = middle
                middle = (left + right) // 2
            elif pattern < text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                right = middle
                middle = (left + right) // 2
            else:
                return sa[middle][1], len(pattern)
        else:
            if pattern[:len(text[sa[middle][1]:])] > text[sa[middle][1]:]:
                left = middle
                middle = (left + right) // 2
            elif pattern[:len(text[sa[middle][1]:])] < text[sa[middle][1]:]:
                right = middle
                middle = (left + right) // 2
            else:
                return None


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
        min_pattern_len = 40
        min_distances = 50 # 1000
        max_distances = 500 # 20000

        seq = str(self.data.seq[0:1000])
        # seq = 'kakoe-to soobcheniekakoe-to s$'
        start_time = time.time()
        # print(seq)
        # print(suffix_array)
        # print(lcp)
        output = []
        for idx in range(len(seq) - (min_pattern_len*2 + min_distances)):
            pattern = seq[idx:idx + min_pattern_len]
            text = seq[idx + min_pattern_len + min_distances:]
            if pattern in text:
                ans = text.index(pattern)
            else:
                ans = None
            if None != ans:
                output.append([idx, ans])
            print(idx, ' из ', len(seq) - (min_pattern_len*2 + min_distances))
        print("--- %s seconds ---" % (time.time() - start_time))
        print(output)
        output = []
        for idx in range(len(seq) - (min_pattern_len*2 + min_distances)):
            pattern = seq[idx:idx + min_pattern_len]
            if idx + min_distances + min_pattern_len < max_distances:
                text = seq[idx + min_pattern_len + min_distances:]
                suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
                lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
                ans = find_pattern(text, pattern, lcp, suffix_array)
            else:
                text = seq[idx + min_pattern_len + min_distances:idx + min_pattern_len + min_distances + max_distances]
                suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
                lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
                ans = find_pattern(text, pattern, lcp, suffix_array)
            if None != ans:
                output.append([idx, ans])
            print(idx, ' из ', len(text) - (min_pattern_len*2 + min_distances))
        print("--- %s seconds ---" % (time.time() - start_time))
        print(output)
        # de_novo_second_step

        # de_novo_last_step

genome = GenomeClass(FILENAME)
genome.run()

# print(genome.__x)
# self.data[0].seq

