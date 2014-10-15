__author__ = 'Eduard Trott'

from Bio import SeqIO
import time
import shelve

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
        # lcp comparison
        if pattern_lcp > lcp[middle]:
            right = middle
            middle = (left + right) // 2
            # print(left, middle, right)
        elif pattern_lcp < lcp[middle]:
            left = middle
            middle = (left + right) // 2
        # pattern comparison
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
        # db = shelve.open('LCP.db', writeback=True)
        # seq = str(self.data.seq)
        # [print(seq[db['young_lcp'][idx][0]:db['young_lcp'][idx][0]+40], seq[db['young_lcp'][idx][1]:db['young_lcp'][idx][1]+40]) for idx in range(15)]
        # print(len(db['young_lcp']))
        # db.close()

    @staticmethod
    def de_novo(self):
        """
        De novo identification of young intact LTR retroelements
        """
        # for one fasta sequence

        # de_novo_first_step(binary searching with LCP array)
        min_pattern_len = 40
        min_distance = 1000
        max_distance = 20000

        seq = str(self.data.seq)
        # start_time = time.time()
        # output = []
        # idx = 0
        # while idx < len(seq) - (min_pattern_len*2 + min_distance):
        #     pattern = seq[idx:idx + min_pattern_len]
        #     if not 'N' in pattern:
        #         text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
        #         ans = None
        #         if pattern in text:
        #             idx += min_pattern_len
        #             output.append([idx, text.index(pattern) + idx + min_pattern_len + min_distance])
        #     idx += 1
        # print("--- %s seconds ---" % (time.time() - start_time))
        # db = shelve.open('LCP.db', writeback=True)
        # # db = [[start_of_pattern, start_of_appropriate_pattern], ...] where pattern has a length equal min_pattern_len
        # # and distance between patterns is in range (min_distance : max_distance)
        # db['young_lcp_parts'] = output
        # db.close()

        ################################################################################################################
        # output = []
        # for idx in range(len(seq) - (min_pattern_len*2 + min_distances)):
        #     pattern = seq[idx:idx + min_pattern_len]
        #     if idx + min_distances + min_pattern_len < max_distances:
        #         text = seq[idx + min_pattern_len + min_distances:]
        #         suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
        #         lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
        #         ans = find_pattern(text, pattern, lcp, suffix_array)
        #     else:
        #         text = seq[idx + min_pattern_len + min_distances:idx + min_pattern_len + min_distances + max_distances]
        #         suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
        #         lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
        #         ans = find_pattern(text, pattern, lcp, suffix_array)
        #     if None != ans:
        #         output.append([idx, ans])
        #     print(idx, ' из ', len(text) - (min_pattern_len*2 + min_distances))
        # print("--- %s seconds ---" % (time.time() - start_time))
        # print(output)
        ################################################################################################################
        # de_novo_second_step
        max_ltr_len = 1000
        db = shelve.open('LCP.db', writeback=True)

        groups_of_ltrs = [[[db['young_lcp_parts'][0][0], db['young_lcp_parts'][0][0] + min_pattern_len], [db['young_lcp_parts'][0][1], db['young_lcp_parts'][0][1] + min_pattern_len]]]

        for lcp_part in db['young_lcp_parts'][1:]:
            if lcp_part[0] + min_distance < groups_of_ltrs[-1][1][0] and ((lcp_part[0] - groups_of_ltrs[-1][0][1] < max_ltr_len) or (lcp_part[1] - groups_of_ltrs[-1][1][1] < max_ltr_len)):
                groups_of_ltrs[-1][0][1] = lcp_part[0] + min_pattern_len
                groups_of_ltrs[-1][1][1] = lcp_part[1] + min_pattern_len
            else:
                groups_of_ltrs.append([[lcp_part[0], lcp_part[0] + min_pattern_len], [lcp_part[1], lcp_part[1] + min_pattern_len]])
        db['young_lcp'] = groups_of_ltrs
        db.close()
        # print(seq[groups_of_ltrs[1][0][0]:groups_of_ltrs[1][0][1]], '\n', seq[groups_of_ltrs[1][1][0]:groups_of_ltrs[1][1][1]])

        # if the distance is less than 1000 then consider this as a duplicates?


        # de_novo_last_step

genome = GenomeClass(FILENAME)
genome.run()

