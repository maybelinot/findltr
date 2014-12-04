from __future__ import generators
__author__ = 'Eduard Trott'

from Bio import SeqIO, Seq, SeqRecord, SeqFeature
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import translate
from BCBio import GFF
import time
import shelve
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

FILENAME = "chrX.fa"


def get_lcp(text, suffix_array):
    lcp = [0] * len(suffix_array)
    rank = [0] * len(suffix_array)
    n = len(suffix_array)
    for idx in range(n):
        rank[suffix_array[idx]] = idx
    h = 0
    for idx in range(n):
        k = rank[idx]
        if k == 0:
            lcp[k] = -1
        else:
            j = suffix_array[k - 1]
            while idx + h < n and j + h < n and text[idx + h] == text[j + h]:
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
    return i + 1


def find_pattern(text, pattern, lcp, sa):
    left = 0
    right = len(text)
    middle = (right + left) // 2
    while True:
        if abs(right - left) <= 1:
            if pattern == text[sa[middle][1]:sa[middle][1] + len(pattern)]:
                return sa[middle][1], len(pattern)
            return None
        pattern_lcp = max([count_lcp(pattern, text[sa[idx][1]:]) for idx in range(left, middle)])
        # lcp comparison
        if pattern_lcp > lcp[middle]:
            right = middle
            middle = (left + right) // 2
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

# Knuth-Morris-Pratt string matching
# David Eppstein, UC Irvine, 1 Mar 2002


def KnuthMorrisPratt(text, pattern):

    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''

    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)

    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift

    # do the actual search
    startPos = 0
    matchLen = 0
    for c in text:
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield startPos

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
        self.de_novo()

    def de_novo(self):
        """
        De novo identification of young intact LTR retroelements
        """

        min_pattern_len = 40
        min_distance = 1000
        max_distance = 20000

        seq = str(self.data.seq).upper()
        # start_time = time.time()
        # output = []
        # idx = 0
        # while idx < len(seq) - (min_pattern_len * 2 + min_distance):
        #     pattern = seq[idx:idx + min_pattern_len]
        #     if not 'N' in pattern:
        #         text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
        #         # ans = [el for el in KnuthMorrisPratt(text, pattern)]
        #         # if ans:
        #         #     output.append([idx, ans[0] + idx + min_pattern_len + min_distance])
        #         #     idx += min_pattern_len
        #         if pattern in text:
        #             output.append([idx, text.index(pattern) + idx + min_pattern_len + min_distance])
        #             idx += min_pattern_len
        #     else:
        #         idx += min_pattern_len
        #     idx += 1
        # print("--- %s seconds ---" % (time.time() - start_time))
        # db = shelve.open('LCP.db', writeback=True)
        # # db = [[start_of_pattern, start_of_appropriate_pattern], ...] where pattern has a length equal min_pattern_len
        # # and distance between patterns is in range (min_distance : max_distance)
        # db['young_lcp_parts'] = output
        # print(output)
        # print(db['young_lcp_parts'])
        # db.close()

        # ###############################################################################################################
        # de_novo_first_step(binary searching with LCP array)
        # start_time = time.time()
        # output = []
        # for idx in range(len(seq) - (min_pattern_len*2 + min_distance)):
        #     pattern = seq[idx:idx + min_pattern_len]
        #     if idx + min_distance + min_pattern_len < min_distance:
        #         text = seq[idx + min_pattern_len + min_distance:]
        #         suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
        #         lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
        #         ans = find_pattern(text, pattern, lcp, suffix_array)
        #     else:
        #         text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + min_distance]
        #         suffix_array = sorted([(text[i:], i) for i in range(0, len(text))])
        #         lcp = get_lcp(text, list(map(lambda x: x[1], suffix_array)))
        #         ans = find_pattern(text, pattern, lcp, suffix_array)
        #     if None != ans:
        #         output.append([idx, ans])
        # print("--- %s seconds ---" % (time.time() - start_time))

        # db['time']['lcpsa'].append(time.time() - start_time)
        # db.close()
        # print(output)
        ################################################################################################################

        # de_novo_second_step
        # !!!!! Add screening by the trailing sequences of LTR
        # !!!!! Add 80% compare to forming LTR groups (search some biopython tool)

        # max_ltr_len = 1000
        # min_ltr_len = 100
        # db = shelve.open('LCP.db', writeback=True)
        #
        # groups_of_ltrs = [[[db['young_lcp_parts'][0][0], db['young_lcp_parts'][0][0] + min_pattern_len],
        #                    [db['young_lcp_parts'][0][1], db['young_lcp_parts'][0][1] + min_pattern_len]]]
        #
        # duplicates = False
        # for lcp_part in db['young_lcp_parts'][1:]:
        #     if lcp_part[0] + min_pattern_len + min_distance > groups_of_ltrs[-1][1][0]:
        #         if lcp_part[0] > groups_of_ltrs[-1][1][1]:
        #             if duplicates\
        #                     or (groups_of_ltrs[-1][0][1] - groups_of_ltrs[-1][0][0]) < min_ltr_len \
        #                     or (groups_of_ltrs[-1][1][1] - groups_of_ltrs[-1][1][0]) < min_ltr_len:
        #                 groups_of_ltrs[-1] = [[lcp_part[0], lcp_part[0] + min_pattern_len],
        #                                       [lcp_part[1], lcp_part[1] + min_pattern_len]]
        #                 duplicates = False
        #             else:
        #                 groups_of_ltrs.append([[lcp_part[0], lcp_part[0] + min_pattern_len],
        #                                        [lcp_part[1], lcp_part[1] + min_pattern_len]])
        #         else:
        #             duplicates = True
        #     elif (lcp_part[0] - groups_of_ltrs[-1][0][0] < max_ltr_len) or \
        #             (lcp_part[1] - groups_of_ltrs[-1][1][0] < max_ltr_len):
        #         groups_of_ltrs[-1][0][1] = lcp_part[0] + min_pattern_len
        #         groups_of_ltrs[-1][1][1] = lcp_part[1] + min_pattern_len
        #     else:
        #         duplicates = True
        # if duplicates or (groups_of_ltrs[-1][0][1] - groups_of_ltrs[-1][0][0]) < min_ltr_len \
        #         or min_ltr_len > (groups_of_ltrs[-1][1][1] - groups_of_ltrs[-1][1][0]):
        #     del groups_of_ltrs[-1]
        #
        # db['young_lcp'] = groups_of_ltrs
        # db.close()

        # creation of sequences records
        db = shelve.open('LCP.db', writeback=True)
        records = []
        # for idx, element in enumerate(db['young_lcp']):
        #     records.append(SeqRecord.SeqRecord(Seq(seq[element [0][1]:element[1][0]], generic_nucleotide).translate(),
        #                                        id=str(idx)))

        gff = SeqRecord.SeqRecord(Seq.Seq(seq), "chrX")
        top_feature = []
        for idx, item in enumerate(db['young_lcp']):
            print()
            print()
            seq1 = SeqRecord.SeqRecord(Seq.Seq(seq[item[0][0]:item[0][1]]), id="seq1")
            seq2 = SeqRecord.SeqRecord(Seq.Seq(seq[item[1][0]:item[1][1]]), id="seq2")
            SeqIO.write(seq1, "seq1.fasta", "fasta")
            SeqIO.write(seq2, "seq2.fasta", "fasta")

            blast_output = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]
            blast_result_record = NCBIXML.read(StringIO(blast_output))
            identity = 0
            for alignment in blast_result_record.alignments:
                for hsp in alignment.hsps:
                    identity = max(hsp.identities/hsp.align_length*100, identity)
            identity = "%.4f" % identity
            identity = identity.rstrip("0")
            identity = identity.rstrip(".")

            sub_qualifiers = {"source": "ltrfind", "ID": "UnknownLTR_"+str(idx+1), "Note": "identity "+identity}
            top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0], item[1][1]),
                            type="SO:0000186", strand=1, qualifiers=sub_qualifiers))
        gff.features = top_feature

        with open('rec.gff', "w") as out_handle:
            GFF.write([gff], out_handle)
        with open('rec.gff', "r+") as out_handle:
            file = out_handle.readlines()
            tmp_gff = [line for line in file if line[0] == '#']
            tmp_gff.extend([line[:-1]+" %\n" for line in file if line[0] != '#'])
        print(tmp_gff)
        with open('rec.gff', "w") as out_handle:
            out_handle.writelines(tmp_gff)

        # if the distance is less than 1000 then consider this seq as a duplicates
        # !!!! Add condition on LTR retroelements inside another LTRs

        # [print(el, '\n') for el in groups_of_ltrs]

        # de_novo_last_step


genome = GenomeClass(FILENAME)
genome.run()

