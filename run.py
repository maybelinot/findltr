__author__ = 'Eduard Trott'

from Bio import SeqIO, Seq, SeqRecord, SeqFeature
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import translate
from BCBio import GFF
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import time, shelve, sys, getopt, os
from io import StringIO
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
        db = shelve.open('LCP.db', writeback=True)
        if db['young_lcp']==[]:
            print('LTRs not found')
            sys.exit(2)
        else:
            unique_name = 'rec_%s.gff' % time.time()
            print('Found LTRs are saved in ' + unique_name)

        records = []

        gff = SeqRecord.SeqRecord(Seq.Seq(seq), "chrX")
        top_feature = []
        for idx, item in enumerate(db['young_lcp']):
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
        db.close()
        gff.features = top_feature
        with open(unique_name, "w") as out_handle:
            GFF.write([gff], out_handle)
        with open(unique_name, "r+") as out_handle:
            file = out_handle.readlines()
            tmp_gff = [line for line in file if line[0] == '#']
            tmp_gff.extend([line[:-1]+" %\n" for line in file if line[0] != '#'])
        with open(unique_name, "w") as out_handle:
            out_handle.writelines(tmp_gff)


if __name__ == "__main__":
    genome = GenomeClass(arguments.arguments(sys.argv[1:]))
    genome.run()