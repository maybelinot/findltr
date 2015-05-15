__author__ = 'Eduard Trott'

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
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

        self.pattern_search = True
        self.grouping = True
        self.gff_writing = True

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

        seq = str(self.data.seq).upper()[16191939:16211070]

        # Pattern searching
        db = shelve.open('LCP.db', writeback=True)
        # db = [[start_of_pattern, start_of_appropriate_pattern], ...] where pattern has a length equal min_pattern_len
        # and distance between patterns is in range (min_distance : max_distance)
        if self.pattern_search:
            db['young_lcp_parts'] = getattr(algorithms, self.algorithm)(seq = seq,\
                                                                        min_pattern_len = self.min_pattern_len,\
                                                                        min_distance = self.min_distance,\
                                                                        max_distance = self.max_distance)
        db.close()

        # de_novo_second_step
        # grouping of patterns to LTRs
        if self.grouping:
            os.system("grouping.py %d %d %d %d" % (self.max_ltr_len, self.min_ltr_len, self.min_pattern_len, self.min_distance))
            
        # db['young_lcp'] - grouped LTR retroelements
        # db = shelve.open('LCP.db', writeback=True)
        db = shelve.open('LCP.db', writeback=True)
        print(db['young_lcp'])
        for idx in range(6):
            seqq = SeqRecord.SeqRecord(Seq(seq[db['young_lcp'][0][0][0]+idx:db['young_lcp'][0][1][1]]).translate(), id="seq")
            with open("example"+str(idx)+".fasta", "w") as output_handle:
                SeqIO.write(seqq, output_handle, "fasta")
        db.close()
                # os.system("C:\Python27\python.exe .\Other\hmmsearch-master\hmmer.py -f %s -d %s -l" % ("example"+str(idx)+".fasta", "./pfam/PF03732_fs.hmm"))
        # separation to certain families:
        ##############################################
        # Ty1/Copia

        # LTR - 100:1300 nt
        # non-coding region
        # PBS - 18 nt
        # gag
        # pol (PR-RT-RH-INT)
        # PPT A/G 10 nt

        ##############################################
        # Ty3/Gypsy

        # LTR - 100:1300 nt 
        # non-coding region
        # PBS - 18 nt
        # gag
        # pol ()
        # PPT A/G 10 nt

        ##############################################        
        # Bel/Pao

        # LTR - 100:900 nt
        # non-coding region
        # PBS - 18 nt
        # gag
        # pol (PR-RT-RH-INT)
        # env-like 4000-10000 nt
        # PPT A/G 10 nt

        ##############################################
        # Retroviridae

        # LTR - 200:1450 nt
        # non-coding region
        # PBS - 18 nt
        # pol (PR-RT-RH-INT)
        # env-like 
        # PPT A/G 10 nt

        ##############################################

        # creation of sequences records
        if self.gff_writing:
            algorithms.gff_writing(seq)



if __name__ == "__main__":
    genome = GenomeClass(arguments.arguments(sys.argv[1:]))
    genome.run()