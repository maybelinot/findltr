import time
import shelve

class GenomeClass:
    """
    Class to represent genome
    """

    def __init__(self):
        db = shelve.open('SEQ.db',writeback = True)
        self.data = db['seq']
        db.close()

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

        seq = self.data
        start_time = time.time()
        output = []
        idx = 0
        while idx < len(seq) - (min_pattern_len * 2 + min_distance):
            pattern = seq[idx:idx + min_pattern_len]
            if not 'N' in pattern:
                text = seq[idx + min_pattern_len + min_distance:idx + min_pattern_len + min_distance + max_distance]
                if pattern in text:
                    output.append([idx, text.index(pattern) + idx + min_pattern_len + min_distance])
                    idx += min_pattern_len
            else:
                idx += min_pattern_len
            idx += 1
        print ("--- %s seconds ---" % (time.time() - start_time))
        db = shelve.open('LCP.db', writeback=True)
        # db = [[start_of_pattern, start_of_appropriate_pattern], ...] where pattern has a length equal min_pattern_len
        # and distance between patterns is in range (min_distance : max_distance)
        db['young_lcp_parts'] = output
        # print(output)
        # print(db['young_lcp_parts'])
        db.close()

genome = GenomeClass()
genome.run()