from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


blast_output = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(blast_output))
print(blast_result_record)