import sys
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

seq1 = sys.argv[1]
seq2 = sys.argv[2]
out = sys.argv[3]

alignment = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -1)
seq1 = alignment[0].seqA
seq2 = alignment[0].seqB

with open(out, "w") as f:
    f.write(seq1 + "\n")
    f.write(seq2 + "\n")