#!/usr/bin/python3

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

def conservation( s1,s2):
    c = ''
    for a,b in zip(s1,s2):
        if a == b:
            c += '*'
        else:
            c += ' '
    return c

def write_clustal( s1, s2, filename):
    c = filename.split('/')[-1][:-4].split('_')
    a = c[0]
    for i in range(1,len(c)-1):
        a += '_' + c[i]
    b = c[-1]
    m = max( len(a), len(b))
    a = a.ljust(m)
    b = b.ljust(m)
    e = ''.ljust(m)
    cons = conservation( s1,s2)
    with open( filename, 'w') as w:
        w.write( 'CLUSTAL W formatted output, created by MutationExplorer\n\n')
        count = 0
        step = 60
        while count < len(s1):
            w.write( a + '\t' + s1[count:count+step] + '\n')
            w.write( b + '\t' + s2[count:count+step] + '\n')
            w.write( e + '\t' + cons[count:count+step] + '\n\n')
            count += step

            
def read_fasta( filename):
    seq = ''
    head = ''
    with open( filename) as r:
        l = r.readline()
        if l[0] != '>':
            print( "WARNING!:",filename, 'does not start with ">" in first line! First line is ignored!!!' )
        head = l[1:].strip()
        for l in r:
            if l[0] == '>':
                return head,seq
            seq += l.strip()
    return head,seq
            


def align( seq1, seq2):
    return pairwise2.align.globalds(seq1, seq2, blosum62, -10, -1)



if __name__ == "__main__":
    import sys

    h1,seq1 = read_fasta( sys.argv[1] )
    h2,seq2 = read_fasta( sys.argv[2] )
    out  = sys.argv[3]

    alignment = align( seq1,seq2)

    seq1 = alignment[0].seqA
    seq2 = alignment[0].seqB

    formate=""
    if len(sys.argv) > 4:
        formate = sys.argv[4]
    
    if formate[0] == "c":
        write_clustal( seq1,seq2, out)
    else:
        with open(out, "w") as f:
            f.write(seq1 + "\n")
            f.write(seq2 + "\n")
