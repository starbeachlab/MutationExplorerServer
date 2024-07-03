
##############################
# author: Rene Staritzbichler
# 10.01.21
##############################

from Bio import PDB
from Bio import AlignIO
from Bio import SeqUtils
import sys


import basic_alignment as bali


aligned = 0.2
grouped = 0.5
factor = 3.0



if len( sys.argv) < 6:
    print( 'script writes sequence conservation into temperature factor of output PDB (optional:ALIGN_SCORE)')
    print( "USAGE:", sys.argv[0], 'PDB CHAIN ALIGNMENT.clw ROW OUTFILE')
    print( "\tALIGNMENT is Clustal formatted multiple sequence alignment")
    print( "\tROW specifies the sequence in the alignment that belongs to PDB")
    print( '\tROW starts with 0')
    print( '\tCHAIN is needed because alignments are single chain only')
    print( '\tALIGN_SCORE is to score aligned positions that are not identical (smaller 1)')
    print( '\tcall multiple times if you want to color more than one chain')
    print()
    help(SeqUtils.IUPACData)
    exit(1)



pdb_file = sys.argv[1]
chain = sys.argv[2]
ali_file = sys.argv[3]
ali_id = int( sys.argv[4] )
out_file = sys.argv[6]

print()
print()
print("outfile", out_file)
print()
print()

if len(sys.argv) > 5:
    aligned = float( sys.argv[5])
    print( 'score aligned positions with:', aligned)

groups = [
    [ 'R', 'H', 'K' ],
    [ 'D', 'E' ],
    [ 'S', 'T', 'N', 'Q'],
    [ 'C', 'U', 'G', 'P'],
    [ 'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V' ]]

def group_id( A):
    for i in range( len(groups)):
        if A in groups[i]:
            return i
    return -1

def same_group( A, B ):
    return group_id( A) == group_id(B)


alignment = bali.ReadAlignment( ali_file)[1]

alignment = [seq for seq in alignment if "*" not in seq]

nr = len( alignment)


if ali_id >= nr:
    print( "ERROR: sequence id too small:", ali_id, nr )
    exit(1)


vals = []
seq = []
for i in range( 0, len( alignment[ ali_id ] )):
    aa = alignment[ ali_id][i]
    if aa == '-':
        print( '-',end=' ' )
        continue
    seq.append( aa )
    count = 0
    for j in range( 0, nr):
        if j == ali_id: continue
        if alignment[j][i] == aa:
            count += 1.0
        elif same_group( alignment[j][i] , aa ):
            count += grouped
        elif alignment[j][i] != '-':
            count += aligned

    vals.append( float(count)/float(nr-1) )



with open(out_file, "w") as w:
    eidi = -1
    previd = -99999
    with open( pdb_file) as r:
        for l in r:
            l = l.strip()
            if "ATOM" == l[0:4] and l[21] == chain:
                new_Id = int(l[22:26])
                if previd != new_Id:
                    previd = new_Id
                    eidi += 1
                aa = SeqUtils.IUPACData.protein_letters_3to1[ l[ 17:20 ].title()]
               # if aa != seq[eidi]:
               #     print( "ERROR: sequence in alignment and pdb do not match:", eidi, seq[eidi],aa)
               #     exit(1)
                w.write( l[:60] + "{:6.2f}".format( factor * vals[eidi] ) + l[66:] + "\n" )
            #elif "ATOM" == l[0:4] or "HETATM" == l[0:6]:
            #    w.write( l[:60] + "{:6.2f}".format( factor ) + l[66:] + "\n" )
            else:
                w.write( l + "\n")
