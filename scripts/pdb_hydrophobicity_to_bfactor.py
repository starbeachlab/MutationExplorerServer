#!/usr/bin/python3

import sys

def residue_name( line ):
    return line[17:20]

def hydrophobicity( residue):
    if residue.upper() == "Ala".upper(): return  1.800  
    if residue.upper() == "Arg".upper(): return -4.500  
    if residue.upper() == "Asn".upper(): return -3.500  
    if residue.upper() == "Asp".upper(): return -3.500  
    if residue.upper() == "Cys".upper(): return  2.500  
    if residue.upper() == "Gln".upper(): return -3.500  
    if residue.upper() == "Glu".upper(): return -3.500  
    if residue.upper() == "Gly".upper(): return -0.400  
    if residue.upper() == "His".upper(): return -3.200  
    if residue.upper() == "Ile".upper(): return  4.500  
    if residue.upper() == "Leu".upper(): return  3.800  
    if residue.upper() == "Lys".upper(): return -3.900  
    if residue.upper() == "Met".upper(): return  1.900  
    if residue.upper() == "Phe".upper(): return  2.800  
    if residue.upper() == "Pro".upper(): return -1.600  
    if residue.upper() == "Ser".upper(): return -0.800  
    if residue.upper() == "Thr".upper(): return -0.700  
    if residue.upper() == "Trp".upper(): return -0.900  
    if residue.upper() == "Tyr".upper(): return -1.300  
    if residue.upper() == "Val".upper(): return  4.200
    print( 'WARNING unknown residue type:', residue)
    return 0.0

if len(sys.argv) < 3:
    print( 'USAGE:', sys.argv[0], 'ORIG.pdb HYDROPH.pdb')
    exit(1)

with open( sys.argv[2], 'w') as w:
    with open( sys.argv[1] ) as r:
        for l in r:
            if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
                x = '{:6.3f}'.format(hydrophobicity( residue_name(l) ) )
                if len(x) > 6:
                    x = x[:6]
                w.write( l[:60] + str(x) + l[66:] )
            else:
                w.write(l)
