#!/usr/bin/python3

import sys

def residue_name( line ):
    return line[17:20]
def resid( line):
    return int( line[22:26].strip() )
def chain( line ):
    return line[21]



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

if len(sys.argv) < 4:
    print('USAGE:',sys.argv[0],'in:ORIG.pdb in:MUTANT.pdb out:MUTANT_DIFF.pdb')
    exit(1)

orig = []
with open( sys.argv[1] ) as r:
    residue = []
    prev_res = -9999
    prev_ch = 'xxxx'
    for l in r:
        if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
            res = resid(l)
            c = chain(l)
            if res != prev_res or c != prev_ch:
                prev_res = res
                prev_ch = c
                if len(residue) > 0:
                    orig.append( residue )
                    residue = []
            residue.append(l)
    if len(residue)>0:
        orig.append(residue)
 
mut = []
with open( sys.argv[2] ) as r:
    residue = []
    prev_res = -9999
    prev_ch = 'xxxx'
    for l in r:
        if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
            res = resid(l)
            c = chain(l)
            if res != prev_res or c != prev_ch:
                prev_res = res
                prev_ch = c
                if len(residue) > 0:
                    mut.append( residue )
                    residue = []
            residue.append(l)
    if len(residue)>0:
        mut.append(residue)


    

if len(orig) != len(mut):
    print( "ERROR, nr of lines do not match:", len(orig), len(mut))
    exit(1)
    
print( len(orig))

with open( sys.argv[3], 'w') as w:
    for i in range( len(orig)):
        a = residue_name(orig[i][0])
        b = residue_name(mut[i][0])
        x = '{:6.3f}'.format( (hydrophobicity(a) - hydrophobicity(b) ) )
        if len(x) > 6:
            x = x[:6]
        for m in mut[i]:
            w.write( m[:60] + x + m[66:] )
