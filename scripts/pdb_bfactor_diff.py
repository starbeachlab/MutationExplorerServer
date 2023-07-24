#!/usr/bin/python3

import sys

def resid( line):
    return int( line[22:26].strip() )
def chain( line ):
    return line[21]
def bfactor( line):
    return float( line[60:66] )
def atom_name( line):
    return line[12:16].strip() 
def residue_name( line ):
    return line[17:20]


if len(sys.argv) < 4:
    print('USAGE:',sys.argv[0],'in:ORIG.pdb in:MUTANT.pdb out:DIFF.pdb')
    print( 'the script works with per residue bfactors')
    exit(1)

orig_bfactors = []

prev_resid = -999999
prev_chain = 'xxx'

with open( sys.argv[1] ) as r:
    for l in r:
        if 'ATOM' == l[0:4] or 'HETATM' == l[0:6]:
            res_id = resid(l)
            chain_name = chain(l)
            if res_id != prev_resid or chain_name != prev_chain:
                orig_bfactors.append( bfactor(l) )
                prev_chain = chain_name
                prev_resid = res_id
                
prev_resid = -999999
prev_chain = 'xxx'
count = -1

with open( sys.argv[4] , 'w' ) as w,  open( sys.argv[3] ) as r:
        for l in r:
            if 'ATOM' == l[0:4] or 'HETATM' == l[0:6]:
                res_id = resid(l)
                chain_name = chain(l)
                if res_id != prev_resid or chain_name != prev_chain:
                    b_factor = bfactor(l)
                    prev_chain = chain_name
                    prev_resid = res_id
                    count += 1
                    x = '{:6.3f}'.format( orig_bfactors[count] - b_factor)
                    if len(x) > 6:
                        x = x[:6]                
                w.write( l[:60] + str(x) + l[66:] )
            else:
                w.write(l) # ?
