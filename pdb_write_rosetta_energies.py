#!/usr/bin/python3


import sys

in_pdb_file = sys.argv[1]
in_score_file = sys.argv[2]
out_pdb_file = sys.argv[3]


# rosetta score file
scores = []
with open( in_score_file ) as r:
    r.readline()
    for l in r:
        if "description" not in l:
            c = l.split()
            scores.append(  float( c[-2] ))
        

# PDB FILES in/out
counter = 0
resid = -999
chain = ''
with open( out_pdb_file, 'w') as w:
    with open( in_pdb_file, 'r') as r:
        for l in r:
            if l[:4] != "ATOM": continue
            
            if resid  == -999 or chain == '':
                resid = int( l[22:26] )
                chain = l[21]
            if l[21] != chain or int(l[22:26]) != resid:
                resid = int(l[22:26])
                chain = l[21]
                counter += 1
            score = scores[counter]
            w.write( l[:60] + ( "%6.2f" % score) + l[66:] )

                
                
    
