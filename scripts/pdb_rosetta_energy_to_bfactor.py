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



lines = []
vals = []
ind = -1
nr = 0
with open( sys.argv[1] ) as r:
    for l in r:
        lines.append( l)
        cols = l.split()
        if "label" == l[:5] and 'total' in l:
            nr = len(cols)
            ind = cols.index( 'total')
        elif ind > -1 and len(cols) == nr  and l[0] != '#' and l[:4] != 'pose' and l[:7] != 'weights':
            vals.append( [ l[:3] , float( l.split()[ind]) ] )
print( len(lines), 'atoms')
print( len(vals), 'residues')
prev_chain = 'XXX'
prev = -99999
mid = -1
with open( sys.argv[2], 'w') as w:
    for l in lines:
        if len(l) > 50 and ('ATOM' == l[:4] or 'HETATM' == l[:6]):
            ch = chain(l)
            res = resid(l)
            if prev != res or prev_chain != ch:
                prev = res
                prev_chain = ch
                mid += 1
                v = vals[mid]
                if v[0] != residue_name(l):
                    print( 'ERROR: resnames do not match!: ', v[0], residue_name(l) , 'at', mid )
                    exit(1)
                x = '{:6.3f}'.format(v[1])
                if len(x) > 6:
                    x = x[:6]
            w.write( l[:60] + str(x) + l[66:] )
        else:
            w.write(l)
