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
    print('USAGE:',sys.argv[0],'in:ORIG.pdb in:MUTANT.pdb out:MUTANT_DIFF.pdb')
    print( 'ORIG and MUTANT must contain Rosetta energies in the same order and have the same number or residues')
    exit(1)


lines_orig = []
vals_orig = []
ind = -1
nr = 0
atom_found = False
with open( sys.argv[1] ) as r:
    for l in r:
        if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
            atom_found = True
        else:
            if atom_found == True:
                lines_orig.append(l)
            if "label" == l[:5] and 'total' in l:
                cols = l.split()
                nr = len(cols)
                ind = cols.index( 'total')
            elif ind > -1 and len(l) > 10 and l[0] != '#' and l[:4] != 'pose' and l[:7] != 'weights':
                cols = l.split()
                vals_orig.append( [ cols[0] , float( cols[ind]) ] )

pdb = []
lines_mut = []
vals_mut = []
ind = -1
nr = 0
atom_found = False
with open( sys.argv[2] ) as r:
    for l in r:
        if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
            atom_found = True
            pdb.append(l)
        else:
            if atom_found == True:
                lines_mut.append(l)
            else:
                pdb.append(l)
                
            if "label" == l[:5] and 'total' in l:
                cols = l.split()
                nr = len(cols)
                ind = cols.index( 'total')
            elif ind > -1 and len(l) > 10 and l[0] != '#' and l[:4] != 'pose' and l[:7] != 'weights':
                cols = l.split()
                vals_mut.append( [ cols[0] , float( cols[ind]) ] )

vals = []
for a,b in zip(vals_mut, vals_orig):
    if a[0] != b[0]:
        print( 'mutated:', a[0], b[0] )
    vals.append( [a[0], a[1] - b[1]] )

            


prev_chain = 'XXX'
prev = -99999
mid = -1
with open( sys.argv[3], 'w') as w:
    for l in pdb:
        if len(l) > 50 and ('ATOM' == l[:4] or 'HETATM' == l[:6]):
            ch = chain(l)
            res = resid(l)
            if prev != res or prev_chain != ch:
                prev = res
                prev_chain = ch
                mid += 1
                v = vals[mid]
                #if residue_name(l) not in v[0]:
                #    print( 'ERROR: resnames do not match!: ', v[0], residue_name(l) , 'at', mid )
                #    exit(1)
                x = '{:6.3f}'.format(v[1])
                if len(x) > 6:
                    x = x[:6]
            w.write( l[:60] + str(x) + l[66:] )
        else:
            w.write(l)
    for l1,l2 in zip( lines_mut, lines_orig):        
        c1 = l1.split()
        c2 = l2.split()
        if len(c1) < 2 or l1[0] == '#' or 'label' == c1[0] or 'weights' == c1[0]:
            w.write(l2)
        else:
            w.write( c1[0] + ' ')
            for s1,s2 in zip( c1[1:], c2[1:]):
                if 'nan' in s1 or 'nan' in s2:
                    w.write( 'nan')
                else:
                    w.write( str( float(s1) - float(s2)) + ' ')
        w.write( '\n' )

                
