#!/usr/bin/python3
from collections import defaultdict
import sys

def resid( line):
    return int( line[22:26].strip() )
def chain( line ):
    return line[21]
def residue_name( line ):
    return line[17:20]
def ReadPDB( fname):
    res = defaultdict(list)
    with open( fname) as r:
        for l in r:
            if 'ATOM' == l[:4] or 'HETATM' == l[:6]:
                res[ chain(l) + str(resid(l)) ].append(l)
    return res


first = ReadPDB( sys.argv[1])
second = ReadPDB( sys.argv[2])

with open( sys.argv[3], 'w') as w:
    for k2,v2 in second.items():
        if residue_name(first[k2][0]) != residue_name(v2[0]):
            for l in v2:
                w.write(l)
