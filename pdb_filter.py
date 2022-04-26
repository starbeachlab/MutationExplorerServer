#!/usr/bin/python3


import sys

with open( sys.argv[2], 'w') as w:
    with open( sys.argv[1], 'r') as r:
        for l in r:  # iterate file line by line
            if l[:6] != "HETATM":
                w.write( l)
    
