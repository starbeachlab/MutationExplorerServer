#!/usr/bin/python3

import sys

with open( sys.argv[1] ) as r, open( sys.argv[2], 'a') as w:
    for l in r:
        if l[:4] == "pose":
            print( l.split()[-1] , file=w)
            exit(0)
