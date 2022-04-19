#!/usr/bin/python3

import sys
from collections import defaultdict

header = []
dic = defaultdict( str )
with open( sys.argv[1]) as r:
    r.readline()
    for l in r:
        l = l.strip()
        if len(l) == 0:
            continue
        c = l.split()

        if len(c) > 1 and len(c[1]) > 0:
            dic[c[0]] +=  c[1] 

for (d,v) in dic.items():
    print( '>' , d )
    print( v )
        
