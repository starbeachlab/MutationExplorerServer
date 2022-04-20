from collections import defaultdict

def ReadAlignment( myfile):
    dic = defaultdict( str )
    with open( myfile) as r:
        r.readline()
        for l in r:
            l = l.strip()
            if len(l) == 0:
                continue
            c = l.split()

            if len(c) > 1 and len(c[1]) > 0:
                dic[c[0]] +=  c[1] 

    alignment = []
    header = []
    for (key,value) in dic.items():
        alignment.append(value)
        header.append( key)
    return [ header, alignment]
