#!/usr/bin/python3


import sys


in_score_file1 = sys.argv[1]
in_score_file2 = sys.argv[2]
in_pdb_file = sys.argv[3]
out_pdb_file = sys.argv[4]


# rosetta score file

def scores_append(file_name): 
    scores = []
    with open(file_name ) as r:
        r.readline()
        for l in r:
            c = l.split()
            scores.append(  float( c[-2] ))
    return scores


#if len(open(in_pdb_file).readlines) == len(open(in_score_file).readlines):
    
score1 = scores_append( in_score_file1 )
score2 = scores_append( in_score_file2 )

if len(score1) != len(score2):
    print("ERROR: score files do not match in length", len(score1), len(score2)) 

difference = []    
for i in range(len(score1)):
    difference.append( score1[i] - score2[i] )
mini = min(difference)
maxi = max(difference)
x = max( abs(mini), abs(maxi) )
mini = -x
maxi = x
    
# PDB FILES in/out
counter = 0
acount = 0
resid = -999
chain = ''
with open( out_pdb_file, 'w') as w:
    with open( in_pdb_file, 'r') as r:
        for l in r:
            if l[:4] != "ATOM": continue
            acount += 1
            if resid  == -999 or chain == '':
                resid = int( l[22:26] )
                chain = l[21]
            if l[21] != chain or int(l[22:26]) != resid:
                resid = int(l[22:26])
                chain = l[21]
                counter += 1
            w.write( l[:60] + ( "%6.2f" % difference[counter] ) + l[66:] )
    acount += 1
    print( 'ATOM  ' + ( "%5d" % acount ) + "  H   DUM X   1       0.000   0.000   0.500      " + ( "{:6.2f}".format(mini) )[:6] + "      X    H ", file=w)
    acount += 1
    print( 'ATOM  ' + ( "%5d" % acount ) + "  H   DUM X   1       0.000   0.000   0.500      " + ( "{:6.2f}".format( maxi) )[:6] + "      X    H ", file=w)
    print( 'RANGE:', mini , maxi, file=w )

                
    



            

            

''''        
mini = min( scores)
maxi = max( scores)
print( 'min:', mini, 'max:', maxi)

maxi = abs( mini)


for i in range( len( scores ) ):
    scores1[i] = min( maxi, scores[i] )


            
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
    print( 'RANGE:', mini , maxi, file=w )
                
                
'''

    
