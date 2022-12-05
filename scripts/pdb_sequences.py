#!/usr/bin/python3


#######################################
###  COPYRIGHT: RENE STARITZBICHLER  ##
###             02.02.2020           ##
#######################################



import sys

filename = "/dev/null"
if len( sys.argv) < 2:
    print( "\nUSAGE: ", sys.argv[0], " PDB  (optional:CHAIN) (optional:OUTFILE)  (optional:'3' for 3-letter AA code, default is 1-letter code)")
    print( "\nNOTE: HETATMs are ignored!\n")
    exit(1)

is_three_letter = False
# use_hetatms = True
chains = []

for arg in sys.argv[2:]:
    if len(arg) > 1:
        filename = arg
    elif arg == '3':
        print( "set to three letter AA code")
        is_three_letter = True
    else:
        chains.append( arg )
        
out = open( filename , 'w')

aa = {}

aa[ "Gly".upper() ] = "G"    								
aa[ "Ala".upper() ] = "A"    								
aa[ "Leu".upper() ] = "L"    								
aa[ "Met".upper() ] = "M"    								
aa[ "Phe".upper() ] = "F"    								
aa[ "Trp".upper() ] = "W"    								
aa[ "Lys".upper() ] = "K"    								
aa[ "Gln".upper() ] = "Q"    								
aa[ "Glu".upper() ] = "E"    								
aa[ "Ser".upper() ] = "S"    						
aa[ "Pro".upper() ] = "P"                           
aa[ "Val".upper() ] = "V"                           
aa[ "Ile".upper() ] = "I"                   
aa[ "Cys".upper() ] = "C"                
aa[ "Tyr".upper() ] = "Y"                
aa[ "His".upper() ] = "H"                
aa[ "Hsd".upper() ] = "H"                
aa[ "Arg".upper() ] = "R"                
aa[ "Asn".upper() ] = "N"                   
aa[ "Asp".upper() ] = "D"     
aa[ "Thr".upper() ] = "T"   

#for k,v in aa.iteritems():
#    print k, v


prev_res = -999999     
prev_chain = "zzz"    
is_first = True       
unknown = []
with open( sys.argv[1]) as f:
    for line in f:

        if line[0:4] != "ATOM" : continue

        resid = int( line[22:27].strip()  )
        chain = line[21]

        if len(chains) > 0 and chain not in chains: continue

        if chain != prev_chain or resid < prev_res:
            if is_first:
                is_first = False
            else:
                print()
                out.write( "\n")
            header = ">" + sys.argv[1].split("/")[-1].split(".")[0] + ":" + chain + " \n"
            out.write( header)
            print( header, end='')
            prev_chain = chain

        if resid != prev_res:
            residue = line[17:20]
            if is_three_letter:
                out.write( residue + " " )
                print( residue, end=' ' )
            else:
                if residue in aa:
                    out.write( aa[ residue ] )
                    print( aa[residue], end='')
                else:
                    out.write( 'X' )
                    print( 'X', end='')
                    if residue not in unknown:
                        unknown.append( residue )
            prev_res = resid
print()
out.write( "\n")
out.close()

if len( unknown) > 0:
    print( "WARNING: UNKNOWN RESIDUE TYPES: ", unknown)

