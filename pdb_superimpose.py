#!/usr/bin/python3


#######################################
###  COPYRIGHT: RENE STARITZBICHLER  ##
###             02.02.2020           ##
#######################################



import Bio.PDB
import Bio.AlignIO
import sys

if len( sys.argv) < 7:
    print( "USAGE_1:",sys.argv[0]," range: ALIGN.pdb  FIRST_RESIDUE_ID  LAST_RESIDUE_ID TEMPLATE.pdb  FIRST_RESIDUE_ID  LAST_RESIDUE_ID  OUT.pdb")
    print ( "USAGE_2:",sys.argv[0]," list:  ALIGN.pdb  ATOM_NAME_1:RESIDUE_ID_1:CHAIN_1 ... TEMPLATE.pdb  ATOM_NAME_1:RESIDUE_ID_1:CHAIN_1 ...  OUT.pdb")
    print ( "USAGE_3:",sys.argv[0]," chain: ALIGN.pdb  CHAIN  TEMPLATE.pdb  CHAIN  OUT.pdb")
    print ( "USAGE_3:",sys.argv[0]," alignment: ALIGN.pdb  CHAIN  TEMPLATE.pdb  CHAIN ALIGNMENT.clw  OUT.pdb")
    exit(1)

first_atoms = []
second_atoms = []

first_file = ""
second_file = ""

mode = sys.argv[1]

first_ids = []
second_ids = []

if mode == "range:":
    first_file = sys.argv[2]    
    id_first_start = int(sys.argv[3])
    id_first_last = int(sys.argv[4])

    second_file = sys.argv[5]
    id_second_start = int(sys.argv[6])
    id_second_last = int(sys.argv[7])

    out_file = sys.argv[8]

    first_ids = range( id_first_start , id_first_last+1 )
    second_ids = range( id_second_start , id_second_last+1 )

elif mode == "chain:":
    first_file = sys.argv[2]
    first_chain = sys.argv[3]
    second_file = sys.argv[4]
    second_chain = sys.argv[5]
    out_file = sys.argv[6]

elif mode == "alignment:":
    first_file = sys.argv[2]
    first_chain = sys.argv[3]
    second_file = sys.argv[4]
    second_chain = sys.argv[5]
    alignment_file = sys.argv[6]
    out_file = sys.argv[7]

elif mode == "list:":
    first_file = sys.argv[2]
    index = 3
    while ".pdb" not in sys.argv[index]:
        first_ids.append( sys.argv[index].split(':') )
        index += 1

    second_file = sys.argv[index]
    index += 1
    while ".pdb" not in sys.argv[index]:
        second_ids.append( sys.argv[index].split(':') )
        index += 1
    
    out_file = sys.argv[index]

    print( first_file, first_ids )
    print( second_file, second_ids)
    print( out_file )
    
pdb_parser = Bio.PDB.PDBParser( QUIET=True )

first_structure = pdb_parser.get_structure( "first", first_file )[0]
second_structure = pdb_parser.get_structure( "second", second_file )[0]

if mode == "range:":
    for chain in first_structure:
        for residue in chain:
            if residue.get_id()[1] in first_ids:
                first_atoms.append( residue['CA'] )

    print ( "first chain, nr atoms:", len(first_atoms))
    for chain in second_structure:
        for residue in chain:
            if residue.get_id()[1] in second_ids:
                #print residue.get_id()[1],
                second_atoms.append( residue['CA'] )
    print( "second chain, nr atoms:", len(second_atoms))

elif mode == "chain:":
    #print "superimpose chains"
    #t = 0
    for chain in first_structure:
        #c = 0
        #for a in chain.get_atoms():
        #    c += 1
        #print( chain.get_id() , c )
        #t += c
        if chain.get_id() == first_chain:
            for residue in chain:
                if 'CA' in residue:
                    first_atoms.append( residue['CA'] )
                else:
                    print( "a residue without CA?:", residue )
    #print "total:", t
    for	chain in second_structure:
        if chain.get_id() == second_chain:
            for	residue	in chain:
                if 'CA' in residue:
                    second_atoms.append( residue['CA'] )
                else:
                    print( "a residue without CA?:", residue )

elif mode == "alignment:":
    alignment = Bio.AlignIO.read( open( alignment_file ), 'clustal')

    for chain in first_structure:
        if chain.get_id() == first_chain:
            residues = list(chain)
            count = 0
            for i in range( 0, len( alignment[0].seq )):
                if alignment[0].seq[i] != '-':
                    if alignment[1].seq[i] != '-':
                        #print( count, residues[count])
                        if count >= len(residues):
                            print( "ERROR: alignment longer than first chain: ", len(alignment[0].seq), len(residues), 'first residues: ', alignment[0].seq[:10], residues[:3])
                        if 'CA' not in residues[count]:
                            print( 'CA not found in residue:', residues[count])
                        first_atoms.append( residues[count]['CA'] )
                    count += 1

    for	chain in second_structure:
        if chain.get_id() == second_chain:
            residues = list(chain)
            count = 0
            for i in range( 0, len( alignment[1].seq )):
                if alignment[1].seq[i] != '-':
                    if alignment[0].seq[i] != '-':
                        second_atoms.append( residues[count]['CA'] )
                    count += 1
                    
elif mode == "list:":
    for chain in first_structure:
        for residue in chain:
            #print( chain.get_id(), residue.get_id() )
            for ids in first_ids:
                if chain.get_id() == ids[2] and residue.get_id()[1] == int( ids[1] ):
                    first_atoms.append( residue[ids[0]] )
                    print( "match", ids)
    for chain in second_structure:
        for residue in chain:
            #print( chain.get_id(), residue.get_id() )
            for ids in second_ids:
                if chain.get_id() == ids[2] and residue.get_id()[1] == int( ids[1] ):
                    second_atoms.append( residue[ids[0]] )
                    print( "match", ids)
                
if len(second_atoms) != len(first_atoms):
    print ("WARNING: number of atoms do not match!", len(first_atoms),len(second_atoms))
            
aligner = Bio.PDB.Superimposer()
aligner.set_atoms( second_atoms, first_atoms ) # place latter on former
print( "RMSD", first_file + ":", aligner.rms)

for chain in first_structure:
    aligner.apply( chain.get_atoms() ) # apply rotation and translation 

#c = 0
#for a in first_structure.get_atoms():
#    c += 1
#print c , "atoms"
               

io = Bio.PDB.PDBIO()
io.set_structure( first_structure )
io.save( out_file)

           
#rotmatrix = aligner.rotran[0]
#transvector = si.rotran[1]
