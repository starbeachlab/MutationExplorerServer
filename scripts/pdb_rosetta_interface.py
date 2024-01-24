"""

pdb_rosetta_interface.py is a script that claculates differences in the energy score between bound and unbound states of the complex. Said scores are then written to the b-factor column of the output .pdb file.

Author: Aleksandra Panfilova

Date of last major changes: 23 Jan 2024

"""

from datetime import datetime

start_import = datetime.now()

#Load packages
from pyrosetta import *
import numpy as np
import argparse
import sys
init()

print('###### Interface Score Script ######')

start = datetime.now()

#Parse the commandline input
parser = argparse.ArgumentParser(description='Calculating per residue interface binding energies')
parser.add_argument('input',
                    type=str,
                    help='Path to the input pdb')
parser.add_argument('output', 
                    type=str,
                    help='Path to the output pdb')
parser.add_argument('-i', '--interface',
                    type=str,
                    default='all',
                    help='Chains composing left side of the interface')

args = parser.parse_args()
input_path = args.input
output_path = args.output
in_face = args.interface

#Function running the InrefaceAnalyzerMover 
def interface_analyzer_func(interface, scorefxn, pose):
    iam = InterfaceAnalyzerMover(interface)
    iam.set_scorefunction(scorefxn)
    iam.set_pack_input(True) #repack the input pose: True
    iam.set_pack_separated(True) #repack the exposed interfaces when calculating binding energy: True
    iam.set_pack_rounds(5) #do 5 rounds of packing (default 1)
    
    #Run IAM 3 times, and get median score for every position
    per_res_dG_sample = np.zeros(shape=(3, len(pose.sequence())))
    for r in range(3):
        iam.apply(pose)
        per_res_dG_sample[r, :] = np.array(iam.get_all_per_residue_data().dG) #gets dG values per residue (see Mover documentation)
    per_res_dG_array = np.median(per_res_dG_sample, axis=0).tolist() #get median scores out of 3
    return per_res_dG_array

#Parse PDB file to Rosetta pose
from pyrosetta.toolbox import cleanATOM
cleanATOM(input_path, out_file='cleaned.pdb', ext="")
pose = pyrosetta.pose_from_pdb('cleaned.pdb')

#loading score function
scorefxn = get_fa_scorefxn()

#Creating a list of chains for unspecified interface
#Get list of chians 
chains = []
for chain_num in np.array(pyrosetta.rosetta.core.pose.get_chains(pose)).tolist():
    chain_id = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(chain_num, pose)
    chains.append(chain_id)
#Getting a list of the last residue in every chain
chain_end_res = pyrosetta.rosetta.core.pose.chain_end_res(pose)
chain_end_res_list = np.array(chain_end_res).tolist()
#Adding a 0 for indexing
chain_end_res_list.insert(0, 0)

if not in_face == 'all':
    check_chains = list(in_face)
    check_chains.remove('_')

#Running the InterfaceAnalyzerMover and assembling a list of dG_binding values
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    
if in_face == 'all':
    print('all vs all mode')
    if len(chains) == 1:
        #To make sure single-chained PDBs do not crush
        per_res_dG = np.zeros(len(pose.sequence()))
    
    elif len(chains) == 2:
        interface = chains[0] + '_' + chains[1] #Interface example: 'A_B'
        per_res_dG = interface_analyzer_func(interface, scorefxn, pose) #extracts dG_binding score and changes type to list. now we have a list of N scores, where N is total residue number

    else:
        #creates a list of zeroes to be later filled with numbers for every chain separately
        #same Mover, but iterates through chains
        per_res_dG = [0] * len(pose.sequence())
        for c in range(0, len(chains)):
            #specifying interface. For chains='ABC' iterations will be: 'A_BC', 'B_AC', 'C_AB',
            #where the left side of the interface is the chain for which scores are saved
            interface = chains[c] + '_' + ''.join(chains[:c]) + ''.join(chains[(c+1):])

            #Get per residue numbers and assembles a list with values from different iterations
            per_res_list = interface_analyzer_func(interface, scorefxn, pose)
            #get indexes for the main (left) chain in this iteration
            chain_start = chain_end_res_list[c]
            chain_end = (chain_end_res_list[c+1])
            #replace zeroes with scores for the residues in the main chain
            per_res_dG[chain_start:chain_end] = per_res_list[chain_start:chain_end]
            
elif sorted(chains) == sorted(check_chains):
    print('interface ok')
    per_res_dG = interface_analyzer_func(in_face, scorefxn, pose)
    
else:
    print("Incomplete interface")
    exit()

#Replace b-factor values in the pose with our dG binding scores
for res in range(1,len(pose.sequence())+1):
    b = per_res_dG[res-1]
    #To avoid negative zeroes after rounding:
    if abs(b) < 0.005:
        b = 0.0
    b = round(b, 2)
    for ii in range(1,(pose.residue(res).natoms()+1)):
        pose.pdb_info().bfactor(res,ii,b)
        
#save pdb file
pose.dump_pdb(output_path)

end = datetime.now()

print('Run time with import: ', end-start_import)
print('Run time without import: ', end-start)