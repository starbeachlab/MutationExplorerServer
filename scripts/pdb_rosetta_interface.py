from datetime import datetime

start_import = datetime.now()

#Load packages
from pyrosetta import *
import numpy as np
import argparse
import sys
init()

start = datetime.now()

#Parse the commandline input
parser = argparse.ArgumentParser(description='Calculating per residue interface binding energies')
parser.add_argument('input', type=str,
                    help='Path to the input pdb')
parser.add_argument('output', type=str,
                    help='Path to the output pdb')
args = parser.parse_args()
input_path = args.input
output_path = args.output

#Parse PDB file to Rosetta pose
pose = pyrosetta.pose_from_pdb(input_path)

#Get list of chians 
chains = []
for chain_num in np.array(pyrosetta.rosetta.core.pose.get_chains(pose)).tolist():
    chain_id = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(chain_num, pose)
    chains.append(chain_id)

#loading score function
scorefxn = get_fa_scorefxn()
    
#Getting a list of the last residue in every chain
chain_end_res = pyrosetta.rosetta.core.pose.chain_end_res(pose)
chain_end_res_list = np.array(chain_end_res).tolist()

#Adding a 0 for indexing
chain_end_res_list.insert(0, 0)

#Running the InterfaceAnalyzerMover and assembling a list of dG_binding values
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
if len(chains) == 2:
    #Interface example: 'A_B'
    interface = chains[0] + '_' + chains[1]
    iam = InterfaceAnalyzerMover(interface)
    iam.set_scorefunction(scorefxn)
    #repack the exposed interfaces when calculating binding energy: True
    iam.set_pack_separated(True)
    iam.apply(pose)
    #adds averaged interface info to the PDB.
    #doesn't make sense for multichain (different for every iteration)
    iam.add_score_info_to_pose(pose)
    #gets all types of scores (see IAMover documentation)
    per_res = iam.get_all_per_residue_data()
    #extracts dG_binding score and changes type to list.
    #now we have a list of N scores, where N is total residue number
    per_res_dG = np.array(per_res.dG).tolist()
else:
    #creates a list of zeroes to be later filled with numbers for every chain separately
    #same Mover, but iterates through chains
    per_res_dG = [0] * len(pose.sequence())
    for c in range(0, len(chains)):
        #specifying interface. For chains='ABC' iterations will be: 'A_BC', 'B_AC', 'C_AB',
        #where the left side of the interface is 
        interface = chains[c] + '_' + ''.join(chains[:c]) + ''.join(chains[(c+1):])

        #loading InterfaceAnalyzer 
        iam = InterfaceAnalyzerMover(interface)
        iam.set_scorefunction(scorefxn)
        iam.set_pack_separated(True)
        iam.apply(pose)

        #Get per residue numbers and assembles a list with values from different iterations
        per_res = iam.get_all_per_residue_data()
        per_res_list = np.array(per_res.dG).tolist()
        #get indexes for the main (left) chain in this iteration
        chain_start = chain_end_res_list[c]
        chain_end = (chain_end_res_list[c+1])
        #replace zeroes with scores for the residues in the main chain
        per_res_dG[chain_start:chain_end] = per_res_list[chain_start:chain_end]
        
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
print('Run time without import: ', end-start)#Load packages


