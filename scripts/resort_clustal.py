import os
from . import io

# amino acid class for reordering sequences
class AminoAcid:
    new_position = None
    def __init__(self, residue_name, chain_id, residue_seq_id, prev_position):
        self.residue_name = residue_name
        self.chain_id = chain_id
        self.residue_seq_id = residue_seq_id
        self.prev_position = prev_position

    def __repr__(self):
        return f'{self.chain_id} - {self.residue_seq_id} - {self.residue_name} - {self.prev_position}/{self.new_position}'

    def __str__(self):
        return f'{self.chain_id} - {self.residue_seq_id} - {self.residue_name} - {self.prev_position}/{self.new_position}'

# checks if line starts with atom or hetatm
def is_atom(line):
    if line.startswith(('ATOM', 'HETATM')):
        return True
    else:
        return False

# checks for each chain the order of the resids
# checks whether chains are mixed up
def is_chain_ordered(structure_lines):
    current_residue_number = -1

    for line in structure_lines:
        if is_atom(line):
            residue_number = int(line[22:26])

            if residue_number < current_residue_number:
                return False
            current_residue_number = residue_number

    return True

# gets all lines for a specifc chain from the atom lines
def get_lines_for_chain(structure_lines, chain):
    chain_lines = []

    for line in structure_lines:
        if is_atom(line):
            tmp_chain = line[21]
            if tmp_chain == chain:
                chain_lines.append(line)
    
    return chain_lines

# create an AminoAcid sequence from all lines within a chain
def create_aa(chain_lines, chain):
    aa = []

    current_res_number = None
    current_position = -1

    for line in chain_lines:
        res_number = int(line[22:26])
        res_name = line[17:20]

        if res_number != current_res_number:
            current_res_number = res_number
            current_position = current_position + 1
            new_aa = AminoAcid(res_name, chain, res_number, current_position)
            aa.append(new_aa)
        
    return aa

# sorts all AminoAcids in a list by their resid
def sort_amino_acids_by_resid(amino_acid_list):
    sorted_amino_acids = sorted(amino_acid_list, key=lambda x: x.residue_seq_id)
    for index, aa in enumerate(sorted_amino_acids):
        aa.new_position = index
    
    return sorted_amino_acids

# create dict mapping new positions of amino acids to their original one
def before_after(amino_acid_list):
    # { before : after }
    positions = { amino_acid.prev_position: amino_acid.new_position for amino_acid in amino_acid_list}
    return positions

def chain_resids_sorted(chain, pdb_file):
    all_lines = io.read_file_lines(pdb_file)

    if all_lines is None:
        print('An error occurred during the reading of the file for the residue sort check.')
        return None

    lines = get_lines_for_chain(all_lines, chain)
    if not is_chain_ordered(lines):
        print(f'Chain {chain} not ordered!')
        return False
    else:
        print(f'Chain {chain} is ordered!')
        return True

def sort_resids(chain, pdb_file):
    all_lines = io.read_file_lines(pdb_file)

    if all_lines is None:
        print('An error occurred during the reading of the file for the residue sorting.')
        return None

    lines = get_lines_for_chain(all_lines, chain)
    amino_acid_list = create_aa(lines, chain)
    amino_acid_list = sort_amino_acids_by_resid(amino_acid_list)
    positions = before_after(amino_acid_list)
    
    return positions    

def reorder_string(original_string, mapping_dict):
    new_order = [None] * len(original_string)

    for original_position, target_position in mapping_dict.items():
        if 0 <= original_position < len(original_string) and 0 <= target_position < len(original_string):
            new_order[target_position] = original_string[original_position]

    for i in range(len(new_order)):
        if new_order[i] is None:
            new_order[i] = original_string[i]

    reordered_string = ''.join(new_order)
    return reordered_string

def reorder_clustal(chain, pdb_file, clustal_file, clustal_file_reordered):
    position_mapping = sort_resids(chain, pdb_file)

    clustal = io.read_clustal(clustal_file)
    conservation = io.read_clustal_conservation(clustal_file)

    chain_names = list(clustal.keys())

    seq1 = clustal[chain_names[0]]
    seq2 = clustal[chain_names[1]]

    seq1_reordered = reorder_string(seq1, position_mapping)
    seq2_reordered = reorder_string(seq2, position_mapping)
    conservation_reordered = reorder_string(conservation, position_mapping)

    io.write_clustal(seq1_reordered, seq2_reordered, chain_names[0], chain_names[1], clustal_file_reordered, conservation_reordered)
    

