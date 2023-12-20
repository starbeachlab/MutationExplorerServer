from . import inout

def check_his_replacement(pdb_file_path):
    all_lines = inout.read_file_lines(pdb_file_path)
    for line in all_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            residue_name = line[17:20].strip()
            
            if residue_name in ['HSD', 'HSE', 'HSP']:
                return True

    return False

def check_chain_id_missing(pdb_file_path):
    all_lines = inout.read_file_lines(pdb_file_path)
    for line in all_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain_id = line[20:22].strip()

            if chain_id == '':
                return True

    return False
