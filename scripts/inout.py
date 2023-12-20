from collections import defaultdict

def read_clustal( clust):
    ali = defaultdict( str)
    with open(clust) as r:
        r.readline()
        for l in r:
            if l.strip() == '' or l[:3] == '   ': continue
            c = l.split()
            if len(c) != 2: continue
            ali[c[0]] += c[1]

    return ali

def get_clustal_seq_start(line):
    c = line.split()
    start_index = line.find(c[1])
    return start_index

def read_clustal_conservation(clustal):
    conservation = ''
    lines = read_file_lines(clustal)
    start_index = -1

    for l in lines:
        if l.strip() == '' or l[:3] == '   ': continue
        c = l.split()
        if len(c) != 2: continue
        index = get_clustal_seq_start(l)
        if index != -1:
            start_index = index
    
    for l in lines:
        if l.strip() == '': continue
        if l[:3] == '   ': 
            conservation += l[start_index: start_index+60]
            continue

    return conservation

def read_file_lines(file_path):
    all_lines = []
    try:
        with open(file_path, 'r') as file:
            all_lines = file.readlines()
        file.close()
    except FileNotFoundError:
        print(f'File {file_path} not found.')
        return None
    except Exception as e:
        print(f'An error ({e}) occurred while trying to read the following file: {file_path}')
        return None
    else:
        return all_lines
    
def write_clustal( s1, s2, c1, c2, filename, conservation_line=None):
    a = c1
    b = c2
    m = max( len(a), len(b))
    a = a.ljust(m)
    b = b.ljust(m)
    e = ''.ljust(m)
    cons = ''
    if conservation_line is None:
        cons = conservation( s1,s2)
    else:
        cons = conservation_line
    try:
        with open( filename, 'w') as w:
            w.write( 'CLUSTAL W formatted output, created by MutationExplorer\n\n')
            count = 0
            step = 60
            while count < len(s1):
                w.write( a + '\t' + s1[count:count+step] + '\n')
                w.write( b + '\t' + s2[count:count+step] + '\n')
                w.write( e + '\t' + cons[count:count+step] + '\n\n')
                count += step
        w.close()
    except IOError as e:
        print(f'Error while writing clustal file. Error: {e}')
    except Exception as e:
        print(f'An unexpected error occured while writing a clustal file: {e}')