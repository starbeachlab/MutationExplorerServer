from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, jsonify
import os, random, subprocess, time
import threading
import requests
import glob
from collections import defaultdict
from werkzeug.utils import secure_filename

app = Flask(__name__)


app.config['USER_DATA_DIR'] = "/disk/user_data/mutation_explorer_gamma/"
app.config['ROSETTA_PATH']  = "/home/hildilab/dev/ext/rosetta/bin/"
app.config['SCRIPTS_PATH']  = "/home/hildilab/app/mutation_explorer_gamma/scripts/"


@app.route('/')
def index():
    return render_template("home.html")


def secure_str( string):
    string = string.replace('.','').replace('/','').replace('*','').replace('?','').replace('!','').replace( ' ', '')
    if len(string) == 0:
        print( "WARNING:", __name__, 'empty string after cleanup!')
    return string

def download_file(url, file_path):
        req = requests.get(url)
        with open(file_path, 'w') as f:
            f.write(req.content)


def start_thread(function, args, name):
    t = threading.Thread(target=function, name=name, args=args)
    t.deamon = True
    t.start()


def build_tree(node, links):
    children = []
    for l in links:
        if l[0] == node:
            children.append(l[1]) 
    print(len(children), 'children')
    
    tree = {}
    for child in children:
        tree[child] = build_tree(child, links)

    return tree


def build_mutation_tree(tag, root):
    info = app.config['USER_DATA_DIR'] + tag + "/info/"
    files = os.listdir(info)

    links = []
    for f in files:
        parent = open(info + f).read().split("\n")[0]
        links.append([parent.split(".")[0], f.split(".")[0]])
    print( len(links), 'links')
    
    # default root: "-"
    return build_tree(root, links)


def name_mutation(base_structure, tag):
    strc = base_structure.split(".")[0]
    mut_tree = build_mutation_tree(tag, strc)

    nums = [0]
    for ke in mut_tree.keys():
        if strc in ke:
            nums.append(int(ke.split("_")[-1]))
    
    return strc + "_" + str(max(nums) + 1) + ".pdb"


def fixbb(tag, structure, resfile, out_file_name, logfile):
    out = app.config['USER_DATA_DIR'] + tag + "/"
    log = open( out + logfile, 'a')
    ### Ich wuerde Funktionen immer so einfach wie moeglich halten und auf eine Aufgabe fokusieren, d.h. hier die Ausfuehrung von fixbb
    

    # generate unique extension (temp file)  ## WARUM?
    while(True):
        ext = "m" + str(random.randint(0, 999))
        if not glob.glob(out + ext + r'.*\.pdb$'):
            break

    # call rosetta
    cmd = "tsp " + app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -in:file:s " + out + structure + " -resfile " + out + resfile + ' -nstruct 1 -linmem_ig 10 -out:pdb  -out:prefix ' + out + ext + '_'
    log.write(cmd+'\n')
    p = subprocess.check_output(cmd.split())
    log.write(p+'\n')

    # rename output file #### WRITE ENERGIES INSTEAD !!!!!
    bfac =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_to_bfactor.py "
    ediff = app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_diff.py "
    hydro = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_to_bfactor.py '
    hdiff = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_diff_to_bfactor.py '
    
    cmd = "tsp mv " + out + ext + "_" + structure.split('.')[0] + "_0001.pdb " + out + out_file_name
    log.write(cmd+'\n')
    p = subprocess.check_output( cmd.split())
    log.write(p+'\n')
    
    cmd = "tsp " + bfac + out + out_file_name + ' ' + out + out_file_name[:-4] + '_absE.pdb'
    log.write(cmd+'\n')
    p = subprocess.check_output(cmd.split())
    log.write(p+'\n')

    if "mut" in structure:
        cmd = "tsp " + ediff + out + structure + ' ' + out + out_file_name + ' ' + out + out_file_name[:-4] + '_diffE.pdb'
        log.write(cmd+'\n')
        p = subprocess.check_output(cmd.split())
        log.write(p+'\n')

    cmd = "tsp " + hydro + out +  out_file_name + ' ' + out +  out_file_name[:-4] + '_HyPh.pdb'
    log.write(cmd+'\n')
    p = subprocess.check_output(cmd.split())
    log.write(p+'\n')

    cmd = "tsp " + hydro + out + structure + ' ' + out + structure[:-4] + '_HyPh.pdb'
    log.write(cmd+'\n')
    p = subprocess.check_output(cmd.split())

    if "mut" in structure:
        cmd = "tsp " + hdiff + out + structure + ' ' + out + out_file_name + ' ' + out + out_file_name[:-4] + '_diffHyPh.pdb'
        log.write(cmd+'\n')
        p = subprocess.check_output(cmd.split())
        log.write(p+'\n')

    log.write("### THREAD FINISHED ###\n")

        


@app.route('/submit', methods=['GET', 'POST'])
def submit(): 
    if request.method == 'GET':
        return render_template("submit.html", error = "")

    # generate tag
    while(True):
        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        if not os.path.exists(outdir):
            break
    os.mkdir(outdir)

    # get form values
    upload = request.files['pdbfile']
    pdb = secure_filename( request.form['pdbid'].strip() )
    af = secure_filename( request.form['alphafoldid'].strip() )

    # save file
    file_path = outdir + "structure.pdb"    
    if upload.filename != "":
        upload.save(file_path)
    elif pdb != "":
        download_file("https://files.rcsb.org/download/" + pdb + ".pdb", file_path)
    elif af != "":
        download_file("https://alphafold.ebi.ac.uk/files/AF-" + af + "-F1-model_v4.pdb", file_path)
    else:
        # no structure -> error
        return render_template("submit.html", error = "Please provide a structure")

    # create log file
    with open(outdir + "log.txt", "w") as f:
        f.write("submit\n")
        # TODO: actually log something

    # create list file
    with open(outdir + "list.txt", "w") as f:
        f.write("structure.pdb\n")

    # create info file
    os.mkdir(outdir + "info/")
    with open(outdir + "info/mut_0.txt", "w") as f:
        f.write("-")

    resfile =  "mut_0_resfile.txt"   
    mutations_to_resfile( [] , outdir + resfile )    

    # relax structure
    start_thread(fixbb, [tag, "structure.pdb", resfile, "mut_0.pdb", "log.txt"], "minimisation")

    return redirect(url_for('mutate', tag = tag))


@app.route('/mutate/<tag>', methods=['GET', 'POST'])
def mutate(tag):
    if request.method == 'GET':
        return render_template("mutate.html", tag = tag, error = "")

    outdir =   app.config['USER_DATA_DIR'] + tag + "/"
    #w = open( outdir + 'log.txt', 'a')
    ###  get form values
    mutations = request.form['mutations'].strip().replace( ' ','').split(',')
    if len(mutations)==1 and mutations[0] == '':
        print(__name__, 'reset')
        mutations = []
    #for m in mutations:
    #    w.write( m + '\n')
    #w.write( '?\n')
    vcf = request.files['vcf']
    #w.write( 'vcf: ' + vcf.filename + '\n')
    # allow multiple for following
    clustal1 = request.files['clustal1']
    clustal2 = request.files['clustal2']
    clustal3 = request.files['clustal3']
    #w.write( 'clw: ' + clustal1.filename + '\n')
    fasta1 = request.files['fasta1']
    fasta2 = request.files['fasta2']
    fasta3 = request.files['fasta3']
    #w.write('fasta1: '+ fasta1.filename + '\n')
    chainF1 = request.form['chainF1']
    chainF2 = request.form['chainF2']
    chainF3 = request.form['chainF3']
    #w.write( 'fasta1 chain: ' + fasta_chain1 + '\n')
    seq_input1 = request.form['sequence1'].strip()
    seq_input2 = request.form['sequence2'].strip()
    seq_input3 = request.form['sequence3'].strip()
    #w.write( "seq1 " + seq_input1+ '\n')
    chainS1 = request.form['chainS1']
    chainS2 = request.form['chainS2']
    chainS3 = request.form['chainS3']
    #w.write( "chain: " + seq_chain1+ '\n')
    uniprot1 = request.form['uniprot1'].strip()
    uniprot2 = request.form['uniprot2'].strip()
    uniprot3 = request.form['uniprot3'].strip()
    #w.write( "uniprot1: " + uniprot1+ '\n')
    chainU1 = request.form['chainU1']
    chainU2 = request.form['chainU2']
    chainU3 = request.form['chainU3']
    
    #w.write( "chain: " +  uniprot_chain1+ '\n')
    #
    #    mail = request.form['email'].strip()
    #    name = request.form['name'].strip() + ".pdb"  
    #    if name == ".pdb":

    name = name_mutation("mut_0.pdb", tag)
    #w.write( name + '\n')
    print(name)

    parent =   "mut_0.pdb"
    resfile =  name[:-4] + "_resfile.txt"
    align =    name[:-4] + ".clw"  # align with parent
    mutfile =  "info/" + name[:-4] + ".txt"  # parent and mutations
    
    ###  return error message if no mutations given
    if len(mutations) == 0 and vcf.filename == '' and  clustal1.filename == '' and  fasta1.filename == '' and seq_input1 == '' and uniprot1 == '':
        print( 'no mutations defined')
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation") # nutzlos, da javascript das gar nicht durchlaesst ohne eingabe 
    else:
        print( 'mutations defined')

    ### wait for parent to exist:
    wait( outdir + parent)
    
    ###  case separation
    if vcf.filename != "":
        vcf_file = os.path.join( outdir,  vcf.filename )
        vcf.save( vcf_file)
        add_mutations_from_vcf( mutations, vcf_file, outdir + parent)        
    if clustal1.filename != "":
        clustal_file = os.path.join( outdir , clustal1.filename )
        clustal1.save( clustal_file)
        add_mutations_from_alignment( mutations, clustal_file, outdir + parent)
    if clustal2.filename != "":
        clustal_file = os.path.join( outdir , clustal2.filename )
        clustal2.save( clustal_file)
        add_mutations_from_alignment( mutations, clustal_file, outdir + parent)
    if clustal3.filename != "":
        clustal_file = os.path.join( outdir , clustal3.filename )
        clustal3.save( clustal_file)
        add_mutations_from_alignment( mutations, clustal_file, outdir + parent)
    if fasta1.filename != "" and chainF1 != "":
        secure_str(chainF1)
        chainF1 = chainF1[0]
        fasta_file =  outdir + fasta1.filename
        fasta1.save(fasta_file)
        target = seq_from_fasta( fasta_file)
        add_mutations_from_sequence( mutations, target, chainF1, outdir+parent)
    if fasta2.filename != "" and chainF2 != "":
        secure_str(chainF2)
        chainF2 = chainF2[0]
        fasta_file =  outdir + fasta2.filename
        fasta2.save(fasta_file)
        target = seq_from_fasta( fasta_file)
        add_mutations_from_sequence( mutations, target, chainF2, outdir+parent)
    if fasta3.filename != "" and chainF3 != "":
        secure_str(chainF3)
        chainF3 = chainF3[0]
        fasta_file =  outdir + fasta3.filename
        fasta3.save(fasta_file)
        target = seq_from_fasta( fasta_file)
        add_mutations_from_sequence( mutations, target, chainF3, outdir+parent)
    if seq_input1 != "" and chainS1 != "":
        secure_str( chainS1)
        chainS1 = chainS1[0]
        secure_str( seq_input1 )
        add_mutations_from_sequence( mutations, seq_input1, chainS1, "seq1", outdir+parent)
    if seq_input2 != "" and chainS2 != "":
        secure_str( seq_input2 )
        chainS2 = chainS2[0]
        add_mutations_from_sequence( mutations, seq_input2, chainS2, "seq2", outdir+parent)
    if seq_input3 != "" and chainS3 != "":
        secure_str( seq_input3 )
        chainS3 = chainS3[0]
        add_mutations_from_sequence( mutations, seq_input3, chainS3, "seq3", outdir+parent)
    if uniprot1 != "" and chainU1 != '':
        uni_file = outdir + uniprot1
        download_uniprot( uniprot1, uni_file)
        target = seq_from_fasta( uni_file)
        add_mutations_from_sequence( mutations, target, chainU1, outdir + parent)
    if uniprot2 != "" and chainU2 != '':
        uni_file = outdir + uniprot2
        download_uniprot( uniprot2, uni_file)
        target = seq_from_fasta( uni_file)
        add_mutations_from_sequence( mutations, target, chainU2, outdir + parent)
    if uniprot3 != "" and chainU3 != '':
        uni_file = outdir + uniprot3
        download_uniprot( uniprot3, uni_file)
        target = seq_from_fasta( uni_file)
        add_mutations_from_sequence( mutations, target, chainU3, outdir + parent)
    print(__name__, 'total number of mutations:', len(mutations))
    if len(mutations) != 0:
        helper_files_from_mutations( mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile)
    else:
        print("no mutations")

    start_thread(fixbb, [tag, parent, resfile, name, "log.txt"], "mutti")
        
    return redirect(url_for('status', tag = tag, filename = name))


@app.route('/vcf', methods=['GET', 'POST'])
def vcf():
    return render_template("vcf.html", error = "")


@app.route('/adjustment', methods=['GET', 'POST'])
def adjustment():
    if request.method == 'GET':
        return render_template("adjustment.html", error = "")

    # generate tag
    while(True):
        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        if not os.path.exists(outdir):
            break
    os.mkdir(outdir)

    # get form values
    seq_upload = request.files['fasta']
    seq_input = request.form['sequence'].strip()
    uniprot = request.form['uniprot'].strip()

    strc_upload = request.files['pdbfile']
    pdb = request.form['pdbid'].strip()
    af = request.form['alphafoldid'].strip()

    # save file
    #    structure
    strc_file_path = outdir + "structure.pdb"    
    if strc_upload.filename != "":
        strc_upload.save(strc_file_path)
    elif pdb != "":
        download_file("https://files.rcsb.org/download/" + pdb + ".pdb", strc_file_path)
    elif af != "":
        download_file("https://alphafold.ebi.ac.uk/files/AF-" + af + "-F1-model_v4.pdb", strc_file_path)
    else:
        # no structure -> error
        return render_template("adjustment.html", error = "Please provide a structure")

    #    sequence
    seq_file_path = outdir + "seq.fasta"    
    if seq_upload.filename != "":
        seq_upload.save(seq_file_path)
    elif seq_input != "":
        with open(seq_file_path, "w") as f:
            f.write("> input seq\n")
            f.write(seq_input)
    elif uniprot != "":
        download_file("https://rest.uniprot.org/uniprotkb/" + uniprot + ".fasta", seq_file_path)
    else:
        # no sequence -> error
        return render_template("adjustment.html", error = "Please provide a sequence")

    # create list file
    with open(outdir + "list.txt", "w") as f:
        f.write("structure.pdb\n")

    # create info file
    os.mkdir(outdir + "info/")
    with open(outdir + "info/structure.txt", "w") as f:
        f.write("-")

    # extract structure sequence
    cmd = app.config['SCRIPTS_PATH'] + "pdb_sequences.py " + strc_file_path
    p = subprocess.check_output(cmd.split())
    with open(outdir + "ex_seq.fasta", "w") as f:
        f.write(p)

    # align sequences
    seq = open(outdir + "seq.fasta").read().split("\n")[1]
    ex_seq = open(outdir + "ex_seq.fasta").read().split("\n")[1]

    cmd = "python3 " + app.config['SCRIPTS_PATH'] + "seq_align.py " + ex_seq + " " + seq + " " + outdir + "aligned.fasta"
    print(cmd)
    p = subprocess.check_output(cmd.split())

    aligned = open(outdir + "aligned.fasta").read().split()
    ex_seq = aligned[0]
    seq = aligned[1]

    # get mutations
    assert len(seq) == len(ex_seq)
    mutations = []
    chain = open(outdir + "ex_seq.fasta").read().split("\n")[0].split(":")[1].strip()
    resid = 1
    for resid in range(1, len(seq) + 1):
        ex_res = ex_seq[resid - 1]
        if ex_res == "-":
            continue
        res = seq[resid - 1]
        if ex_res != res and res != "-":
            mutations.append(chain + ":" + str(resid) + res)

    # mutate structure
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    start_thread(fixbb, [tag, "structure.pdb", mutations, "mut_structure.pdb", outdir + "log.txt"], "minimisation")

    return redirect(url_for('status', tag = tag, filename = "mut_structure.pdb"))
    


@app.route('/get_status/<tag>/<filename>')
def get_status(tag, filename):
    dirname = os.path.join( app.config['USER_DATA_DIR'], tag + "/" + filename )
    done = os.path.isfile(dirname)
    status = "done"
    msg = ""
    #print( 'get_status', dirname, str(done))
    return jsonify({'done': done, 'status': status, 'message': msg})


@app.route('/status/<tag>/<filename>')
def status(tag, filename):
    return render_template("status.html", tag = tag, filename = filename)


def build_list(d):
    s = ""
    for ke in d:
        link = "<a id='" + ke + "' href='' class='structures'>" + ke + "</a><br>"
        if d[ke] != {}:
            s += "<li><i class='toggle material-icons tiny rot'>play_arrow</i>" + link + "</li>"
            s += "<ul class='browser-default' style='list-style-type: none'>" + build_list(d[ke]) + "</ul>"
        else:
            s += "<li><i class='material-icons tiny grey-text'>play_arrow</i>" + link + "</li>"

    return s


@app.route('/explore/<tag>/<filename>', methods=['GET', 'POST'])
def explore(tag,filename):
    if request.method == 'GET':
        mut_tree = build_mutation_tree(tag, "-")
        print('explore::tree', mut_tree)
        structures = "<ul>" + build_list(mut_tree) + "</ul>"
        print('explore::tree:', structures)
        parent = ""
        mutations = ""
        with open( app.config['USER_DATA_DIR'] + tag + "/info/" + filename[:-4] + ".txt") as r:
            parent = r.readline().strip()
            mutations += r.readline().strip()
            for l in r:
                mutations += ',' + l.strip()
        print( __name__, filename , tag)
        return render_template("explore.html", tag = tag, structures = structures, parent=parent, mutations = mutations, filename=filename)

    # get form values
    mutations = request.form['mutations'].strip().replace(' ', '').split(',')
    strc = request.form['structure'].strip() 
    name = request.form['name'].strip() + ".pdb"
    if name == ".pdb":
        name = name_mutation(strc, tag)
    
    print('novel mutant:', name)
    # OUTDIR ???
        
    # mutate structure
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    helper_files_from_mutations( mutations, outdir + strc, outdir + name[:-4] + '_resfile.txt', outdir + name[:-4] + '.clw', outdir + "info/" + name[:-4] + '.txt' ) 
    start_thread(fixbb, [tag, strc,  name[:-4] + '_resfile.txt', name, "log.txt"], "remutate")
    
    return redirect(url_for('status', tag = tag, filename = name))


# Testing molstar / mdsrv
@app.route('/molstar/<tag>')
def molstar(tag):
    mut_tree = build_mutation_tree(tag, "-")
    structures = "<ul>" + build_list(mut_tree) + "</ul>"
    return render_template("explore4.html", tag = tag, structures = structures)
    
    
@app.route('/examples')
def examples():
    return render_template("examples.html")


@app.route('/faq')
def faq():
    return render_template("faq.html")


@app.route('/downloads/<tag>/<filename>')
def download(tag, filename):
    path = app.config['USER_DATA_DIR'] + tag + "/"
    if filename == "results" + tag + ".zip":
        if not os.path.isfile(path + filename):
            # create zip file
            os.chdir(path)
            cmd = ['zip','results' + tag + '.zip','list.txt']
            for f in os.listdir('.'):
                if ".pdb" in f:
                    cmd.append(f)
            p = subprocess.check_output(cmd)

    return send_from_directory(path,filename)


@app.route('/fixbbtest')
def fbt():
    fixbb("/disk/user_data/mutation_explorer_gamma/863040/", "structure_1.pdb", ['A:1A'], "", "log")
    return render_template("faq.html")


if __name__ == "__main__":
    app.run()

    
def mutations_to_resfile( mutations, resfile):
    with open(resfile, 'w') as f:
        f.write('NATRO\n')
        f.write('start\n')
        prev_chains = []
        for mut in mutations:
            if len(mut) == 0: continue
            print(__name__, mut)
            chain = mut[0]
            resid = mut[2:-1]
            res = mut[-1]
            
            if chain not in prev_chains:
                f.write('* ' + chain + ' NATAA \n')
                prev_chains.append(chain)
                
            f.write(resid + ' ' + chain + ' PIKAA ' + res + '\n')
                

            
def resid( line):
    return int( line[22:26].strip() )


def residue_name( line ):
    return line[17:20]


def chain( line ):
    return line[21]


def single_letter( residue ):
    if residue.upper() == "Gly".upper(): return "G"
    if residue.upper() == "Ala".upper(): return "A"
    if residue.upper() == "Leu".upper(): return "L"
    if residue.upper() == "Met".upper(): return "M"
    if residue.upper() == "Phe".upper(): return "F"
    if residue.upper() == "Trp".upper(): return "W"
    if residue.upper() == "Lys".upper(): return "K"
    if residue.upper() == "Gln".upper(): return "Q"
    if residue.upper() == "Glu".upper(): return "E"
    if residue.upper() == "Ser".upper(): return "S"
    if residue.upper() == "Pro".upper(): return "P"                           
    if residue.upper() == "Val".upper(): return "V"                           
    if residue.upper() == "Ile".upper(): return "I"                   
    if residue.upper() == "Cys".upper(): return "C"                
    if residue.upper() == "Tyr".upper(): return "Y"                
    if residue.upper() == "His".upper(): return "H"                
    if residue.upper() == "Hsd".upper(): return "H"                
    if residue.upper() == "Arg".upper(): return "R"                
    if residue.upper() == "Asn".upper(): return "N"                   
    if residue.upper() == "Asp".upper(): return "D"     
    if residue.upper() == "Thr".upper(): return "T"   
    return "X"



def sequence_chain_resids( parent):
    chains = defaultdict(list)
    with open( parent) as r:
        prev_chain = 'xxx'
        prev_id = -9999
        for l in r:
            if l[:4] != "ATOM" and l[:6] != "HETATM": continue
            c = chain(l)
            res = resid(l)
            if prev_chain != c or prev_id != res:
                chains[c].append( [ single_letter( residue_name(l)) , res ] )
                prev_chain = c
                prev_id = res
    return chains



def alignment_from_mutations( mutations, parent, align,mutant_file):
    chains = sequence_chain_resids( parent)
    print( len(chains))
    mutseq = ""
    parent_str = parent.split('/')[-1][:-4]
    mutant_str = mutant_file.split('/')[-1][:-4]
    for c, l in chains.items():
        w = open( align[:-4] + '_' + c + '.clw', 'w') 
        w.write( 'CLUSTAL W formatted output, created by MutationExplorer\n\n')
        curr_par = parent_str #+ ':' + c
        curr_mut = mutant_str #+ ':' + c
        l1 = len(curr_par)
        l2 = len(curr_mut)
        count = 0
        if l1 < l2:
            curr_par = curr_par.ljust(l2)
            maxl = l2
        elif l2 < l1:
            curr_mut = curr_mut.ljust(l1)
            maxl = l1
        else:
            maxl = l1
        curr_par += '\t'
        curr_mut += '\t'
        curr_match = ''.rjust(maxl) + '\t'
        astr = ''
        bstr = ''
        cstr = ''
        for aa,rid in l:
            astr += aa
            #print( c,i,s)
            mstr = c + ':' + str(rid)
            #print( mstr)
            mutti = ''
            for m in mutations:
                #print( 'm: ', m)
                if mstr == m[:-1]:
                    #print( 'found')
                    mutti = m[-1]
                    break
            if mutti != '':
                bstr += mutti
                cstr += ' '
            else:
                bstr += aa
                cstr += '*'
            if len(astr) == 60:
                w.write( curr_par + astr + '\n')
                w.write( curr_mut + bstr + '\n')
                w.write( curr_match + cstr + '\n\n')
                astr = ''
                bstr = ''
                cstr = ''
        if len(astr) > 0:
            w.write( curr_par + astr + '\n')
            w.write( curr_mut + bstr + '\n')
            w.write( curr_match + cstr + '\n\n')
        w.close()



            
            
def mutation_parent_file( mutations, parent, mutfile):   
    with open( mutfile, 'w') as f:
        if mutations:
            f.write(parent.split('/')[-1] + "\n")
            for mut in mutations:
                f.write(mut + "\n")
        else:
            f.write("-\n")

            
def helper_files_from_mutations( mutations, parent, resfile, align, mutfile):
    mutations_to_resfile( mutations, resfile)
    mutation_parent_file( mutations, parent, mutfile)
    alignment_from_mutations( mutations, parent, align, mutfile)

    
    
def add_mutations_from_vcf( mutations, vcf_file, parent):
    print('chris...')


def write_clustal( s1, s2, filename):
    c = filename.split('/')[-1][:-4].split('_')
    a = c[0]
    for i in range(1,len(c)-1):
        a += '_' + c[i]
    b = c[i]
    m = max( len(a), len(b))
    a = a.ljust(m)
    b = b.ljust(m)
    e = ''.ljust(m)
    cons = conservation( s1,s2)
    with open( filename, 'w') as w:
        w.write( 'CLUSTAL W formatted output, created by MutationExplorer\n\n')
        count = 0
        step = 60
        while count < len(s1):
            w.write( a + '\t' + s1[count:count+step] + '\n')
            w.write( b + '\t' + s2[count:count+step] + '\n')
            w.write( e + '\t' + cons[count:count+step] + '\n\n')
            count += step
    

def add_mutations_from_sequence( mutations, target, chain, idy, parent):
    h1 = idy
    h2 = os.path.basename(parent)[:-4]
    outdir = os.path.dirname( parent) + '/'
    f1 =  outdir + h1 + '.fa'
    f2 =  outdir + h2 + '.fa'
    f3 =  outdir + h2 + '_' + h1 + '.clw'
    pdb_seq = pdb2seq( parent)
    write_fasta(f1,target,h1)
    write_fasta(f2,pdb_seq[chain][0],h2)
    cmd = "python3  " + app.config['SCRIPTS_PATH'] + "seq_align.py " + f1 + ' ' + f2 + ' ' + f3 + ' clw'
    print(cmd)
    p = subprocess.check_output( cmd.split())
    print(p)
    mutations.extend( mutations_from_alignment( f3, parent) )

    
def seq_from_fasta( filename):
    seq = ''
    head = ''
    with open( filename) as r:
        l = r.readline()
        if l[0] != '>':
            print( "WARNING!:",filename, 'does not start with ">" in first line! First line is ignored!!!' )
        head = l[1:].strip()
        for l in r:
            if l[0] == '>':
                return head,seq
            seq += l.strip()
    return head,seq

def write_fasta( filename, seq, header=""):
    with open(filename,'w') as w:
        w.write('>' + header + '\n')
        w.write(seq + '\n')
    
def add_mutations_from_alignment( mutations, clustal_file, parent):
    mutations.extend( mutations_from_alignment( clustal_file, parent ) )


        
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



def sequences( ali):
    seqs = []
    for k,v in ali.items():
        seqs.append( v.replace('-','') )
    return seqs



# nearly same as sequence_chain_resids()
def pdb2seq( pdb):
    chains = defaultdict( lambda : ['',[]] )
    with open( pdb) as r:
        prev_chain = 'xxx'
        prev_id = -9999
        for l in r:
            if l[:4] != "ATOM" and l[:6] != "HETATM": continue
            c = chain(l)
            res = resid(l)
            if prev_chain != c or prev_id != res:
                chains[c][0] += single_letter( residue_name(l))
                chains[c][1].append( res )
                prev_chain = c
                prev_id = res
    return chains


    
def mutations_from_alignment( clustal, parent ):
    mutations = []
    ali = read_clustal( clustal)
    #seqs = sequences( ali)
    chains = pdb2seq( parent)
    #print( ali)
    #print( chains)
    chain_in_pdb = ''
    ali_id = ''
    count = 0
    for c,v in chains.items():
        for a,b in ali.items():
            if v[0] == b:
                count += 1
                chain_in_pdb = c 
                ali_id = a
                print( 'mutations(): matching sequences found:', a,c)
    if count != 1:
        print( 'ERROR: mutations() found ',count,'identical matches')
    if count == 0:
        print('bye')
        exit(1)
    ### this will not work for MSA!!
    for a,b in ali.items():
        if a != ali_id:
            aligned = b
            break
    pdbseq = chains[chain_in_pdb][0]
    pdbres = chains[chain_in_pdb][1]
    #print( __name__)
    #print( pdbseq)
    #print( aligned)
    #print( pdbres[:3] )
    count = 0
    for i in range( len( pdbseq )):
        while aligned[count] == '-':
            count += 1
        if aligned[count] != pdbseq[i]:
            mutations.append( chain_in_pdb + ':' + str(pdbres[i]) + aligned[count] )
        count += 1
    print( __name__, 'found', len(mutations), 'mutations')
    return mutations

                    
            
def wait( filename):
    while not os.path.isfile(filename):
        time.sleep(1)


        
def download_uniprot( unid, filename):
    link = "https://www.uniprot.org/uniprot/" + secure_filename(unid) + '.fasta'
    req = requests.get(link)
    with open(filename, "w") as f:
        f.write(req.content)

