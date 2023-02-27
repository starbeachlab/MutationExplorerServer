from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, jsonify
import os, random, subprocess, time
import threading
import requests
import glob
from collections import defaultdict
from werkzeug.utils import secure_filename

app = Flask(__name__)


app.config['USER_DATA_DIR'] = "/disk/user_data/mutation_explorer_gamma/"
app.config['EXAMPLE_DIR'] = "/disk/data/mutation_explorer_gamma/examples/"
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


def get_alphafold_id(af):
    # get_alphafold_id tries to match various fractions of "AF-<UniProt>-F<X>-model_v<y>.pdb" to the full ID
    if af[:3] != "AF-":
        af = "AF-" + af
    if af[-3:] != "pdb":
        if af[-9:-1] != "-model_v":
            if af[-3:-1] != "-F":
                af += "-F1"
            af += "-model_v4"
        af += ".pdb"
    return af



def start_thread(function, args, name):
    t = threading.Thread(target=function, name=name, args=args)
    t.deamon = True
    t.start()


def build_tree(node, links):
    children = []
    for l in links:
        if l[0] == node:
            children.append(l[1]) 
    #print(len(children), 'children')
    
    tree = {}
    for child in children:
        tree[child] = build_tree(child, links)

    return tree


def build_mutation_tree(out, tag, root):
    info = out + tag + "/info/"
    files = os.listdir(info)

    links = []
    for f in files:
        parent = open(info + f).read().split("\n")[0]
        links.append([parent.split(".")[0], f.split(".")[0]])
    #print( len(links), 'links')
    
    # default root: "-"
    return build_tree(root, links)


def name_mutation(out, base_structure, tag):
    strc = base_structure.split(".")[0]
    mut_tree = build_mutation_tree(out, tag, strc)

    nums = [0]
    for ke in mut_tree.keys():
        if strc in ke:
            nums.append(int(ke.split("_")[-1]))
    
    return strc + "_" + str(max(nums) + 1)  + ".pdb"


def bash_cmd(cmd, log):
    log.write(cmd + '\n')
    p = subprocess.check_output(cmd.split())
    log.write(p + '\n')


def fixbb(tag, structure, resfile, out_file_name, logfile):
    #print("######## fixbb  #######")
    #print('infile:', structure, "out_file_name: ", out_file_name)
    out = app.config['USER_DATA_DIR'] + tag + "/"
    log = open( out + logfile, 'a')
    
    
    if out_file_name[-4:] == '.pdb':
        out_file_name = out_file_name[:-4]

    ext = out_file_name

    # call rosetta
    cmd = "tsp " + app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -use_input_sc -in:file:s " + out + structure + " -resfile " + out + resfile + ' -nstruct 1 -out:pdb -out:prefix ' + out
    if structure == "structure.pdb":
        cmd += " -linmem_ig 10  -ex1 -ex2  "
    bash_cmd(cmd, log)

    # rename output file #### WRITE ENERGIES INSTEAD !!!!!
    bfac =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_to_bfactor.py "
    ediff = app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_diff.py "
    hydro = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_to_bfactor.py '
    hdiff = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_diff_to_bfactor.py '
    mutti = app.config['SCRIPTS_PATH'] + 'pdb_mutated_aa.py '
    
    cmd = "tsp mv " + out + structure[:-4] + "_0001.pdb " + out + out_file_name + ".pdb"
    bash_cmd(cmd, log)
    
    cmd = "tsp " + bfac + out + out_file_name + '.pdb ' + out + out_file_name + '_absE.pdb'
    bash_cmd(cmd, log)

    if "mut" in structure:
        cmd = "tsp " + ediff + out + structure + ' ' + out + out_file_name + '.pdb ' + out + out_file_name + '_diffE.pdb' 
        bash_cmd(cmd, log)

    cmd = "tsp " + hydro + out +  out_file_name + '.pdb ' + out +  out_file_name + '_HyPh.pdb'
    bash_cmd(cmd, log)

    cmd = "tsp " + hydro + out + structure + ' ' + out + structure[:-4] + '_HyPh.pdb' # TODO: wieso structure[:-4] statt out_file_name?
    bash_cmd(cmd, log)
    
    if "mut" in structure:
        cmd = "tsp " + hdiff + out + structure + ' ' + out + out_file_name + '.pdb ' + out + out_file_name + '_diffHyPh.pdb'
        bash_cmd(cmd, log)

        cmd = "tsp " + mutti + out + structure + " " + out + out_file_name + '.pdb ' + out + out_file_name + '_aa.pdb'
        bash_cmd(cmd, log)

    
    if not os.path.exists(out + "fin/"):
        os.mkdir(out + "fin/")

    cmd = "tsp touch " + out + "fin/" + out_file_name + ".pdb"
    bash_cmd(cmd, log)
        
    log.write("### THREAD FINISHED ###\n")
    log.close()


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
    chain_filter = request.form['chain_filter'].strip().upper()
    hetatom_filter = request.form.get('hetatom_filter')  # .get() needed for checkbox

    email = request.form['email'].strip()
    min_type = request.form['min-selector'] # = short | long

    if email:
        results_link = "https://proteinformatics.uni-leipzig.de" + url_for('explore', tag = tag, filename = "")
        write_email(outdir + "mail.txt", email, results_link)

    # save file
    file_path = outdir + "structure.pdb"    
    if upload.filename != "":
        upload.save(file_path)
    elif pdb != "":
        if pdb[-4:] == ".pdb":
            pdb = pdb[:-4]
        download_file("https://files.rcsb.org/download/" + pdb + ".pdb", file_path)
    elif af != "":
        af = get_alphafold_id(af)
        print( 'alphafold:', af)
        download_file("https://alphafold.ebi.ac.uk/files/" + af, file_path)
    else:
        # no structure -> error
        return render_template("submit.html", error = "Please provide a structure")

    if not is_pdb( file_path):
        return render_template("submit.html", error = "It was not possible to upload the structure you provided.")

    # filter structure
    remove_hets = (hetatom_filter is not None)
    if chain_filter != "" or remove_hets:
        chains = chain_filter.split()

        with open(file_path, "r") as f_in:
            with open(outdir + "structure2.pdb", "w") as f_out:
                for line in f_in:
                    if line[:4] == "ATOM":
                        if chain(line) in chains:
                            continue
                    if line[:6] == "HETATM":
                        if remove_hets:
                            continue
                    f_out.write(line)

        bash_cmd("mv " + outdir + "structure2.pdb " + file_path, open(outdir + "temp.log", "a")) # TODO: log file

    # create log file
    print("submit\n")
    if os.path.isfile( file_path):
        print( file_path + ' is downloaded\n')
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
    start_thread(fixbb, [tag, "structure.pdb", resfile, "mut_0", "log.txt"], "minimisation")
    print( 'fixbb started for initial upload\n')
    return redirect(url_for('mutate', tag = tag))


def add_mutations(tag, mutant, inputs):

    outdir =   app.config['USER_DATA_DIR'] + tag + "/"
    
    mutations = inputs["mutations"]
    if len(mutations) == 1 and mutations[0] == '':
        print(__name__, 'reset')
        mutations = []

    #mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    print( 'mutate: ' + mutant + '\n')
    print(mutant)

    parent =   "mut_0.pdb"
    resfile =  mutant[:-4] + "_resfile.txt"
    align =    mutant[:-4] + ".clw"  # align with parent
    mutfile =  "info/" + mutant[:-4] + ".txt"  # parent and mutations
    
    """
    ###  return error message if no mutations given
    if len(mutations) == 0 and   clustal1.filename == '' and  fasta1.filename == '' and seq_input1 == '' and uniprot1 == '':
        print( 'no mutations defined')
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation") # nutzlos, da javascript das gar nicht durchlaesst ohne eingabe 
    else:
        print( 'mutations defined')
    """
    
    print( 'wait for parent to exist\n')
    ### wait for parent to exist:
    if wait( outdir + parent, 1, 900) == False:
        #return render_template("mutate.html", tag = tag, error = "Your structure could not be uploaded.")
        return

    
    # get all mutations
    i = 0
    for clustal in inputs["clustals"]:
        if clustal.filename != "":
            clustal_file = os.path.join( outdir , clustal.filename )
            clustal.save( clustal_file)
            add_mutations_from_alignment( mutations, clustal_file, outdir + parent)
    for fasta, chainF in zip(inputs["fastas"], inputs["chainFs"]):
        i += 1
        if fasta.filename != "" and chainF != "":
            secure_str(chainF)
            chainF = chainF[0]
            fasta_file =  outdir + fasta.filename
            fasta.save(fasta_file)
            head, target = seq_from_fasta( fasta_file) 
            add_mutations_from_sequence( mutations, target, chainF, "fa" + (i % 3), outdir+parent)
    for seq_input, chainS in zip(inputs["seq_inputs"], inputs["chainSs"]):
        i += 1
        if seq_input != "" and chainS != "":
            secure_str(chainS)
            chainS = chainS[0]
            secure_str(seq_input)
            add_mutations_from_sequence( mutations, seq_input, chainS, "seq" + (i % 3), outdir+parent)
    for uniprot, chainU in zip(inputs["uniprots"], inputs["chainUs"]):
        i += 1
        if uniprot != "" and chainU != '':
            uni_file = outdir + uniprot
            download_uniprot( uniprot, uni_file)
            target = seq_from_fasta( uni_file)
            add_mutations_from_sequence( mutations, target, chainU, "uni" + (i % 3), outdir + parent)

    print(__name__, 'total number of mutations:', len(mutations))
    if len(mutations) != 0:
        helper_files_from_mutations( mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile)
    else:
        print("no mutations")
        #return render_template("mutate.html", tag = tag, error = "Please provide a mutation")
        return
    
    #start_thread(fixbb, [tag, parent, resfile, mutant, "log.txt"], "mutti") 
    fixbb(tag, parent, resfile, mutant, "log.txt")

    send_email(outdir + "mail.txt")


@app.route('/mutate/<tag>', methods=['GET', 'POST'])
def mutate(tag):

    outdir =   app.config['USER_DATA_DIR'] + tag + "/" # TR
    mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    
    if request.method == 'GET':
        chains = get_chains( outdir + "structure.pdb")
        return render_template("mutate.html", tag = tag, chains=chains, error = "")

    ###  get form values
    mutations = request.form['mutations'].strip().replace(' ','').split(',')

    """ # TR
    if len(mutations) == 1 and mutations[0] == '':
        print(__name__, 'reset')
        mutations = []
    """

    #for m in mutations:
    #    w.write( m + '\n')
    #w.write( '?\n')
    #vcf = request.files['vcf']
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

    inputs = {}
    inputs["mutations"] = mutations
    inputs["clustals"] = [clustal1, clustal2, clustal3]
    inputs["fastas"] = [fasta1, fasta2, fasta3]
    inputs["chainFs"] = [chainF1, chainF2, chainF3]
    inputs["seq_inputs"] = [seq_input1, seq_input2, seq_input3]
    inputs["chainSs"] = [chainS1, chainS2, chainS3]
    inputs["uniprots"] = [uniprot1, uniprot2, uniprot3]
    inputs["chainUs"] = [chainU1, chainU2, chainU3]

    start_thread(add_mutations, [tag, mutant, inputs], "add_muts")
    #add_mutations(tag, mutant, inputs)
    
    #w.write( "chain: " +  uniprot_chain1+ '\n')
    #
    #    mail = request.form['email'].strip()
    #    name = request.form['name'].strip() + ".pdb"  
    #    if name == ".pdb":

    """ # TR
    mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    print( 'mutate: ' + mutant + '\n')
    print(mutant)

    parent =   "mut_0.pdb"
    resfile =  mutant[:-4] + "_resfile.txt"
    align =    mutant[:-4] + ".clw"  # align with parent
    mutfile =  "info/" + mutant[:-4] + ".txt"  # parent and mutations
    
    ###  return error message if no mutations given
    if len(mutations) == 0 and   clustal1.filename == '' and  fasta1.filename == '' and seq_input1 == '' and uniprot1 == '':
        print( 'no mutations defined')
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation") # nutzlos, da javascript das gar nicht durchlaesst ohne eingabe 
    else:
        print( 'mutations defined')

    print( 'wait for parent to exist\n')
    ### wait for parent to exist:
    if wait( outdir + parent, 1, 900) == False:
        return render_template("mutate.html", tag = tag, error = "Your structure could not be uploaded.")
    
    ###  case separation
    #    if vcf.filename != "":
    #        vcf_file = os.path.join( outdir,  vcf.filename )
    #        vcf.save( vcf_file)
    #        add_mutations_from_vcf( mutations, vcf_file, outdir + parent)        
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
        head,target = seq_from_fasta( fasta_file) # TODO: sollte target nicht ein String sein? Ist ein Error bei mir
        add_mutations_from_sequence( mutations, target, chainF1, "fa1", outdir+parent)
    if fasta2.filename != "" and chainF2 != "":
        secure_str(chainF2)
        chainF2 = chainF2[0]
        fasta_file =  outdir + fasta2.filename
        fasta2.save(fasta_file)
        head,target = seq_from_fasta( fasta_file)
        add_mutations_from_sequence( mutations, target, chainF2, "fa2", outdir+parent)
    if fasta3.filename != "" and chainF3 != "":
        secure_str(chainF3)
        chainF3 = chainF3[0]
        fasta_file =  outdir + fasta3.filename
        fasta3.save(fasta_file)
        head,target = seq_from_fasta( fasta_file)
        add_mutations_from_sequence( mutations, target, chainF3, "fa3", outdir+parent)
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
        add_mutations_from_sequence( mutations, target, chainU1, "uni1", outdir + parent)
    if uniprot2 != "" and chainU2 != '':
        uni_file = outdir + uniprot2
        download_uniprot( uniprot2, uni_file)
        target = seq_from_fasta( uni_file)
        add_mutations_from_sequence( mutations, target, chainU2, "uni2", outdir + parent)
    if uniprot3 != "" and chainU3 != '':
        uni_file = outdir + uniprot3
        download_uniprot( uniprot3, uni_file)
        target = seq_from_fasta( uni_file)
        add_mutations_from_sequence( mutations, target, chainU3, "uni3", outdir + parent)
    print(__name__, 'total number of mutations:', len(mutations))
    if len(mutations) != 0:
        helper_files_from_mutations( mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile)
    else:
        print("no mutations")
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation")
    
    start_thread(fixbb, [tag, parent, resfile, mutant, "log.txt"], "mutti")

    print( 'started with nr mutations: ' + str( len(mutations) ) + '\n')
    """

    return redirect(url_for('status', tag = tag, filename = mutant))


@app.route('/vcf', methods=['GET', 'POST'])
def vcf():
    if request.method == 'GET':
        return render_template("vcf.html", error = "")
  
    # generate tag
    while(True):
        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        if not os.path.exists(outdir):
            break
    os.mkdir(outdir)

    vcf = request.files['vcf']
    vcf_file = vcf.filename
    if vcf_file == "":
        return render_template("vcf.html", error = "no filename was given")
    vcf.save( outdir + vcf_file)
    cmd = "tsp " + app.config['SCRIPTS_PATH'] + "run_vcf.sh " + outdir + ' ' + vcf_file
    print(cmd)
    p = subprocess.check_output(cmd.split())
    print(p)

    if wait( outdir +  vcf_file[:-4] + '_missense.csv', 1, 900) == False:
        return render_template("vcf.html", error = "No missense was found.")
    
    # create info file
    os.mkdir(outdir + "info/")
    with open(outdir + "info/mut_0.txt", "w") as f:
        f.write("-")
    with open( outdir + "mut_0_resfile.txt", 'w') as w:
        w.write('NATAA\nstart\n')

    alphafold,mutations = mutations_from_vcf( outdir + vcf_file[:-4] + '_missense.csv')
    print( 'alphafold:', alphafold)
    download_file("https://alphafold.ebi.ac.uk/files/" + alphafold.strip() + "-model_v4.pdb", outdir + 'structure.pdb' )
    if not is_pdb( outdir + 'structure.pdb'):
        return render_template("vcf.html", error = "It was not possible to upload the AlphaFold model: " + alphafold + "<br>Currently only the first candidate can be uploaded")

            
    # relax structure
    start_thread(fixbb, [tag, "structure.pdb", "mut_0_resfile.txt", "mut_0", "log.txt"], "minimisation")
    print( 'fixbb started for initial upload\n')

    if wait( outdir + 'mut_0.pdb', 1, 920) == False:
        return render_template("vcf.html", error = "Your structure could not be minimized.")

    helper_files_from_mutations( mutations,  outdir + 'mut_0.pdb',  outdir + 'mut_0_1_resfile.txt',  outdir + 'mut_0_1.clw',  outdir + 'info/mut_0_1.txt')

    start_thread(fixbb, [tag, 'mut_0.pdb', 'mut_0_1_resfile.txt', 'mut_0_1.pdb', "log.txt"], "mutti")


    
    return redirect(url_for('status', tag = tag, filename = "mut_0_1.pdb"))



@app.route('/interface', methods=["GET", "POST"])
def interface():
    if request.method == 'POST':

        # get data from form

        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        while os.path.exists(outdir):
            tag = str(random.randint(0, 999999))
            outdir = app.config['USER_DATA_DIR'] + tag + "/"

        print("############## DELTA INTERFACE ##############")
        print(tag)
        os.mkdir(outdir)

        #conv_file_upload = True
        conv_filename = ""

        file_conv = request.files['file_conv']
        if file_conv.filename != "":
            file_conv.save(outdir + secure_filename(file_conv.filename))
            conv_filename = secure_filename(file_conv.filename)
        else:
            conv_file_upload = False
            pdb_conv = request.form.getlist('pdb_conv')[0]
            file_conv_link = "https://files.rcsb.org/download/" + pdb_conv + ".pdb"
            req = requests.get(file_conv_link)
            with open(outdir + pdb_conv + ".pdb", "w") as f:
                f.write(req.content)
            conv_filename = pdb_conv + ".pdb"

        superimpose = False
        #super_file_upload = True
        super_filename = ""

        file_super = request.files['file_super']
        if file_super.filename != "":
            superimpose = True
            file_super.save(outdir + secure_filename(file_super.filename))
            super_filename = secure_filename(file_super.filename)
        else:
            pdb_super = request.form.getlist('pdb_super')[0]
            if pdb_super != "":
                file_super_link = "https://files.rcsb.org/download/" + pdb_super + ".pdb"
                superimpose = True
                #super_file_upload = False
                req = requests.get(file_super_link)
                with open(outdir + pdb_super + ".pdb", "w") as f:
                    f.write(req.content)
                super_filename = pdb_super + ".pdb"


        alignment_link = request.form.getlist('alignment_link')[0]
        #alignment_link = "https://www.bioinfo.mpg.de/AlignMeBeta/work/" + alignment_link.split("work/")[1]
        #alignment = outdir + "alignment.aln"
        print("########### alignment link")
        print(alignment_link)
        
        """
        #alignment_link = app.config['USER_DATA_DIR'] + "alignment.aln"

        req = requests.get(alignment_link)
        with open(alignment, "w") as f:
            f.write(req.content)
        """
        return "<h1>Hallo Welt</h1>"

    

@app.route('/get_status/<tag>/<filename>')
def get_status(tag, filename):
    status = ""
    msg = ""
    dirname = os.path.join( app.config['USER_DATA_DIR'], tag + "/fin/" + filename)
    done = os.path.isfile(dirname)
    if done:
        dirname = os.path.join( app.config['USER_DATA_DIR'], tag + "/" + filename )
        successful = os.path.isfile(dirname)
        if successful:
            status = "done"
            msg = ""
        else:
            status = "error"
            msg = "Rosetta calculation encountered an error"
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


def load_explore_page(out, tag, filename):
    mut_tree = build_mutation_tree(out, tag, "-")
    print('explore::tree', mut_tree)
    structures = "<ul>" + build_list(mut_tree) + "</ul>"
    print('explore::tree:', structures)
    parent = ""
    mutations = ""
    if filename == "":
        filename = "mut_0.pdb"
    with open( out + tag + "/info/" + filename[:-4] + ".txt") as r:
        parent = r.readline().strip()
        mutations += r.readline().strip()
        for l in r:
            mutations += ',' + l.strip()
    outdir = out + tag + "/"
    if parent == "-":
        chains = get_chains(outdir + filename)
    else:
        chains = get_chains( outdir + parent)
    energy = get_energy (outdir + filename)
    print( __name__, filename , tag, chains)
    return render_template("explore.html", tag = tag, structures = structures, parent=parent, mutations = mutations, filename=filename , chains = chains, energy=energy)


@app.route('/explore/<tag>/<filename>', methods=['GET', 'POST'])
@app.route('/explore/<tag>/', methods=['GET', 'POST'])
def explore(tag, filename = ""):
    if request.method == 'GET':
        return load_explore_page(app.config['USER_DATA_DIR'], tag, filename)

    # get form values
    mutations = request.form['mutations'].strip().replace(' ', '').split(',')
    parent = request.form['fname'].strip()
    if parent[-4:] != ".pdb":
        parent += '.pdb'
    mutant = name_mutation(app.config['USER_DATA_DIR'], parent, tag)
    
    print(__name__, 'original:', parent, 'novel mutant:', mutant)
    # OUTDIR ???
        
    # mutate structure
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    helper_files_from_mutations( mutations, outdir + parent, outdir + mutant[:-4] + '_resfile.txt', outdir + mutant[:-4] + '.clw', outdir + "info/" + mutant[:-4] + '.txt' ) 

    start_thread(fixbb, [tag, parent,  mutant[:-4] + '_resfile.txt', mutant, "log.txt"], "remutate")
    
    return redirect(url_for('status', tag = tag, filename = mutant))



# Testing molstar / mdsrv
@app.route('/molstar/<tag>/')
def molstar(tag):
    return render_template("molstar.html")
    
    
@app.route('/examples')
def examples():
    return render_template("examples.html")


@app.route('/examples/<tag>/')
def load_example(tag):
    return load_explore_page(app.config['EXAMPLE_DIR'], tag, "mut_0_1.pdb")


@app.route('/faq')
def faq():
    return render_template("faq.html")


@app.route('/downloads/<tag>/<filename>')
def download(tag, filename):
    if tag.isdigit():
        path = app.config['USER_DATA_DIR'] + tag + "/"
    else:
        path = app.config['EXAMPLE_DIR'] + tag + "/"
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
        f.write('NATAA\n')
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
"""                
def mutations_from_resfile( resfile):
    mutations = []
    with open( resfile) as r:
        r.readline()
        r.readline()
        for l in r:
            if 'NATAA' in 
"""

def mutations_from_vcf( fname):
    mutations = []
    alphafold = ""
    print( 'get mutations from:', fname)
    with open( fname) as r:
        r.readline()
        for l in r:
            c = l.split(',')
            if alphafold == "":
                alphafold = c[-1]
                mutations.append( 'A:' + c[-2][1:] )
            elif c[-1] == alphafold:
                mutations.append( 'A:' + c[-2][1:] )
    return alphafold,mutations

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
            if l[:4] != "ATOM": continue
            c = chain(l)
            res = resid(l)
            if prev_chain != c or prev_id != res:
                chains[c].append( [ single_letter( residue_name(l)) , res ] )
                prev_chain = c
                prev_id = res
    return chains



def alignment_from_mutations(mutations, parent, align, mutant_file):
    chains = sequence_chain_resids( parent)
    print( len(chains))
    mutseq = ""
    parent_str = parent.split('/')[-1][:-4]
    mutant_str = mutant_file.split('/')[-1][:-4]
    for c, l in chains.items():
        w = open( align[:-4] + '_' + c + '.clw', 'w') 
        w.write( 'CLUSTAL W formatted output, created by MutationExplorer\n\n')
        curr_par = parent_str + ':' + c
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


'''    
def helper_files_from_resfile( resfile, parent, mutfile, align):
    mutations = mutations_from_resfile( resfile)
    mutation_parent_file( mutations, parent, mutfile)
    alignment_from_mutations(mutations, parent, align, mutfile)
'''    


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
    print("add_mutations_from_sequence")
    print("idy: ", idy)
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
                print("######### check")
                print(prev_chain)
                print(prev_id)
                print(chains)
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
    at_end = False
    for i in range( len( pdbseq )):
        while aligned[count] == '-':
            count += 1
            if count == len(aligned) - 1:
                at_end = True
                break
        if at_end:
            break
        if aligned[count] != pdbseq[i]:
            mutations.append( chain_in_pdb + ':' + str(pdbres[i]) + aligned[count] )
        count += 1
    print( __name__, 'found', len(mutations), 'mutations')
    return mutations

                    
            
def wait( filename, step, maxw):
    for i in range(0, maxw):
        time.sleep(step)
        if os.path.isfile(filename):
            return True
    return False
        
def download_uniprot( unid, filename):
    link = "https://www.uniprot.org/uniprot/" + secure_filename(unid) + '.fasta'
    req = requests.get(link)
    with open(filename, "w") as f:
        f.write(req.content)

def get_chains( fname):
    chains = ""
    print( __name__, 'get', fname)
    with open( fname) as r:
        for l in r:
            if l[:4] == "ATOM" or l[:6] == "HETATM":
                c = chain(l)
                if c not in chains:
                    chains += c
    return chains

def get_energy( fname):
    print( __name__, 'get', fname)
    with open( fname) as r:
        for l in r:
            if l[:4] == "pose":
                return float( l.split()[-1])


def is_pdb( fname):
    if not os.path.isfile(fname):
        return False
    with open( fname) as r:
        for l in r:
            if l[:4] == "ATOM" or l[:6] == "HETATM":
                return True
    return False


def write_email(fil, user, link):
    with open(fil, 'w') as out:
        out.write('From: mutex@proteinformatics.de\n')
        out.write('Subject: MutationExplorer results \n')
        out.write('To: ' + user + '\n\n')
        out.write('Hello  ' + user + '!\n\n')
        out.write('Your MutationExplorer calculation is done. \n')
        out.write('You can view the results here: \n\n' + link + '\n\n')
        out.write('Thanks for using MutEx.\n\n')
        out.write('Have a nice day!\n\n')

def send_email(fil):
    if(os.path.exists(fil)):
        os.system('sendmail -t < ' + fil)


"""
def send_email(user, link):
    #log = open( 'calc.log', 'w' )
    #log.write( 'send mail\n')
    #log.write( os.getcwd() + '\n')
    with open( '/home/hildilab/app/mutation_explorer_delta/mail/mail.txt', 'w') as out:
        out.write('From: voronoia@proteinformatics.de\n')
        out.write('Subject: MutationExplorer results \n')
        out.write('To: ' + user + '\n\n')
        out.write('hello  ' + user + '!\n\n')
        out.write('your voronoia calculation is done. \n')
        out.write('you can view the results here: \n\n' + link + '\n\n')
        out.write('thanks for using voronoia.\n\n')
        out.write('have a nice day!\n\n')
    os.system( 'sendmail -t < /home/hildilab/app/mutation_explorer_delta/mail/mail.txt' )
    #log.close()
"""