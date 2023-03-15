from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, jsonify
import os, random, subprocess, time
import threading
import requests
import glob
from collections import defaultdict
from werkzeug.utils import secure_filename
import shutil
import datetime

cfg_file = 'app.cfg'

app = Flask( __name__ )
app.config.from_pyfile( cfg_file )





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
        with open(file_path, 'wb') as f:
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
        fstr = f.split(".")[0]
        # avoid failed mutations to be listed in tree:
        if os.path.isfile( out + tag + '/' + fstr + '.pdb' ):
            links.append([parent.split(".")[0], fstr])
    #print( len(links), 'links')
    #print( links)
    
    # default root: "none"
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
    log.write(cmd + "\n")
    p = subprocess.check_output(cmd.split())
    log.write(str(p) + "\n")


def fixbb(tag, structure, resfile, out_file_name, logfile, longmin=False, path_to_store=""):
    print("######## fixbb  #######")
    print(tag)
    print('infile:', structure, "out_file_name: ", out_file_name)
    out = app.config['USER_DATA_DIR'] + tag + "/"
    log = open( out + logfile, 'a')

    if out_file_name[-4:] == '.pdb':
        out_file_name = out_file_name[:-4]

    #ext = out_file_name
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    status = "Start+fixbb+for+" + out_file_name
    cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
    bash_cmd(cmd, log)



    # call rosetta
    cmd = "tsp " + app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -use_input_sc -in:file:s " + out + structure + " -resfile " + out + resfile + ' -nstruct 1  -linmem_ig 10 -out:pdb -out:prefix ' + out
    if longmin == True:
        cmd += " -ex1 -ex2  "
    bash_cmd(cmd, log)
    print("rosetta done")	


    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    status = "fixbb+for+" + out_file_name + "+done"
    cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
    bash_cmd(cmd, log)

    if(path_to_store != ""):
        cmd = "tsp cp " + out + structure[:-4] + "_0001.pdb " + path_to_store
        bash_cmd(cmd, log)

    cmd = "tsp mv " + out + structure[:-4] + "_0001.pdb " + out + out_file_name + ".pdb"
    bash_cmd(cmd, log)

    calc_rasp(tag, structure,out_file_name,logfile, path_to_store)
    print("MV Done")
    file_processing( tag, structure, out_file_name, logfile)
    print("File processing")
    log.close()



def calc_rasp(tag, structure, out_file_name, logfile, path_to_store=""):
    out = app.config['USER_DATA_DIR'] + tag + "/"

    if(len(glob.glob(out+"rasp.status")) > 0):
        log = open( out + logfile, 'a')
        chains = get_chains(out + "structure.pdb")
        chain_list = list(chains)

        for chain in chain_list:

            if(len(glob.glob( path_to_store + "_" + chain + ".csv" )) > 0 ):
                print("exists")
                listig =  glob.glob(  path_to_store + "_" + chain + ".csv" )

                if len(listig) == 1:
                    print( listig[0], out + "cavity_pred_" + out_file_name + "_" + chain + ".csv")
                    shutil.copyfile(listig[0],out + "cavity_pred_" + out_file_name + "_" + chain + ".csv")
            else:
                status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
                status = "Start+RaSP+calculation+for+" + out_file_name + "+with+chain+" + chain
                cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
                bash_cmd(cmd, log)


                cmd = "tsp " + "bash -i " + app.config['RASP_PATH'] + "calc-rasp.sh " + out + out_file_name + ".pdb " + chain + " " + out_file_name + " " + out
                print(cmd)
                bash_cmd(cmd, log)
           

                if(path_to_store != ""):
                    cmd = "tsp cp -d " + out + "cavity_pred_" + out_file_name + "_" + chain + ".csv " + path_to_store + "_" + chain + ".csv"
                    print(cmd)
                    bash_cmd(cmd, log)
        


                status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
                status = "RaSP+calculation+for+" + out_file_name + "+with+chain+" + chain + "+done"
                cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path + " " + "check_rasp " + os.path.join( app.config['USER_DATA_DIR'], tag + "/cavity_pred_" + out_file_name + "_" + chain + ".csv")
                print(cmd)
                bash_cmd(cmd, log)


    



def file_processing( tag, structure, out_file_name, logfile):
    # rename output file #### WRITE ENERGIES INSTEAD !!!!!
    bfac =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_to_bfactor.py "
    ediff = app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_diff.py "
    hydro = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_to_bfactor.py '
    hdiff = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_diff_to_bfactor.py '
    mutti = app.config['SCRIPTS_PATH'] + 'pdb_mutated_aa.py '
    
    out = app.config['USER_DATA_DIR'] + tag + "/"
    log = open( out + logfile, 'a')
    
    if out_file_name[-4:] == '.pdb':
        out_file_name = out_file_name[:-4]
    
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

    cmd = "tsp touch " + out + "fin/" + out_file_name + ".pdb"  ### TODO: what was this for??
    bash_cmd(cmd, log)

    cmd = "tsp " +  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_append.py " + out + out_file_name + ".pdb " + out + "info/" + out_file_name + ".txt"
    bash_cmd(cmd, log)
    
    log.close()



def create_user_dir():
    # generate tag
    while(True):
        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        if not os.path.exists(outdir):
            break
    os.mkdir(outdir)

    return outdir, tag



def save_pdb_file(file_path, upload, pdb, af):
    original_name = ""
    error = False
    error_message = ""
    msg = "x"

    if upload.filename != "":
        original_name = upload.filename
        upload.save(file_path)
    elif pdb != "":
        original_name = pdb
        if pdb[-4:] == ".pdb":
            pdb = pdb[:-4]
        if is_in_db( pdb):
            msg="found"
            cp_from_db(pdb,file_path)
        else:
            msg="notfound"
            download_file("https://files.rcsb.org/download/" + pdb + ".pdb", file_path)
    elif af != "":
        af_id = get_alphafold_id(af)
        original_name = af_id
        print( 'alphafold:', af_id)
        if is_in_db( af): # @Daniel: sollte das nicht auch af_id sein?
            msg="found"
            cp_from_db(af,file_path)
        else:
            msg="notfound"
            download_file("https://alphafold.ebi.ac.uk/files/" + af_id, file_path)
            
    else:
        # no structure -> error
        error = True
        error_message = "Please provide a structure" 
        #return render_template("submit.html", error = "Please provide a structure")

    if not is_pdb( file_path):
        error = True
        error_message = "It was not possible to upload the structure you provided."
        #return render_template("submit.html", error = "It was not possible to upload the structure you provided.")

    # extract first model in PDB (NMR structures)  ### REPLACE THIS WITH CLEANUP SCRIPT! 
    with open( file_path) as r, open( file_path[:-4] + '.tmp', 'w') as w:
        model_count = 0
        for l in r:
            if l[:5] == "MODEL":
                model_count += 1
            if l[:6] == 'ENDMDL' or model_count > 1:
                break
            w.write(l)
        os.rename( file_path[:-4] + '.tmp', file_path )

    return original_name, error, error_message, msg


def filter_structure(outdir, file_path, chain_filter, remove_hets):
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

        print("mv " + outdir + "structure2.pdb " + file_path, open(outdir + "temp.log", "a")) 
        bash_cmd("mv " + outdir + "structure2.pdb " + file_path, open(outdir + "temp.log", "a")) # TODO: log file
        return True
    return False



def relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af):
    resfile =  "mut_0_resfile.txt"   
    structure_file = "structure.pdb"
    log_file = "log.txt"

    mutations_to_resfile( [] , outdir + resfile )    

    path = ""

    if msg != "found":

        if (not filtered) and longmin == True:
            if(pdb != ""):
                rose = app.config['ROSEMINT_PATH']
                path = rose + "pdb/" + pdb.upper() + ".pdb"
            if(af != ""):    
                rose = app.config['ROSEMINT_PATH']
                path = rose + "alphafold/" + af.upper() + ".pdb"

        print(path)
        start_thread(fixbb, [tag, structure_file, resfile, "mut_0", log_file, longmin, path ], "minimisation")
    else:

        if (not filtered) and longmin == True:
            shutil.copyfile( outdir + structure_file, outdir + "mut_0.pdb")
            if(pdb != ""):
                rose = app.config['ROSEMINT_PATH']
                path = rose + "pdb/" + pdb.upper() + ".pdb"
            if(af != ""):    
                rose = app.config['ROSEMINT_PATH']
                path = rose + "alphafold/" + af.upper() + ".pdb"

            print("path rasp: " + path)
            calc_rasp(tag, structure_file, "mut_0", log_file, path ) # TODO

            file_processing( tag, structure_file, "mut_0",  log_file)

        else:
            print("Fixbb since filtering")
            start_thread(fixbb, [tag, structure_file, resfile, "mut_0", log_file, longmin], "minimisation")
            msg =  "notfound"




@app.route('/submit', methods=['GET', 'POST'])
def submit(): 
    if request.method == 'GET':
        return render_template("submit.html", error = "")


    ### get form values

    upload = request.files['pdbfile']
    pdb = secure_filename( request.form['pdbid'].strip() )
    af = secure_filename( request.form['alphafoldid'].strip() )

    chain_filter = request.form['chain_filter'].strip().upper()
    hetatom_filter = request.form.get('hetatom_filter')  # .get() needed for checkbox
    remove_hets = (hetatom_filter is not None)

    email = request.form['email'].strip()
    min_type = request.form['min-selector'] # = short | long
    longmin = (min_type == 'long')



    ### processing

    outdir, tag = create_user_dir()

    rasp_calculation = False
    rasp_checkbox = request.form.get('rasp-checkbox') # on none

    if(rasp_checkbox == "on"):
        rasp_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/rasp.status")
        with open(rasp_path, "w") as f:
            f.write(rasp_checkbox)
    
    print(rasp_calculation)



    # prewrite email (is sent seperately)
    if email:
        results_link = "https://proteinformatics.uni-leipzig.de" + url_for('explore', tag = tag, filename = "")
        write_email(outdir + "mail.txt", email, results_link)

    # save file
    file_path = outdir + "structure.pdb"    
    original_name, unsuccessful, error_message, msg = save_pdb_file(file_path, upload, pdb, af)
    if unsuccessful:
        return render_template("submit.html", error = error_message)

    # filter (remove) chains, heteroatoms
    filtered = filter_structure(outdir, file_path, chain_filter, remove_hets)

    print("submit\n")
    if os.path.isfile( file_path):
        print( file_path + ' is downloaded\n')

    # create info file for mut_0
    os.mkdir(outdir + "info/")
    with open(outdir + "info/mut_0.txt", "w") as f:
        f.write("none\n")

    #create status file
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name.log")
    # print("status")
    # print(status_path)
    with open(status_path, "w") as f:
        f.write(get_current_time()+"+Start+Calculation\n")

    with open(name_path, "w") as f:
        f.write(original_name)
        
    # relax structure
    relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af)
    print('fixbb started for initial upload\n')

    return redirect(url_for('mutate', tag = tag, msg=msg))


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
        #return render_template("mutate.html", tag = tag, error = "Your structure could not be uploaded.")    ### @Nikola: warum ist das auskommentiert? Rene
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
            add_mutations_from_sequence( mutations, target, chainF, "fa" + str(i % 3), outdir+parent)
    for seq_input, chainS in zip(inputs["seq_inputs"], inputs["chainSs"]):
        i += 1
        if seq_input != "" and chainS != "":
            secure_str(chainS)
            chainS = chainS[0]
            secure_str(seq_input)
            print( "mutations from sequence:", chainS, seq_input)
            add_mutations_from_sequence( mutations, seq_input, chainS, "seq" + str(i % 3), outdir+parent)
    for uniprot, chainU in zip(inputs["uniprots"], inputs["chainUs"]):
        i += 1
        if uniprot != "" and chainU != '':
            uni_file = outdir + uniprot
            download_uniprot( uniprot, uni_file)
            target = seq_from_fasta( uni_file)
            add_mutations_from_sequence( mutations, target, chainU, "uni" + str(i % 3), outdir + parent)

    print(__name__, 'total number of mutations:', len(mutations))
    if len(mutations) != 0:
        helper_files_from_mutations( mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile)
    else:
        print("no mutations")
        #return render_template("mutate.html", tag = tag, error = "Please provide a mutation")
        status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
        status = "no+mutations+defined+" + outdir + parent
        cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
        print(cmd)
        print(outdir +  "log.txt")
        log = open( outdir + "log.txt", 'a')
        bash_cmd(cmd,  log)
        return
    
    #start_thread(fixbb, [tag, parent, resfile, mutant, "log.txt"], "mutti") # fixbb sends all cmds to threads 
    fixbb(tag, parent, resfile, mutant, "log.txt")

    #add_energy(  outdir + mutant, outdir +  mutfile ) # part of fixbb now !!
   #TODO Define Chain    
    


    #out = app.config['USER_DATA_DIR'] + tag + "/"
    #log = open( out + "log.txt", 'a')
    #out_file_name = mutant
    #if out_file_name[-4:] == '.pdb':
    #    out_file_name = out_file_name[:-4]

    #cmd = "tsp " + "bash -i " + app.config['RASP_PATH'] + "calc-rasp.sh " + out + out_file_name + ".pdb " +  "A " + out_file_name
    #bash_cmd(cmd, log)
    #print("bash -i " + app.config['RASP_PATH'] + "calc-rasp.sh " + out + out_file_name + ".pdb " +  "A " + out_file_name)
    #os.system("bash -i " + app.config['RASP_PATH'] + "calc-rasp.sh " + out + out_file_name + ".pdb " +  "A " + out_file_name)



    send_email(outdir + "mail.txt")




@app.route('/mutate/<tag>', methods=['GET', 'POST'])
@app.route('/mutate/<tag>/<msg>', methods=['GET', 'POST'])
def mutate(tag,msg=""):

    if request.method == 'GET':
        outdir = app.config['USER_DATA_DIR'] + tag + "/"

        # get chains, resid-ranges from uploaded structure
        chains_range = get_chains_and_range( outdir + "structure.pdb")
        chains = ''
        for w in chains_range.split(",")[0:-1]:
            w = w.strip()
            if len(w) > 0:
                chains += w[0]
        print( 'chains: ', chains)

        # check if pdb in DB
        status = ""
        if msg=="found":
            status="Your PDB was found in our DB, no minimization will be performed."
        elif msg == "notfound":
            status = "Your PDB was not found in our DB, minimization will be performed."
            
        return render_template("mutate.html", tag = tag, chains=chains, chains_range=chains_range, status=status, error = "")


    ###  get form values

    inputs = {}

    inputs["mutations"] = request.form['mutations'].strip().replace(' ','').split(',')

    inputs["clustals"] = [
        request.files['clustal1'],
        request.files['clustal2'],
        request.files['clustal3'],
    ]

    inputs["fastas"] = [
        request.files['fasta1'],
        request.files['fasta2'],
        request.files['fasta3'],
    ]

    inputs["chainFs"] = [
        request.form['chainF1'],
        request.form['chainF2'],
        request.form['chainF3'],
    ]

    inputs["seq_inputs"] = [
        request.form['sequence1'].strip(),
        request.form['sequence2'].strip(),
        request.form['sequence3'].strip(),
    ]

    inputs["chainSs"] = [
        request.form['chainS1'],
        request.form['chainS2'],
        request.form['chainS3'],
    ]

    inputs["uniprots"] = [
        request.form['uniprot1'].strip(),
        request.form['uniprot2'].strip(),
        request.form['uniprot3'].strip(),
    ]

    inputs["chainUs"] = [
        request.form['chainU1'],
        request.form['chainU2'],
        request.form['chainU3'],
    ]

    # get all mutations
    i = 0
   # print(len(request.form['mutations'].strip()))
    if len(request.form['mutations'].strip()) != 0:
        i = i+1

    for clustal in inputs["clustals"]:
        if clustal.filename != "":
           # print("clustal")
            i = i+1

    for fasta, chainF in zip(inputs["fastas"], inputs["chainFs"]):
        if fasta.filename != "" and chainF != "":
            #print("fasta")
            i += 1
    for seq_input, chainS in zip(inputs["seq_inputs"], inputs["chainSs"]):
        if seq_input != "" and chainS != "":
           # print("seq")
            i += 1
    for uniprot, chainU in zip(inputs["uniprots"], inputs["chainUs"]):
        if uniprot != "" and chainU != '':
           # print("uniprot")
            i += 1

   # print("i " + str(i))
    if i == 0:
        print("no mutations")
        #
        status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
        status = "no+mutations+defined+" + os.path.join( app.config['USER_DATA_DIR'], tag) 
        cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
        print(cmd)
        print(os.path.join( app.config['USER_DATA_DIR'], tag)  +  "log.txt")
        log = open( os.path.join( app.config['USER_DATA_DIR'], tag)  + "log.txt", 'a')
        bash_cmd(cmd,  log)
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation")


    ### start calculations

    mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    start_thread(add_mutations, [tag, mutant, inputs], "add_muts") 
    

    return redirect(url_for('status', tag = tag, filename = mutant, msg="-"))


@app.route('/vcf', methods=['GET', 'POST'])
def vcf():
    if request.method == 'GET':
        return render_template("vcf.html", error = "")
  

    ### get form values

    vcf = request.files['vcf']
    vcf_file = vcf.filename
    if vcf_file == "":
        return render_template("vcf.html", error = "no filename was given")

    min_type = request.form['min-selector'] # = short | long
    longmin = (min_type == 'long')


    ### processing

    outdir, tag = create_user_dir()


    rasp_calculation = False
    rasp_checkbox = request.form.get('rasp-checkbox') # on none

    if(rasp_checkbox == "on"):
        rasp_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/rasp.status")
        with open(rasp_path, "w") as f:
            f.write(rasp_checkbox)
    
    print(rasp_calculation)

    vcf.save( outdir + vcf_file)

    # call run_vcf 
    cmd = "tsp " + app.config['SCRIPTS_PATH'] + "run_vcf.sh " + outdir + ' ' + vcf_file
    print(cmd)
    p = subprocess.check_output(cmd.split())
    print(p)
    if wait( outdir +  vcf_file[:-4] + '_missense.csv', 1, 900) == False:
        return render_template("vcf.html", error = "No missense was found.")
    
    # get mutations
    alphafold,mutations = mutations_from_vcf( outdir + vcf_file[:-4] + '_missense.csv')
    print( 'alphafold:', alphafold)

    # retrieve structure
    msg = ""
    file_path = outdir + "structure.pdb"
    if is_in_db( alphafold.strip().upper()):
        msg="found"
        cp_from_db(alphafold.strip().upper(),file_path)
    else:
        msg="notfound"
        download_file("https://alphafold.ebi.ac.uk/files/" + alphafold.strip() + "-model_v4.pdb", outdir + 'structure.pdb' )
    print(msg)
    if not is_pdb( outdir + 'structure.pdb'):
        error_message = "It was not possible to upload the AlphaFold model: " + alphafold + "<br>Currently only the first candidate can be uploaded"
        return render_template("vcf.html", error = error_message)

    # create info file
    os.mkdir(outdir + "info/")
    with open(outdir + "info/mut_0.txt", "w") as f:
        f.write("none\n")
    with open( outdir + "mut_0_resfile.txt", 'w') as w:
        w.write('NATAA\nstart\n')

    # relax structure
    if msg != "found":

        path = ""
        rose = app.config['ROSEMINT_PATH']
        path = rose + "alphafold/" + alphafold.strip().upper() + ".pdb"

        print("Store " + path)
        start_thread(fixbb, [tag, "structure.pdb", "mut_0_resfile.txt", "mut_0", "log.txt",longmin, path ], "minimisation")
    else:
        shutil.copyfile( outdir + "structure.pdb", outdir + "mut_0.pdb")
        rose = app.config['ROSEMINT_PATH']
        path = rose + "alphafold/" + alphafold.strip().upper() + ".pdb"

        print("path rasp: " + path)
        calc_rasp(tag, "structure.pdb", "mut_0", "log.txt", path )
        file_processing( tag, "structure.pdb", "mut_0", "log.txt" )
        file_processing( tag, "structure.pdb", "mut_0", "log.txt" )
    print( 'fixbb started for initial upload\n')

    if wait( outdir + 'mut_0.pdb', 1, 920) == False:
        return render_template("vcf.html", error = "Your structure could not be minimized.")

    # mutate
    helper_files_from_mutations( mutations,  outdir + 'mut_0.pdb',  outdir + 'mut_0_1_resfile.txt',  outdir + 'mut_0_1.clw',  outdir + 'info/mut_0_1.txt')
    start_thread(fixbb, [tag, 'mut_0.pdb', 'mut_0_1_resfile.txt', 'mut_0_1.pdb', "log.txt"], "mutti")

    
    return redirect(url_for('status', tag = tag, filename = "mut_0_1.pdb"))



def interface_calculation(outdir, tag, msg, filtered, pdb, af, mutant, clustal):
    # TODO: log

    # relax provided structure
    print("####### INT: relax structure")
    longmin = True 
    relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af)


    ### get mutations 

    print("####### INT: get mutations")
    mutations = []

    parent = "mut_0.pdb"
    resfile = mutant[:-4] + "_resfile.txt"
    align = mutant[:-4] + ".clw"
    mutfile = "info/" + mutant[:-4] + ".txt"

    # wait for parent file
    if wait(outdir + parent, 1, 900) == False:
        return

    add_mutations_from_alignment(mutations, clustal, outdir + parent)

    if len(mutations) == 0:
        return


    ### calculation

    # helpers
    print("####### INT: generate helper files")
    helper_files_from_mutations(mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile)

    # start mutation calculation
    print("####### INT: start fixbb")
    fixbb(tag, parent, resfile, mutant, "log.txt")

    # send email when done
    send_email(outdir + "mail.txt")

    print("####### INT: done")



@app.route('/interface', methods=["GET", "POST"])
def interface():
    if request.method == 'GET':
        return "<h1>Interface</h1>"


    ### get form values

    print("################# delta interface")

    # first file - "conv"
    upload = request.files['file_conv']
    pdb = secure_filename( request.form['pdb_conv'].strip() )
    af = "" # secure_filename( request.form['af_conv'].strip() )

    # second file - "super"
    upload_super = request.files['file_super']
    pdb_super = secure_filename( request.form['pdb_super'].strip() )
    af_super = "" # secure_filename( request.form['af_super'].strip() )

    # alignment TODO: generalize
    alignment_link = request.form.getlist('alignment_link')[0]
    alignment_link = "https://www.bioinfo.mpg.de/AlignMeBeta/work/" + alignment_link.split("work/")[1]



    print("######## got form values")


    ### processing

    outdir, tag = create_user_dir()

    # save alignment file
    clustal = outdir + "alignment.aln"
    req = requests.get(alignment_link) 
    with open(clustal, "w") as f:
        f.write(req.content)

    # save file
    file_path = outdir + "structure.pdb"    
    original_name, unsuccessful, error_message, msg = save_pdb_file(file_path, upload, pdb, af)
    if unsuccessful:
        # TODO: return error_message
        return "there was an error\n" + error_message

    # filter (remove) chains, heteroatoms
    chain_filter = ""
    remove_hets = False
    filtered = filter_structure(outdir, file_path, chain_filter, remove_hets)

    # create info file for mut_0
    os.mkdir(outdir + "info/")
    with open(outdir + "info/mut_0.txt", "w") as f:
        f.write("none\n")

    #create status file
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name.log")
    with open(status_path, "w") as f:
        f.write(get_current_time()+"+Start+Calculation\n")

    with open(name_path, "w") as f:
        f.write(original_name)

    # start calculation 
    mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    start_thread(interface_calculation, [outdir, tag, msg, filtered, pdb, af, mutant, clustal], "interface calc")

    return redirect(url_for('status', tag = tag, filename = mutant, msg="-"))


    

@app.route('/get_status/<tag>/<filename>')
def get_status(tag, filename):
    status = ""
    msg = ""
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")

    check_status = os.path.isfile(status_path)
    if check_status:
        status_file = open(os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log"), "r")
        msg = status_file.read()
        print(msg)
        msg = msg.replace("+", " ")
        msg = msg.replace("\n", "<br>")
        if("no mutations defined" in msg):
            print("exit status") 
            return jsonify({'done': True, 'status': "skip", 'message': "No mutation was defined"})




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
@app.route('/status/<tag>/<filename>/<msg>')
def status(tag, filename, msg=""):
    return render_template("status.html", tag = tag, filename = filename, msg=msg)

@app.route('/info/<tag>/<filename>')
def info(tag, filename):
    if tag.isdigit():
        path = app.config['USER_DATA_DIR'] + tag + "/info/"
    else:
        path = app.config['EXAMPLE_DIR'] + tag + "/info/"
    print( 'info', tag, filename)
    mutations = ""
    with open( path + filename + ".txt") as r:
        lines = r.readlines()
        parent = lines[0].strip()
        energy = lines[-1].strip()
        if len(lines) > 2:
            mutations += lines[1].strip()
            for i in range(2,len(lines)-1):
                mutations += ',' + lines[i].strip()
    if tag.isdigit():
        path = app.config['USER_DATA_DIR'] + tag + "/name.log"
    else:
        path = app.config['EXAMPLE_DIR'] + tag + "/name.log"
    name_file = open(path, "r")
    name = name_file.read()
    print(name)
    return render_template("info.html", tag = tag, parent=parent, mutations = mutations,  energy=energy, name = name)


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
    mut_tree = build_mutation_tree(out, tag, "none")
    print('explore::tree', mut_tree)
    structures = "<ul>" + build_list(mut_tree) + "</ul>"
    print('explore::tree:', structures)
    parent = ""
    mutations = ""
    if filename == "":
        filename = "mut_0.pdb"
    with open( out + tag + "/info/" + filename[:-4] + ".txt") as r:
        parent = r.readline().strip()
        energy = r.readline().strip()
        mutations += r.readline().strip()
        for l in r:
            mutations += ',' + l.strip()
    outdir = out + tag + "/"
    if parent == "none":
        chains = get_chains(outdir + filename)
    else:
        chains = get_chains( outdir + parent)
    #energy = get_energy (outdir + filename)
    print( __name__, filename , tag, chains)


    return render_template("explore.html", tag = tag, structures = structures, parent=parent, mutations = mutations, filename=filename , chains = chains, energy=energy)


@app.route('/explore/<tag>/<filename>', methods=['GET', 'POST'])
@app.route('/explore/<tag>/', methods=['GET', 'POST'])
def explore(tag, filename = ""):
    if request.method == 'GET':
        return load_explore_page(app.config['USER_DATA_DIR'], tag, filename)


    ### get form values

    mutations = request.form['mutations'].strip().replace(' ', '').split(',')
    parent = request.form['fname'].strip()
    if parent[-4:] != ".pdb":
        parent += '.pdb'
    mutant = name_mutation(app.config['USER_DATA_DIR'], parent, tag)
    
    print(__name__, 'original:', parent, 'novel mutant:', mutant) # TODO: OUTDIR ???


    ### mutate structure

    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    helper_files_from_mutations( mutations, outdir + parent, outdir + mutant[:-4] + '_resfile.txt', outdir + mutant[:-4] + '.clw', outdir + "info/" + mutant[:-4] + '.txt' ) 
    start_thread(fixbb, [tag, parent,  mutant[:-4] + '_resfile.txt', mutant, "log.txt"], "remutate")
    

    return redirect(url_for('status', tag = tag, filename = mutant))



# Testing molstar / mdsrv
@app.route('/molstar/<tag>')
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

@app.route('/tutorial')
def tutorial():
    return render_template("tutorial.html")


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
            cmd = ['zip','results' + tag + '.zip']
            for f in os.listdir('.'):
                if ".pdb" in f:
                    cmd.append(f)
                if ".clw" in f:
                    cmd.append(f)
                if ".csv" in f:
                    cmd.append(f)
            p = subprocess.check_output(cmd)

    return send_from_directory(path,filename)


@app.route('/fixbbtest')
def fbt():
    fixbb("/disk/user_data/mutation_explorer_gamma/863040/", "structure_1.pdb", ['A:1A'], "", "log")
    return render_template("faq.html")

"""
@app.route('/rasp/<tag>/<model>')
@app.route('/rasp/<tag>/<model>/<chain>')
def rasp( tag, model, chain = ""):
    outname =  app.config['USER_DATA_DIR'] + tag + 'rasp_' + chain + '.csv'
    if model[:-4] == '.pdb':
        model = model[:-4]
    cmd = app.config['RASP_PATH'] + 'rasp.py ' + app.config['USER_DATA_DIR'] + tag + '/' + model + '.pdb ' + chain + ' ' + outname
    bash_cmd( cmd, log)
    return jsonfy( ) # oder file schreiben?
    #while not file_exists( outname ):   needed?
    #    sleep(5s)
"""        


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
            f.write("none\n")

            
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


    
def mutations_from_alignment( clustal, parent, target_seq="", preselected_chain="" ):
    mutations = []


    ### read alignment, pdb sequence
    ali = read_clustal( clustal)
    #seqs = sequences( ali)
    chains = pdb2seq( parent)
    chain_in_pdb = preselected_chain
    ali_id = ''
    count = 0

    # find pdb sequence in alignment
    if chain_in_pdb:
        for a,b in ali.items():
            if chains[chain_in_pdb][0] == b:
                count += 1
                ali_id = a
                print( 'mutations(): matching sequences found:', a, chain_in_pdb)
    else:
        # if no chain was provided, last chain w/ matching sequence is chosen
        for c,v in chains.items():
            for a,b in ali.items():
                if v[0] == b:
                    count += 1
                    chain_in_pdb = c 
                    ali_id = a
                    print( 'mutations(): matching sequences found:', a, c)

    if count != 1:
        print( 'ERROR: mutations() found ',count,'identical matches')
    if count == 0:
        print('bye')
        exit(1)


    ### select target sequence
    aligned = target_seq
    if not aligned:
        # if no target sequence was provided, first non-identical sequence is chosen instead
        for a,b in ali.items():
            if a != ali_id:
                aligned = b
                break
    pdbseq = chains[chain_in_pdb][0]
    pdbres = chains[chain_in_pdb][1]


    ### find mutations
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
    # TODO: wait sollte nur in seperaten Threads aufgerufen werden, sonst muss der user auch waiten
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

def get_chains_and_range( fname):
    chains = {}
    with open( fname) as r:
        for l in r:
            if l[:4] == "ATOM":  # does not make sense: or l[:6] == "HETATM":
                c = chain(l)
                r = resid(l)
                if c not in chains.keys():
                    chains[c] = [r,r]
                else:
                    chains[c][1] = r
    chainstr = ""
    for c,v in chains.items():
        chainstr += c + ': ' + str(v[0]) + '-' + str(v[1]) + ', '
    #print( chainstr)
    return chainstr
    

def get_energy( fname):
    print( __name__, 'get', fname)
    with open( fname) as r:
        for l in r:
            if l[:4] == "pose":
                return float( l.split()[-1])
    print( 'WARNING: no energy found in', fname)
    return -0

def add_energy( fromfile, tofile):
    with open( tofile, 'a') as w:
        w.write( str( get_energy(fromfile)) + '\n')
        
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


def is_in_db( pdb):
    rose = app.config['ROSEMINT_PATH']
    return len(glob.glob( rose + 'pdb/' + pdb.upper() + '*.pdb')) > 0 or len(glob.glob( rose + 'fixbb/' + pdb.upper() + '*.pdb')) > 0 or len(glob.glob( rose + 'alphafold/' + pdb.upper() + '*.pdb')) > 0

def cp_from_db( pdb, outfile):
    rose = app.config['ROSEMINT_PATH']
    listig =  glob.glob( rose + 'fixbb/' + pdb.upper() + '*.pdb')
    listig.extend( glob.glob( rose + 'pdb/' + pdb.upper() + '*.pdb'))
    listig.extend( glob.glob( rose + 'alphafold/' + pdb.upper() + '*.pdb'))
    if len(listig) == 1:
        print( listig[0], outfile)
        shutil.copyfile(listig[0],outfile)
    else:
        mini = 1e9
        best = ''
        for l in listig:
            energy = get_energy( l)
            if energy < mini:
                best = l.strip()
                mini = energy
        shutil.copyfile(best,outfile)

  
def get_current_time():
   dt = datetime.datetime.now()
   x = dt.strftime("%Y-%m-%d+%H:%M:%S")

   return x
