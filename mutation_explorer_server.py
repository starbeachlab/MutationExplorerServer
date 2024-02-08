from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, jsonify
import os, random, subprocess, time
import threading
import requests
import glob
from collections import defaultdict
from werkzeug.utils import secure_filename
import shutil
import datetime

from scripts import resort_clustal
from scripts import check_pdb

# Mail Stuff
import smtplib, ssl



# fatal error messages
NO_MUTATIONS = "No valid mutations were defined"
RELAXATION_FAILED = "Relaxation of initial structure failed. Check your input PDB. Potentially unknown HETATMs (if this is the case, try using the 'remove hetatms' option in the 'Filter structure' button on the 'upload structure or model' page). Calpha only traces can not be handled. RaSP requires side-chains to be present."
MUTATION_FAILED = "Mutation failed"
STRUCTURE_NOT_IN_ALIGNMENT = "No sequence in the alignment matches the sequence of the provided structure"

READ_FAILED = 'Something has gone wrong reading a file. '
WRITE_FAILED = 'Something went wrong when writing to a file. '
READ_WRITE_FAILED = 'While reading or writing files, something went wrong. '
FILE_NOT_FOUND = 'File not found. '
UNEXPECTED = 'An unexpected error occured. '

STATUS_UPDATE_FAILED = 'Something went wrong during the status update.'
FIXBB_FAILED = 'Error in fixbb\n'
COPY_FAILED = 'Something went wrong while copying or moving a file.'
INTERFACE_SCORE_FAILED = 'There was an error during the interface score calculation.'
RASP_FAILED = 'There was an error during the RaSP calculation.'
RASP_CP_FAILED = 'Something went wrong while copying the RaSP files.'
SUPERIMOPSE_FAILED = 'There was an error during the superimposing of the structures.'
CONSERVATION_FAILED = 'There was an error during the conservation calculation.'
ENERGY_FILE_FAILED = 'There was an error while writing the energy to the b factor.'
ENERGY_DIFF_FAILED = 'There was an error while calculating the energy difference between two files.'
HYDRO_FILE_FAILED = 'There was an error while writing the hydrophobicity to the b factor.'
HYDRO_DIFF_FAILED = 'There was an error while calculating the hydrophobicity difference between two files.'
INTERFACE_SCORE_DIFF_FAILED = 'There wasn an error while calculating the difference of the interface score between two files.'
ENERGY_APPEND_FAILED = 'There was an error while handling the energy files.'
SCORE_FAILED = 'There was an error while scoring the structure.'
VCF_FAILED = 'There was an error during the vcf calculation.'


# wait times before error is returned (in seconds)
# could be useful to adjust for minimization option
WAIT_RELAXATION = 900
WAIT_MUTATION = 900
WAIT_VCF = 900
WAIT_SUPERIMPOSE = 900
TASKS = 2



class PrefixMiddleware(object):

    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):

        if environ['PATH_INFO'].startswith(self.prefix):
            environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
            return self.app(environ, start_response)
        else:
            start_response('404', [('Content-Type', 'text/plain')])
            return ["This url does not belong to the app.".encode()]



cfg_file = 'app.cfg'

app = Flask( __name__ )
app.config.from_pyfile( cfg_file )
app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/mutation_explorer')

## Setting the slots for tsp
cmd = "tsp -S " + str(TASKS)
print(cmd)
subprocess.run(cmd.split(), check=True, capture_output=True, text=True).stdout






def bash_cmd(cmd, tag):
    out = app.config['USER_DATA_DIR'] + tag + "/"

    logfile = "log.txt"
    idfile = "id.txt"
    log = open( out + logfile, 'a')
    wid = open( out + idfile, 'a')


    line = ""
    # Read line by line.
    try:
        with open(out + idfile, "r") as file:
            for line in file:
                pass
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND)
    except IOError:
        fatal_error(tag, READ_FAILED)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e)

    print("ID: " + line.strip())
    #p = subprocess.check_output(cmd.split())
    #log.write(str(p) + "\n")

    if any(char.isdigit() for char in str(line)):
        cmd = "tsp -D " + line.strip() + " " + cmd
    else:
        cmd = "tsp " + cmd
    print(cmd)
    log.write(cmd + "\n")
    pid = subprocess.run(cmd.split(), check=True, capture_output=True, text=True).stdout
    print( 'new process id:', pid)
    wid.write(pid)
    log.close()
    wid.close()
    return pid


def fatal_error(tag, msg):
    outdir = app.config['USER_DATA_DIR'] + tag + "/"

    with open(outdir + "fatal.log", "a") as f:
        f.write(msg)
    send_error_mail(tag)

    exit(1)


def status_update(tag, status, check_rasp=""):
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    status_path = outdir + "status.log"

    cmd = "bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
    if check_rasp != "":
        cmd = cmd + " check_rasp " + check_rasp

    bash_cmd(cmd, tag) 


def dev_status(tag, status):
    # for development purposes; to be removed (visible for user)
    status_update(tag, status)



@app.route('/')
def index():
    return render_template("home.html")


def secure_str( string):
    string = string.replace('.','').replace('/','').replace('*','').replace('?','').replace('!','').replace( ' ', '')
    if len(string) == 0:
        print( "WARNING:", __name__, 'empty string after cleanup!')
    return string

def download_file(url, file_path, tag):
    req = requests.get(url)

    try:
        with open(file_path, 'wb') as f:
            f.write(req.content)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND)
    except IOError:
        fatal_error(tag, WRITE_FAILED)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e)

    return req.status_code

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


def fixbb(tag, structure, resfile, out_file_name, logfile, longmin=False, path_to_store="", ifscore=""):
    # TODO: remove logfile


    print("######## fixbb  #######")
    print(tag)
    print('infile:', structure, "out_file_name: ", out_file_name)
    out = app.config['USER_DATA_DIR'] + tag + "/"

    if out_file_name[-4:] == '.pdb':
        out_file_name = out_file_name[:-4]

    #ext = out_file_name
    status = "Start+fixbb+for+" + out_file_name
    status_update(tag, status)



    # call rosetta
    cmd = "" + app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -use_input_sc -in:file:s " + out + structure + " -resfile " + out + resfile + ' -nstruct 1  -linmem_ig 10 -out:pdb -out:prefix ' + out
    if longmin == True:
        cmd += " -ex1 -ex2  "
    pid = bash_cmd(cmd, tag)



    print("rosetta done")       

    status = "fixbb+for+" + out_file_name + "+done"
    status_update(tag, status)

    sendTestTsp(tag, FIXBB_FAILED)
    print("give the threads some time to terminate")
    time.sleep(10)

    if(path_to_store != ""):
        cmd = "cp " + out + structure[:-4] + "_0001.pdb " + path_to_store
        pid = bash_cmd(cmd, tag)

    cmd = "mv " + out + structure[:-4] + "_0001.pdb " + out + out_file_name + ".pdb"
    bash_cmd(cmd, tag)
    print("MV Done")



    sendTestTsp(tag, COPY_FAILED)

    # TODO:: move one level up, should not be within fixbb or rename to fixbb_rasp
    calc_rasp(tag, structure,out_file_name,logfile, path_to_store)
    if ifscore != "":
        calc_interface( tag, out + out_file_name + ".pdb" , out+out_file_name + "_IF.pdb", ifscore)
    file_processing( tag, structure, out_file_name, logfile, ifscore)
    print("File processing")


def calc_interface( tag, in_file, out_file, parameter):
    cmd = "bash -i " + app.config['SCRIPTS_PATH'] + "pdb_rosetta_interface.sh " + in_file + " " + out_file + " " + parameter
    print(cmd)
    bash_cmd(cmd, tag)


    sendTestTsp(tag, INTERFACE_SCORE_FAILED)

    # should only be true once it is done!
    status = "Interface+calculation+for+" + in_file.split('/')[-1] + "+out+" + out_file.split('/')[-1] + "+done"
    status_update(tag, status)
    
def calc_interface_initial_structure(tag, outdir, minimize, ifscore=""):
    if ifscore != "":
        if minimize == "True":
            calc_interface(tag, outdir + "structure.pdb", outdir + "mut_0_IF.pdb", ifscore)
        else:
            calc_interface(tag, outdir + "mut_0.pdb", outdir + "mut_0_IF.pdb", ifscore)

def calc_rasp(tag, structure, out_file_name, logfile, path_to_store=""):
    out = app.config['USER_DATA_DIR'] + tag + "/"

    if(len(glob.glob(out+"rasp.status")) > 0):
        chains = get_chains(out + "structure.pdb", tag)
        chain_list = list(chains)

        for chain in chain_list:

            if(len(glob.glob( path_to_store + "_" + chain + ".csv" )) > 0 ):
                print("exists")
                listig =  glob.glob(  path_to_store + "_" + chain + ".csv" )

                if len(listig) == 1:
                    print( listig[0], out + "cavity_pred_" + out_file_name + "_" + chain + ".csv")
                    shutil.copyfile(listig[0],out + "cavity_pred_" + out_file_name + "_" + chain + ".csv")
            else:
                status = "Start+RaSP+calculation+for+" + out_file_name + "+with+chain+" + chain
                status_update(tag, status)


                cmd =  "bash -i " + app.config['RASP_PATH'] + "calc-rasp.sh " + out + out_file_name + ".pdb " + chain + " " + out_file_name + " " + out + " " + out + "rasp-error_" + chain + ".log" 
                print(cmd)
                bash_cmd(cmd, tag)

                sendTestTsp(tag, RASP_FAILED)
           

                if(path_to_store != ""):
                    cmd = "cp -d " + out + "cavity_pred_" + out_file_name + "_" + chain + ".csv " + path_to_store + "_" + chain + ".csv"
                    print(cmd)
                    bash_cmd(cmd, tag)
                    sendTestTsp(tag, RASP_CP_FAILED)
        
 

                status = "RaSP+calculation+for+" + out_file_name + "+with+chain+" + chain + "+done"
                status_update(tag, status)



def superimpose(tag, align_structure, align_chain, template_structure, template_chain, alignment):
    si = "python3 " + app.config['SCRIPTS_PATH'] + 'pdb_superimpose.py alignment: '

    tmp = align_structure[:-4] + '_tmp.pdb'
    cmd =  si + align_structure + ' ' + align_chain + ' ' + template_structure + ' ' + template_chain + ' ' + alignment + ' ' + tmp
    pid = bash_cmd(cmd, tag)

    # wait for superimposed structure
   # if not wait(tmp, 1, WAIT_SUPERIMPOSE):
    #return
    if(not waitID(pid)):
        error_message = "Superimposing failed!."
        fatal_error(tag, error_message)
        return


    os.rename(tmp, align_structure)



def calc_conservation(tag, structure, alignment, pdb_chain, seq_id, logfile):
    cons = "python3 " + app.config['SCRIPTS_PATH'] + 'pdb_conservation.py '

    cmd =  cons + structure + ' ' + pdb_chain + ' ' + alignment + ' ' + str(seq_id) + ' 0.2 ' + structure[:-4] + '_cons.pdb'
    pid = bash_cmd(cmd, tag)
    return pid
    

def mutant_calc_conservation(tag, structure, logfile, my_id):
    try:
        with open( logfile, 'a') as w:
            print('mutant_calc_conservation:', tag, structure, my_id, file=w)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + logfile)
    except IOError:
        fatal_error(tag, WRITE_FAILED + logfile)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + logfile)
        
    # wait for mutant structure  ### why in here? 
    print("ID Again " + my_id)
    if(not waitID(my_id)):
        error_message = "Conservation calculation failed due to previous errors. Most likely the RaSP or Rosetta interface score calculation failed. RaSP would fail on PDBs containing no side-chain atoms."
        fatal_error(tag, error_message)
        return


   # if(not wait(structure, 1, WAIT_MUTATION)):
   #     return


    # find mutated chain
    chains = pdb2seq(structure, tag)

    print( "mutant_calc_conservation, chains: ", chains.keys())

    cons = structure[:-4] + "_cons.pdb"
    # cp structure into tmp
    tmp = structure[:-4] + ".tmp"
    shutil.copy( structure, tmp)
    
    for c in chains:
        chain_alignment = structure[:-4] + '_' + c + '.clw'

        if not os.path.isfile(chain_alignment):
            continue

        cid = structure[:-4].split("/")[-1]

        sid = get_seq_id(chain_alignment, cid, tag)
        if sid is None:
            print("seq id not found")
            continue

        print( "mutant_calc_conservation:", c, sid )
        pid = calc_conservation(tag, structure, chain_alignment, c, sid, logfile)
        # cp structure + "cons.pdb" into structure
        if waitID(pid):
            print("start waiting for conserved file")
            shutil.copy(cons, structure)
        sendTestTsp(tag, CONSERVATION_FAILED)
     #   if wait( cons, 1, WAIT_SUPERIMPOSE):
     #       shutil.copy( cons, structure )
            
    # mv tmp back to structure, restore original state
    os.rename( tmp, structure) 


def file_processing( tag, structure, out_file_name, logfile, ifscore=""):
    # rename output file #### WRITE ENERGIES INSTEAD !!!!!
    bfac =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_to_bfactor.py "
    ediff = app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_diff.py "
    hydro = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_to_bfactor.py '
    hdiff = app.config['SCRIPTS_PATH'] + 'pdb_hydrophobicity_diff_to_bfactor.py '
    mutti = app.config['SCRIPTS_PATH'] + 'pdb_mutated_aa.py '
    
    out = app.config['USER_DATA_DIR'] + tag + "/"
    
    if out_file_name[-4:] == '.pdb':
        out_file_name = out_file_name[:-4]
    
    cmd = bfac + out + out_file_name + '.pdb ' + out + out_file_name + '_absE.pdb'
    bash_cmd(cmd, tag)

    if "mut" in structure:
        cmd = ediff + out + structure + ' ' + out + out_file_name + '.pdb ' + out + out_file_name + '_diffE.pdb'
        bash_cmd(cmd, tag)

    cmd = hydro + out +  out_file_name + '.pdb ' + out +  out_file_name + '_HyPh.pdb'
    bash_cmd(cmd, tag)

    cmd = hydro + out + structure + ' ' + out + structure[:-4] + '_HyPh.pdb' # TODO: wieso structure[:-4] statt out_file_name? # rene: gute frage, 
    bash_cmd(cmd, tag)
    
    if "mut" in structure:
        cmd = hdiff + out + structure + ' ' + out + out_file_name + '.pdb ' + out + out_file_name + '_diffHyPh.pdb'
        bash_cmd(cmd, tag)

        cmd = mutti + out + structure + " " + out + out_file_name + '.pdb ' + out + out_file_name + '_aa.pdb'
        my_id = bash_cmd(cmd, tag)
        print("Last id: " + my_id)

        start_thread(mutant_calc_conservation, [tag, out + out_file_name + '.pdb', logfile, my_id], "mutant conservation")

        if ifscore != "":
            #if wait( out + structure[:-4] + '_IF.pdb', 1, WAIT_MUTATION) and wait( out + out_file_name + "_IF.pdb", 1, WAIT_MUTATION):
            cmd =  app.config['SCRIPTS_PATH'] + "pdb_bfactor_diff.py " + out + structure[:-4] + "_IF.pdb " + out + out_file_name + "_IF.pdb " + out + out_file_name + "_diffIF.pdb"
            print(cmd)
            bash_cmd(cmd, tag)
    
    if not os.path.exists(out + "fin/"):
        os.mkdir(out + "fin/")

    cmd = "touch " + out + "fin/" + out_file_name + ".pdb"
    bash_cmd(cmd, tag)

    cmd =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_append.py " + out + out_file_name + ".pdb " + out + "info/" + out_file_name + ".txt"
    bash_cmd(cmd, tag)

    cmd =  app.config['SCRIPTS_PATH'] + "pdb_rosetta_energy_append.py " + out + out_file_name + ".pdb " + out + "info/" + out_file_name + ".txt"

    sendTestTsp(tag, ENERGY_DIFF_FAILED)


def sendTestTsp(tag, error):
    cmd = "bash"
    my_id = bash_cmd(cmd, tag)
    print("Check id: " + my_id)
    if(not waitID(my_id)):
        error_message = error
        fatal_error(tag, error_message)
        return

def create_user_dir():
    # generate tag
    while(True):
        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        if not os.path.exists(outdir):
            break
    os.mkdir(outdir)

    return outdir, tag

def save_pdb_file(file_path, upload, pdb, af, tag):
    original_name = ""
    error = False
    error_message = ""
    msg = "x"

    print( 'save pdb file:', file_path, pdb, af)
    if upload.filename != "":
        original_name = upload.filename
        upload.save(file_path)
    elif pdb != "":
        original_name = pdb
        if pdb[-4:] == ".pdb":
            pdb = pdb[:-4]
        if is_in_db( pdb):
            msg="found"
            cp_from_db(pdb,file_path, tag)
        else:
            msg="notfound"
            status = download_file("https://files.rcsb.org/download/" + pdb + ".pdb", file_path, tag)

            if status == 404:
                error = True
                error_message = 'The PDB id you provided was not found in the PDB database.'
                return original_name, error, error_message, msg

        print('pdb in db:', msg)
    elif af != "":
        af_id = get_alphafold_id(af)
        original_name = af_id
        print( 'alphafold:', af_id)
        if is_in_db( af): # @Daniel: sollte das nicht auch af_id sein?
            msg="found"
            cp_from_db(af,file_path, tag)
        else:
            msg="notfound"
            status = download_file("https://alphafold.ebi.ac.uk/files/" + af_id, file_path, tag)

            if status == 404:
                error = True
                error_message = 'The AlphaFold id you provided was not found in the AlphaFold database.'
                return original_name, error, error_message, msg
            
    else:
        # no structure -> error
        error = True
        error_message = "Please provide a structure" 

    if not is_pdb(file_path, tag):
        error = True
        error_message = "It was not possible to upload the structure you provided, or the PDB was empty. This may also indicate that it was a C-alpha only modell, which cannot be handled by MutationExplorer currently."
        #return render_template("submit.html", error = "It was not possible to upload the structure you provided.")

    # extract first model in PDB (NMR structures)  ### REPLACE THIS WITH CLEANUP SCRIPT! 
    if error == False:
        try:
            with open( file_path) as r, open( file_path[:-4] + '.tmp', 'w') as w:
                model_count = 0
                for l in r:
                    if l[:5] == "MODEL":
                        model_count += 1
                    if l[:6] == 'ENDMDL' or model_count > 1:
                        break
                    w.write(l)
                os.rename( file_path[:-4] + '.tmp', file_path )
                #order_pdb_resids( file_path[:-4] + '.tmp', file_path )
        except FileNotFoundError:
            fatal_error(tag, READ_WRITE_FAILED + FILE_NOT_FOUND + file_path)
        except IOError:
            fatal_error(tag, READ_WRITE_FAILED + file_path)
        except Exception as e:
            fatal_error(tag, READ_WRITE_FAILED + UNEXPECTED + e + ' ' + file_path)

    return original_name, error, error_message, msg


def filter_structure(tag, outdir, file_path, chain_filter, remove_hets):
    if chain_filter != "" or remove_hets:
        chains = chain_filter.split()

        try:
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
        except FileNotFoundError:
            fatal_error(tag, READ_WRITE_FAILED + FILE_NOT_FOUND + file_path + '/' + outdir + 'structure2.pdb')
        except IOError:
            fatal_error(tag, READ_WRITE_FAILED + file_path + '/' + outdir + 'structure2.pdb')
        except Exception as e:
            fatal_error(tag, READ_WRITE_FAILED + UNEXPECTED + e + ' ' + file_path + '/' + outdir + 'structure2.pdb')

        bash_cmd("mv " + outdir + "structure2.pdb " + file_path, tag) # why in queue? os.rename
        return True
    return False

def filter_chain( inpdb, fchain, outpdb, tag):
    atoms = []
    energies = []

    try:
        with open(inpdb, "r") as f_in:
            count = 0
            prev_resid = -999999
            min_id = 9999999
            max_id = -9999999
            good_boy = False
            for l in f_in:
                if l[:3] == "TER" and good_boy:
                    atoms.append(l)
                    good_boy  = False
                elif l[:4] == "ATOM" or l[:6] == "HETATM":
                    curr_resid = resid(l)
                    if curr_resid != prev_resid:
                        prev_resid = curr_resid
                        count += 1
                    if chain(l) == fchain:
                        atoms.append( l)
                        min_id = min( min_id, count)
                        max_id = max( max_id, count)
                        good_boy = True
                    else:
                        good_boy = False
                elif l[0] != '#' and len(l) > 20:
                    good_boy = False
                    c = l.split()
                    if 'label' in l and 'total' in l:
                        energies.append(l)
                    elif c[0].find('_') != -1:
                        energies.append(l)
                else:
                    good_boy = False
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + f_in)
    except IOError:
        fatal_error(tag, READ_FAILED + f_in)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + f_in)
    
    try:
        with open(outpdb, "w") as f_out:
            for l in atoms:
                f_out.write(l)
            count = 1
            for l in energies:
                if 'label' in l and 'total' in l:
                    f_out.write(l)
                    continue
                c = l.split()
                cid = c[0].rfind( '_' )
                myid = int(c[0][cid+1:])
                if myid >= min_id and myid <= max_id:
                    f_out.write( c[0][:cid+1] + str(count))
                    for i in range(1,len(c)):
                        f_out.write( ' ' + c[i])
                    f_out.write('\n')
                    count += 1
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + f_out)
    except IOError:
        fatal_error(tag, WRITE_FAILED + f_out)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + f_out)
            
def relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af, name, structure, ifscore=""):
    # name = "mut_0"
    # structure = "structure.pdb"
    
    resfile =  name + "_resfile.txt"   
    log_file = "log.txt"

    mutations_to_resfile( [] , outdir + resfile, tag)    

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
        start_thread(fixbb, [tag, structure, resfile, name, log_file, longmin, path], "minimisation")
    else:

        if (not filtered) and longmin == True:
            shutil.copyfile( outdir + structure, outdir + name + ".pdb")
            if(pdb != ""):
                rose = app.config['ROSEMINT_PATH']
                path = rose + "pdb/" + pdb.upper() + ".pdb"
            if(af != ""):    
                rose = app.config['ROSEMINT_PATH']
                path = rose + "alphafold/" + af.upper() + ".pdb"

            print("path rasp: " + path)
            calc_rasp(tag, structure, name, log_file, path ) # TODO

            # interface score for initial structure is going to be calculated later
            # because the selection for the chains was not done yet
            if ifscore != '':
                print( "calc interface from relax_initial_structure")
                calc_interface( tag, outdir + structure, outdir + name + "_IF.pdb", ifscore)

            file_processing( tag, structure, name,  log_file)

        else:
            print("Fixbb since filtering")
            start_thread(fixbb, [tag, structure, resfile, name, log_file, longmin], "minimisation")
            msg =  "notfound"

def score_structure(tag, outdir, name, structure, ifscore=''):
    logfile = "log.txt"

    status_update(tag, "Start+scoring+for+" + name)

    cmd =  app.config['ROSETTA_PATH'] + 'score_jd2.static.linuxgccrelease -in:file:s ' + outdir + structure + ' -nstruct 1 -out:pdb -out:prefix ' + outdir
    bash_cmd(cmd, tag)

    status_update(tag, "scoring+for+" + name + "+done")

    cmd = "mv " + outdir + structure[:-4] + '_0001.pdb ' + outdir + name + '.pdb'
    bash_cmd(cmd, tag)

    # TODO: rasp 
    # calc_rasp(tag, structure, name, logfile, path_to_store)
    # see fixbb / relax_initial_structure

    if ifscore != '':
        cmd = "bash -i " + app.config['SCRIPTS_PATH'] + 'pdb_rosetta_interface.sh ' + outdir + name + '.pdb ' + outdir + name + '_IF.pdb'
        bash_cmd( cmd, tag)
        status_update( tag, "interface+calculated+" + name + "_IF.pdb" )
    
    file_processing(tag, structure, name, logfile, ifscore)

@app.route('/submit', methods=['GET', 'POST'])
def submit(): 
    if request.method == 'GET':
        return render_template("submit.html", error = "")


    ### get form values

    upload = request.files['pdbfile']
    if upload and not allowed_file(upload.filename, {'pdb'}):
        error_message = 'You can only upload PDB files. Please try again with the right format.'
        return render_template('submit.html', error = error_message)

    pdb = secure_filename( request.form['pdbid'].strip() )
    af = secure_filename( request.form['alphafoldid'].strip() )

    chain_filter = request.form['chain_filter'].strip().upper()
    hetatom_filter = request.form.get('hetatom_filter')  # .get() needed for checkbox
    remove_hets = (hetatom_filter is not None)

    #email = False
    email = request.form['email'].strip()

    min_type = request.form['min-selector'] # = long | short | none
    minimize = (min_type != 'none')
    longmin = (min_type == 'long')


    ### processing

    outdir, tag = create_user_dir()

    rasp_calculation = False
    rasp_checkbox = request.form.get('rasp-checkbox') # on none

    if(rasp_checkbox == "on"):
        rasp_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/rasp.status")

        try:
            with open(rasp_path, "w") as f:
                f.write(rasp_checkbox)
        except FileNotFoundError:
            fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + rasp_path)
        except IOError:
            fatal_error(tag, WRITE_FAILED + rasp_path)
        except Exception as e:
            fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + rasp_path)
    
    print(rasp_calculation)

    # prewrite email (is sent seperately)
    if email:
        results_link = app.config["SERVER_URL"]+ url_for('explore', tag = tag, filename = "mut_0_1.pdb") 
        write_email(outdir + "mail.txt", email, results_link)

    # save file
    file_path = outdir + "structure.pdb"    
    original_name, unsuccessful, error_message, msg = save_pdb_file(file_path, upload, pdb, af, tag)
    if unsuccessful:
        print( error_message)
        return render_template("submit.html", error = error_message)
    
    # check pdb for HIS protonation
    if check_pdb.check_his_replacement(file_path):
        error_message = 'The PDB you have provided seems to contain HSE/HSD/HSP instead of HIS. Please replace HSE/HSD/HSP with HIS and try again.'
        return render_template('submit.html', error = error_message)
    # check pdb for missing chain id
    if check_pdb.check_chain_id_missing(file_path):
        error_message ='The pDB you have provided seems to be missing some chain ids. Please make sure all entries have a chain id.'
        return render_template('submit.html', error = error_message)    

    protype = protein_type(file_path, tag)
    print( 'protein type:',protype)
    if protype == "void":
        error_message = "The PDB you have provided seems to contain no protein."
        return render_template( "submit.html", error = error_message)
    elif protype =="calpha":
        error_message = "The PDB you have provided contains C-alpha atoms only. We currently cannot handle these."
        return render_template( "submit.html", error = error_message)
    elif protype == "backbone" and rasp_checkbox == "on":
        error_message = "The PDB you have provided contains the protein backbone only. RaSP requires side-chains to be present. Either upload a full atom modell or do not use RaSP."
        return render_template( "submit.html", error =  error_message)
    
    # filter (remove) chains, heteroatoms
    filtered = filter_structure(tag, outdir, file_path, chain_filter, remove_hets)

    print("submit\n")
    if os.path.isfile( file_path):
        print( file_path + ' is downloaded\n')

    # create info file for mut_0
    os.mkdir(outdir + "info/")

    try:
        with open(outdir + "info/mut_0.txt", "w") as f:
            f.write("none\n")
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + "info/mut_0.txt")
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + "info/mut_0.txt")
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + "info/mut_0.txt")

    #create status file
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name.log")
    # print("status")
    # print(status_path)

    try:
        with open(status_path, "w") as f:
            f.write(get_current_time()+"+Start+Calculation\n")
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + status_path)
    except IOError:
        fatal_error(tag, WRITE_FAILED + status_path)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + status_path)

    try:
        with open(name_path, "w") as f:
            f.write(original_name)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + status_path)
    except IOError:
        fatal_error(tag, WRITE_FAILED + status_path)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + status_path)
        
    # relax structure
    if minimize:
        relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af, "mut_0", "structure.pdb")
    else:
        score_structure(tag, outdir, "mut_0", "structure.pdb")

    print('fixbb started for initial upload\n')

    return redirect(url_for('mutate', tag = tag, msg=msg, minimize=minimize))




def add_mutations(tag, mutant, inputs, ifscore=""):
    outdir = app.config['USER_DATA_DIR'] + tag + "/"

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
    ### wait for mut_0.pdb to exist

    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, RELAXATION_FAILED + " from add_mutations()")

    # get all mutations
    i = 0
    for clustal_file, chainC in zip(inputs["clustal_files"], inputs["chainCs"]):
        if clustal_file != "":
            secure_str(chainC)
            chainC = chainC[0]
            add_mutations_from_alignment( mutations, clustal_file, outdir + parent, tag, base_chain=chainC)
    for fasta_file, chainF in zip(inputs["fasta_files"], inputs["chainFs"]):
        i += 1
        if fasta_file != "" and chainF != "":
            secure_str(chainF)
            chainF = chainF[0]
            head, target = seq_from_fasta( fasta_file, tag)
            add_mutations_from_sequence( mutations, target, chainF, "fa" + str(i % 3), outdir+parent, tag)
    for seq_input, chainS in zip(inputs["seq_inputs"], inputs["chainSs"]):
        i += 1
        if seq_input != "" and chainS != "":
            secure_str(chainS)
            chainS = chainS[0]
            secure_str(seq_input)
            print( "mutations from sequence:", chainS, seq_input)
            add_mutations_from_sequence( mutations, seq_input, chainS, "seq" + str(i % 3), outdir+parent, tag)
    for uniprot, chainU in zip(inputs["uniprots"], inputs["chainUs"]):
        i += 1
        if uniprot != "" and chainU != '':
            uni_file = outdir + uniprot
            print( 'get uniprot:', uniprot, uni_file, i)
            download_uniprot( uniprot, uni_file, tag)
            target = seq_from_fasta( uni_file, tag)
            print( 'mutations before:', mutations, 'target', target, chainU, outdir + parent)
            add_mutations_from_sequence( mutations, target, chainU, "uni" + str(i % 3), outdir + parent, tag)
            print( 'mutations after:', mutations)
            
    print(__name__, 'total number of mutations:', len(mutations))
    if len(mutations) != 0:
        helper_files_from_mutations( mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile, tag)
    else:
        fatal_error(tag, NO_MUTATIONS)  # second layer of check for no mutation

    fixbb(tag, parent, resfile, mutant, "log.txt", ifscore=ifscore)

    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, MUTATION_FAILED + " (fixbb) ")
    # wait for mutation
    #if wait(mutant, 1, WAIT_MUTATION) == False:
    #    fatal_error(tag, MUTATION_FAILED + " (1)")

@app.route('/mutate/<tag>', methods=['GET', 'POST'])
@app.route('/mutate/<tag>/<msg>/<minimize>', methods=['GET', 'POST'])
def mutate(tag,msg="",minimize="True"):

    outdir = app.config['USER_DATA_DIR'] + tag + "/"

    # get chains, resid-ranges from uploaded structure
    chains_range = get_chains_and_range( outdir + "structure.pdb", tag)

    try: 
        with open( outdir + 'chains.txt', 'w') as w:
            w.write( chains_range + '\n')
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + 'chains.txt')
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + 'chains.txt')
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + 'chains.txt')
        
    chains = ''
    for w in chains_range.split(",")[0:-1]:
        w = w.strip()
        if len(w) > 0:
            chains += w[0]
    print( 'chains: ', chains)

    structure = ""

    try:
        with open( outdir + 'name.log') as r:
            structure = r.readline().strip()
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + outdir + 'name.log')
    except IOError:
        fatal_error(tag, READ_FAILED + outdir + 'name.log')
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + outdir + 'name.log')

    # check if pdb in DB
    status = ""
    if msg=="found":
        status="PDB entry " + structure + " has already been minimized in our database. No additional minimization will be performed."
    elif msg == "notfound":
        status = "PDB entry " + structure + " was not previously minimized in our database. Minimization will be performed."

    if request.method == 'GET':
        return render_template("mutate.html", tag = tag, chains=chains, chains_range=chains_range, status=status, error = "")

    ###  get form values
    inputs = {}

    # manual mutations
    inputs["mutations"] = request.form['mutations'].strip().replace(' ','').split(',')

    # alignment
    inputs["clustals"] = [
        request.files['clustal1'],
        request.files['clustal2'],
        request.files['clustal3'],
    ]

    if len(inputs['clustals']) != 0:
        for clustal in inputs['clustals']:
            if len(clustal.filename) != 0 and not allowed_file(clustal.filename, {'aln', 'clw'}):
                error_message = 'You can only upload .aln or .clw files for the ClustalW files. Please try again with the right format.'
                return render_template("mutate.html", tag = tag, error = error_message)

    inputs["chainCs"] = [
        request.form.get('chainC1'),
        request.form.get('chainC2'),
        request.form.get('chainC3'),
    ]

    # target sequence
    inputs["fastas"] = [
        request.files['fasta1'],
        request.files['fasta2'],
        request.files['fasta3'],
    ]

    if len(inputs['fastas']) != 0:
        for fasta in inputs['fastas']:
            if len(fasta.filename) != 0 and not allowed_file(fasta.filename, {'fasta', 'fas', 'fa', 'fna', 'ffn', 'faa', 'mpfa', 'frn'}):
                error_message = 'You can only upload .fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, or .frn files for the FASTA files. Please try again with the right format.'
                return render_template("mutate.html", tag = tag, error = error_message)

    inputs["chainFs"] = [
        request.form.get('chainF1'),
        request.form.get('chainF2'),
        request.form.get('chainF3'),
    ]

    inputs["seq_inputs"] = [
        request.form['sequence1'].strip(),
        request.form['sequence2'].strip(),
        request.form['sequence3'].strip(),
    ]

    inputs["chainSs"] = [
        request.form.get('chainS1'),
        request.form.get('chainS2'),
        request.form.get('chainS3'),
    ]

    inputs["uniprots"] = [
        request.form['uniprot1'].strip(),
        request.form['uniprot2'].strip(),
        request.form['uniprot3'].strip(),
    ]

    inputs["chainUs"] = [
        request.form.get('chainU1'),
        request.form.get('chainU2'),
        request.form.get('chainU3'),
    ]

    # get type of interface score calculation
    if_score_option = request.form["ifscore"] # = none | manual | all
    if_score = ""

    if if_score_option == "none":
        if_score = ""
    elif if_score_option == "all":
        if_score = "all"
    else:
        # manual calculation
        # get information which chains are on the left and right side
        left = ''
        right = ''
        chains_range = open_chains_range_file(outdir, tag)
        chains = get_chain_letters(chains_range)

        for c in chains:
            chain = "chain_"+c
            chain_choice = request.form["chain_"+c] # left | right
            if chain_choice == "left":
                left += c
            else:
                right += c
        if_score = left + "_" + right
        print('manual interface definition:', if_score)


    email = request.form['email'].strip() 
    if email:
        results_link = app.config["SERVER_URL"]+ url_for('explore', tag = tag, filename = "mut_0_1.pdb") 
        write_email(outdir + "mail.txt", email, results_link)

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
        #status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
        #status = "no+mutations+defined+" + os.path.join( app.config['USER_DATA_DIR'], tag) 
        #cmd = "tsp bash " + app.config['SCRIPTS_PATH'] + "write-status.sh " + status + " " + status_path
        #print(cmd)
        #print(os.path.join( app.config['USER_DATA_DIR'], tag)  +  "log.txt")
        #log = open( os.path.join( app.config['USER_DATA_DIR'], tag)  + "log.txt", 'a')
        #bash_cmd(cmd,  log)
        return render_template("mutate.html", tag = tag, chains=chains, chains_range=chains_range, status=status, error = "Please provide a mutation")

    # files have to be saved in this function, out of reasons beyond my understanding
    inputs["clustal_files"] = []
    for clustal in inputs["clustals"]:
        if clustal.filename != "":
            clustal_file = os.path.join( outdir , clustal.filename )
            clustal.save( clustal_file)
            inputs["clustal_files"].append(clustal_file)
        else:
            inputs["clustal_files"].append("")

    inputs["fasta_files"] = []
    for fasta, chainF in zip(inputs["fastas"], inputs["chainFs"]):
        i += 1
        if fasta.filename != "" and chainF != "":
            fasta_file =  outdir + fasta.filename
            fasta.save(fasta_file) 
            inputs["fasta_files"].append(fasta_file)
        else:
            inputs["fasta_files"].append("")

        
    # calculate the interface score for the initial structure 
    # could not be done earlier because above selection was not known

    
    start_thread(calc_interface_initial_structure, [tag, outdir, minimize, if_score], "init_interface")    
    #calc_interface_initial_structure(tag, outdir, minimize, ifscore=if_score)

    ### start calculations
    mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
    start_thread(add_mutations, [tag, mutant, inputs, if_score], "add_muts")     
    return redirect(url_for('status', tag = tag, filename = mutant, msg="-"))

def vcf_calculation(tag, inputs):

    vcf = inputs["vcf"]
    print(vcf)
    vcf_file = vcf.split("/").pop()
    print(vcf_file)
    minimize = inputs["minimize"]
    longmin = inputs["longmin"]
    outdir = os.path.join( app.config['USER_DATA_DIR'], tag + "/")
    rasp_calculation = inputs["rasp_calculation"]

    ifscore_calculation = inputs['ifscore_calculation']
    ifscore = ''
    if ifscore_calculation:
        ifscore = 'all'

    # rasp.status
    if rasp_calculation:
        rasp_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/rasp.status")

        try:
            with open(rasp_path, "w") as f:
                #f.write(rasp_checkbox)
                f.write(str(rasp_calculation))
        except FileNotFoundError:
            fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + rasp_path)
        except IOError:
            fatal_error(tag, WRITE_FAILED + rasp_path)
        except Exception as e:
            fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + rasp_path)
    
    print(rasp_calculation)

    # call run_vcf 
    cmd =  app.config['SCRIPTS_PATH'] + "run_vcf.sh " + outdir +" " + vcf_file
    pid = bash_cmd(cmd, tag)
    print(outdir +  vcf_file[:-4] + '_missense.csv')
    # wait for missense file
    if(waitID(pid)) == False:
        return render_template("vcf.html", error="No missense was found.")
    # if wait( outdir +  vcf_file[:-4] + '_missense.csv', 1, WAIT_VCF) == False:
    #     return render_template("vcf.html", error = "No missense was found.")
    # print("test")
    # get mutations
    alphafold,mutations = mutations_from_vcf( outdir + vcf_file[:-4] + '_missense.csv', tag)
    print( 'alphafold:', alphafold)

    # retrieve structure
    msg = ""
    file_path = outdir + "structure.pdb"
    if is_in_db( alphafold.strip().upper()):
        msg="found"
        cp_from_db(alphafold.strip().upper(),file_path, tag)
    else:
        msg="notfound"
        download_file("https://alphafold.ebi.ac.uk/files/" + alphafold.strip() + "-model_v4.pdb", outdir + 'structure.pdb', tag)

    if not is_pdb(outdir + 'structure.pdb', tag):
        error_message = "It was not possible to upload the AlphaFold model: " + alphafold + "<br>Currently only the first candidate can be uploaded"
        fatal_error(tag, error_message)

    # create info file
    os.mkdir(outdir + "info/")

    try:
        with open(outdir + "info/mut_0.txt", "w") as f:
            f.write("none\n")
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + "info/mut_0.txt")
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + "info/mut_0.txt")
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + "info/mut_0.txt")
    
    try:
        with open( outdir + "mut_0_resfile.txt", 'w') as w:
            w.write('NATAA\nstart\n')
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + "mut_0_resfile.txt")
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + "mut_0_resfile.txt")
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + "mut_0_resfile.txt")


    # save file name
    name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name.log")

    try:
        with open(name_path, "w") as f:
            f.write(alphafold)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + name_path)
    except IOError:
        fatal_error(tag, WRITE_FAILED + name_path)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + name_path)

    status_update(tag, "Start+Calculation")

    # relax structure
    if minimize:
        if msg != "found":

            path = ""
            rose = app.config['ROSEMINT_PATH']
            path = rose + "alphafold/" + alphafold.strip().upper() + ".pdb"

            print("Store " + path)
            start_thread(fixbb, [tag, "structure.pdb", "mut_0_resfile.txt", "mut_0", "log.txt", longmin, path, ifscore], "minimisation")
            error_msg = "Minimization of initial structure using Rosettas fixbb failed. Check your input PDB."
            print( 'fixbb started for initial upload\n')
        else:
            shutil.copyfile( outdir + "structure.pdb", outdir + "mut_0.pdb")
            rose = app.config['ROSEMINT_PATH']
            path = rose + "alphafold/" + alphafold.strip().upper() + ".pdb"

            print("path rasp: " + path)
            calc_rasp(tag, "structure.pdb", "mut_0", "log.txt", path )
            #file_processing( tag, "structure.pdb", "mut_0", "log.txt" ) #rene: warum war das 2x hier?
            if ifscore_calculation:
                calc_interface(tag, outdir + "structure.pdb", outdir + "mut_0" + "_IF.pdb", 'all')
            file_processing( tag, "structure.pdb", "mut_0", "log.txt", 'all') # TODO: fehler ?
            error_msg = "Deriving energies and writing them into a PDB for visualization failed."
    else:
        score_structure(tag, outdir, "mut_0", "structure.pdb", ifscore)

    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, RELAXATION_FAILED)

    # wait for relaxation
    #if wait( outdir + 'mut_0.pdb', 1, WAIT_RELAXATION) == False:
    #    fatal_error(tag, RELAXATION_FAILED)

    # mutate
    helper_files_from_mutations( mutations,  outdir + 'mut_0.pdb',  outdir + 'mut_0_1_resfile.txt',  outdir + 'mut_0_1.clw',  outdir + 'info/mut_0_1.txt', tag)
    start_thread(fixbb, [tag, 'mut_0.pdb', 'mut_0_1_resfile.txt', 'mut_0_1.pdb', "log.txt", False, '', ifscore], "mutti") # thread isnt needed anymore

    # check if mutation successful
    # if wait( outdir + 'mut_0_1.pdb', 1, WAIT_MUTATION) == False:
    #     fatal_error(tag, MUTATION_FAILED + " (2)")
    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, MUTATION_FAILED + " (fixbb in vcf pipeline)")

@app.route('/vcf', methods=['GET', 'POST'])
def vcf():
    if request.method == 'GET':
        return render_template("vcf.html", error = "")
  
    ### get form values

    outdir, tag = create_user_dir()

    inputs = {}

    vcf = request.files['vcf']

    if vcf and not allowed_file(vcf.filename, {'vcf'}):
        error_message = 'You can only upload .vcf files. Please try again with the right format.'
        return render_template('vcf.html', error = error_message)

    vcf_file = vcf.filename
    if vcf_file == "":
        return render_template("vcf.html", error = "no filename was given")

    vcf.save( outdir + vcf_file)
    inputs["vcf"] = outdir + vcf_file


    min_type = request.form['min-selector'] # = long | short | none
    inputs["minimize"] = (min_type != 'none')
    inputs["longmin"] = (min_type == 'long')

    email = request.form['email'].strip()

    # prewrite email (is sent seperately)
    if email:
        results_link = app.config["SERVER_URL"]+ url_for('explore', tag = tag, filename = "mut_0_1.pdb") 
        write_email(outdir + "mail.txt", email, results_link)

    rasp_calculation = False
    rasp_checkbox = request.form.get('rasp-checkbox') # on none
    inputs["rasp_calculation"] = (rasp_checkbox == 'on')

    ifscore_checkbox = request.form.get('ifscore_checkbox')
    inputs['ifscore_calculation'] = (ifscore_checkbox == 'on')

    ### start calculation

    start_thread(vcf_calculation, [tag, inputs], "vcf calc")

    return redirect(url_for('status', tag = tag, filename = "mut_0_1.pdb"))

def interface_one_structure(tag, mutant, inputs):
    # TODO: rasp_calculation (option to not do rasp calculation)

    outdir = app.config['USER_DATA_DIR'] + tag + '/'

    clustal = inputs["clustal"]
    base_clustal_id = inputs["base_clustal_id"]
    target_clustal_id = inputs["target_clustal_id"]

    pdb = inputs["base_pdb"]
    af = inputs["base_af"]
    chain = inputs["base_chain"]

    filtered = inputs["base_filtered"]
    msg = inputs["base_msg"]

    minimize = inputs["minimize"]
    longmin = inputs["longmin"]

    ifscore_calculation = inputs['ifscore_calculation']
    ifscore = ''
    if ifscore_calculation:
        ifscore = 'all'

    if minimize:
        # relax provided structure
        relax_initial_structure(outdir, tag, msg, filtered, longmin, pdb, af, "mut_0", "structure.pdb", ifscore)
        error_msg="Relaxation of initial structure failed. Check PDB."
    else:
        score_structure(tag, outdir, "mut_0", "structure.pdb", ifscore)
        error_msg="Scoring of initial structure failed. Check PDB."

    ### get mutations 

    parent = "mut_0.pdb"
    resfile = mutant[:-4] + "_resfile.txt"
    align = mutant[:-4] + ".clw"
    mutfile = "info/" + mutant[:-4] + ".txt"

    #inputs["connector_string"] = parent + ":" + align + "," + base_clustal_id + "," + base_chain + ";" + mutant + ":" + align + "," + target_clustal_id + "," + target_chain

    # wait for mut_0.pdb
    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, error_msg)

    # wait for mut_0.pdb
    #if wait(outdir + parent, 1, WAIT_RELAXATION) == False:
    #    fatal_error(tag, RELAXATION_FAILED)


    mutations, noncanonical_residues = mutations_from_alignment(clustal, outdir + parent, tag, base_clustal_id=base_clustal_id, target_clustal_id=target_clustal_id, base_chain=chain)

    if len(mutations) == 0:
        fatal_error(tag, "NO MUTATIONS could be extracted from alignment. Identical sequences or no aligned positions?")

    if noncanonical_residues:
        status_update(tag, "mut_0.pdb containes noncanonical residues, which are ignored")


    ### calculation

    # helpers
    helper_files_from_mutations(mutations, outdir + parent, outdir + resfile, outdir + align, outdir + mutfile, tag)

    # start mutation calculation
    fixbb(tag, parent, resfile, mutant, "log.txt", ifscore=ifscore)

    # check mutation success
    if wait(mutant, 1, WAIT_MUTATION) == False:
        fatal_error(tag, MUTATION_FAILED + " (fixbb in AlignMe interface using single PDB)")

    # check mutation success
    pid = getLastID(outdir, tag)

    if not waitID(pid):
        fatal_error(tag, MUTATION_FAILED + " (3)")
        
def interface_two_structures(tag, inputs):
    # TODO: errors

    outdir = app.config['USER_DATA_DIR'] + tag + '/'

    clustal = inputs["clustal"]
    base_clustal_id = inputs["base_clustal_id"]
    target_clustal_id = inputs["target_clustal_id"]

    base_pdb = inputs["base_pdb"]
    base_af = inputs["base_af"]
    base_chain = inputs["base_chain"]

    base_filtered = inputs["base_filtered"]
    base_msg = inputs["base_msg"]

    target_pdb = inputs["target_pdb"]
    target_af = inputs["target_af"]
    target_chain = inputs["target_chain"]

    target_filtered = inputs["target_filtered"]
    target_msg = inputs["target_msg"]

    minimize = inputs["minimize"]
    longmin = inputs["longmin"]

    ifscore_calculation = inputs['ifscore_calculation']
    ifscore = ''
    if ifscore_calculation:
        ifscore = 'all'
        status_update( tag,  'calc+interface+score')

    #inputs["connector_string"] = "mut_0.pdb:mut_0_1.clw," + base_clustal_id + "," + base_chain + ";mut_1.pdb:mut_0_1.clw," + target_clustal_id + "," + target_chain
    
    if base_chain != '':
        status_update( tag, "filter+base+chain")
        filter_chain( outdir + "structure.pdb", base_chain, outdir + "tmp.pdb", tag)
        os.rename( outdir+ "tmp.pdb", outdir + "structure.pdb" )
    if target_chain != '':
        status_update( tag, "filter+target+chain")
        filter_chain( outdir + "structure2.pdb", target_chain, outdir + "tmp2.pdb", tag)
        os.rename( outdir+ "tmp2.pdb", outdir + "structure2.pdb" )
    
    if not minimize:
        status_update( tag, "score+base+structure")
        score_structure(tag, outdir, "mut_0", "structure.pdb", ifscore)
        status_update( tag, "score+target+structure")
        score_structure(tag, outdir, "mut_1", "structure2.pdb", ifscore)
    else:
        # relax provided structure
        status_update( tag, "relax+base+structure")  
        # minimize structures
        relax_initial_structure(outdir, tag, base_msg, base_filtered, longmin, base_pdb, base_af, "mut_0", "structure.pdb", ifscore)

        pid = getLastID(outdir, tag)

        if not waitID(pid):
            fatal_error(tag, "Relaxation of base structure failed. Does the selected chain has ATOMs in the PDB?")

        #if(not wait(outdir + "mut_0.pdb", 1, WAIT_RELAXATION)):
        #    fatal_error(tag, RELAXATION_FAILED)
        status_update( tag, "relax+target+structure")  

        relax_initial_structure(outdir, tag, target_msg, target_filtered, longmin, target_pdb, target_af, "mut_1", "structure2.pdb", ifscore)

        pid = getLastID(outdir, tag)

        if not waitID(pid):
            fatal_error(tag, "Relaxation of target structure failed. Does the selected chain has ATOMs in the PDB?")

        #if(not wait(outdir + "mut_1.pdb", 1, WAIT_RELAXATION)):
        #    fatal_error(tag, RELAXATION_FAILED)

    base_strc = outdir + "mut_0.pdb"
    target_strc = outdir + "mut_1.pdb"

    dev_status( tag, "find base pdb in alignment")

    # get clustal ids, chains, sequence ids
    pdb_match, base_noncanonical_residues = find_pdb_in_alignment(clustal, base_strc, tag, chain=base_chain, clustal_id=base_clustal_id)
    if pdb_match == ["",""]:
        fatal_error(tag, STRUCTURE_NOT_IN_ALIGNMENT)
    base_clustal_id, base_chain = pdb_match

    if base_noncanonical_residues:
        status_update(tag, "mut_0.pdb containes noncanonical residues, which are ignored")

    dev_status( tag, "get base seq id in clustal")    
    base_seq_id = get_seq_id(clustal, base_clustal_id, tag)
    if base_seq_id is None:
        status_update( tag, "no+base+seq+in+alignment")
        return

    status_update( tag, 'base+ids+matched')
    
    dev_status( tag, "find target pdb in alignment")

    pdb_match, target_noncanonical_residues = find_pdb_in_alignment(clustal, target_strc, tag, chain=target_chain, clustal_id=target_clustal_id)
    if pdb_match == ["",""]:
        status_update( tag, "structure not in alignment")
        fatal_error(tag, STRUCTURE_NOT_IN_ALIGNMENT)
    target_clustal_id, target_chain = pdb_match

    match_string = base_clustal_id + ',' + base_chain + ',' + target_clustal_id + ',' + target_chain
    write_alignment_match(outdir, match_string, tag)

    if target_noncanonical_residues:
        status_update(tag, "mut_1.pdb containes noncanonical residues, which are ignored")

    dev_status( tag, "get target seq id in clustal")    
    target_seq_id = get_seq_id(clustal, target_clustal_id, tag)
    if target_seq_id is None:
        status_update( tag, "no+target+seq+in+alignment")
        return

    if base_seq_id == target_seq_id:
        fatal_error(tag, NO_MUTATIONS)

    status_update( tag, 'target+ids+matched')
    status_update( tag, 'superimpose')
    dev_status(tag, "superimpose")

    # superimpose
    if base_seq_id < target_seq_id:
        superimpose(tag, base_strc, base_chain, target_strc, target_chain, clustal)
    else:
        superimpose(tag, target_strc, target_chain, base_strc, base_chain, clustal)

    dev_status(tag, "calc conservation")
    status_update(tag, 'calc+conservation')

    # calculate conservation
    calc_conservation(tag, base_strc, clustal, base_chain, base_seq_id, "log.txt")

    calc_conservation(tag, target_strc, clustal, target_chain, target_seq_id, "log.txt")

    dev_status(tag, "finished calculation")

    
@app.route('/interface_post', methods=["GET", "POST"])
def interface_post():
    if request.method == 'GET':
        redirect(url_for('index'))


    # alignment   TODO: generalize
    alignment_link = request.form.getlist('alignment_link')[0]
    alignment_link = "https://www.bioinfo.mpg.de/AlignMeBeta/work/" + alignment_link.split("work/")[1]

    outdir, tag = create_user_dir()

    # save alignment file
    clustal = outdir + "alignment.aln"
    req = requests.get(alignment_link)
    print( "ali:\n", req.content)
    content = req.content.decode('ASCII')
    print( "with conservation:\n", add_conservation(content))

    try:
        with open(clustal, "w") as f:
            f.write(add_conservation(content)) # add conservation because it is needed for spheres in mol*
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + clustal)
    except IOError:
        fatal_error(tag, WRITE_FAILED + clustal)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + clustal)

    return redirect(url_for('interface', tag = tag))


@app.route('/interface/<tag>', methods=["GET", "POST"])
def interface(tag):
    if request.method == 'GET':

        # get sequences from alignment
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        ali = read_clustal(outdir + "alignment.aln", tag)
        seqs = ",".join(ali.keys())

        return render_template("interface.html", seqs = seqs, error = "")

    outdir = app.config['USER_DATA_DIR'] + tag + "/"


    ### get form values

    inputs = {}

    inputs["clustal"] = outdir + "alignment.aln"

    # clustal id selection
    inputs["base_clustal_id"] = request.form.get("base_seq")
    inputs["target_clustal_id"] = request.form.get("target_seq")

    # base file
    inputs["base_upload"] = request.files['base_pdbfile']
    inputs["base_pdb"] = secure_filename( request.form['base_pdbid'].strip() )
    inputs["base_af"] = secure_filename( request.form['base_alphafoldid'].strip() )

    inputs["base_chain"] = secure_filename( request.form['base_chain'].strip() )

    # target file
    inputs["target_upload"] = request.files['target_pdbfile']
    inputs["target_pdb"] = secure_filename( request.form['target_pdbid'].strip() )
    inputs["target_af"] = secure_filename( request.form['target_alphafoldid'].strip() )

    inputs["target_chain"] = secure_filename( request.form['target_chain'].strip() )

    # options
    min_type = request.form['min-selector'] # = long | short | none
    inputs["minimize"] = (min_type != 'none')
    inputs["longmin"] = (min_type == 'long')

    rasp_checkbox = request.form.get('rasp-checkbox') # = on | none
    inputs["rasp_calculation"] = (rasp_checkbox == 'on')

    ifscore_checkbox = request.form.get('ifscore_checkbox')
    inputs['ifscore_calculation'] = (ifscore_checkbox == 'on')

    ### processing

    # save base file 
    base_file_path = outdir + "structure.pdb"    
    base_original_name, unsuccessful, error_message, base_msg = save_pdb_file(base_file_path, inputs["base_upload"], inputs["base_pdb"], inputs["base_af"], tag)
    if unsuccessful:
        ali = read_clustal(outdir + "alignment.aln", tag)
        seqs = ",".join(ali.keys())

        return render_template("interface.html", seqs = seqs, error = error_message)

    inputs["base_original_name"] = base_original_name
    inputs["base_msg"] = base_msg
    inputs["base_filtered"] = False

    # save target file 
    target_file_path = outdir + "structure2.pdb"    
    target_original_name, unsuccessful, error_message, target_msg = save_pdb_file(target_file_path, inputs["target_upload"], inputs["target_pdb"], inputs["target_af"], tag)
    target_given = not unsuccessful

    inputs["target_original_name"] = target_original_name
    inputs["target_msg"] = target_msg
    inputs["target_filtered"] = False


    print("status update")
    status_update(tag, "Start+Calculation")

    # save original filename
    name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name.log")

    try:
        with open(name_path, "w") as f:
            f.write(base_original_name)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + name_path)
    except IOError:
        fatal_error(tag, WRITE_FAILED + name_path)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + name_path)

    if target_given:
        name_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/name2.log")

        try:
            with open(name_path, "w") as f:
                f.write(target_original_name)
        except FileNotFoundError:
            fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + name_path)
        except IOError:
            fatal_error(tag, WRITE_FAILED + name_path)
        except Exception as e:
            fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + name_path)

    print("info files")

    # create info file for mut_0, (mut_1)
    if not os.path.isdir( outdir + 'info/'):
        os.mkdir(outdir + "info/")

    try:
        with open(outdir + "info/mut_0.txt", "w") as f:
            f.write("none\n")
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + "info/mut_0.txt")
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + "info/mut_0.txt")
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + "info/mut_0.txt")

    if target_given:
        #print("target given")
        dev_status(tag, "Target+given")

        try:
            with open(outdir + "info/mut_1.txt", "w") as f:
                f.write("none\n")
        except FileNotFoundError:
            fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + "info/mut_1.txt")
        except IOError:
            fatal_error(tag, WRITE_FAILED + outdir + "info/mut_1.txt")
        except Exception as e:
            fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + "info/mut_1.txt")

 
    # start calculation 
    if not target_given:
        mutant = name_mutation(app.config['USER_DATA_DIR'], "mut_0", tag)
        start_thread(interface_one_structure, [tag, mutant, inputs], "interface calc one structure")
        #connector_string = inputs["connector_string"]  # needed for load pdb in explore.html
        #print( 'connector string:', connector_string)
        return redirect(url_for('status', tag = tag, filename = mutant, msg="")) #, connector_string = connector_string ))
    

    start_thread(interface_two_structures, [tag, inputs], "interface calc two structures")
    #connector_string = inputs["connector_string"]  # needed for load pdb in explore.html
    #print( 'connector string:', connector_string)
    return redirect(url_for('status', tag = tag, filename = "mut_1.pdb", msg="" )) #, connector_string = connector_string ))
    


    

@app.route('/get_status/<tag>/<filename>')
def get_status(tag, filename):
    status = ""
    msg = ""
    status_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/status.log")
    fatal_path = os.path.join( app.config['USER_DATA_DIR'], tag + "/fatal.log")

    check_status = os.path.isfile(status_path)
    if check_status:
        status_file = open(status_path)
        msg = status_file.read()
        # print(msg)
        msg = msg.replace("+", " ")
        msg = msg.replace("\n", "<br>")
        if("no mutations defined" in msg):
            print("exit status") 
            return jsonify({'done': True, 'status': "skip", 'message': "No mutation was defined"})


    if os.path.isfile(fatal_path):
        msg = open(fatal_path, "r").read()
        return jsonify({'done': True, 'status': 'skip', 'message': msg})



    dirname = os.path.join( app.config['USER_DATA_DIR'], tag + "/fin/" + filename)
    done = os.path.isfile(dirname)
    if done:
        dirname = os.path.join( app.config['USER_DATA_DIR'], tag + "/" + filename )
        successful = os.path.isfile(dirname)
        if successful:
            status = "done"
            msg = ""
            print("send mail?")
            send_email(os.path.join( app.config['USER_DATA_DIR'], tag + "/mail.txt"))
        else:
            status = "error"
            msg = "Rosetta calculation encountered an error"
    #print( 'get_status', dirname, str(done))
    return jsonify({'done': done, 'status': status, 'message': msg})


@app.route('/status/<tag>/<filename>')
@app.route('/status/<tag>/<filename>/')
@app.route('/status/<tag>/<filename>/<msg>')
#@app.route('/status/<tag>/<filename>/<msg>')
#@app.route('/status/<tag>/<filename>/<msg>/<connector_string>')
def status(tag, filename, msg=""): #,connector_string=""):
    return render_template("status.html", tag = tag, filename = filename, msg=msg) #, connector_string=connector_string)

@app.route('/info/<tag>/<filename>')
@app.route('/info/<tag>/<filename>/<two_structures>')
def info(tag, filename, two_structures=""):
    if tag.isdigit():
        path = app.config['USER_DATA_DIR'] + tag + "/info/"
    else:
        path = app.config['EXAMPLE_DIR'] + tag + "/info/"
    print( 'info', tag, filename)
    mutations = ""

    try:
        with open( path + filename + ".txt") as r:
            lines = r.readlines()
            parent = lines[0].strip()
            energy = str(round(float(lines[-1].strip()),2))
            if len(lines) > 2:
                mutations += lines[1].strip()
                for i in range(2,len(lines)-1):
                    mutations += ',' + lines[i].strip()
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + path + filename + ".txt")
    except IOError:
        fatal_error(tag, READ_FAILED + path + filename + ".txt")
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + path + filename + ".txt")

    ediff = '-'
    if filename.count('_') > 1:

        try:
            with open( path + filename[:-2] + ".txt") as r:
                lines = r.readlines()
                ediff = str( round( float(lines[-1]) - float(energy), 2) )
        except FileNotFoundError:
            fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + path + filename[:-2] + ".txt")
        except IOError:
            fatal_error(tag, READ_FAILED + path + filename[:-2] + ".txt")
        except Exception as e:
            fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + path + filename[:-2] + ".txt")
        
                
    name_file = "/name.log"
    if two_structures != "":
        name_file = "/name2.log"
    if tag.isdigit():
        path = app.config['USER_DATA_DIR'] + tag + name_file
    else:
        path = app.config['EXAMPLE_DIR'] + tag + name_file
    name_file = open(path, "r")
    name = name_file.read()
    print(name)

    return render_template("info.html", tag = tag, parent=parent, mutations = mutations,  energy=energy, name = name, diff = ediff)


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

def load_explore_page(out, tag, filename):  #, connector_string = ""):
    # TODO: rename out (out should contain tag)
    print('load_explore_page', tag, filename)
    mut_tree = build_mutation_tree(out, tag, "none")
    print('explore::tree', mut_tree)
    structures = "<ul>" + build_list(mut_tree) + "</ul>"
    print('explore::tree:', structures)
    parent = ""
    mutations = ""
    if filename == "":
        filename = "mut_0.pdb"
        print( 'filename now:', filename)

    try:
        with open( out + tag + "/info/" + filename[:-4] + ".txt") as r:
            parent = r.readline().strip()
            energy = r.readline().strip()
            mutations += r.readline().strip()
            for l in r:
                mutations += ',' + l.strip()
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + out + tag + "/info/" + filename[:-4] + ".txt")
    except IOError:
        fatal_error(tag, READ_FAILED + out + tag + "/info/" + filename[:-4] + ".txt")
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + out + tag + "/info/" + filename[:-4] + ".txt")

    outdir = out + tag + "/"

    if parent == "none":
        chains = get_chains(outdir + filename, tag)
    else:
        chains = get_chains( outdir + parent, tag)
    
    chains_range = open_chains_range_file(outdir, tag)

    #energy = get_energy (outdir + filename)
    print( __name__, filename , tag, chains)

    two_structures = os.path.isfile(out + tag + "/mut_1.pdb")
    match_string = ''
    # assuming that this is only true when coming from AlignMe ...
    if two_structures:
        chains += get_chains( out + tag + "/mut_1.pdb", tag)
        two_structures = get_alignment_ids( out + tag + "/alignment.aln", tag)
        match_string = read_txt_file(outdir, 'alignment.txt', tag)
    else:
        two_structures = ''

    print("###############")
    print("###############")
    print(out + "mut_1.pdb")
    return render_template("explore.html", tag = tag, structures = structures, parent=parent, mutations = mutations, filename=filename , chains = chains, energy=energy, two_structures = two_structures, match_string = match_string, chains_range = chains_range) #, connector_string = connector_string)

@app.route('/chain_resids_sorted', methods=['POST'])
def chain_resids_sorted():
    out = app.config['USER_DATA_DIR']

    data = request.get_json()
    chain = data['chain']
    pdb = data['pdbFile']
    tag = data['tag']
    # extraxt clustal file name without suffix
    clustal = data['clustalUrl'].split('/')[-1].split('.')[0]

    pdb_file = f'{out}/{tag}/{pdb}.pdb'
    clustal_file_old = f'{out}/{tag}/{clustal}.clw'
    clustal_file_new = f'{out}/{tag}/{clustal}_reordered.clw'

    # checks if for this chain a reordered clustal file already exists
    if os.path.exists(clustal_file_new):
        print(f'Reordered clustal file already exists for chain {chain}')
        return jsonify({'done': 'done'})

    # checks if the resids for this chain are sorted
    sorted = resort_clustal.chain_resids_sorted(chain, pdb_file)

    if sorted is None:
        # error during check
        print(f'An error occurred during the pdb file sort check and subsequent file reordering.')
        print('Copying the original file')
        shutil.copyfile(clustal_file_old, clustal_file_new)
    elif sorted:
        # sorted
        print(f'Chain {chain} is sorted.')
        print('Copying clustal file')
        shutil.copyfile(clustal_file_old, clustal_file_new)
    else:
        # not sorted
        print(f'Chain {chain} is not sorted.')
        print('Reordering clustal file')
        resort_clustal.reorder_clustal(chain, pdb_file, clustal_file_old, clustal_file_new)

    return jsonify({'done': 'done'})

#@app.route('/explore/<tag>/<filename>/<connector_string>', methods=['GET', 'POST'])
@app.route('/explore/<tag>/<filename>/', methods=['GET', 'POST'])
@app.route('/explore/<tag>/<filename>', methods=['GET', 'POST'])
@app.route('/explore/<tag>/', methods=['GET', 'POST'])
def explore(tag, filename = ""):  #, connector_string = ""):
    if request.method == 'GET':
        #if connector_string == '':
        #    connector_string = "mut_0_1_diffE.pdb:mut_0_1_A.clw,mut_0_1,A;mut_0_1_diffE.pdb:mut_0_1_B.clw,mut_0_1,B;"
        return load_explore_page(app.config['USER_DATA_DIR'], tag, filename)#, connector_string)


    ### get form values

    mutations = request.form['mutations'].strip().replace(' ', '').split(',')
    parent = request.form['fname'].strip()
    if parent[-4:] != ".pdb":
        parent += '.pdb'
    mutant = name_mutation(app.config['USER_DATA_DIR'], parent, tag)
    
    print(__name__, 'original:', parent, 'novel mutant:', mutant) # TODO: OUTDIR ???


    ### mutate structure

    outdir = app.config['USER_DATA_DIR'] + tag + "/"

    email = request.form['email'].strip()
    # prewrite email (is sent seperately)
    if email:
        results_link = app.config["SERVER_URL"]+ url_for('explore', tag = tag, filename = mutant) 
        write_email(outdir + "mail.txt", email, results_link)

    helper_files_from_mutations( mutations, outdir + parent, outdir + mutant[:-4] + '_resfile.txt', outdir + mutant[:-4] + '.clw', outdir + "info/" + mutant[:-4] + '.txt', tag) 
    start_thread(fixbb, [tag, parent,  mutant[:-4] + '_resfile.txt', mutant, "log.txt"], "remutate")
    
    return redirect(url_for('status', tag = tag, filename = mutant))#, connector_string = connector_string))



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

@app.route('/documentation')
def documentation():
    return render_template('documentation.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

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
                elif ".clw" in f:
                    cmd.append(f)
                elif ".csv" in f:
                    cmd.append(f)
                elif ".aln" in f:
                    cmd.append(f)
                #elif "log.txt" in f:
                    #cmd.append(f)
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

    
def mutations_to_resfile( mutations, resfile, tag):
    try:
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
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + resfile)
    except IOError:
        fatal_error(tag, WRITE_FAILED + resfile)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + resfile)
"""                
def mutations_from_resfile( resfile):
    mutations = []
    with open( resfile) as r:
        r.readline()
        r.readline()
        for l in r:
            if 'NATAA' in 
"""

def mutations_from_vcf( fname, tag):
    mutations = []
    alphafold = ""
    print( 'get mutations from:', fname)

    try:
        with open( fname) as r:
            r.readline()
            for l in r:
                c = l.split(',')
                if alphafold == "":
                    alphafold = c[-1]
                    mutations.append( 'A:' + c[-2][1:] )
                elif c[-1] == alphafold:
                    mutations.append( 'A:' + c[-2][1:] )
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + fname)
    except IOError:
        fatal_error(tag, READ_FAILED + fname)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + fname)
    
    return alphafold,mutations

def resid( line):
    return int( line[22:26].strip() )


def residue_name( line ):
    return line[17:20]

def atom_name( line):
    return line[12:16].strip()

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



def sequence_chain_resids(parent, tag):
    chains = defaultdict(list)

    try:
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
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + parent)
    except IOError:
        fatal_error(tag, READ_FAILED + parent)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + parent)
    
    return chains



def alignment_from_mutations(mutations, parent, align, mutant_file, tag):
    chains = sequence_chain_resids(parent, tag)
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



            
            
def mutation_parent_file( mutations, parent, mutfile, tag):   
    try:
        with open( mutfile, 'w') as f:
            if mutations:
                f.write(parent.split('/')[-1] + "\n")
                for mut in mutations:
                    f.write(mut + "\n")
            else:
                f.write("none\n")
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + mutfile)
    except IOError:
        fatal_error(tag, WRITE_FAILED + mutfile)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + mutfile)

            
def helper_files_from_mutations( mutations, parent, resfile, align, mutfile, tag):
    mutations_to_resfile( mutations, resfile, tag)
    print("give the threads some time to terminate")
    time.sleep(5)
    mutation_parent_file( mutations, parent, mutfile, tag)
    print("give the threads some time to terminate")
    time.sleep(5)
    alignment_from_mutations( mutations, parent, align, mutfile, tag)


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
    

def add_mutations_from_sequence( mutations, target, chain, idy, parent, tag):
    print("add_mutations_from_sequence")
    print("idy: ", idy, 'chain:', chain, 'target:', target, 'parent:', parent, 'mutations:', mutations)
    h1 = idy
    h2 = os.path.basename(parent)[:-4]
    outdir = os.path.dirname( parent) + '/'
    f1 =  outdir + h1 + '.fa'
    f2 =  outdir + h2 + '.fa'
    f3 =  outdir + h2 + '_' + h1 + '.clw'
    pdb_seq = pdb2seq(parent, tag)
    write_fasta(f1 ,target ,tag, h1)
    write_fasta(f2, pdb_seq[chain][0], tag, h2)
    cmd = "python3  " + app.config['SCRIPTS_PATH'] + "seq_align.py " + f1 + ' ' + f2 + ' ' + f3 + ' clw'
    print(cmd)
    p = subprocess.check_output( cmd.split())
    print(p)
    add_mutations_from_alignment(mutations, f3, parent, tag)

    
def seq_from_fasta(filename, tag):
    seq = ''
    head = ''

    try:
        with open( filename) as r:
            l = r.readline()
            if l[0] != '>':
                print( "WARNING!:",filename, 'does not start with ">" in first line! First line is ignored!!!' )
            head = l[1:].strip()
            for l in r:
                if l[0] == '>':
                    return head,seq
                seq += l.strip()
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + filename)
    except IOError:
        fatal_error(tag, READ_FAILED + filename)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + filename)

    return head,seq

def write_fasta(filename, seq, tag, header=""):
    print( 'write_fasta:', filename, seq, header)

    try:
        with open(filename,'w') as w:
            w.write('>' + header + '\n')
            w.write(seq + '\n')
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + filename)
    except IOError:
        fatal_error(tag, WRITE_FAILED + filename)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + filename)
    
def add_mutations_from_alignment( mutations, clustal_file, parent, tag, base_chain=""):
    mutations.extend( mutations_from_alignment( clustal_file, parent, tag, base_chain=base_chain)[0] )


        
def read_clustal(clust, tag):
    ali = defaultdict( str)

    try:
        with open(clust) as r:
            r.readline()
            for l in r:
                if l.strip() == '' or l[:3] == '   ': continue
                c = l.split()
                if len(c) != 2: continue
                ali[c[0]] += c[1]
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + clust)
    except IOError:
        fatal_error(tag, READ_FAILED + clust)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + clust)

    print( 'read_clustal:', ali)
    return ali



def sequences( ali):
    seqs = []
    for k,v in ali.items():
        seqs.append( v.replace('-','') )
    return seqs



# nearly same as sequence_chain_resids()
def pdb2seq(pdb, tag, noncanon_warning=False):
    noncanonical_residues = False
    chains = defaultdict( lambda : ['',[]] )

    try:
        with open( pdb) as r:
            prev_chain = 'xxx'
            prev_id = -9999
            for l in r:
                if l[:4] != "ATOM" and l[:6] != "HETATM": continue
                c = chain(l)
                res = resid(l)
                if prev_chain != c or prev_id != res:
                    #print(prev_chain)
                    #print(prev_id)
                    #print(chains)
                    aa = single_letter( residue_name(l))
                    if aa == 'X':
                        noncanonical_residues = True
                        continue
                    chains[c][0] += aa
                    chains[c][1].append( res )
                    prev_chain = c
                    prev_id = res
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + pdb)
    except IOError:
        fatal_error(tag, READ_FAILED + pdb)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + pdb)

    if noncanon_warning:
        return chains, noncanonical_residues
    return chains





def get_seq_id(clustal, cid, tag):
    # returns position of a sequence in an alignment
    # very similar to read_clustal

    i = 0
    try:
        with open(clustal) as r:
            r.readline()
            for l in r:
                if l.strip() == '' or l[:3] == '   ': continue
                c = l.split()
                if len(c) != 2: continue

                if c[0] == cid:
                    return i
                i += 1
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + clustal)
    except IOError:
        fatal_error(tag, READ_FAILED + clustal)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + clustal)

    return None


def find_pdb_in_alignment(clustal, structure, tag, chain="", clustal_id=""):
    # returns a (clustal id, chain) pair with matching sequences, if possible with provided chain and/or clustal id
    # if no such pair exists, None is returned

    print("################")
    print("find pdb in alignment")
    print(clustal)
    if os.path.isfile(structure):
        print(structure)
    else:
        print( structure, 'does not yet exist, waiting ..')
        count = 0
        while not os.path.isfile(structure):
            print( 'wait')
            time.sleep(5)
            count += 1
            if count >= 18:
                print('giving up')
                exit(1)
    print(chain)
    print(clustal_id)


    # warum doppelt?
    alignment = read_clustal(clustal, tag) # dict; key: clustal id, value: seq
    pdb_seq, noncanonical_residues = pdb2seq(structure, tag, noncanon_warning = True)
    pdb_chains = {k: v[0] for (k, v) in pdb_seq.items()} # dict; key: chain, value: seq

    print( "alignment:")
    print(alignment)
    print( "chains:")
    print(pdb_chains)

    # find all (clustal id, chain) pairs with matching sequences
    matches = []
    for (cid, clustal_seq) in alignment.items():
        for (c, chain_seq) in pdb_chains.items():
            if chain_seq == clustal_seq.replace("-", ""):
                matches.append([cid, c])

    if len(matches) == 0:
        print("nothing found")
        # no matching pair exists
        return ["",""], noncanonical_residues


    # select (base clustal id, base chain) pair based on provided parameters
    selected_match = ["", ""] # [clustal id, chain] 
    
    if len(matches) == 1:
        selected_match = matches[0]

    else:
        # multiple matches

        available_chains = [m[1] for m in matches]
        available_clustal_ids = [m[0] for m in matches]

        chain_given = (chain != "" and chain in available_chains)
        cid_given = (clustal_id != "" and clustal_id in available_clustal_ids)

        if chain_given and not cid_given:
            # if only the base chain is given, select that chain; cid randomly
            selected_match[1] = chain
            available_clustal_ids = [m[0] for m in matches if m[1] == selected_match[1]] 
            selected_match[0] = available_clustal_ids[0]

        else:
            # else select clustal id first

            # select clustal id
            if cid_given:
                # if possible, take provided base clustal id
                selected_match[0] = clustal_id
            else:
                # otherwise, take random (first) clustal id
                selected_match[0] = available_clustal_ids[0]

            print("############")
            print("point of no return")
            print(matches)
            print(clustal_id)

            # select chain
            available_chains = [m[1] for m in matches if m[0] == selected_match[0]]
            print(available_chains)
            if chain != "" and chain in available_chains:
                # if possible, take provided base chain
                selected_match[1] = chain
            else:
                # otherwise, take random (first) chain
                selected_match[1] = available_chains[0]

    return selected_match, noncanonical_residues

    

def mutations_from_alignment(clustal, base_structure, tag, base_clustal_id="", target_clustal_id="", base_chain="" ):

    print( "mutations_from_alignment")
    alignment = read_clustal(clustal, tag) # dict; key: clustal id, value: seq
    pdb_chains = {k: v[0] for (k, v) in pdb2seq(base_structure, tag).items()} # dict; key: chain, value: seq

    resids = {k: v[1] for (k, v) in pdb2seq(base_structure, tag).items()} # dict; key: chain, value: resids of seq


    # find matching pdb chain, clustal id pair
    pdb_match, noncanonical_residues = find_pdb_in_alignment(clustal, base_structure, tag, chain=base_chain, clustal_id=base_clustal_id)
    if pdb_match is None:
        #exit(1) # TODO: error?
        return []
    pdb_cid, pdb_chain = pdb_match


    # select target clustal id
    target_cid = ""
    available_cids = [cid for cid in alignment.keys() if alignment[cid] != alignment[pdb_cid]]
    if len(available_cids) == 0:
        return []

    if target_clustal_id != "" and target_clustal_id in available_cids:
        target_cid = target_clustal_id
    else:
        target_cid = available_cids[0]


    # get mutations
    #base_seq = pdb_chains[pdb_chain]
    base_seq = alignment[pdb_cid]
    target_seq = alignment[target_cid]

    base_resids = resids[pdb_chain]


    if len(base_seq) != len(target_seq):
        print( "seqs differ in length:",  len(base_seq) , len(target_seq) )

    print("################")
    print("get mutations from alignment")
    print('base:', base_seq)
    print('target:', target_seq)

    mutations = []

    count = 0
    for i in range( len( base_seq )):
        # i = index for sequences in alignment; might include gaps
        # cound = index for base sequence exluding all gaps

        # no mutations from gaps
        if base_seq[i] == '-':
            continue
        if target_seq[i] == '-':
            count += 1
            continue

        # sequences contain different residues (no gaps) => mutation
        if target_seq[i] != base_seq[i]:
            mutations.append( pdb_chain + ':' + str(base_resids[count]) + target_seq[i] )
        count += 1

    print( __name__, 'found', len(mutations), 'mutations')

    return mutations, noncanonical_residues

def waitID(myid):
    print('wait for process: ', myid)
    state = True

    while(state):
        cmd = 'tsp -s ' + str(myid)
        pid = subprocess.run(cmd.split(), check=True, capture_output=True, text=True).stdout
        #print(pid)
        if("finished" in pid):
            return state
        if(not ("queued" in pid or "running" in pid)):
            return False








                    
            
def wait( filename, step, maxw):
    # TODO: wait sollte nur in seperaten Threads aufgerufen werden, sonst muss der user auch waiten
    print( 'wait: ' , filename, step, maxw)
    for i in range(0, maxw):
        time.sleep(step)
        if os.path.isfile(filename):
            print('found',i)
            return True
    return False
        
def download_uniprot( unid, filename, tag):
    link = "https://rest.uniprot.org/uniprotkb/" + secure_filename(unid) + '.fasta'
    req = requests.get(link)

    try:
        with open(filename, "w") as f:
            f.write(req.text)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + filename)
    except IOError:
        fatal_error(tag, WRITE_FAILED + filename)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + filename)
        

def get_chains(fname, tag):
    chains = {}
    print( __name__, 'get', fname)

    try:
        with open( fname) as r:
            for l in r:
                #if l[:4] == "ATOM" or l[:6] == "HETATM":
                    #c = chain(l)
                    #if c not in chains:
                        #chains += c


                if l[:4] == "ATOM":  # does not make sense: or l[:6] == "HETATM":
                    c = chain(l)
                    r = resid(l)
                    if c not in chains.keys():
                        chains[c] = [[r,r]]
                    else:
                        prev = chains[c][-1][1]
                        # case: numbering jumps backward or there is gap
                        if prev > r or r > prev + 1:
                            chains[c].append( [r,r] )
                        elif r == prev + 1:
                            chains[c][-1][1] = r
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + fname)
    except IOError:
        fatal_error(tag, READ_FAILED + fname)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + fname)
                
    chainstr = ""
    for c,v in chains.items():
        for x in v:
            if not c in chainstr:
                chainstr += c 
    print( chainstr)
    return chainstr

def get_chains_and_range(fname, tag):
    chains = {}

    try:
        with open( fname) as r:
            for l in r:
                if l[:4] == "ATOM":  # does not make sense: or l[:6] == "HETATM":
                    c = chain(l)
                    r = resid(l)
                    if c not in chains.keys():
                        chains[c] = [[r,r]]
                    else:
                        prev = chains[c][-1][1]
                        # case: numbering jumps backward or there is gap
                        if prev > r or r > prev + 1:
                            chains[c].append( [r,r] )
                        elif r == prev + 1:
                            chains[c][-1][1] = r
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + fname)
    except IOError:
        fatal_error(tag, READ_FAILED + fname)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + fname)

    chainstr = ""
    for c,v in chains.items():
        for x in v:
            chainstr += c + ': ' + str(x[0]) + '-' + str(x[1]) + ', '
    #print( chainstr)
    return chainstr

def open_chains_range_file(outdir, tag):
    chains_range = ''
    try:
        with open( outdir + "chains.txt") as r:
            chains_range = r.read()
            chains_range = chains_range[:-1]
            print("chain_range ", chains_range)
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + outdir + "chains.txt")
    except IOError:
        fatal_error(tag, READ_FAILED + outdir + "chains.txt")
    except Exception as e:
        print("An error occurred while trying to read the chain ranges file.")
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + outdir + "chains.txt")
        
    return chains_range

def read_txt_file(outdir, filename, tag):
    string = ''
    try:
        with open(outdir + filename) as r:
            string = r.read()
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + outdir + filename)
    except IOError:
        fatal_error(tag, READ_FAILED + outdir + filename)
    except Exception as e:
        print("An error occurred while trying to read a file " + filename)
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + outdir + filename)
        
    return string

def write_alignment_match(outdir, string, tag):
    try:
        with open( outdir + 'alignment.txt', 'w') as w:
            w.write(string)
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + outdir + 'alignment.txt')
    except IOError:
        fatal_error(tag, WRITE_FAILED + outdir + 'alignment.txt')
    except Exception as e:
        print("An error occurred while trying to write the alignment match string into a file.")
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + outdir + 'alignment.txt')

def get_chain_letters(raw_ranges):
    ranges = raw_ranges.split(",")
    letters = []
    for r in ranges:
        if len(r.strip()) == 0:
            continue
        letter = r.strip()[0]
        if letter not in letters:
            letters.append(letter)
    return letters

def allowed_file(filename, extensions):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in extensions

def get_energy(fname, tag):
    print( __name__, 'get', fname)

    try:
        with open( fname) as r:
            for l in r:
                if l[:4] == "pose":
                    return float( l.split()[-1])
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + fname)
    except IOError:
        fatal_error(tag, READ_FAILED + fname)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + fname)

    print( 'WARNING: no energy found in', fname)
    return -0

def add_energy(fromfile, tofile, tag):
    try:
        with open( tofile, 'a') as w:
            w.write( str( get_energy(fromfile, tag)) + '\n')
    except FileNotFoundError:
        fatal_error(tag, WRITE_FAILED + FILE_NOT_FOUND + tofile)
    except IOError:
        fatal_error(tag, WRITE_FAILED + tofile)
    except Exception as e:
        fatal_error(tag, WRITE_FAILED + UNEXPECTED + e + ' ' + tofile)
        
def is_pdb(fname, tag):
    if not os.path.isfile(fname):
        return False
    
    try:
        with open( fname) as r:
            for l in r:
                if l[:4] == "ATOM" or l[:6] == "HETATM":
                    return True
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + fname)
    except IOError:
        fatal_error(tag, READ_FAILED + fname)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + fname)

    return False


def write_email(fil, user, link):
    with open(fil, 'w') as out:
        out.write(user + "\n")
        out.write('Hello  ' + user + '!\n\n')
        out.write('Your MutationExplorer calculation is done. \n')
        out.write('You can view the results here: \n\n' + link + '\n\n')
        out.write('Thanks for using MutationExplorer.\n\n')
        out.write('Have a nice day!\n\n')

def send_email(fil):
   if(os.path.exists(fil)):
   #     os.system('sendmail -t < ' + fil)
        smtp_server = app.config['SMTP_SERVER']
        port = 587
        sender_user = app.config['SENDER_LOGIN']
        sender_email = app.config['SENDER_MAIL']
        password = app.config['MAIL_PASSWORD']
        receiver = ""
        message = "Subject: MutationExplorer results \n\n"

        with open(fil, 'r') as mail:
           receiver = mail.readline()
           line = mail.readline()
           while(line != ''):

                message = message + line
                line = mail.readline()


        print(receiver)
        #print(message)


        context = ssl.create_default_context()
        receiver_email = receiver
        try:
            server = smtplib.SMTP(smtp_server,port)
            server.ehlo() 
            server.starttls(context=context) 
            server.ehlo()
            server.login(sender_user, password)
            server.sendmail(sender_email, receiver_email, message)
            os.remove(fil)
        except Exception as e:
            print(e)
        finally:
            server.quit() 

def send_error_mail(tag):
  
        smtp_server = app.config['SMTP_SERVER']
        port = 587
        sender_user = app.config['SENDER_LOGIN']
        sender_email = app.config['SENDER_MAIL']
        password = app.config['MAIL_PASSWORD']
        receiver = ""
        message = "Subject: MutationExplorer failed! \n\n"
        message = message + tag


        
        context = ssl.create_default_context()
        receiver_email = "daniel@informatik.uni-leipzig.de"
        try:
            server = smtplib.SMTP(smtp_server,port)
            server.ehlo() 
            server.starttls(context=context) 
            server.ehlo()
            server.login(sender_user, password)
            server.sendmail(sender_email, receiver_email, message)
            os.remove(fil)
        except Exception as e:
            print(e)
        finally:
            server.quit() 





def is_in_db( pdb):
    rose = app.config['ROSEMINT_PATH']
    return len(glob.glob( rose + 'pdb/' + pdb.upper() + '*.pdb')) > 0 or len(glob.glob( rose + 'fixbb/' + pdb.upper() + '*.pdb')) > 0 or len(glob.glob( rose + 'alphafold/' + pdb.upper() + '*.pdb')) > 0

def cp_from_db(pdb, outfile, tag):
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
            energy = get_energy(l, tag)
            if energy < mini:
                best = l.strip()
                mini = energy
        shutil.copyfile(best,outfile)

  
def get_current_time():
   dt = datetime.datetime.now()
   x = dt.strftime("%Y-%m-%d+%H:%M:%S")

   return x



def add_conservation( ali ):
    s = defaultdict( list )
    rows = ali.split('\n')
    newali = rows[0] + '\n\n'
    for l in rows[2:]:
        if l.strip() == '' or l[:4] == '    ': continue # ignore empty lines or existing conservation lines
        c = l.split()
        if len(c) != 2: continue
        s[c[0]].append( c[1])
    k = list(s.keys())
    m = max( len(k[0]), len(k[1]) )
    a = k[0].ljust(m)
    b = k[1].ljust(m)
    e = ''.ljust( m)
    
    if len(k) != 2 or ( len( s[k[0]] ) != len( s[k[1]] ) ):
        print('ERROR in alignment: input:', content, '\nnow:', s )
        return ''
    for i in range( len( s[k[0]] ) ):
        s1 = s[k[0]][i]
        s2 = s[k[1]][i]
        newali +=  a + '\t' + s1 + '\n'
        newali +=  b + '\t' + s2 + '\n'
        newali +=  e + '\t'
        for j in range( len(s1)):
            if s1[j] == s2[j]:
                newali += '*'
            else:
                newali += ' '
        
        newali += '\n\n'
    return newali


def getLastID(outdir, tag):
    try:
        with open(outdir+"id.txt") as f:
            for line in f:
                pass
            pid = line
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + outdir+"id.txt")
    except IOError:
        fatal_error(tag, READ_FAILED + outdir+"id.txt")
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + outdir+"id.txt")

    return pid



def protein_type(file_path, tag):
    bb = ['C','CA','N','O']
    atom_names = []

    try:
        with open( file_path) as r:
            for l in r:
                if l[:4] != "ATOM":continue
                name = atom_name( l)
                if name not in atom_names and ('C' in name or 'N' in name or 'O' in name):
                    atom_names.append(name) # ignore H, S, ...
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + file_path)
    except IOError:
        fatal_error(tag, READ_FAILED + file_path)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + file_path)

    print( 'heavy atoms:', atom_names)
    nr = len(atom_names)
    if nr == 0:
        return "void"
    if nr == 1 and atom_names[0] == 'CA':
        return "calpha"
    if nr == len(bb) and sorted( bb) == sorted( atom_names):
        return "backbone"
    return "fullatom"

def get_alignment_ids(filename, tag):
    idstr = ""

    try:
        with open( filename) as r:
            r.readline()
            ids = []
            for l in r:
                l = l.strip()
                if len(l) == 0:
                    continue
                c = l.split()
                if len(c) == 2 and c[0] not in ids:
                    ids.append( c[0] )
            if len(ids) > 0:
                idstr = ids[0]
            for i in range(1,len(ids)):
                idstr += ';' + ids[i]
    except FileNotFoundError:
        fatal_error(tag, READ_FAILED + FILE_NOT_FOUND + filename)
    except IOError:
        fatal_error(tag, READ_FAILED + filename)
    except Exception as e:
        fatal_error(tag, READ_FAILED + UNEXPECTED + e + ' ' + filename)
            
    return idstr

