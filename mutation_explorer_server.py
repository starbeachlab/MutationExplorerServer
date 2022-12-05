from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, jsonify
import os, random, subprocess, time
import threading
import requests
import glob


app = Flask(__name__)


app.config['USER_DATA_DIR'] = "/disk/user_data/mutation_explorer_gamma/"
app.config['ROSETTA_PATH']  = "/home/hildilab/dev/ext/rosetta/bin/"
app.config['SCRIPTS_PATH']  = "/home/hildilab/app/mutation_explorer_gamma/scripts/"


@app.route('/')
def index():
    return render_template("home.html")


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


def fixbb(tag, structure, mutations, name, log):
    out = app.config['USER_DATA_DIR'] + tag + "/"

    # create resfile
    if not mutations:
        resfile = app.config['SCRIPTS_PATH'] + "resfile.txt"
    if mutations:
        resfile = out + "resfile.txt"
        with open(resfile, 'w') as f:
            f.write('NATRO\n')
            f.write('start\n')
            prev_chains = []
            for mut in mutations:
                chain = mut[0]
                resid = mut[2:-1]
                res = mut[-1]

                if chain not in prev_chains:
                    f.write('* ' + chain + ' NATAA \n')
                    prev_chains.append(chain)

                f.write(resid + ' ' + chain + ' PIKAA ' + res + '\n')

    # wait for previous structure
    while not os.path.isfile(out + structure):
        time.sleep(1)

    # generate unique extension (temp file)
    while(True):
        ext = "m" + str(random.randint(0, 999))
        if not glob.glob(out + ext + r'.*\.pdb$'):
            break

    # call rosetta
    cmd = "tsp " + app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -in:file:s " + out + structure + " -resfile " + resfile + ' -nstruct 1 -linmem_ig 10 -out:pdb  -out:prefix ' + out + ext + '_'
    print(cmd)
    p = subprocess.check_output(cmd.split())

    # rename output file
    cmd = "tsp mv " + out + ext + "_" + structure.split('.')[0] + "_0001.pdb " + out + name
    print(cmd)
    p = subprocess.check_output(cmd.split())

    # create info file
    with open(out + "info/" + name.split(".")[0] + ".txt", "w") as f:
        if mutations:
            f.write(structure + "\n")
            for mut in mutations:
                f.write(mut + " ")
        else:
            f.write("-")

    # add to file listing
    with open(out + "list.txt", "a") as f:
        f.write(name + "\n")

    # delete resfile
    if mutations:
        cmd = "tsp rm " + out + "resfile.txt"
        p = subprocess.check_output(cmd.split())

    print("### THREAD FINISHED ###")

        


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
    upload = request.files['file']
    pdb = request.form['pdb'].strip()
    af = request.form['alphafold'].strip()

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
        f.write("hi")
        # TODO: actually log something

    # create list file
    with open(outdir + "list.txt", "w") as f:
        f.write("structure.pdb\n")

    # create info file
    os.mkdir(outdir + "info/")
    with open(outdir + "info/structure.txt", "w") as f:
        f.write("-")

    # relax structure
    start_thread(fixbb, [tag, "structure.pdb", [], "relaxed_structure.pdb", outdir + "log.txt"], "minimisation")

    return redirect(url_for('mutate', tag = tag))


@app.route('/mutate/<tag>', methods=['GET', 'POST'])
def mutate(tag):
    if request.method == 'GET':
        return render_template("mutate.html", tag = tag, error = "")

    # get form values
    mutations = request.form['mutations'].strip().split('; ')
    mail = request.form['email'].strip()

    name = request.form['name'].strip() + ".pdb"
    if name == ".pdb":
        name = name_mutation("relaxed_structure.pdb", tag)

    # return error message if no mutations given
    if not mutations:
        return render_template("mutate.html", tag = tag, error = "Please provide a mutation")

    # mutate structure
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    start_thread(fixbb, [tag, "relaxed_structure.pdb", mutations, name, outdir + "log.txt"], "minimisation")
    
    return redirect(url_for('status', tag = tag, file = name))


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

    strc_upload = request.files['file']
    pdb = request.form['pdb'].strip()
    af = request.form['alphafold'].strip()

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

    return redirect(url_for('status', tag = tag, file = "mut_structure.pdb"))
    


@app.route('/get_status/<tag>/<file>')
def get_status(tag, file):
    dir = app.config['USER_DATA_DIR'] + tag + "/"
    done = os.path.isfile(dir + file)
    status = "done"
    msg = ""
    return jsonify({'done': done, 'status': status, 'message': msg})


@app.route('/status/<tag>/<file>')
def status(tag, file):
    return render_template("status.html", tag = tag, file = file)


def build_list(d):
    s = ""
    for ke in d:
        link = "<a id='" + ke + "' href='' class='structures'>" + ke + "</a><br>"
        if d[ke] != {}:
            s += "<li><i class='toggle material-icons tiny rot'>play_arrow</i>" + link + "</li>"
            s += "<ul class='browser-default non' style='list-style-type: none'>" + build_list(d[ke]) + "</ul>"
        else:
            s += "<li><i class='material-icons tiny grey-text'>play_arrow</i>" + link + "</li>"

    return s


@app.route('/explore/<tag>', methods=['GET', 'POST'])
def explore(tag):
    if request.method == 'GET':
        mut_tree = build_mutation_tree(tag, "-")
        structures = "<ul>" + build_list(mut_tree) + "</ul>"

        return render_template("explore.html", tag = tag, structures = structures)

    # get form values
    mutations = request.form['mutations'].strip().split('; ')
    strc = request.form['structure'].strip()
    name = request.form['name'].strip() + ".pdb"
    if name == ".pdb":
        name = name_mutation(strc, tag)

    # mutate structure
    outdir = app.config['USER_DATA_DIR'] + tag + "/"
    start_thread(fixbb, [tag, strc, mutations, name, outdir + "log.txt"], "minimisation")
    
    return redirect(url_for('status', tag = tag, file = name))
    

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
