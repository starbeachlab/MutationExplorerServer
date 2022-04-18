from flask import Flask, render_template, request, url_for, send_from_directory, jsonify
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, FileField
from flask_bootstrap import Bootstrap
from werkzeug import secure_filename

import os, datetime, random, string
import requests, subprocess, sys


app = Flask(__name__)
app.config['USER_DATA_DIR'] = "/disk/user_data/mutation_explorer/"
app.config['APP_PATH']      = "/home/hildilab/app/mutation_explorer/"
app.config['SECRET_KEY']    = 'klsdf23*&%..'

bootstrap = Bootstrap( app )



@app.route('/')
def index():
    return render_template('home.html')


@app.route('/submit', methods=['GET', 'POST'])
def submit():
    if request.method == 'POST':

        # get data from form

        tag = str(random.randint(0, 999999))
        outdir = app.config['USER_DATA_DIR'] + tag + "/"
        while os.path.exists(outdir):
            tag = str(random.randint(0, 999999))
            outdir = app.config['USER_DATA_DIR'] + tag + "/"
        os.mkdir(outdir)

        file_conv = request.files['file_conv']
        file_conv_text = request.form['file_conv_text'].strip()
        if file_conv.filename != "":
            file_conv_fn = outdir + secure_filename(file_conv.filename)
            file_conv.save(file_conv_fn)
        elif file_conv_text != "":
            file_conv_text = secure_filename(file_conv_text).lower()
            url = 'https://files.rcsb.org/view/{0}.pdb'.format(file_conv_text.upper())
            r = requests.get(url, allow_redirects=True)
            file_conv_fn = outdir + file_conv_text + '.pdb'
            open(fn, 'wb').write(r.content)

        superimpose = False

        file_super = request.files['file_super']
        file_super_text = request.form['file_super_text'].strip()
        if file_super.filename != "":
            superimpose = True
            file_super_fn = outdir + secure_filename(file_super.filename)
            file_super.save(file_super_fn)
        elif file_super_text != "":
            superimpose = True
            file_super_text = secure_filename(file_super_text).lower()
            url = 'https://files.rcsb.org/view/{0}.pdb'.format(file_super_text.upper())
            r = requests.get(url, allow_redirects=True)
            file_super_fn = outdir + file_super_text + '.pdb'
            open(fn, 'wb').write(r.content)

        alignment_link = request.form.getlist('alignment_link')[0]
        alignment_link = "https://www.bioinfo.mpg.de/AlignMeBeta/work/" + alignment_link.split("work/")[1]
        alignment = outdir + "alignment.aln"
        
        #alignment_link = app.config['USER_DATA_DIR'] + "alignment.aln"

        req = requests.get(alignment_link)
        with open(alignment, "w") as f:
            f.write(req.content)


        ###  determine chain and seqid
        chain, seqid = ChainSeqidFromPDBandAlignment( file_conv_fn, alignment )        

        print( 'check', chain, seqid)

        error = ""
        if chain == '' or seqid == '':
            error = "The program was not able to find any sequence of the alignment in the PDB. Please note, that one chain has to match EXACTLY to one of the sequences within the alignment!"

        if not superimpose:
        
            pdb_name = file_conv_fn
            aliname = alignment

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py ' + pdb_name + ' ' + chain + ' ' + aliname + ' ' + seqid + ' 0.2  > ' + outdir + 'colored.pdb'
            print( cmd )
            os.system( cmd )


            return render_template('overview.html', tag=tag, file_conv=url_for('download', tag=tag, filename="colored.pdb"), errmsg=error) # Maybe needs changing
            
        
        chain_1 = chain
        seqid_1 = seqid
        pdb_name_1 = file_conv_fn
        ali_name = alignment

        pdb_name_2 = file_super_fn

        chain_2, seqid_2 = ChainSeqidFromPDBandAlignment( file_super_fn, alignment )        

        if chain_2 == '' or seqid_2 == '':
            error = "The program was not able to find any sequence of the alignment in the PDB. Please note, that one chain has to match EXACTLY to one of the sequences within the alignment!"
        
        if seqid_1 == '0' and seqid_2 == '1':
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_superimpose.py alignment: ' + pdb_name_1 + ' ' + chain_1 + ' ' + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + outdir + 'super.pdb'
            
            print( cmd )
            sys.stdout.flush()
            os.system( cmd )

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + outdir + 'super.pdb' + ' ' + chain_1 + ' ' + ali_name + ' ' + seqid_1 + ' 0.2  > ' + outdir + 'super_colored.pdb'
            print( cmd )
            os.system( cmd )

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + seqid_2 + ' 0.2  > ' + outdir + 'template_colored.pdb'
            print( cmd )
            os.system( cmd )

        elif seqid_1 == '1' and seqid_2 == '0':
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_superimpose.py alignment: ' + pdb_name_2 + ' ' + chain_2 + ' ' + pdb_name_1 + ' ' + chain_1 + ' ' + ali_name + ' ' + outdir + 'super.pdb'
            
            print( cmd )
            sys.stdout.flush()
            os.system( cmd )

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + outdir + 'super.pdb' + ' ' + chain_2 + ' ' + ali_name + ' ' + seqid_2 + ' 0.2  > ' + outdir + 'super_colored.pdb'
            print( cmd )
            os.system( cmd )

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + pdb_name_1 + ' ' + chain_1 + ' ' + ali_name + ' ' + seqid_1 + ' 0.2  > ' + outdir + 'template_colored.pdb'
            print( cmd )
            os.system( cmd )

        
        # return render_template( 'fullmenu.html', path=tag, pdb="super_colored.pdb", chain=chain_1, pdb2="template_colored.pdb", chain2 = chain_2 )  # use same html ?
        return render_template('overviewSuper.html', tag=tag, file_conv=url_for('download', tag=tag, filename="super_colored.pdb"), file_super=url_for('download', tag=tag, filename="template_colored.pdb"), errmsg=error) # Maybe needs changing

    return render_template('home.html')


@app.route('/menu/<tag>')
def menu(tag):
    if os.path.exists(app.config['USER_DATA_DIR'] + tag + "/super_colored.pdb"):
        return render_template('menuSuper.html', tag=tag, file_conv=url_for('download', tag=tag, filename="super_colored.pdb"), file_super=url_for('download', tag=tag, filename="template_colored.pdb"))
    return render_template('menu.html', tag=tag, file_conv=url_for('download', tag=tag, filename="colored.pdb"))


@app.route('/downloads/<tag>/<filename>')
def download(tag, filename):
    path = app.config['USER_DATA_DIR'] + '/' + tag 
    return send_from_directory(path, filename)


#####################       GENERAL     ##################################

def ChainSeqidFromPDBandAlignment( pdb, ali):
    #print( 'ChainSeqid from PDB & Alignment')
    cmd = app.config['APP_PATH'] + 'pdb_sequences.py'
    #print( cmd)
    #pdb_fastas = subprocess.run( [ 'python3', cmd , alignment_link ], stdout=subprocess.PIPE).stdout.decode('utf-8') # from python 3.6 on
    pdb_fastas = subprocess.check_output( [ 'python3', cmd , pdb ]) 
    #print( pdb_fastas)
    #print( len(pdb_fastas))
    s = pdb_fastas.split('\n')
    pdb_seqs = s[1::2]
    pdb_chains = [ f.split(':')[-1] for f in s[0::2] ]
    #print( pdb_chains)
    #print( pdb_seqs)
    
    cmd = app.config['APP_PATH'] + 'aln2fasta.py'
    #print( cmd)
    aln_seqs =  subprocess.check_output( [ 'python3', cmd , ali ] ).split('\n')
    #print(aln_seqs)
    aln_seqs = [ f.replace('-','') for f in aln_seqs[1::2] ]
    #print( aln_seqs)
    chain = ''
    seqid = ''
    for i in range( len(pdb_seqs)):
        for j in range( len( aln_seqs)):
            if pdb_seqs[i] == aln_seqs[j]:
                chain = pdb_chains[i]
                seqid = str(j)
                break
        # needed?
        if chain != '' and seqid != '':
            break
    if chain == '' or seqid == '':
        # print error message in form!
        print('WARNING, no match between PDB and alignment found, no identical sequence!')
    print( 'chain:', chain, 'seqid:', seqid)
    return chain, seqid



@app.route( '/get/<mydir>/<pdb>' )
def get(mydir,pdb):
    out = app.config['USER_DATA_DIR'] + mydir # + "/colored.pdb"
    print( 'get:', out, pdb)
    return send_from_directory( out, pdb)



def randstr( nr ):
    return ''.join(random.choice(string.ascii_lowercase) for i in range(nr))



def check_ending( form, field ):
    if field.data:
        ext = os.path.splitext( field.data.filename)[1].strip('.').lower()
        if ext != '.pdb' and ext != '.aln' and ext != '.clw':
            raise validators.ValidationError( 'Allowed file endings are ".pdb", ".aln", or ".clw"')
    else:
        raise validators.ValidationError( 'Please provide file' )

    

#@app.route( '/json', methods=['GET','POST'])
#def json():
#    content = request.get_json(silent=True)
#    print("json:", content)
#    nr = int( content['nr'])
#    print( nr, 'entries')
#    for i in range( nr):
#        print( content[ str(i) ] )
#    


    
#####################   CONSERVATION    ##################################
    
class ConservationForm(FlaskForm):
    pdb = FileField('Upload a pdb file:' , validators = [check_ending] )
    ali = FileField('Upload an alignment:', validators = [check_ending] )
    seqID = StringField( 'Select sequence identifier matching PDB:')
    chain = StringField( 'Chain to be used from PDB' )
    submit = SubmitField( 'Visualize')

    

@app.route('/conservation', methods=['GET', 'POST'])
def conservation():
    form = ConservationForm()

    if request.method == 'GET':
        return render_template( 'conservation.html', form=form )

    
@app.route( '/conserve', methods=['POST'] )
def conserve():
    form = ConservationForm()
    
    pdb   = request.files[ 'pdb' ]
    ali   = request.files[ 'ali' ]
    chain = request.form[ 'chain' ]
    seqid = request.form[ 'seqID' ]

    check = True
    today = str(datetime.date.today()).replace('-', '')
    job = ""
    while check:
        job = today + '_' + randstr(5)
        outdir = app.config['USER_DATA_DIR'] + job + '/'
        #print( outdir )
        if not os.path.exists( outdir ):
            os.mkdir( outdir )
            check = False
            
    pdb_name = outdir + os.path.basename( pdb.filename )
    pdb.save( pdb_name )

    ali_name = outdir + os.path.basename( ali.filename )
    ali.save( ali_name )

    cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py ' + pdb_name + ' ' + chain + ' ' + ali_name + ' ' + seqid + ' > ' + outdir + 'colored.pdb'
    print( cmd )
    os.system( cmd )

    return render_template( 'display.html', path=job, pdb="colored.pdb" )



#####################  SUPERIMPOSITION  ###################################
    
class SuperForm(FlaskForm):
    pdb_1 = FileField('Upload first PDB file:' , validators = [check_ending] )
    pdb_2 = FileField('Upload second PDB file:' , validators = [check_ending] )
    ali = FileField('Upload an alignment:', validators = [check_ending] )
    seqID_1 = StringField( 'Select sequence identifier matching first PDB:')
    seqID_2 = StringField( 'Select sequence identifier matching second PDB:')
    chain_1 = StringField( 'Chain to be used from first PDB' )
    chain_2 = StringField( 'Chain to be used from second PDB' )
    submit = SubmitField( 'Visualize')

   

@app.route('/superimposition', methods=['GET', 'POST'])
def superimposition():
    form = SuperForm()

    if request.method == 'GET':
        return render_template( 'superimposition.html', form=form )



    
@app.route( '/superx', methods=['POST']  )
def superx():
    print('super')
    form = SuperForm()
    
    pdb_1   = request.files[ 'pdb_1' ]
    pdb_2   = request.files[ 'pdb_2' ]
    ali   = request.files[ 'ali' ]
    chain_1 = request.form[ 'chain_1' ]
    chain_2 = request.form[ 'chain_2' ]
    seqid_1 = request.form[ 'seqid_1' ]
    seqid_2 = request.form[ 'seqid_2' ]

    check = True
    today = str(datetime.date.today()).replace('-', '')
    job = ""
    while check:
        job = today + '_' + randstr(5)
        outdir = app.config['USER_DATA_DIR'] + job + '/'
        print( outdir )
        if not os.path.exists( outdir ):
            os.mkdir( outdir )
            check = False
            
    pdb_name_1 = outdir + os.path.basename( pdb_1.filename )
    pdb_1.save( pdb_name_1 )
    pdb_name_2 = outdir + os.path.basename( pdb_2.filename )
    pdb_2.save( pdb_name_2 )

    ali_name = outdir + os.path.basename( ali.filename )
    ali.save( ali_name )
    

    cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_superimpose.py alignment: ' + pdb_name_1 + ' ' + chain_1 + ' ' + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + outdir + 'super.pdb'
    print( cmd )
    sys.stdout.flush()
    os.system( cmd )

    cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + outdir + 'super.pdb' + ' ' + chain_1 + ' ' + ali_name + ' ' + seqid_1 + ' 0.2  > ' + outdir + 'super_colored.pdb'
    print( cmd )
    os.system( cmd )

    cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + seqid_2 + ' 0.2  > ' + outdir + 'template_colored.pdb'
    print( cmd )
    os.system( cmd )

    
    return render_template( 'fullmenu.html', path=job, pdb="super_colored.pdb", chain=chain_1, pdb2="template_colored.pdb", chain2 = chain_2 )  # use same html ?


