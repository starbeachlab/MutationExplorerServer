from flask import Flask, render_template, request, url_for, send_from_directory, jsonify, session,redirect
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, FileField, RadioField, HiddenField
from flask_bootstrap import Bootstrap
from werkzeug import secure_filename


import os, datetime, random, string
import requests, subprocess, sys
import basic_alignment as bali


app = Flask(__name__)
app.config['USER_DATA_DIR'] = "/disk/user_data/mutation_explorer/"
app.config['APP_PATH']      = "/home/hildilab/app/mutation_explorer_beta/"    ### ADJUST !!!!
app.config['ROSETTA_PATH']  = "/home/hildilab/dev/ext/rosetta/bin/"
app.config['SECRET_KEY']    = 'klsdf23*&%..'

bootstrap = Bootstrap( app )




@app.route('/' , methods=['GET','POST'] )
def index():
    return render_template( 'home.html')


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
        alignment_link = "https://www.bioinfo.mpg.de/AlignMeBeta/work/" + alignment_link.split("work/")[1]
        alignment = outdir + "alignment.aln"
        
        #alignment_link = app.config['USER_DATA_DIR'] + "alignment.aln"

        req = requests.get(alignment_link)
        with open(alignment, "w") as f:
            f.write(req.content)


        ###  determine chain and seqid
        #chain, seqid = ChainSeqidFromPDBandAlignment( outdir + secure_filename( file_conv.filename), alignment )        
        chain, seqid = ChainSeqidFromPDBandAlignment( outdir + conv_filename, alignment )        

        print( 'check', chain, seqid)

        error = ""
        if chain == '' or seqid == '':
            error = "The program was not able to find any sequence of the alignment in the PDB. Please note, that one chain has to match EXACTLY to one of the sequences within the alignment!"

        if not superimpose:
        
            #pdb_name = outdir + secure_filename(file_conv.filename)
            pdb_name = outdir + conv_filename
            aliname = alignment
            #new_name = secure_filename( file_conv.filename)[:-4] + '_consrv.pdb'
            new_name = conv_filename[:-4] + '_consrv.pdb'
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py ' + pdb_name + ' ' + chain + ' ' + aliname + ' ' + seqid + ' 0.2  > ' + outdir + new_name
            print( cmd )
            os.system( cmd )


            return render_template('overview.html', tag=tag, pdb=new_name, errmsg=error)
            
        
        chain_1 = chain
        seqid_1 = seqid
        #pdb_name_1 = outdir + secure_filename(file_conv.filename)
        pdb_name_1 = outdir + conv_filename
        ali_name = alignment

        #pdb_name_2 = outdir + secure_filename(file_super.filename)
        pdb_name_2 = outdir + super_filename

        chain_2, seqid_2 = ChainSeqidFromPDBandAlignment( outdir + super_filename, alignment )        

        if chain_2 == '' or seqid_2 == '':
            error = "The program was not able to find any sequence of the alignment in the PDB. Please note, that one chain has to match EXACTLY to one of the sequences within the alignment!"
            return render_template('overview.html', errmsg=error);
        
        if seqid_1 == '0' and seqid_2 == '1':
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_superimpose.py alignment: ' + pdb_name_1 + ' ' + chain_1 + ' ' + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + outdir + 'super.pdb'
            
            print( cmd )
            sys.stdout.flush()
            os.system( cmd )

            #new_name_1 = secure_filename(file_conv.filename)[:-4] + '_consvr_suprmp.pdb'
            new_name_1 = conv_filename[:-4] + '_consvr_suprmp.pdb'
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + outdir + 'super.pdb' + ' ' + chain_1 + ' ' + ali_name + ' ' + seqid_1 + ' 0.2  > ' + outdir + new_name_1
            print( cmd )
            os.system( cmd )

            #new_name_2 = secure_filename( file_super.filename)[:-4] + "_consvr.pdb" 
            new_name_2 = super_filename[:-4] + "_consvr.pdb" 
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + pdb_name_2 + ' ' + chain_2 + ' ' + ali_name + ' ' + seqid_2 + ' 0.2  > ' + outdir + new_name_2
            print( cmd )
            os.system( cmd )

        elif seqid_1 == '1' and seqid_2 == '0':
            
            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_superimpose.py alignment: ' + pdb_name_2 + ' ' + chain_2 + ' ' + pdb_name_1 + ' ' + chain_1 + ' ' + ali_name + ' ' + outdir + 'super.pdb'
            
            print( cmd )
            sys.stdout.flush()
            os.system( cmd )

            #new_name_2 = secure_filename(file_super.filename)[:-4] + '_consvr_suprmp.pdb'
            new_name_2 = super_filename[:-4] + '_consvr_suprmp.pdb'

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + outdir + 'super.pdb' + ' ' + chain_2 + ' ' + ali_name + ' ' + seqid_2 + ' 0.2  > ' + outdir + new_name_2
            print( cmd )
            os.system( cmd )

            #new_name_1 = secure_filename(file_conv.filename)[:-4] + '_consvr.pdb'
            new_name_1 = conv_filename[:-4] + '_consvr.pdb'

            cmd =  'python3 ' + app.config['APP_PATH'] + 'pdb_conservation.py '  + pdb_name_1 + ' ' + chain_1 + ' ' + ali_name + ' ' + seqid_1 + ' 0.2  > ' + outdir + new_name_1
            print( cmd )
            os.system( cmd )

        
        return render_template('overviewSuper.html', tag=tag, pdb1=new_name_1, pdb2=new_name_2, errmsg=error)



    return render_template('home.html')

# this can handle both cases of single and pair of PDBs: same could be done for overview
@app.route('/menu/<tag>/<pdb1>')
@app.route('/menu/<tag>/<pdb1>/<pdb2>')
def menu(tag, pdb1, pdb2 = ""):
    return render_template( 'fullmenu.html', path=tag, pdb=pdb1,  pdb2=pdb2)

# same as get() below...
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
    
    #print( cmd)
    aln_seqs =  bali.ReadAlignment( ali )[1]
    #print(aln_seqs)
    aln_seqs = [ f.replace('-','') for f in aln_seqs ]
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
        print( pdb_seqs)
        print( aln_seqs)
    #print( 'chain:', chain, 'seqid:', seqid)
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




###  MUTATION EXPLORER ITSELF  ##########

class MutExForm( FlaskForm):
    pdb = FileField( 'Upload PDB from your machine' )
    idstr = StringField( 'OR download from server by identifier:' , render_kw={'placeholder':'PDB or alpha Fold ID'})
    source = SelectField( 'Select server:', choices=[('A', 'RCSB (PDB)'), ('B','OPM - membrane proteins'), ('C', 'alpha Fold')], default='A')
    opensession = StringField( 'OR open session by ID')
    mutations = StringField( 'Define mutations for W:T')
    display = SelectField( 'Display PDB colored by energy')
    download = SelectField( 'Download files from this session', choices=[('A','(none)'), ('B', 'ZIP'),('C', 'TAR GZIPPED')] , default='A')
    submit = SubmitField( 'Go')


    
"""
session / sammlung von cookies / damit folgeseiten auf fruehere werte zugreifen
- sessionid
- history 
- protein
- nr
"""


@app.route('/mutantX' , methods=['GET','POST'] )
def mutant():
    form = MutExForm()
    print('mutex')
    if request.method == 'POST':
        #if form.validate_on_submit():
        if 'sessionid' in session:
            print( 'sessionid', session.get('sessionid'))
            sessionid = session.get( 'sessionid')
            outdir = app.config['USER_DATA_DIR'] + sessionid + '/'
            tmpdir = outdir + 'tmp/'
            
        ##### pdb upload #####
        pdb_file = secure_filename( form.pdb.data.filename)  # local file upload
        idstr = form.idstr.data   # upload via RCSB or OPM or alpha Fold
        if pdb_file != "" or idstr != "":  # if new session, create session id and directory, add to history
            sessionid = str( random.randint(0, 9999999))
            sessionid = "66827"
            outdir = app.config['USER_DATA_DIR'] + sessionid + "/"
            tmpdir = outdir + 'tmp/'
            #while os.path.exists(outdir):
            #    sessionid = str(random.randint(0, 9999999))
            #    outdir = app.config['USER_DATA_DIR'] + sessionid + "/"
            #os.mkdir(outdir)
            session['nr'] = 0
            if 'history' in session and not sessionid in session.get( 'history'):
                session['history'] = session.get('history') + ' ' + sessionid
            else:
                session['history'] = sessionid
        if pdb_file != "":   # upload local pdb file
            print( 'upload local file', pdb_file)
            pdb.save( os.path.join( tmpdir , pdb_file))
            session['protein'] = pdb_file[:-4]
            form.pdb.data = ''
        elif idstr != "":   # upload via server
            print( 'upload from server:', idstr)
            if idstr[:-4] == ".pdb":
                idstr = idstr[:-4]
            session['protein'] = idstr
            idstr += ".pdb"
            if form.source.data == 'A':
                url = "https://files.rcsb.org/download/" + idstr
            elif form.source.data == 'B':
                url = "https://opm-assets.storage.googleapis.com/pdb/" + idstr
            elif form.source.data == 'C':
                url = "https://alphafold.ebi.ac.uk/files/" + idstr 
            print( url)
            req = requests.get( url )
            with open( os.path.join( tmpdir, idstr) , 'w') as f:
                f.write( req.content )
            form.idstr.data = ''
        elif form.opensession.data != "":  # open previous session
            print( 'open previous:', form.opensession.data)
            sessionid = form.opensession.data
            session['sessionid'] = sessionid
            outdir = app.config['USER_DATA_DIR'] + sessionid + '/'
            tmpdir = outdir + 'tmp'
            maxi = 0
            if os.path.exists( outdir):
                # add message if ID was not found
                for f in os.listdir( outdir ):
                    if '_m' in f:
                        print( 'mut: ', f, f.index('_m'))
                        kid = int( f[ f.index('_m') + 2 : -4])
                        maxi = max( maxi, kid)
                    elif '_wt' in f:
                        print( 'wt:', f)
                        session['protein'] = f[:-7]
            session['nr'] = maxi
            form.opensession.data = ''
            print ( 'found', session.get( 'nr'), 'mutations of ', session.get('protein'))
            
        # if file was uploaded (new sesssion), filter pdb and write energies into PDBID_wt.pdb
        if form.pdb.data.filename != "" or form.idstr.data != "":
            # filter (ONLY HETATMs at the moment!!!!) TODO: undef RESTYPES !!!
            #print(' CWD: ', os.getcwd())
            print( 'form:',  form.pdb.data.filename ,  form.idstr.data)
            pdb_file = session.get( 'protein')
            if not os.path.exists( tmpdir ):
                os.mkdir( tmpdir )
            filt = tmpdir + pdb_file + "_filtered.pdb"
            rose = filt[:-4] + '.out'
            if not os.path.exists( filt ):
                cmd =  app.config['APP_PATH'] + "pdb_filter.py " + tmpdir + pdb_file + '.pdb ' + filt
                print( cmd )
                os.system( cmd )
            # calc energies (rosetta)
            if not os.path.exists( rose ):
                cmd = app.config['ROSETTA_PATH'] + "per_residue_energies.static.linuxgccrelease -in:file:s " + filt + " -out:file:silent " + rose
                print( cmd )
                os.system(cmd)
            # write energies to pdb
            if not os.path.exists( tmpdir + pdb_file + "_wt.pdb" ):
                cmd =  app.config['APP_PATH'] + "pdb_write_rosetta_energies.py " + filt + " " + rose + " "  + outdir + pdb_file + "_wt.pdb"
                print( cmd)
                os.system(cmd)

            

        #### MUTATION SECTION
        if form.mutations.data != "":
            nr = session.get('nr')
            nr += 1
            session['nr'] = nr
            mstr = 'm' + str(nr)
            print( 'mutations:', form.mutations.data)
            muts = form.mutations.data.split(',')
            prev_chains = []
            with open( tmpdir + 'resfile' + str(nr) + '.txt', 'w' ) as w:
                w.write( 'NATRO\n' )
                w.write( 'start\n')
                for m in muts:
                    c = m.split(':')
                    res1 = c[-1][0]
                    res2 = c[-1][-1]
                    chain = ''
                    if len(c) > 1:
                        chain = c[0]
                    if chain not in prev_chains:
                        w.write( '* ' + chain + ' NATAA \n' )
                        prev_chains.append( chain )
                        
                    resid = c[-1][1:-1]
                    print( 'mutate ' + res1 + ' ' + chain + ' ' + resid + ' to ' + res2 )
                    #print( resid, chain, 'PIKAA', res2, file=w)
                    w.write( resid + ' ' + chain + ' PIKAA ' + res2 + '\n')
                    
            # call rosetta design / mutate
            protein = session.get('protein')
            wt = outdir + protein + "_wt.pdb"
            cmd = app.config['ROSETTA_PATH'] + "fixbb.static.linuxgccrelease -in:file:s " + wt + " -resfile " + tmpdir + 'resfile' + str(nr) + '.txt -nstruct 1 -linmem_ig 10 -out:pdb  -out:prefix ' + tmpdir + mstr + '_  >& ' + outdir + 'tmp/log.txt'
            print( cmd )
            os.system(cmd)

            # find mutant
            mutant = ''
            for f in os.listdir( tmpdir ):
                if mstr + '_' in f and '_0001.pdb' in f:
                    mutant = f

            # per residue energy
            if mutant != '':
                rose = tmpdir  + protein + '_' + mstr + '.out'
                cmd = app.config['ROSETTA_PATH'] + "per_residue_energies.static.linuxgccrelease -in:file:s " + tmpdir + mutant + " -out:file:silent " + rose
                print( cmd  )
                os.system(cmd)

                # write energies into bfactor of pdb
                cmd =  app.config['APP_PATH'] + "pdb_write_rosetta_energies.py " + tmpdir + mutant + " " + rose + " "  + outdir + protein + '_' + mstr + ".pdb"
                print( cmd)
                os.system(cmd)
                form.mutations.data = ''
            ### END of mutation block ###

            
            
        ###  DISPLAY SECTION  #####
        choices = []
        file_names = []
        protein = session.get('protein')
        c = 0
        for f in os.listdir( outdir ):
            if f[-4:] == ".pdb":
                choices.append(( f[:-4], f[:-4]) )
                file_names.append( f)
                c+=1
        current_choice = form.display.data
        current = ''
        form.display.choices = choices
        if current_choice != '' and current_choice != 'A':
            current = current_choice + '.pdb'
        else:
            # display latest mutation by default
            default = ''
            if session.get( 'nr') > 0:
                default = protein + '_m' + str( session.get('nr') ) 
            else:
                default = protein + '_wt'
                form.display.default = default
                print( 'form.display.default' ,  form.display.default , default )

                #current = pdb_file[:-4] + "_wt.pdb"
                current = ''
                if default != '':
                    current =  default + ".pdb"

        print( 'current pdb to display:' , current)


        ### DOWNLOAD SECTION ###
        if form.download.data != 'A':
            protein = session.get( 'protein')
            ### create thread / def to avoid issue with return statement at the end ??
            print('Download Section')
            if form.download.data == 'B':
                out_name = protein + '.zip'
                cmd = "cd " + outdir + "; zip " + out_name
                for y in file_names:
                    cmd += ' ' + y
                print(cmd)
                os.system(cmd)
                return get( sessionid , out_name )
            elif form.download.data == 'C':
                out_name = protein + '.tgz'
                cmd = "cd " + outdir + "; tar -czPf " + out_name
                for y in file_names:
                    cmd += ' ' + y
                print(cmd)
                os.system(cmd)
                return redirect( url_for( 'get', mydir=sessionid , pdb=out_name ))

            
        form.download.data = 'A'
        

        ### HERE???
        if sessionid != '':
            session['sessionid'] = sessionid
        else:
            sessionid = session.get('sessionid')
        if 'history' not in session:
            session['history'] = sessionid
        print( 'rebuild', outdir, sessionid, pdb_file)
        print( 'session', session)
        return render_template( 'mutations.html', form=form, mydir=sessionid, pdb=current, history=session['history']) 

    # GET
    form.display.choices = [('A','Nothing to display')]
    session.clear()
    print( 'request', request )
    return render_template('mutations.html', form=form)




#https://alphafold.ebi.ac.uk/files/AF-B8JHS2-F1-model_v2.pdb
 
