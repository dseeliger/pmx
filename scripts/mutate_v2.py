#!/usr/bin/env python
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

__doc__="""
Program to insert mutated residues in structure files for
free energy simulations.
"""

import sys,os
from pmx import *
from pmx.parser import *
from pmx import library
from pmx.mutdb import *
from pmx.geometry import *

class UnknownResidueError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)
class RangeCheckError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)
class mtpError(Exception):
    def __init__(self,s):
        self.s = s
    def __str__(self):
        return repr(self.s)



ext_one_letter = {
    'ALA':'A',
    'ARG':'R',
    'ASN':'N',
    'ASP':'D',
    'ASPH':'B',
    'ASPP':'B',
    'ASH':'B',
    'CYS':'C',
    'CYS2':'C',
    'CYN':'C',
    'CYX':'C',
    'CYM':'C',
    'CYSH':'C',
    'GLU':'E',
    'GLUH':'J',
    'GLUP':'J',
    'GLH':'J',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H',
    'HIE':'X',
    'HISE':'X',
    'HSE':'X',
    'HIP':'Z',
    'HSP':'Z',
    'HISH':'Z',
    'HID':'H',
    'HSD':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K',
    'LYSH':'K',
    'LYP':'K',
    'LYN':'O',
    'LSN':'O',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V',
}

def check_residue_name( res ):
    if res.resname == 'LYS':
        if res.has_atom( 'HZ3'):
            res.set_resname('LYP')
    elif res.resname == 'HIS':
        if res.has_atom('HD1') and \
           res.has_atom('HE2'):
            res.set_resname('HIP')
        elif res.has_atom('HD1') and not \
                 res.has_atom('HE2'):
            res.set_resname( 'HID' )
        elif not res.has_atom('HD1') and  \
                 res.has_atom('HE2'):
            res.set_resname( 'HIE' )
    elif res.resname == 'ASP':
        if res.has_atom('HD2'):
            res.set_resname('ASH')
    elif res.resname == 'GLU':
        if res.has_atom('HE2'):
            res.set_resname('GLH')
    elif res.resname == 'CYS':
        if not res.has_atom('HG'):
            print >>sys.stderr,' Cannot mutate SS-bonded Cys %d' % res.id
    
def check_OPLS_LYS( res ):
    if res.has_atom( 'HZ3'):
        return('K')
    else:
	return('O')

def get_restype(r):
    if r.resname in ['DA','DT','DC','DG']:
        return 'DNA'
    elif r.resname in ['RA','RU','RC','RG']:
        return 'RNA'
    else: return 'PEPTIDE'

def read_script(fn):
    return read_and_format(fn,"is")

def int_input():
    inp = raw_input()
    try:
        inp = int(inp)
        return inp
    except:
        print 'You entered "%s" -> Try again' % inp
        return None

def check_residue_range(m, idx):
    valid_ids = range(1, len(m.residues)+1)
    if idx not in valid_ids: return False
    return True

def select_residue(m):
    valid_ids = range(1, len(m.residues)+1)
    print '\nSelect residue to mutate:'
    for i,r in enumerate(m.residues):
        if r.resname not in library._ions+library._water: 
            sys.stdout.write('%6d-%s-%s' % (r.id,r.resname,r.chain_id))
            if r.id % 6 == 0: print
    print
    selected_residue_id = None
    while not selected_residue_id:
        sys.stdout.write('Enter residue number: ')
        selected_residue_id = int_input()
        if selected_residue_id is not None and selected_residue_id not in valid_ids:
            print 'Residue id %d not in range %d-%d -> Try again' % (selected_residue_id,1,len(residues))
            selected_residue_id = None
    return selected_residue_id

def select_mutation(m, selected_residue_id, ffpath):

    residue = m.residues[selected_residue_id - 1]
    if get_restype(residue) == 'PEPTIDE':
        return select_aa_mutation(residue,ffpath)
    elif get_restype(residue) in ['DNA','RNA']:
        return select_nuc_mutation(residue)

def select_nuc_mutation(residue):
    aa = None
    print '\nSelect new base for %s-%s: ' % (residue.id,residue.resname) 
    sys.stdout.write('One-letter code: ')
    while aa is None:
        aa = raw_input().upper()
        if get_restype(residue) == 'DNA' and aa not in ['A','C','G','T']:
            sys.stdout.write('Unknown DNA residue "%s"!\nOne-letter code: ' % aa)
            aa = None
        elif get_restype(residue) == 'RNA' and aa not in ['A','C','G','U']:
            sys.stdout.write('Unknown RNA residue "%s"!\nOne-letter code: ' % aa)
            aa = None
        if aa:
            print 'Will apply mutation %s->%s on residue %s-%d' % (residue.resname[1],aa,residue.resname,residue.id)
        return aa

def select_aa_mutation(residue,ffpath):
    check_residue_name( residue )
    print '\nSelect new amino acid for %s-%s: ' % (residue.id,residue.resname) 
    sys.stdout.write('Three- or one-letter code (or four-letter for ff specific residues): ')
    if residue.resname in ['HIE','HISE','HSE']: rol = 'X'
    elif residue.resname in ['HIP','HISH','HSP']: rol = 'Z'
    elif residue.resname in ['GLH','GLUH','GLUP']: rol = 'J'
    elif residue.resname in ['ASH','ASPH','ASPP']: rol = 'B'
    elif residue.resname in ['LYN','LYS','LSN']: rol = 'O'
    else:
        rol = library._one_letter[residue.resname]
    aa = None
    ol = library._aacids_dic.keys()
    tl = library._aacids_dic.values()
    if('amber' in ffpath):
	    ol = library._aacids_ext_amber.keys()
	    tl = library._aacids_ext_amber.values()
    if('opls' in ffpath):
            ol = library._aacids_ext_oplsaa.keys()
            tl = library._aacids_ext_oplsaa.values()
    if('charmm' in ffpath):
            ol = library._aacids_ext_charmm.keys()
            tl = library._aacids_ext_charmm.values()

    while aa is None:
        aa = raw_input().upper()
        if len(aa) != 1 and len(aa)!=3 and len(aa)!=4:
            sys.stdout.write('Nope!\nThree- or one-letter code (or four-letter for ff specific residues): ')
            aa = None
        elif (len(aa) == 1 and aa not in ol+['B','J','O','X','Z']) or (len(aa)==3 and aa not in tl) or (len(aa)==4 and aa not in tl):
            sys.stdout.write('Unknown aa "%s"!\nThree- or one-letter code (or four-letter for ff specific residues): ' % aa)
            aa = None
        if aa and (len(aa)==3 or len(aa)==4): aa = ext_one_letter[aa]
    print 'Will apply mutation %s->%s on residue %s-%d' % (rol,aa,residue.resname,residue.id)
    return aa


def interactive_selection(m,ffpath):
    residue_id = select_residue(m)
    mutation = select_mutation(m, residue_id, ffpath )
    return residue_id, mutation
    
def ask_next():
    sys.stdout.write('\nApply another mutation [y/n]? ')
    res = raw_input().lower()
    if res == 'y': return True
    elif res == 'n': return False
    else: return ask_next()

def convert_aa_name( aa ):
    if len(aa) == 1: return aa.upper()
    elif len(aa) == 3: return ext_one_letter[aa.upper()]
    elif len(aa) == 4: return ext_one_letter[aa.upper()]
    else: raise UnkownResidueError(aa)
    
def rename_to_match_library(res):
    name_hash = {}
    atoms = res.atoms
    for atom in atoms:
	foo = atom.name
	## for serine
	if (atom.resname == 'SER') and (atom.name == 'HG1'):
	    atom.name = 'HG'
        if ('S2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
	## for cysteine
        if (atom.resname == 'CYS') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('C2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
#	print atom.resname,atom.name
	name_hash[atom.name] = foo
    return name_hash

def rename_back( res, name_hash ):
    for atom in res.atoms:
        atom.name = name_hash[atom.name]

def set_conformation(old_res, new_res, rotdic):
    old_res.get_real_resname()
    dihedrals = library._aa_dihedrals[old_res.real_resname]
    for key, lst in rotdic.items():
        new = new_res.fetchm(lst)
        rotdic[key] = new
    chis = []
    for key in rotdic.keys():
        at1,at2 = key.split('-')
        for d in dihedrals:
            if d[1] == at1 and d[2] == at2 \
                   and d[-1] != -1:
                chis.append(d)
    for d in chis:
        atoms = old_res.fetchm(d[:4])
        phi = atoms[0].dihedral(atoms[1], atoms[2], atoms[3])
        atoms2 = new_res.fetchm(d[:4])
        phi2 = atoms2[0].dihedral(atoms2[1], atoms2[2], atoms2[3])
        diff = phi-phi2
        a1,a2 = new_res.fetchm(d[1:3])
        key= a1.name+'-'+a2.name
        atoms = rotdic[key]
        rot = Rotation(a1.x,a2.x)
        for atom in atoms:
            atom.x = rot.apply(atom.x,diff)
    for atom in new_res.atoms:
        if atom.name[0] != 'D':
            atom.x = old_res[atom.name].x


def apply_nuc_mutation(m, residue, new_nuc_name, mtp_file):

    hybrid_residue_name = residue.resname+new_nuc_name
    print 'log_> Residue to mutate: %d | %s | %s ' % ( residue.id, residue.resname, residue.chain_id)
    print 'log_> Mutation to apply: %s->%s' % (residue.resname[1], new_nuc_name)
    print 'log_> Hybrid residue name: %s' % hybrid_residue_name
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)
    hybrid_res.nm2a()
    nuc_super( residue, hybrid_res )
    for atom in hybrid_res.atoms:
        if atom.name[0] != 'D':
            atom.x = residue[atom.name].x
    m.replace_residue( residue, hybrid_res)
    print 'log_> Inserted hybrid residue %s at position %d (chain %s)' %\
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id)

    
def apply_aa_mutation(m, residue, new_aa_name, mtp_file, bStrB, infileB):

    if residue.resname == 'ILE': rename_ile( residue )
    olkey = convert_aa_name( residue.resname )

    # olkey should contain the correct one letter name of the WT residue
    # however, due to the different namings of the residues in the FFs
    # Lys needs to be checked once again: in OPLS Lys is non-protonated, while in the other FFs it is protonated
    if ('opls' in mtp_file) and ('LYS' in residue.resname):
        olkey = check_OPLS_LYS( residue )

    hybrid_residue_name = olkey+'2'+new_aa_name
    print 'log_> Residue to mutate: %d | %s | %s ' % ( residue.id, residue.resname, residue.chain_id)
    print 'log_> Mutation to apply: %s->%s' % (olkey, new_aa_name)
    print 'log_> Hybrid residue name: %s' % hybrid_residue_name
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)
    #hybrid_res.nm2a()
    bb_super(residue, hybrid_res )

    ## VG rename residue atoms
    hash1 = rename_to_match_library(residue)
    hash2 = rename_to_match_library(hybrid_res)
    set_conformation(residue, hybrid_res, rotdic)
    if bStrB:
	print "log_> Set Bstate geometry according to the provided structure"
   	mB = Model(infileB)
   	rename_atoms_to_gromacs( mB )
	mB.nm2a()
	residueB = mB.residues[residue.id-1] 
    	bb_super(residue, residueB )
	for atom in hybrid_res.atoms:
            if atom.name[0] == 'D':
	        for atomB in residueB.atoms:
		    if atomB.name == hybrid_res.morphes[atom.name]['n1']:
	 	        atom.x = atomB.x
    rename_back(residue,hash1)
    rename_back(hybrid_res,hash2)
    ## VG rename residue atoms back

    m.replace_residue( residue, hybrid_res)
    print 'log_> Inserted hybrid residue %s at position %d (chain %s)' %\
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id)


def apply_mutation(m, mut, mtp_file, bStrB, infileB):
    residue_id = mut[0]
    if not check_residue_range(m, residue_id):
        raise RangeCheckError(residue_id)
    residue = m.residues[residue_id - 1]
    if get_restype(residue) == 'PEPTIDE':
        new_aa_name = convert_aa_name( mut[1] )
        apply_aa_mutation(m, residue, new_aa_name, mtp_file, bStrB, infileB)
    elif get_restype(residue) in ['DNA','RNA']:
        new_nuc_name = mut[1].upper()
        apply_nuc_mutation(m, residue, new_nuc_name, mtp_file)
        
    
def get_hybrid_residue(residue_name, mtp_file = 'ffamber99sb.mtp'):
    print 'log_> Scanning database for %s ' % residue_name
    resi, bonds, imps, diheds, rotdic = read_mtp_entry(residue_name, filename = mtp_file, version = 'new')
    if len(resi.atoms) == 0:
        raise mtpError("Hybrid residue %s not found in %s" % (residue_name, mtp_file) )
    return resi, bonds, imps, diheds, rotdic


    
def rename_ile(residue):
    dic = {'CD':'CD1',
           'HD1':'HD11',
           'HD2':'HD12',
           'HD3':'HD13'
           }           
    for key, value in dic.items():
        try:
            atom = residue[key]
            atom.name = value
        except:
            pass
       
def rename_atoms_to_gromacs( m ):
    for atom in m.atoms:
        if atom.name[0].isdigit():
            atom.name =  atom.name[1:]+atom.name[0]
    



def get_restype(r):
    if r.resname in ['DA','DT','DC','DG']:
        return 'DNA'
    elif r.resname in ['RA','RU','RC','RG']:
        return 'RNA'
    else: return 'PEPTIDE'


def get_ff_path( ff ):
    ff_path = None
    if not os.path.isdir(ff):
	### VG ###
#        gmxlib = os.environ.get('GMXDATA')
	gmxlib = os.environ.get('GMXLIB')
	p = os.path.join(gmxlib,ff)
	### VG ###
        if not os.path.isdir(p):
            print >>sys.stderr,' Error: forcefield path "%s" not found' % ff
        else:
            ff_path = p
    else:
        ff_path = ff
    print 'Opening forcefield: %s' % ff_path
    return ff_path

    
def main(argv):

   options = [
        Option( "-resinfo", "bool", False, "print a 3-letter -> 1-letter residue list"),
##         Option( "-r", "rvec", [1,2,3], "some string"),
##         Option( "-b", "bool", True, "bool"),
##         Option( "-r2", "rvec", [1,2,3], "some vector that does wonderful things and returns always segfaults")
        ]
    
   files = [
       FileOption("-f", "r",["pdb","gro"],"protein.pdb", "input structure file"),
       FileOption("-fB", "r",["pdb","gro"],"proteinB.pdb", "input structure file of the Bstate (optional)"),
       FileOption("-o", "w",["pdb","gro"],"out.pdb", "output structure file"),
       FileOption("-ff", "dir",["ff"],"amber99sbmut", "path to mutation forcefield"),
       FileOption("-script", "r",["txt"],"mutations.txt", "text file with mutations to insert"),
       ]
    
   help_text = ('This script applies mutations of residues in a structure file ',
                'for subsequent free energy calculations like FEP, TI, etc.',
                'The mutation information and dummy placements are taken from',
                'the hybrid residue database "mutres.mtp". The best way to use',
                'this script is to take a pdb/gro file that has been written with pdb2gmx',
                'with all hydrogen atoms present.'
                'The program can either be executed interactively or via script.',
                'The script file simply has to consist of "resi_number target_residue." pairs.',
                'The script uses an extended one-letter code for amino acids to account for',
                'different protonation states. Use the -resinfo flag to print the dictionary.',
                'Currently available force fields:',
                '    - amber99sbmut (Hornak et al, 2006)',
                '    - amber99sb-star-ildn-mut (Best & Hummer, 2009; Lindorff-Larsen et al, 2010)',
                '    - charmm22starmut.ff (Piana et al, 2011)',
                '    - charmm36mut (Best et al, 2012)',
                '    - oplsaamut (Jorgensen et al, 1996; Kaminski et al, 2001)',
                '',
                '',
                'Please cite:',
		'Vytautas Gapsys, Servaas Michielssens, Daniel Seeliger and Bert L. de Groot.',
		'Automated Protein Structure and Topology Generation for Alchemical Perturbations.',
		'J. Comput. Chem. 2015, 36, 348-354. DOI: 10.1002/jcc.23804',
		'',
		'Old pmx (pymacs) version:',
                'Daniel Seeliger and Bert L. de Groot. Protein Thermostability Calculations Using',
                'Alchemical Free Energy Simulations, Biophysical Journal, 98(10):2309-2316 (2010)',
                '',
                '',
                '',
                )

    
   cmdl = Commandline( argv, options = options,
                       fileoptions = files,
                       program_desc = help_text,
                       check_for_existing_files = False )
    

   if cmdl['-resinfo']:
       print 'Residue dictionary:'
       lst = ext_one_letter.items()
       lst.sort(lambda a,b: cmp(a,b))
       for key, val in lst:
           print "%5s %4s" % (key, val)
       sys.exit(0)

   bStrB = False
   infile = ''
   if cmdl['-fB']:
	bStrB = True
	infileB = cmdl['-fB']

   ffpath = get_ff_path(cmdl['-ff'])
   mtp_file = os.path.join( ffpath,'mutres.mtp')
   infile = cmdl['-f']
   m = Model(infile)
   rename_atoms_to_gromacs( m )
#   m.write('ll.pdb')
   m.nm2a()
#   m.rename_atoms()
   mutation_list = []
   if cmdl.opt['-script'].is_set:
       mutations_to_make = read_script( cmdl['-script'] )
       for mut in mutations_to_make:
	   check_residue_name( m.residues[ mut[0]-1 ] )
           apply_mutation( m, mut, mtp_file, bStrB, infileB )
   else:
       do_more = True
       while do_more:
           mutation = interactive_selection(m,ffpath)
           apply_mutation( m, mutation, mtp_file, bStrB, infileB )
           if not ask_next(): do_more = False
       

   m.write(cmdl['-o'])
   print
   print 'mutations done...........'
   print
if __name__=='__main__':
    main(sys.argv)

