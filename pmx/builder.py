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
This module comtains some functions to build protein
structures from scratch.

Most important:
- build_chain

Usage:
    >>> from pmx.builder import *
    >>> ch = build_chain('AAAAAA')    # builds an ordinary polyalanine
    >>> ch.write('polyala.pdb')       # store file
    >>> ch = build_chain('FGHRTCV',hydrogens = False) # build chain without H
    >>> ch = build_chain('FGHRTCV',ss='HHHHHHH') build as helix
    >>> ch = build_chain('FGHRTCV',dihedrals = ((phi1,psi1,omega1),(...)))
    build chain with defined dihedral angles
    
"""
import sys, os
from library import pmx_data_file
from geometry import *
from atom import Atom
from model import Model


def cross(x,y):
    return array([x[1]*y[2]
                  -x[2]*y[1],
                  x[2]*y[0]
                  -x[0]*y[2],
                  x[0]*y[1]
                  -x[1]*y[0]])


    
def add_bp(m, strand = None):
    if strand:
        N = len(strand)/2 
    else:
        N = 1
#    print N, len(strand)
    r = Rotation([0,0,0],[0,0,1])
    phi = N*36*pi/180.
    for atom in m.atoms:
        atom.x = r.apply(atom.x, phi)
    for atom in m.atoms:
        atom.x[2]+=N*.34

def make_3ter(r):
    # we add a proton at O5'

    c3, o3 = r.fetchm(['C3\'', 'O3\''])
    v = array(o3.x)-array(c3.x)
    normed = v*1./linalg.norm(v)
    newpos = array(o3.x) + normed*.1
    a = Atom(name='H3T', x = newpos)
    r.append(a)

def make_5ter(r):
    del r['O1P']
    del r['O2P']
    a = r.fetch_atoms('P')[0]
    a.name = 'H5T'
    h5t, o5 = r.fetchm(['H5T', 'O5\''])
    v = array(h5t.x)-array(o5.x)
    normed = v*1./linalg.norm(v)
    newpos = array(o5.x) + normed*.1
    h5t.x = newpos


def build_dna_strand(seq):

    dic = pmx_data_file('bp.pkl')

    seq = seq.lower()
    ss = []
    for x in seq:
        ss.append(x)
    new = []
    new.extend(dic[ss[0]].residues)
    ss.pop(0)
    while ss:
        newbp = ss.pop(0)
        newm = dic[newbp].copy()
        add_bp(newm, new)
        new.extend(newm.residues)
    chA = []
    chB = []
    a = 1
    b = 1
    for r in new:
        if r.chain_id == 'B':
            r.set_resid(a)
            a+=1
            chA.append(r)
        else:
            chB.append(r)
    chB.reverse()
    for r in chB:
        r.set_resid(b)
        b+=1
    mm = Model(residues=chA+chB)
    mm.chains[1].set_chain_id('C')
    mm.chains[0].set_chain_id('A')
    mm.chains[1].set_chain_id('B')
    for chain in mm.chains:
        r = chain.residues[0]
        r.set_resname(r.resname+'5')
        make_5ter(r)
        r = chain.residues[-1]
        r.set_resname(r.resname+'3')
        make_3ter(r)
    return mm


def get_fragments():
    dic = pmx_data_file('fragments.pkl')
    n = len(dic.keys())
    print >>sys.stderr,"pmx__> # Fragments loaded: %d" % n
    return dic


def read_pdb_with_connect(f):
    m = Model().read(f)
    m.nm2a()
    mol = m.residues[0]
    l = open(f).readlines()
    conn = []
    for line in l:
        if line[0:6] == 'CONECT':
            conn.append( [int(x) for x in line[6:].split()] )
    for lst in conn:
        if len(lst) > 1:
            atoms = mol.get_by_id(lst)
            a0 = atoms[0]
            for atom in atoms[1:]:
                a0.bonds.append(atom)
    return mol

def write_pdb_with_connect(mol, f, n = 1):
    if not hasattr(f,"write"):
        fp = open(f,"w")
    else:
        fp = f
    print >>fp, "MODEL%5d" % n
    for atom in mol.atoms:
        print  >>fp, atom
    for atom in mol.atoms:
        s= "CONECT%5d" % atom.id
        for a in atom.bonds:
            s+='%5d' % a.id
        print >>fp, s


def attach_group(atom, mol):
    master = atom.molecule
    
    bb = atom.bonds[0]
    R = mol.fetch_atoms('R#')[0]
    bR = R.bonds[0]
    diff = array(R.x)-array(bb.x)
    if bR.name[0] == 'C':
        l = 1.54
    elif bR.name[0] == 'S':
        l = 1.8
    elif bR.name[0] == 'O':
        l = 1.4
    else:
        l = 1.5
    for a in mol.atoms:
        a.x = array(a.x)-diff
    ang = R.angle(bR, atom)
    v1 = array(atom.x)-array(bb.x)
    v2 = array(bR.x)-array(R.x)
    cr = cross(v1, v2)
    x = 1./linalg.norm(cr)
    cr = cr*x
    cr = R.x + cr
    r = Rotation(cr, R.x)
    for a in mol.atoms:
        a.x = r.apply(a.x, ang)
    vec = array(bb.x) - array(bR.x)
    x = 1./linalg.norm(vec)
    for a in mol.atoms:
        a.x += vec*(1-x*l)
    for a in master.atoms:
        if atom in a.bonds:
            a.bonds.remove(atom)
    master.remove_atom(atom)
    bR.bonds.remove(R)
    mol.remove_atom(R)
    bR.bonds.append(bb)
    bb.bonds.append(bR)
    for a in mol.atoms:
        master.append(a)


