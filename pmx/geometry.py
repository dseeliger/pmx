# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2016 by Daniel Seeliger
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
This file contains the tools to do basic operations on molecule geometries
like rotating around abound or computing deviation from planarity from a set of atoms.
A Rotation instance is built with two vectors.
Usage:
     >>> v1 = [1,2,3]
     >>> v2 = [2,3,4]
     >>> r = Rotation(v1,v2)  # create rotation object around v2-v1
     >>> v3 = [4,5,6]
     >>> v3 = r.apply(v3)    # rotate v3 around v2-v1
     
     

"""
from __future__ import print_function
from numpy import *
from atom import Atom
import _pmx as _p

import logging

logger = logging.getLogger()

class Rotation:

    def __init__(self,v1,v2):
        """ creates a rotation object
        around the vector v2-v1.

        :param v1: list or numpy array of length 3
        :param v2: list or numpy array of length 3

        """
        
        self.v1 = array(v1)
#        self.v2 = array(v2)
        self.v2 = [v2[0], v2[1], v2[2]] #array(v2)
        tmp = array(v2)
        self.vec = tmp-self.v1
        self.m1 = None
        self.m2 = None
        x = 1./linalg.norm(self.vec)
        self.norm_vec = self.vec*x
        self.__rm1()
        self.__rm2()

    def __rm1(self):

        a = self.norm_vec
        self.m1 =  [
            [ a[0]*a[0], a[0]*a[1], a[0]*a[2]], 
            [ a[1]*a[0], a[1]*a[1], a[1]*a[2]], 
            [ a[2]*a[0], a[2]*a[1], a[2]*a[2]] 
            ] 
        
    def __rm2(self):

        a = self.norm_vec
        self.m2 =  [
            [ 0.0, -a[2], a[1]],
            [ a[2], 0.0, -a[0]],
            [ -a[1], a[0], 0.0]
            ] 

    def apply(self,v, phi):
        return _p.apply_rotation( self, [v[0], v[1], v[2]], phi)



    
def vec_ang(v1,v2):
    """ compute the angle between to vectors """
    x1 = linalg.norm(v1)
    x2 = linalg.norm(v2)
    return arccos(inner(v1,v2)/(x1*x2))

def bb_super(mol1,mol2, use_orig_mc_coords = True):
    """ Superpose a protein residue on the backbone atoms of an existing one.

    :param mol1: molecule instance to fit on
    :param mol2: molecule instance to be fitted on mol1
    :param use_orig_mc_coords: mol2 inherits the mainchain coordinates of molecule 1

    """

    N1,CA1,C1 = mol1.fetchm(['N','CA','C'])
    N2,CA2,C2 = mol2.fetchm(['N','CA','C'])

    if( (mol1.resname=='GLY') or (mol2.resname=='GLY') or (mol2.resname[:2]=='G2') ):
        fit_atoms( [N1,CA1,C1], [N2,CA2,C2], mol2.atoms )
    else:
        N1,CA1,C1,CB1 = mol1.fetchm(['N','CA','C','CB'])
        N2,CA2,C2,CB2 = mol2.fetchm(['N','CA','C','CB'])
        fit_atoms( [N1,CA1,C1,CB1], [N2,CA2,C2,CB2], mol2.atoms )

#    fit_atoms( [N1,CA1,C1], [N2,CA2,C2], mol2.atoms )

    if use_orig_mc_coords:
        atom_set = ['N','CA','C','H','O','HA','HN']
        gly_atom_set = ['N','CA','C','H','O','HA1','HN']
        if mol1.resname == 'GLY':
            atoms1 = mol1.fetchm(gly_atom_set)
        else:
            atoms1 = mol1.fetchm(atom_set)
        if mol2.resname == 'GLY' or mol2.resname[:2] == 'G2':
            atoms2 = mol2.fetchm(gly_atom_set)
        else:
            atoms2 = mol2.fetchm(atom_set)
        assert len(atoms1) == len(atoms2), "%s -> %s" % ( '-'.join( map(lambda a: a.name, atoms1)),'-'.join( map(lambda a: a.name, atoms2)) ) 
        for atom1, atom2 in zip(atoms1, atoms2):
            atom2.x = atom1.x

def nuc_super(mol1,mol2):
    """ Superpose a nucleic acid molecule on another
    
    :param mol1: molecule to fit on
    :param mol2: molecule to be fitted

    """

    if mol1.resname[:2] in ['DT','DC','RC','RU']:
        fit1_atoms = ['C1\'', 'C6','N1','C2','C5','N3']
    else:
        fit1_atoms = ['C1\'', 'C8','N9','C4','N7','C5']

    if mol2.resname[:2] in ['DT','DC','RC','RU']:
        fit2_atoms = ['C1\'', 'C6','N1','C2','C5','N3']
    else:
        fit2_atoms = ['C1\'', 'C8','N9','C4','N7','C5']

    atoms1 = mol1.fetchm(fit1_atoms)
    atoms2 = mol2.fetchm(fit2_atoms)

    fit_atoms( atoms1, atoms2, mol2.atoms)



def planarity(atom_list):
    """ Computes the deviation from planarity for a set of atoms.
    It computes the optimal plain and subsequently computes the average deviation
    from the atoms to this plain.

    :param atom_list: list of atom instances
    :returns: average deviation from planarity
    """

    coords = map(lambda a: a.x, atom_list)
    plan = _p.planarity(coords)
    return plan

def apply_fit_R( atoms, R):
    """ Applies a fit matrix to a set of atoms

    :param atoms: list with atom instances
    :param R: fit matrix
    """

    for atom in atoms:
        x_old = map(lambda x: x, atom.x)
        for r in range(3):
            atom.x[r] = 0
            for c in range(3):
                atom.x[r]+=R[r][c]*x_old[c]

def fit(model1, model2, atom_names = []):
    """ Fit a model onto another. model2 is fitted on model1.
    
    :param model1: an instance containing an atomlist (model, chain or molecule)
    :param model2: an instance containing an atomlist (mode, chain or molecule)
    :param atom_names: optionally a list with atom names can be provides which are used for fitting.
    """


    if atom_names:
        subset1 = model1.fetch_atoms( atom_names )
        subset2 = model2.fetch_atoms( atom_names )
        cs1 = map(lambda a: a.x, subset1)
        cs2 = map(lambda a: a.x, subset2)
    else:
        cs1 = model1.coords()
        cs2 = model2.coords()

    assert( len(cs1) == len(cs2) )
    m = map(lambda x: 1., cs1) # dummy array
    v = _p.center_vec( cs1 )
    v2 = _p.center_vec( cs2 )
    R = _p.calc_fit_R(cs1, cs2, m)
    model2.translate( [-v2[0], -v2[1], -v2[2] ] )
    apply_fit_R( model2.atoms, R)
    model2.translate( v )


def fit_atoms( fit_atoms1, fit_atoms2, rot_atoms2 ):
    """ Computes a fit matrix for a list list of atoms to be fitted. Then the 
    fit matrix is applied to second list of atoms.

    :param fit_atoms1: list of atoms to fit on
    :param fit_atoms2: list of atoms to fit on fit_atoms2
    :param rot_atoms2: list of atoms to which transformation will be applied

    """
    cs1 = map(lambda a: a.x, fit_atoms1)
    cs2 = map(lambda a: a.x, fit_atoms2)
    assert len(cs1) == len(cs2)
    m = map(lambda x: 1., cs1) # dummy array
    v = _p.center_vec( cs1 )
    v2 = _p.center_vec( cs2 )
    R = _p.calc_fit_R(cs1, cs2, m)
    for atom in rot_atoms2:
        atom.x[0]-=v2[0]
        atom.x[1]-=v2[1]
        atom.x[2]-=v2[2]
    apply_fit_R( rot_atoms2, R)
    for atom in rot_atoms2:
        atom.x[0]+=v[0]
        atom.x[1]+=v[1]
        atom.x[2]+=v[2]
        



