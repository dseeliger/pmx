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

import sys, os

from ctypes import cdll, c_int, c_float, c_double
from ctypes import c_char_p,POINTER,c_int,byref

c_real = c_float
rvec=c_real*3
matrix=c_real*3*3


class Trajectory:

    def __init__(self, filename):
        self.filename = filename
        self.libgmx_path = self.__check_gmxlib()
        self.libgmx = cdll.LoadLibrary(self.libgmx_path)
        self.fp = None
        self.natoms = c_int()
        self.step = c_int()
        self.time = c_real()
        self.prec = c_real()
        self.bOK = c_int()
        self.box =  matrix()
        self.x = POINTER(rvec)()
        self.__have_first_frame = False
        self.open_xtc(filename)
        self.read_first_xtc()
        
    def __check_gmxlib(self):
        p = os.environ.get('GMX_DLL')
        if not p:
            print >> sys.stderr, "pmx_Error_> Path to Gromacs shared libraries is not set (GMX_DLL)"
            print >> sys.stderr, "pmx_Error_> Cannot load \"libgmx.so\""
            sys.exit(1)
        return os.path.join(p,'libgmx.so')
    
      
    def open_xtc(self, filename):
	self.libgmx.open_xtc.restype = POINTER(c_char_p)
        self.fp = self.libgmx.open_xtc(filename, "r")
        
    def close_xtc(self):
        self.libgmx.close_xtc(self.fp)
        
    def read_first_xtc(self):

        ret = self.libgmx.read_first_xtc(self.fp,byref(self.natoms),byref(self.step),byref(self.time),self.box,byref(self.x),byref(self.prec),byref(self.bOK))
        if not ret:
            print >> sys.stderr, "pmx_Error_> in Trajectory.read_first_xtc()"
            sys.exit(1)
        self.__have_first_frame = True

    def read_next_xtc(self):
        if not self.__have_first_frame:
            print >> sys.stderr, "pmx_Error_> First frame not read"
            sys.exit(1)
        ret = self.libgmx.read_next_xtc(self.fp,self.natoms,byref(self.step),byref(self.time),self.box,self.x,byref(self.prec),byref(self.bOK))
        return ret

    def update_box( self, box ):
        for i in range(3):
            for k in range(3):
                box[i][k] = self.box[i][k]

    def update_atoms( self, atom_sel ):
        assert len(atom_sel.atoms) == self.natoms.value
        for i, atom in enumerate(atom_sel.atoms):
            if atom_sel.unity == 'A':
                atom.x[0] = self.x[i][0]*10
                atom.x[1] = self.x[i][1]*10
                atom.x[2] = self.x[i][2]*10
            else:
                atom.x[0] = self.x[i][0]
                atom.x[1] = self.x[i][1]
                atom.x[2] = self.x[i][2]

    def update( self, atom_sel ):
        self.update_atoms(atom_sel )
        self.update_box( atom_sel.box )

    def get_box(self):
        return [
            [self.box[0][0], self.box[0][1], self.box[0][1]],
            [self.box[1][0], self.box[1][1], self.box[1][1]],
            [self.box[2][0], self.box[2][1], self.box[2][1]],
            ]

    def get_natoms(self):
        return self.natoms.value
    def get_time(self):
        return self.time.value
    def get_step(self):
        return self.step.value
    def get_prec(self):
        return self.prec.value
    def get_x(self, unity='nm'):
        x = []
        for i in range(self.get_natoms):
            if unity == 'A':
                x.append( [self.x[i][0]*10, self.x[i][1]*10, self.x[i][2]*10 ])
            else:                
                x.append( [self.x[i][0], self.x[i][1], self.x[i][2] ])
        return x
    
    def __str__(self):

        s = '< Trajectory: %s | natoms = %d | step = %d | time = %8.2f >' % (self.filename, self.get_natoms(), self.get_step(), self.get_time())
        return s
    

    def __iter__(self):
        return self

    def next(self):

        ok = self.read_next_xtc()
        if not ok:
            raise StopIteration
        return self
        


