from numpy import *
import os, re
from DataObj import DataObj
from SelObj import SelObj, Model, Chain, Residue, Atom #maybe not load those!
from options import Commandline, Option, FileOption
import library as lib
import geometry as geo


M = Model('KcsA_System.pdb')
M.load_Trajectory('md.part0001.xtc', stop=10)

W = M.water

W0 = W.T[0].N[:12]
W9 = W.T[9].N[:12]

P0 = M.protein.T[0]
P1 = M.protein.T[3:6]
P2 = M.protein.T[6:9]

#from NGMX_geometry import linAssignment

W0.fetch('OW').write("W0.pdb")
W9.fetch('OW').write("WU.pdb")

W9.perm_reduction(W0)
W9.fetch('OW').write("WC.pdb")