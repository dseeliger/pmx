import sys, os
from pmx import *

from pmx.geometry import Rotation, Rotation2

import time

R = Rotation([0,0,0],[1,0,0])

m = Model("cdk2.pdb")

t1 = time.clock()

for r in m.residues:
    r.set_phi(-139, True)
    r.set_psi(135, True)
## for i in range(100):
##     for atom in m.atoms:
##         atom.x = R.apply( atom.x, 60*pi/180.)
    

t2 = time.clock()
print t2-t1
m.write("out2.pdb")

