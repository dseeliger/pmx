import sys, os
from pmx import *
from pmx.ndx import *

m = Model( "gmx.pdb" )

ndx = IndexFile("index.ndx")

print ndx
print ndx['Backbone']

atoms = ndx['Backbone'].select_atoms( m )
del ndx['Backbone']
grp = IndexGroup( "Backbone", atoms = atoms )

ndx.add_group( grp )



print ndx
#for atom in atoms:
#    print atom

