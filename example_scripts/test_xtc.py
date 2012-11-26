from pmx import *
from pmx.xtc import *
import sys, os

m = Model("out.pdb")


atom1 = m.residues[0]['CA']    # CA-atom of first residues 
atom2 = m.residues[-1]['CA']  # CA-atom of last residue 
distance_12 = lambda a,b: a-b   # function that returns the distance between atom1 and atom2 
fp = open("analysis.dat","w")  # open output file 


trj = Trajectory("concoord.xtc")          # open xtc file and read first frame
for frame in trj:                        # go over each frame 
    print frame                          # print some info 
    trj.update( m )    # update coords in model 
    d = distance_12( atom1, atom2 )     # calculate observable 
    print >>fp, "%8.3f %8.3f" % (frame.time.value, d ) # store data 
fp.close()       # close output file 
trj.close_xtc() # close xtc 
from pylab import *
from pmx.parser import *
data = read_and_format("analysis.dat","ff")
time = map( lambda a: a[0], data)
dist = map(lambda a:a[1], data)
plot( time, dist, "r-", lw = 2)
xlabel( "time [ps]")
ylabel(r"end to end distance [$\AA$]")
savefig("plot.png")

