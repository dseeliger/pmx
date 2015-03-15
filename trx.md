MD simulations usually provide a huge amount of data which makes analysis quite time consuming and, hence, Python-based scripts might indeed be too slow for many cases.
However, if the amount of data is manageable you can write quite fancy tools with Python, particularly in combination with powerful plotting libraries like matplotlib.

So what I usually do is not to take trajectories directly written by mdrun but preprocessed trajectories with only the atoms of interest inside( usually protein+ligands ) and only a limited number of frames.
It is also important to note that coordinates are simply read from the trajectory, which means that no periodic boundary corrections or whatever are applied. So do these things before with trjconv.


Ok, reading trajectories is qickly explained:
Since xtc files only store coordinates, you usually want to update coordinates in a model instance.

```
>>> m = Model("protein.pdb")
```
Let's define a simple oberservable

```
>>> atom1 = m.residues[0]['CA']    # CA-atom of first residues
>>> atom2 = m.residues[-1]['CA']  # CA-atom of last residue
>>> ee_distance = lambda a,b: a-b # function that returns the distance between atom1 
                                  # and atom2
>>> fp = open("analysis.dat","w")  # open output file
```

Then we initalize a Trajectory instance with an xtc file.

```
>>> from pmx.xtc import *
>>> trj = Trajectory("traj.xtc")          # open xtc file and read first frame
>>> for frame in trj:                        # go over each frame
>>>     print frame                          # print some info
>>>     trj.update( m )    # update coords and box in model 
>>>     d = ee_distance( atom1, atom2 )     # calculate observable
>>>     print >>fp, "%8.3f %8.3f" % (frame.time.value, d ) # store data
>>> fp.close()       # close output file
>>> trj.close_xtc() # close xtc
```

Finally we can make a nice plot

```
>>> from pylab import *
>>> from pmx.parser import *
>>> data = read_and_format("analysis.dat","ff")
>>> time = map( lambda a: a[0], data)
>>> dist = map(lambda a:a[1], data)
>>> plot( time, dist, "r-", lw = 2)
>>> xlabel( "time [ps]")
>>> ylabel("end to end distance [A]")
>>> savefig("plot.png")
```

<img src='http://pmx.googlecode.com/files/plot.png' alt='plot' />