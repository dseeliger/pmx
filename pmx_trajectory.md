This module contains the Trajectory class which ( so far ) can read xtc files

## class xtc.Trajectory() ##

### methods: ###

```
Trajectory.open_xtc( xtc_file )          ( open xtc file for reading )
Trajectory.read_first_xtc()                ( read first frame )
Trajectory.read_next_xtc()                ( read next frame )
Trajectory.update_atoms( model )                ( update coords in the model )
Trajectory.close_xtc()                ( close file )
```

### Constructors: ###
```
>>> trj = Trajectory( "traj.xtc")        # open xtc file for reading
```