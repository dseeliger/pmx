This module contains classes and functions for geometric operations.

## class geometry.Rotation(): ##

### Constructor: ###

```
>>> r = Rotation(vec1, vec2)       # creates rotation matrices for rotating vectors around the vector vec2-vec1
```

Example:

```
>>> r = Rotation( [0,0,0], [1,0,0])  # rotate around x
>>> for atom in model.atoms: atom.x = r.apply( atom.x, 20)    # rotate all atoms by 20 degrees
```

### Functions: ###

```
def fit( model1, model2, atom_names = [] ):  # fit model2 onto model1. use optionally only subset of atoms for fitting
```

```
Example:
>>> fit( model1, model2, atom_names = ['CA','N','C'])      # backbone fit

def fit_atoms( atoms1, atoms2, rot_atoms2 ):  # atoms1 and atoms2 are lists of atoms,
                                                                   # the transformation matrix is applied to rot_atoms2

def bb_super( mol1, mol2)               # superpose the backbone of mol2 onto mol1
```