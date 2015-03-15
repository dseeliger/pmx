Functions to deal with Dunbrack rotamer libraries

```
def load_bbdep()             # load library from database

def get_rotamers( bbdep, resname, phi, psi ):
>>> bbdep = load_bbdep()                                    # load rotamer library
>>> mol = model.fetch_residues("VAL")[0]            # get first valine in model
>>> phi = mol.get_phi( degree = True )                 # get phi
>>> psi = mol.get_psi( degree = True )                  # get psi
>>> rot = get_rotamers( bbdep, "VAL", phi, psi )     # get rotamers for valine for a given phi, psi pair
>>> mol.set_conformation( rot[0] )                        # set sidechain dihedrals to the most likely rotamer conformation
```