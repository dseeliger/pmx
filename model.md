## class model.Model( Atomselection ) ##


Some attributes:

  * model.title        ( name )
  * model.atoms    ( list of atoms)
  * model.residues ( list of residues )
  * model.chains    ( list of chains )
  * model.chdic      ( dictionary to retrieve chains by id )
  * model.box        ( unit cell or simulation box )

Some methods:

see also [Atomselection](pmx_atomselection.md) methods

  * model.fetch\_residues(["SER","TRP"])       ( return list of all serine and tryptophane residues )
  * del model['A']                      ( remove chain A )
  * model.remove\_residue( molecule )     ( remove a residue )
  * model.remove\_atom( atom )              ( remove atom )
  * model.make\_chains()                       ( make chain list )
  * model.make\_residues()                    ( make residue list )
  * model.write("out.pdb")                      ( write structure file )

Construction methods:
```
>>> m = Model()                   # empty model
>>> m = Model("protein.pdb")  # create model from pdb or gro file
>>> m = Model( atoms = list_of_atoms )     # create model from list of atoms
>>> m = Model( residues = list_of_molecules )     # create model from list of residues
>>> m = Model( chains = list_of_chains )     # create model from list of chains
>>> m = Model( pdbline = line )                 # from pdbfile as single string ( m = Model( open('in.pdb').read() ) )
```