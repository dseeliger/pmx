# molecule.Molecule() #

## class molecule.Molecule(Atomselection) ##

Some attributes:

  * molecule.id            ( molecule or residue id)
  * molecule.resname  ( name )
  * molecule.atoms      ( list of atoms )


Some methods:

see also [Atomselection](pmx_atomselection.md) methods

There are a couple of methods that are protein specific

  * molecule.get\_phi( degree = True)        ( return phi angle )
  * molecule.get\_psi()
  * molecule.get\_omega()
  * molecule.get\_chi( 1 )
  * molecule.set\_phi( new\_phi, propagate = True )    ( set new phi angle and propagate change through the chain )
  * molecule.set\_psi(new\_psi,...)
  * molecule.set\_omega( new\_omega )
  * molecule.set\_chi(1, new\_chi )                             ( new value for chi 1 )


Other methods:

  * molecule.set\_resname( new\_resname )        ( new residue name, atom residue names are renamed, too )
  * molecule.insert\_atom( position, Atom)          ( insert a new atom at a certain position)
  * molecule.append( Atom )                             ( append new atom at the end of the atom list )
  * molecule.remove(Atom)                               ( remove this atom )
  * del molecule['CA']                                       ( remove atom with name CA )
  * molecule.write( "out.pdb")                            ( write molecule as pdb )
  * molecule.copy()                                          ( make a hard copy of this molecule)
  * molecule.new\_aa( "A")                                ( create new alanine residue )


Construction options:
```
>>> mol = Molecule()                                       # empty molecule
>>> mol = Molecule( atoms = list_of_atoms )     # create molecule from a list of atoms
>>> mol = Molecule().new_aa( "SER")               # create serine
```