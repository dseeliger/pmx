## class atomselection.Atomselection() ##

This class is a base class for other classes that contain lists of atoms ([Molecule](pmx_molecule.md), [Chain](pmx_chain.md), [Model](pmx_model.md)).
These classes inherit a bunch of methods from the Atomselection class.

# Inherited methods: #

  * atomselection.write("out.pdb")      ( write atoms in pdb format )
  * atomselection.com()                    ( put atoms to center of mass )
  * atomselection.a2nm()                  ( convert from Angstroem to nm )
  * atomselection.nm2a()                  ( convert from nm to Angstroem )
  * atomselection.get\_symbol()         ( determine element type for each atom )
  * atomselection.make\_long\_name()  ( make 4 charactor atom name )
  * atomselection.coords()                  ( return coordinates of all atoms as matrix )
  * atomselection.fetch\_atoms(["CA","C","N"])      ( return atoms with these names )
  * atomselection.translate( [1,2,3] )       ( move all atoms by this vector )
  * atomselection.rotate( axis, degree )   ( rotate all atoms around axis by this angle )
  * atomselection.random\_rotation()       ( rotate all atoms around a random vector )
  * atomselection.search\_neighbors()     ( do grid based neighbor search )
  * atomselection.get\_by\_id( list\_of\_atom\_ids )      ( return list of atoms )