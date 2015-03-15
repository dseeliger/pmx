## pmx.atom ##

### class atom.Atom() ###

This class stores all properties (attributes) related to an atom. E.g., when pmx reads a pdb file, each line starting with "ATOM" or "HETATM" is converted to an Atom instance.
Hence, and atom object constructed from a pdb line contains:

  * atom id (atom.id)
  * atom name (atom.name)
  * residue name (atom.resname)
  * residue id (atom.resnr)
  * chain id (atom.chain\_id)
  * coordinates (atom.x)
  * occupancy (atom.occ)
  * b-factor (atom.bfac)


In addition to storage properties, the atom class also has a number of methods:

  * atom.dist( atom2 ) ( distance between two atoms)
  * atom - atom2       ( same as above with overloaded "-" operator )
  * atom.dist2( atom2) (squared distance, faster than distance)
  * atom.angle( atom2, atom3, degree=True) angle calculation, returns rad by default )
  * atom.dihedral( atom2, atom2, atom4, degree=True) (dihedral calculation)
  * print atom   ( writes atom as pdb line)


for full documentation use help(Atom) from the interpreter.