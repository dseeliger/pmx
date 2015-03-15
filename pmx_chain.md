## class chain.Chain(Atomselection) ##


Some arttributes:

  * chain.id                ( chain identifier )
  * chain.atoms          ( list of atoms )
  * chain.residues       ( list of residues )



Protein specific methods:

  * chain.sequence()               ( return peptide sequence as string )
  * chain.get\_sequence()         ( return sequence in fasta format )
  * chain.create( "AFGCFRT")  ( create chain instance from sequence )
  * chain.nbuild("R")                ( attach arginine at the n-terminus )
  * chain.cbuild("R")                ( attach arginine at the c-terminus )
  * chain.fuse( other\_chain )    ( attach a second chain at the c-terminus )
  * chain.add\_nterm\_cap()       ( attach an ACE molecule at the n-terminus )
  * chain.add\_cterm\_cap()       ( attach NME molecule at the c-terminus )



Other methods:


  * chain.fetch\_residues(["ALA","GLY","TRP"])        ( return all ALA, GLY and TRP residues)
  * del chain[5](5.md)                        ( remove residue 5 from the chain )
  * chain.insert\_residue(pos, mol )  ( insert a new residue into the chain )
  * chain.renumber\_residues()
  * chain.set\_chain\_id("Z")          ( set new chain id, change is propagated to all atoms in the chain )
  * chain.copy()                         ( hard copy of the chain )


Construction options:
```
>>> c = Chain()                                      # empty chain
>>> c = Chain().create("AFGTRL")           # create peptide chain
>>> c = Chain( atoms = list_of_atoms )    # create a chain from a list of atoms
>>> c = Chain( residues = list_of_molecules )   # create a chain from a list of molecules
```