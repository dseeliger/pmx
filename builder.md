This section describes how to build peptides and dna strands. There are several possibilities to do this.

## Proteins ##
### Constructing amino acids ###
Let's start by building a single aminoacid.
```
>>> from pmx import *
>>> mol = Molecule().new_aa( "ALA")          # build alanine
```

if we wish to build a peptid chain starting with this residue we first have to initialize the peptide chain
```
>>> c = Chain( residues = [mol] )               # build peptide chain
```
now we can extend this chain either at the n-terminus or at the c-terminus
```
>>> c.nbuild("T")                    # extend at n-terminus
>>> c.cbuild("R")                   # extend at c-terminus
```
if this is all we want to, we simply write the chain as pdb file
```
>>> c.write("TAR.pdb")
```

### Building peptide chains from sequence ###

usually we have a peptide sequence to be transformed to 3d coords. This can be done from the chain level
```
>>> sequence = "FRTLKNCWQ"
>>> c = Chain().create( sequence )
>>> c.write("extended_pep.pdb")
```

This builds a peptide chain in an extended conformation.

<img src='http://pmx.googlecode.com/files/extended_pep.png' alt='extended peptide' />



### Changing secondary structure ###

If we now want to change the peptide conformation
we can change the backbone dihedral angles
for each residue. E.g., the phi/psi angles for
an ideal alpha-helix are -57 and -47 degree.

```
>>> for resi in c.residues:
>>>      resi.set_phi( -57, True )
>>>      resi.set_psi( -47, True )
>>> c.write("helix_pep.pdb")
```

transforms the peptide into a helical conformation.
<img src='http://pmx.googlecode.com/files/helix_pep.png' alt='helical peptide' />

### Chain fusion ###

We can also attach a new series of peptides to
the helix peptides, either by using the nbuild() or
cbuild() functions or by creating a new chain
from scratch that should be fused with the
first one

So let's first build a new chain

```
>>> c2 = Chain().create("LMNFRTS")
```
and attach this peptide to the helix peptide

```
>>> c.fuse( c2 )
>>> c.write( "fused_pep.pdb")
```

Chain c2 is now merged into chain c.
Creating chains of repeating units

Changing secondary structure
If we would like to construct a chain of
repeating units, this can also be done
in a straightforward manner. We
simply construct the repeating unit
and extend by copies of itself.

```
>>> rep_seq = "HEATLLFT"
>>> c = Chain().create(rep_seq)
```

now we create a copy of this chain

```
>>> new_chain = c.copy()
```

and attach the new chain 5 times

```
>>> for i in range(5): 
>>>     c.fuse( new_chain.copy() )
>>> c.write("repeat_pep.pdb")
>>> print c.sequence()
HEATLLFTHEATLLFTHEATLLFTHEATLLFTHEATLLFTHEATLLFT
```

and finally set some residues in helical conformation

```
>>> helical_residues = c.fetch_residues(["HIS","GLU","ALA","THR","LEU"])
>>> for resi in helical_residues:
>>>    resi.set_phi( -57, True )
>>>    resi.set_psi( -47, True )
>>> c.write("repeat_helix.pdb")
```

<img src='http://pmx.googlecode.com/files/rep_helix.png' alt='helical peptide' />


## DNA ##

There is not much functionality regarding DNA up to now.
Just a simple DNA builder.

### Building DNA strands ###

This requires the pmx.builder module
```
>>> from pmx.builder import *
>>> model = build_dna_strand("ACGTGTCA")
>>> model.write("dna.pdb")
```

<img src='http://pmx.googlecode.com/files/dna_strand.png' alt='helical peptide' />