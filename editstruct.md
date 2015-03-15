# Reading and writing structure files #

pmx reads structure files into the "Model" class.
This is simply done by
```
>>> m = Model("protein.pdb")
```
or if you want to load a gro file
```
>>> m = Model("protein.gro")
```

an imporant difference is that in gro files coordinates are in nanometer
whereas in pdb files we have Angstroems.
To check the coordinate set of your model use
```
>>> print m.unity
```

which is either 'A' or 'nm'

Conversion is simply done by
```
>>> m.nm2a()
```
or
```
>>> m.a2nm()
```

For writing a structure file you don't have to care about this, it is checked.
```
>>> m.write("out.gro")
```
will write nanometers, whereas
```
>>> m.write("out.pdb")
```
will write Angstroems


## Iterating over atoms, residues and chains ##

A structure in the "Model" class is organized in a Model->Chain->Molecule->Atom hirarchy as indicated below


The Model contains list of atoms, residues and chains which allows to navigate on whatever level is appropriate.
Let's first iterate over all atoms in a model

```
>>> for atom in m.atoms:
>>>    print atom.id, atom.name, atom.resname
>>>    if atom.name == 'CA':
>>>        atom.name = 'CAX'               # rename some atoms
```

Same works for residues:

```
>>> for resi in m.residues:
>>>     if resi.resname in ['ALA','SER']:     # select some residues
>>>        print resi
>>>        for atom in resi.atoms:              # loop over atoms in these residues
>>>            print atom.bfac                     # print some properties
```

and also for chains

```
>>> for c in m.chains:
>>>     print c.id, len(c.residues), len(c.atoms)        # print chain id, number of
                                                         # residues and number of atoms
```

chains can also be accessed by their identifier

```
>>> chainB = m.chdic['B']
>>> for atom in chainB.atoms:
>>>     print atom.occ
```

From each level you can also reach the next upper level. Hence, an atom knows its residue, chain and model.

```
>>> atom = m.atoms[0]          # first atom in model
>>> print atom.molecule         # the residue of this atom
>>> print atom.chain              # the chain of this residue
>>> print atom.model             # the model
```

### Selecting subsets of a structure ###

The model.chdic usage shown above is already an example of how to select a subset of a structure.
There are some other functions for this task, e.g. for selecting atoms.

```
>>> subset = m.fetch_atoms("CA")       # returns all ca atoms
>>> subset = m.fetch_atoms(["CA","C","N"])  # returns all backbone atoms
>>> subset = m.fetch_atoms("C", how="byelem") # returns all carbon atoms
>>> subset = m.fetch_atoms("H", how = "byelem", inv = True)   # return all non-hydrogen atoms
```

This function also works on the chain and the molecule level

```
>>> resi = m.chdic['B'].residues[-1]            # last residue of chain B
>>> bb_atoms = resi.fetch_atoms(["N","CA","C"])  
```

There is a very similar function called fetchm
```
>>> bb_atoms = resi.fetchm(["N","CA","C"])  
```

The difference between both is that
```
>>> atoms = resi.fetch_atoms(["C","N"])
```

does not care about which atom comes first. Hence, since N usually comes before C the first atom returned will be N, the second C.
If you use the fetchm function
```
>>> atoms = fetchm(["C","N"])
```

atoms will be strictly returned in the order of the list entries.

On the residue level you have also a shortcut for selecting single atoms using
```
>>> CA = resi['CA']
```
Selecting residues works similarly (works of course only from the model or chain level )

```
resl = m.fetch_residues(["ALA","SER"])       # returns all alanine and serine residues
resl = m.fetch_residues(["ALA","SER"], inv = True)       # returns all non-alanine and non-serine residues

reslB = m.chdic['B'].fetch_residues('PRO')       # proline residues in chain B
```

### Adding and deleting ###

Deleting atoms, molecules and chains

Deleting atoms is straightforward
On each level that contains an atomlist (Molecule, Chain, Model) you can remove an atom by
```
>>> atom = m.atoms[0]        # first atom
>>> m.remove_atom(atom)
```

on the molecule level you can also use
```
>>> del resi['CA']       # which remove the CA atom
```

deleting residues is similar

```
>>> mol = m.residues[0]
>>> m.remove_residue( mol )
```

whereas deleting a chain works like
```
>>> m.remove_chain("A")   # to remove chain A
>>> del m['A']                   # does the same
```


In general it should be noted that deleting atoms or residues is not very fast since the entire internal structure is rebuilt each time you delete something. Hence, if you want to delete a large number of atoms it might be faster to construct a new model
with a subset of atoms of the first one. E. g., deleting all hydrogen atoms woul be faster in the following way

```
>>> non_h_atoms = m.fetch_atoms("H",how="byelem", inv = True) # select all non-hydrogen atoms
>>> new = Model( atoms = non_h_atoms )            # construct new model
>>> new.write("heavy_atoms_only.pdb")             # write pdb
```

Adding atoms, molecules and chains
Atoms can be added to molecules using the "insert\_atom" or "append" function.

```
>>> a = Atom( name = "X", x = [0,0,0] )      # construct atom
>>> resi = m.residues[-1]                # last residue
>>> resi.insert_atom(0, a )              # insert atom at position 0
>>> resi.append( atom )                 # add atom at the end of the atomlist
```

Molecules can be inserted into chains

```
>>> mol = Molecule().new_aa( "TRP")      # construct residue
>>> chB = m.chdic['B']
>>> chB.insert_residue( 5, mol)           # insert new residue at position 5
>>> chB.append( mol )                    # add residue at the end of the residue list
```

Molecules can also be inserted into models

```
>>> m.insert_residue( 5, mol, "B")          # insert residue at position 5 of chain B
```
Inserting chains is just the same

```
>>> c = Chain().create("AAA")             # new chain
>>> c.set_chain_id('D')
>>> m.insert_chain(3, c)                     # add chain at position 3
>>> m.append( c )                             # add chain at the end
```

you can also insert a chain into an existing chain

```
>>> chB = m.chdic['B']
>>> c = Chain().create("AAA")             # new chain
>>> chB.insert_chain( 6, c )                # insert chain
```

In this case the new chain is not independent anymore but merged into chB.

IMPORTANT NOTE 1:
THE INSERT AND APPEND FUNCTIONS DO NOT CHANGE THE COORDINATES OF ANY ATOM.
HENCE, NEW ATOMS, RESIDUES AND CHAINS ARE JUST PUT AT THE CORRECT POSITION IN
THE STRUCTURE TREE BUT ARE NOT MODELLED INTO THE STRUCTURE WHICH IS SOME
ORDERS OF MAGNITUDE MORE COMPLEX.
IF YOU WANT EXTEND CHAINS USE THE FUNCTIONS OUTLINED IN THIS SECTION.

IMPORTANT NOTE 2:
IT IS STRONGLY RECOMMENDED TO USE THESE FUNCTIONS TO ADD OR DELETE ATOMS,
MOLECULES OR CHAINS, SINCE CALLING THIS FUNCTIONS TRIGGERS AN UPDATE
OF THE DATA STRUCTURE. IF YOU JUST USE THE DELETE FUNCTIONALITY OF A LIST
YOU WILL MESS THIS UP COMPLETELY

Example:

```
>>> c = Chain().create("ALT")    
>>> print c.atoms[0].name
"N"
>>> del c.atoms[0]                    # delete list item, first atom is gone
>>> print c.atoms[0].name                  # first atom is a different one now
"H"
>>> atom = c.residues[0].atoms[0]     # select first atom from the first residue
>>> print atom.name                        # should be the same as above, but it is not
"N"
```

### Useful shortcuts ###

It is sometimes necessary to select atoms based on their coordinates, e.g. choosing atoms above a membrane or something like that. To do this you can use the Python "filter" function.

```
>>> atoms = filter( lambda a: a.x[2] > 10 and a.x[2] < 20, m.atoms)
```
this will return a list of all atoms with a z-coordinate between 10 and 20.