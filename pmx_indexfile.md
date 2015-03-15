Index files are frequently needed to select subsets of a simulation system.
In pmx you can read, change and write index files quite easy.
To use this functionality you need the pmx.ndx module

```
>>> from pmx.ndx import *
```
Reading and writing index files

Reading is done by

```
>>> ndx_file = IndexFile( "index.ndx")
```

to check what you read use
```
>>> print ndx_file
```
Accessing an index group works in dictionary style
```
>>> bb = ndx_file['Backbone']
>>> print bb
```

If you have a structure file corresponding to your index file you can
select atom subsets from the index file groups

```
>>> m = Model("gmx.pdb")
>>> bb_atoms = ndx_file['Backbone'].select_atoms( m )  # list with backbone atoms
>>> for atom in bb_atoms:            
>>>     print atom
```

Writing index files is done by

```
>>> ndx_file.write("new_index.ndx")
```

Constructing index groups

Index groups can be constructed with lists atom integers ( atom ids ) or lists of atoms

```
>>> ndx_group = IndexGroup( "backbone", atoms = bb_atoms ) # with list of atoms
>>> atom_ids = map( lambda atom: atom.id, m.atoms[:100] ) # atom ids of first 100 atoms
>>> ndx_group = IndexGroup( "first_100", ids = atom_ids )  # construct with ids
```

Adding and deleting index groups

```
>>> ndx_file.add_group( ndx_group )               # adds the group to the index file
>>> del ndx_file['Backbone']                         # deletes group Backbone
```