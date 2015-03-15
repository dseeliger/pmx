The ndx module contains classes to handle Gromacs index files. It's pretty straightforward. An index group is handled by the class "IndexGroup" and an entire index file is handled by the class "IndexFile" which contains a list of IndexGroups.

## class IndexGroup: ##

### attributes: ###

  * IndexGroup.name         ( name of the index group )
  * IndexGroup.ids             ( list of atom ids )

### methods: ###

print IndexGroup         ( writes group in Gromacs .ndx format

### Construction methods: ###

```
>>> grp = IndexGroup( ids = list_of_atom_ids, name = "backbone")         # from a list of atom ids
>>> grp = IndexGroup().fromFile( "[ CA ]", open("index.ndx").readlines() )  # from file
```


## class IndexFile: ##

### attributes: ###

  * IndexFile.groups         ( list of index groups )
  * IndexFile.names             ( names of index groups )
  * IndexFile.dic                ( access index groups by name )


### methods: ###

  * IndexFile.write( "index.ndx" )

### Construction methods: ###

```
>>> ndx = IndexFile( groups = list_of_index_groups)         # from a list of index groups
>>> ndx = IndexFile("index.ndx")                    # from file
```