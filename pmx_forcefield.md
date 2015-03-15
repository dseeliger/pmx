Classes and functions to handle topologies and itp file. Only tested for amber forcefields

## class forcefield.Topology(): ##

### Constructors: ###

```
>>> top = Topology("topol.top")
>>> for atom in top.atoms: print atom.q, atom.atomtype         # print charge and atomtype
>>> top.write_top("new.top")                 # output new topology
```