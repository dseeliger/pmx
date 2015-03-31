# pmx
## Introduction ##

pmx (formerly pymacs) has started as a small bunch of classes to read structure files such as pdb or gro and trajectory data in gromacs xtc format. Over the years it has been extended
towards a versatile (bio-) molecular structure manipulation package with some additional functionalities, e.g. gromacs file parsers and scripts for setup and analysis of free energy calculations.

## Citations: ##

D. Seeliger and Bert L. de Groot, Biophys. J. 98(10):2309-2316 (2010)

V. Gapsys, S. Michielssens, D. Seeliger, B. L. de Groot. J. Comput. Chem. 2014, DOI: 10.1002/jcc.23804


## Purpose ##

I mostly use pmx to write short scripts which perform some changes in pdb files, e.g. changing atom or residue names, applying some geometric transformations or doing some kind of analysis. The critical issue for these things is usually not calculation time but straightforward selection of some atoms/residues of interest, quick file parsing and data visualization, which renders a well-organized data structure and easy programming style more important than computation performance. Hence, and ideal task for Python.

## Installation ##

Checkout the source code and run the usual python installation
```
git clone https://code.google.com/p/pmx/
cd pmx 
sudo python setup.py install
```

## Software Requirements ##

  * [numpy](http://numpy.scipy.org/)
  * [scipy](http://www.scipy.org/)
  * [matplotlib](http://matplotlib.org/) ( for analysis scripts )

## Getting Started ##

pmx stores structure data in Python classes. The "Model" class is the uppermost class
which contains severals lists of Atoms, Molecules and Chains.
The following script reads a pdb file, prints some atom properties and writes the structure in gro format.

<img src='http://pmx.googlecode.com/files/pmx_data.jpg' alt='pmx data structure' />

The figure above shows the most important data structure. A "Model" instance contains list of chains, residues and atoms. A "Chain" instance of residues and atoms and a "Molecule" instance of a list of atoms only. Check the example scripts for how to navigate trough particular storage classes.

## pmx Modules ##
pmx currently contains the following modules. Click on the links below for a (short) documentation

  * [atom](pmx_atom.md)
  * [molecule](pmx_molecule.md)
  * [chain](pmx_chain.md)
  * [model](model.md)
  * [atomselection](pmx_atomselection.md)
  * [options](pmx_options.md)
  * [library](pmx_library.md)
  * [ndx](pmx_ndx.md)
  * [geometry](pmx_geometry.md)
  * [parser](pmx_parser.md)
  * [xtc](pmx_xtc.md)
  * [builder](pmx_builder.md)
  * [forcefield](pmx_forcefield.md)
  * [rotamer](pmx_rotamer.md)
  * [ffparser](pmx_ffparser.md)


## pmx Classes (the most important ones) ##

  * [Atom](pmx_atom.md)
  * [Molecule](pmx_molecule.md)
  * [Chain](pmx_chain.md)
  * [Model](model.md)
  * [Trajectory](pmx_trajectory.md)
  * [Topology](pmx_topology.md)
  * [IndexFile](pmx_indexfile.md)

## Using pmx ##

  * [Editing structure files, selecting atoms, residues, chain, etc.](editstruct.md)
  * [The pmx commandline parser, documenting programs and scripts](pmx_options.md)
  * [Building peptides and DNA strands from scratch](builder.md)
  * [Writing MD analysis tools, access xtc data from Python](trx.md)
  * [Modifying Gromacs topologies](topol.md)
  * [Generating Gromacs index files](pmx_indexfile.md)
