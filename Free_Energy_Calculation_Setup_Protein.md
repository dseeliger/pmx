## Howto set up a free energy calculation for point mutations in a protein ##

This tutorial is also available in the Download section. (TI\_examples.tgz)

some general things:
  * the whole mutation stuff requires an extended forcefield which is packed in ffamber99sb\_mut.tgz.
  * whenever you start doing free energy calculations, copy this tgz to you working directory and unpack it
  * the setup requires gromacs version 4.0.x, it is not ported to 4.5 yet ( you can do the mdrun with 4.5, but not the setup)

Ok, let's go:

  * the input pdb is barnase.pdb, prepared from pdb 1bni
  * the individual steps are the same as in the DNA example, but there is a script
> > that can be used as a shortcut to set up simulations and also the free energy runs

Start with the setup of a plain MD run
```
python ~/software/pmx/scripts/mdsetup_407.py -f barnase.pdb -min_mdp mdp/em.mdp
```

This does the usual MD setup and yields md\_in.gro, which contains the solvated protein

For a crooks runs you can use the same script. And you need a file that defines the mutations you want to make. It looks like:

```
2|V|_ 9|G|_ # two mutations in one rune
2|V|_       # single mutation
9|G|_       # single mutation
6|A|_       # D->A mutation, charged mutation
17|A|_      # K->A mutation, charged mutation
```

The command to set up these runs is:

```
python ~/software/pmx/scripts/mdsetup_407.py -f md_in.gro -min_mdp mdp/em.mdp 
-skip_md_setup -fe.crooks -crooks_run_time 50 -crooks_mdp 
mdp/crooks_equilibration_stateA.mdp  -n_crooks_runs 1
```

This generates a directory for each mutation, prepares the topologies and the input files.
For charged mutations I keep the box neutral by scaling charges of 5 randomly picked Na+ ions.
It generated one run per state ( -n\_crooks\_runs 1), of 50ns length ( -crooks\_run\_time ) when runs are finished you can use a script to set up the non-eq runs

```
python ~/software/pmx/scripts/prepare_crooks_runs.py -d <run_dir> -mdp 
mdp/crooks_non_eq.mdp -sw_time 50 
```


submit all the stuff
```
foreach S( runA.0/morphes/frame* runB.0/morphes/frame* )
cd $S
qsub ....
cd ../../../
end
```

and finally analyze the dgdl.xvg files ( or dhdl.xvg of you used gmx 4.5 )

```
python ~/software/pmx/scripts/analyze_crooks.py -pa runA/morphes/frame*/dgdl.xvg
 -pb runB/morphes/frame*/dgdl.xvg
```

<img src='http://pmx.googlecode.com/files/W_over_t.png' alt='crooks analysis' />

writes some plots and results.dat with dG and error estimation

good luck! ;)