# Data for simulations of bothropstoxin-I (BthTx-I) and mutants in water and membrane

This repository presents the parameters, input files and raw data corresponding to the study:

**Conformational plasticity and membrane binding govern functional dimer states of lys49-pla₂**

 -Input files for MD simulations of the systems (mdp files) using the GROMACS software;
 
 -Topology files;
 
 -Final structures for all simulations;
 
 -Raw data of the analysis presented in the paper;
 
 -Scripts used for analysis


```bash

Summary of content

protein_in_water
.
├── compact_model_dimer
│   ├── w77h-l10w
│   │   ├── input
│   │   │   ├── em_solv.gro
│   │   │   ├── posre_Protein_chain_A.itp
│   │   │   ├── posre_Protein_chain_B.itp
│   │   │   ├── topol_Protein_chain_A.itp
│   │   │   ├── topol_Protein_chain_B.itp
│   │   │   └── topol.top
│   │   ├── repl1
│   │   │   ├── final.gro
│   │   │   ├── fret.dat
│   │   │   ├── numcont.dat
│   │   │   ├── Rg.dat
│   │   │   └── rmsd.dat
│   │   ├── repl2
│   │   └── repl3
│   ├── w77h-v31w
│   └── wt
├── extended_model_dimer
└── mdp
    ├── minim.mdp
    ├── npt-1.mdp
    ├── npt-2.mdp
    ├── npt-3.mdp
    ├── npt-4.mdp
    ├── npt-5.mdp
    └── nvt.mdp

```
