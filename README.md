# Data for simulations of bothropstoxin-I (BthTx-I) and mutants in water and membrane

This repository presents the parameters, input files and raw data corresponding to the study:
    Conformational plasticity and membrane binding govern functional dimer states of lys49-pla₂
    DOI: **************

-Input files for MD simulations of the systems (mdp files) using the GROMACS software;

-Topology files;

-Final structures for all simulations;

-Raw data of the analysis presented in the paper;

-Scripts used for analysis

```bash
	Summary of content

    |
    |   
    └── protein_in_water

    	"Atomistic simulations of wild-type BthTx-I and double the double mutants W77H/L10W and W77H/V31W. Three independet replica simulations were performed for each protein (repl1, repl2, repl3)"
	
        │
        ├── mdp
        │   ├── 1. minim.mdp
        │   ├── 1. nvt.mdp
        │   ├── 2. npt-1.mdp
        │   ├── 3. npt-2.mdp
        │   ├── 4. npt-3.mdp
        │   ├── 5. npt-4.mdp    
        │   └── 6. npt-5.mdp
        │
        ├── compact_model_dimer
        │   │
        │   ├── wt
        │   │   ├── input
        │   │   │   ├── 1. em_solv.gro
        │   │   │   ├── 2. posre_Protein_chain_A.itp
        │   │   │   ├── 3. posre_Protein_chain_B.itp
        │   │   │   ├── 4. topol_Protein_chain_A.itp
        │   │   │   ├── 5. topol_Protein_chain_B.itp
        │   │   │   └── 6. topol.top
        │   │   │
        │   │   ├── repl1
        │   │   │   ├──	1. final.gro
        │   │   │   ├──	2. fret.dat
        │   │   │   ├──	3. numcont.dat
        │   │   │   ├──	4. Rg.dat
        │   │   │   └──	5. rmsd.dat
        │   │   │
        │   │   ├── repl2
        │   │   └── repl3
        │   │   
        │   ├── w77h-l10w
        │   └── w77h-v31w
        │
        └── extended_model_dimer
```




