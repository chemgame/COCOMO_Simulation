# cocomo
The COCOMO force field defines a single-bead-per-residue coarse grained model for molecular dynamics of proteins and RNA. This helps accelerate simulations of processes such as biocondensate formation via liquid-liquid phase separation. Please read the papers in the reference for details. This repository contains codes to run MD simulations using the COCOMO/COCOMO2 force field developed by Feig and coworkers. Although the COCOMO2 model supports simulations of multidomain proteins, this functionality is not yet included in this repository. The current version supports IDR/RNA simulations only.
Functionality of all codes are available via ```code_name.py -h```. I **do not** own rights to this force field. I have developed some codes that helps the user run simulations using the COCOMO force field in a streamlined fashion.

## Table of contents
- [References](##References)
- [seq2len.py](##seq2len.py)
- [gen_slab.py](##gen_slab.py)

## References
1) [J. Chem. Theory Comput. 2023, 19, 2, 669–678.](https://doi.org/10.1021/acs.jctc.2c00856)
2) [J. Chem. Theory Comput. 2025, 21, 4, 2095–2107.](https://doi.org/10.1021/acs.jctc.4c01460)

The following codes are included in this repository:

## seq2len.py
This script determines the length of a protein/RNA chain if all beads were arranged linearly, based on the sigma (diameter) values listed in the force field parameters. This helps in approximating the minimum box dimensions if all chains were initially generated as linear. It requires the sequence string (1-letter names) as an input.
```
usage: seq2len.py [-h] -s

Estimate minimum box length from a protein sequence

options:
  -h, --help        show this help message and exit
  -s , --sequence   Protein sequence (str)
```
## ```gen_slab.py```
This code generates the intitial simulation box. It takes a sequence file as an input. This file must have three coulmns: PRO/RNA number-of-chains sequence. PRO/RNA tells the code the type of biomolecule to consider. Initially chains can be placed near the box center or randomly inthe box. The initial chain configuration can be coiled or straight.
```
usage: gen_slab.py [-h] -s  -o  -box    [-d] [-g] [-p]

Generate a PDB file with protein/RNA chains arranged in a 3D simulation box

options:
  -h, --help            show this help message and exit
  -s , --seq_file       Path to the sequence file. Each line: [type PRO/RNA] [chain_count] [sequence].
  -o , --out_file       Output PDB file name.
  -box   , --box_size   
                        Box dimensions in nm (Lx Ly Lz).
  -d , --chain_distance 
                        Distance (nm) between 2 chains if exactly two. Default: 0.1 nm.
  -g , --geometry       Geometry for chain generation: 'straight' or 'coil'. Default: 'straight'.
  -p , --position       Positioning of chains in the box: 'center' or 'random'. Default: 'center'.
```
## ```genpdb_cocomo.py```
This code generates a pdb file that is compatible with the way the cocomo scripts deermine system topology.
```
Generate input PDB for COCOCMO simulation

options:
  -h, --help        show this help message and exit

I/O:
  -i , --infile     Input PDB file (str)
  -o , --outfile    Output PDB file (str) [default: cocomo_input.pdb]

Elastic network:
  -el , --elastic   Chains for elastic bonds (str) [default: " "]
```
## ```cocomo.py```/```cocomo2.py```
These are the main force field and simulation scripts for COCOMO and COCOMO2.
```
usage: cocomo.py [-h] [-prm]

Run COCOCMO MD simulation

options:
  -h, --help            show this help message and exit
  -prm , --parameters   Input parameter file (str) [default: params.dat]
```
## ```dcd2pdbFrame.py```
```
usage: dcd2pdbFrame.py [-h] -i  -top  -o  [-dump]

Extract a time frame from a .dcd trajectory as a pdb file

options:
  -h, --help            show this help message and exit

Required arguments:
  -i , --infile         Input DCD file (str)
  -top , --intop        Input topology file (str)
  -o , --outfile        Output PDB file (str)
  -dump , --dumpframe   Time frame to extract in ps (float) [default: -1]
```
This scripts generates PDB file from a a DCD trajectory. The time frame can be given as input. By default it generates the PDB for the last frame.
## ```simulation.py```
This scripts streamliunes the simulation process. It takes a configuration (YAML format) as input.
## ```gen_config.py```
Generates a YAML configuration file for the simulation. Be default it generates the default configurations. These can be tuned using user-input flags.
## ```gyrate.py```
Calculate time-trace of radius of gyration of the system.
## ```contacts.py```
Calculate number of intra- and inter-chain contacts in the simulation trajectory (normalized with respect to the maximum possible contacts).
## ```center.py```
Recenter selected atoms in the trajectory.

All the programs can be used as stand-alone codes or may be used as modules imported in external codes.

## Usage
```python simulation.py -c config.yaml```

## Contact
saumyak.mukherjee@biophys.mpg.de
