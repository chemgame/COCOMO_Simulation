# cocomo
This repository contains codes to run MD simulations using the COCOMO force field developed by Feig and coworkers. 
I do not own rights to this force field. I have developed some codes that helps the user run simulations using the COCOMO force field in a streamlined fashion.

# References: 
1) J. Chem. Theory Comput. 2023, 19, 2, 669–678.
2) 2) J. Chem. Theory Comput. 2025, 21, 4, 2095–2107

N.B. Although the COCOMO2 model supports simulations of multidomain proteins, this functionality is not yet included in this repository. The current version supports IDR/RNA simulations only.
Functionality of all codes are available via code_name.py -h

# Description:
The COCOMO force field defines a single-bead-per-residue coarse grained model for molecular dynamics of proteins and RNA. This helps accelerate simulations of processes such as biocondensate formation via liquid-liquid phase separation. Please read the papers in the reference for details.
This repositiory contains the following codes:

1. seq2len.py: This script determines the length of a protein/RNA chain if all beads were arranged linearly, based on the sigma (diameter) values listed in the force field parameters. This helps in approximating the minimum box dimensions if all chains were initially generated as linear. It requires the sequence string (1-letter names) as an input.
2. gen_slab.py: This code generates the intitial simulation box. It takes a sequence file as an input. This file must have three coulmns: PRO/RNA number-of-chains sequence. PRO/RNA tells the code the type of biomolecule to consider. Initially chains can be placed near the box center or randomly inthe box. The initial chain configuration can be coiled or straight.
3. genpdb_cocomo.py: This code generates a pdb file that is compatible with the way the cocomo scripts deermine system topology. 
4. cocomo.py/cocomo2.py: These are the main force field and simulation scripts for COCOMO and COCOMO2. 
5. dcd2pdbFrame.py: This scripts generates PDB file from a a DCD trajectory. The time frame can be given as input. By default it generates the PDB for the last frame.
6. simulation.py: This scripts streamliunes the simulation process. It takes a configuration (YAML format) as input.
7. gen_config.py: Generates a YAML configuration file for the simulation. Be default it generates the default configurations. These can be tuned using user-input flags.
8. gyrate.py: Calculate time-trace of radius of gyration of the system.
9. contacts.py: Calculate number of intra- and inter-chain contacts in the simulation trajectory (normalized with respect to the maximum possible contacts).
10. center.py: Recenter selected atoms in the trajectory.

All the programs can be used as stand-alone codes or may be used as modules imported in external codes.

# Usage:
python simulation.py -c config.yaml

# Contact: 
saumyak.mukherjee@biophys.mpg.de
