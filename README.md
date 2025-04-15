# cocomo
The COCOMO force field defines a single-bead-per-residue coarse grained model for molecular dynamics of proteins and RNA. This helps accelerate simulations of processes such as biocondensate formation via liquid-liquid phase separation. Please read the papers in the reference for details. This repository contains codes to run MD simulations using the COCOMO/COCOMO2 force field developed by Feig and coworkers. Although the COCOMO2 model supports simulations of multidomain proteins, this functionality is not yet included in this repository. The current version supports IDR/RNA simulations only. Functionality of all codes are available via ```code_name.py -h```. I **do not** own rights to this force field. I have developed some codes that helps the user run simulations using the COCOMO force field in a streamlined fashion.

### References
1) [J. Chem. Theory Comput. 2023, 19, 2, 669–678.](https://doi.org/10.1021/acs.jctc.2c00856)
2) [J. Chem. Theory Comput. 2025, 21, 4, 2095–2107.](https://doi.org/10.1021/acs.jctc.4c01460)

The following codes are included in this repository:

## Table of contents
- [seq2len.py](#seq2len.py)
- [gen_slab.py](#gen_slab.py)
- [genpdb_cocomo.py](#genpdb_cocomo,py)
- [cocomo.py/cocomo2.py](#cocomo.py/cocomo2.py)
- [dcd2pdbFrame.py](#dcd2pdbFrame.py)
- [simulation.py](#simulation.py)
- [gen_config.py](#gen_config.py)
- [gyrate.py](#gyrate.py)
- [contacts.py](#contacts.py)
- [center.py](#center.py)

### ```seq2len.py```
This script determines the length of a protein/RNA chain if all beads were arranged linearly, based on the sigma (diameter) values listed in the force field parameters. This helps in approximating the minimum box dimensions if all chains were initially generated as linear. It requires the sequence string (1-letter names) as an input.
```
usage: seq2len.py [-h] -s

Estimate minimum box length from a protein sequence

options:
  -h, --help        show this help message and exit
  -s , --sequence   Protein sequence (str)
```
### ```gen_slab.py```
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
### ```genpdb_cocomo.py```
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
### ```cocomo.py```/```cocomo2.py```
These are the main force field and simulation scripts for COCOMO and COCOMO2.
```
usage: cocomo.py [-h] [-prm]

Run COCOCMO MD simulation

options:
  -h, --help            show this help message and exit
  -prm , --parameters   Input parameter file (str) [default: params.dat]
```
### ```dcd2pdbFrame.py```
This scripts generates PDB file from a a DCD trajectory. The time frame can be given as input. By default it generates the PDB for the last frame.
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
### ```simulation.py```
This scripts streamliunes the simulation process. It takes a configuration (YAML format) as input.
```
usage: simulation.py [-h] [-c]

    COCOMO Simulation Pipeline
    (A sample configuration file can be generated using the gen_config.py script)
    
options:
  -h, --help      show this help message and exit
  -c , --config   Path to YAML configuration file (default: config.yaml)
```
### ```gen_config.py```
Generates a YAML configuration file for the simulation. Be default it generates the default configurations. These can be tuned using user-input flags.
```
usage: gen_config.py [-h] [-s] [-c] [-r] [-f] [-T] [-p] [-t] [-F] [-em] [-en] [-eos] [-eod] [-es] [-ep] [-epi] [-eax] [-eay] [-eaz] [-egg] [-echk] [-pm] [-pn] [-pos] [-pod] [-ps] [-pp] [-ppi] [-pax] [-pay] [-paz] [-pg] [-pchk] [-b] [-g] [-d] [-ca] [-sa] [-cd] [-cp] [-tp] [-ip]
                     [-pf] [-ed] [-el] [-epdb] [-pd] [-pl] [-ppdb] [-epar] [-ppar] [-l] [-o]

Generate a default configuration file for the COCOMO simulation pipeline.

options:
  -h, --help            show this help message and exit
  -s , --seq_file       Path to the sequences file (default: "sequences.dat")
  -c , --codes_path     Directory containing simulation modules (default: "/home/tb/samukher/DATA/.bin/COCOCMO/MODULES")
  -r , --resources      Compute resource (default: "CUDA")
  -f , --forcefield     Forcefield version, options: "cocomo", "cocomo2" (default: "cocomo2")
  -T , --temp           Simulation temperature in Kelvin (default: 300.0)
  -p , --press          Pressure in bar (default: 1.0)
  -t , --time_step      Time step in ps (default: 0.02)
  -F , --friction       Friction coefficient (default: 0.01)
  -em , --eq_mini       Minimization steps for equilibration (default: 5000000)
  -en , --eq_nstep      Equilibration simulation steps (default: 500000)
  -eos , --eq_nstout    Log output interval for equilibration (default: 5000)
  -eod , --eq_nstdcd    Trajectory output interval for equilibration (default: 5000)
  -es, --eq_switching   Enable switching for equilibration (default: false)
  -ep, --eq_press_couple
                        Enable pressure coupling for equilibration (default: false)
  -epi , --eq_p_iso     Enable isotropic pressure for equilibration (default: true)
  -eax , --eq_p_aniso_x 
                        Enable anisotropic pressure along x for equilibration (default: false)
  -eay , --eq_p_aniso_y 
                        Enable anisotropic pressure along y for equilibration (default: false)
  -eaz , --eq_p_aniso_z 
                        Enable anisotropic pressure along z for equilibration (default: false)
  -egg , --eq_genpsf    Generate PSF file during equilibration (default: true)
  -echk , --eq_checkpoint 
                        Intervals for saving checkpoint files during equilibration (default: 100000)
  -pm , --prod_mini     Minimization steps for production (default: 5000000)
  -pn , --prod_nstep    Production simulation steps (default: 5000000)
  -pos , --prod_nstout 
                        Log output interval for production (default: 5000)
  -pod , --prod_nstdcd 
                        Trajectory output interval for production (default: 5000)
  -ps, --prod_switching
                        Enable production switching (default: false)
  -pp, --prod_press_couple
                        Enable production pressure coupling (default: false)
  -ppi , --prod_p_iso   Enable isotropic pressure for production (default: true)
  -pax , --prod_p_aniso_x 
                        Enable anisotropic pressure along x for production (default: false)
  -pay , --prod_p_aniso_y 
                        Enable anisotropic pressure along y for production (default: false)
  -paz , --prod_p_aniso_z 
                        Enable anisotropic pressure along z for production (default: false)
  -pg , --prod_genpsf   Generate PSF file during production (default: false)
  -pchk , --prod_checkpoint 
                        Intervals for saving checkpoint files during production (default: 100000)
  -b , --box_type       Box type: "cube" or "slab" (default: "cube")
  -g , --geometry       Initial chain geometry: "coil" or "straight" (default: "coil")
  -d , --box_dim        Custom box dimensions in nm as comma-separated values (e.g., "15.0,15.0,15.0"). Optional; if omitted, dimensions are computed based on box_type.
  -ca , --cube_add      Extra nm for cube box dimensions (default: 0)
  -sa , --slab_add      Extra nm added to Z dimension for slab box (default: 1000)
  -cd , --chain_dist    Chain separation in nm (default: 0.5)
  -cp , --chain_pos     Chain position (default: "center")
  -tp , --temp_pdb      Temporary PDB filename (default: "temp.pdb")
  -ip , --initial_pdb   Simulation-ready initial PDB filename (default: "initial.pdb")
  -pf , --psf_file      PSF filename (default: "topol.psf")
  -ed , --eq_dcd        Equilibration trajectory filename (default: "eqm.dcd")
  -el , --eq_log        Equilibration log filename (default: "eqm.log")
  -epdb , --eq_pdb      Equilibration PDB filename (default: "eqm.pdb")
  -pd , --prod_dcd      Production trajectory filename (default: "prod.dcd")
  -pl , --prod_log      Production log filename (default: "prod.log")
  -ppdb , --prod_pdb    Production PDB filename (default: "prod.pdb")
  -epar , --eq_param    Equilibration parameter file (default: "eqm.dat")
  -ppar , --prod_param 
                        Production parameter file (default: "prod.dat")
  -l , --log            Main simulation log file (default: "simulation.log")
  -o , --output         Output YAML file name (default: "config.yaml")
```
### ```gyrate.py```
Calculate time-trace of radius of gyration of the system.
```
usage: gyrate.py [-h] -f  -s  [-o] [-b] [-e] [-dt]

Calculate the mass-weighted total radius of gyration (Rg) using MDAnalysis.

options:
  -h, --help          show this help message and exit
  -f , --trajectory   Trajectory file (e.g., DCD) (default: None)
  -s , --topology     Topology file (PDB or PSF) containing mass information. (default: None)
  -o , --output       Output file for Rg vs time (default: gyrate.dat) (default: gyrate.dat)
  -b , --begin        Begin time in ps (default: 0.0)
  -e , --end          End time in ps. If not specified, process until the end of trajectory. (default: None)
  -dt , --skip        Process every nth frame (default: 1) (default: 1)
```
### ```contacts.py```
Calculate number of intra- and inter-chain contacts in the simulation trajectory (normalized with respect to the maximum possible contacts).
```
usage: contacts.py [-h] -f  -s  [-b] [-e] [-dt] [-o [...]] [-c]

Residue contacts with MDTraj.

options:
  -h, --help            show this help message and exit
  -f , --trajectory     Trajectory file (e.g., DCD) (default: None)
  -s , --topology       Topology file (PDB/PSF) (default: None)
  -b , --begin          Begin time (ps) (default: 0.0)
  -e , --end            End time (ps) (default: None)
  -dt , --skip          Process every nth frame (default: 1)
  -o [ ...], --output [ ...]
                        Chain pair specifiers (e.g., AA, AB) (default: None)
  -c , --cutoff         Cutoff in Å (default: 5.0)
```
### ```center.py```
Recenter selected atoms in the trajectory.
```
usage: center.py [-h] [-f] -s  -sel  -o  [-b] [-e] [-sk]

Center a molecular system by shifting the selected atoms' center to the box center.

options:
  -h, --help            show this help message and exit
  -f , --trajectory     Trajectory file (e.g., dcd or xtc). Omit for a static structure
  -s , --topology       Topology file (e.g., psf or tpr)
  -sel , --atom-selection 
                        MDAnalysis atom selection string
  -o , --output         Output file name. Format will match the input trajectory if provided
  -b , --begin          Index of the starting frame for processing
  -e , --end            Index of the ending frame for processing (-1 for last frame)
  -sk , --skip          Skip interval between frames
```
All the programs can be used as stand-alone codes or may be imported as modules in external codes.

## Contact
saumyak.mukherjee@biophys.mpg.de
