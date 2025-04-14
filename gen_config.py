#!/usr/bin/env python
"""
Generate a configuration YAML file for the COCOMO simulation pipeline.
Each configuration option has a short and a long flag.
If not specified, default values are used.
The output YAML file name can be specified (default: config.yaml).
"""

import argparse
import yaml
from datetime import datetime
from pathlib import Path
import os

def main():
    parser = argparse.ArgumentParser(
        description="Generate a default configuration file for the COCOMO simulation pipeline."
    )
    
    # ---------------------------
    # Input files and codes
    # ---------------------------
    parser.add_argument("-s", "--seq_file", type=str, default="sequences.dat",
                        help='Path to the sequences file (default: "sequences.dat")', metavar="")
    parser.add_argument("-c", "--codes_path", type=str,
                        default="/home/tb/samukher/DATA/.bin/COCOCMO/MODULES",
                        help='Directory containing simulation modules (default: "/home/tb/samukher/DATA/.bin/COCOCMO/MODULES")',
                        metavar="")

    # ---------------------------
    # General simulation options
    # ---------------------------
    parser.add_argument("-r", "--resources", type=str, default="CUDA",
                        help='Compute resource (default: "CUDA")', metavar="")
    parser.add_argument("-f", "--forcefield", type=str, default="cocomo2",
                        help='Forcefield version, options: "cocomo", "cocomo2" (default: "cocomo2")', metavar="")
    parser.add_argument("-T", "--temp", type=float, default=300.0,
                        help='Simulation temperature in Kelvin (default: 300.0)', metavar="")
    parser.add_argument("-p", "--press", type=float, default=1.0,
                        help='Pressure in bar (default: 1.0)', metavar="")
    parser.add_argument("-t", "--time_step", type=float, default=0.02,
                        help='Time step in ps (default: 0.02)', metavar="")
    parser.add_argument("-F", "--friction", type=float, default=0.01,
                        help='Friction coefficient (default: 0.01)', metavar="")

    # ---------------------------
    # Equilibration simulation parameters
    # ---------------------------
    parser.add_argument("-em", "--eq_mini", type=int, default=5000000,
                        help='Minimization steps for equilibration (default: 5000000)', metavar="")
    parser.add_argument("-en", "--eq_nstep", type=int, default=500000,
                        help='Equilibration simulation steps (default: 500000)', metavar="")
    parser.add_argument("-eos", "--eq_nstout", type=int, default=5000,
                        help='Log output interval for equilibration (default: 5000)', metavar="")
    parser.add_argument("-eod", "--eq_nstdcd", type=int, default=5000,
                        help='Trajectory output interval for equilibration (default: 5000)', metavar="")
    parser.add_argument("-es", "--eq_switching", action="store_true",
                        help='Enable switching for equilibration (default: false)')
    parser.add_argument("-ep", "--eq_press_couple", action="store_true",
                        help='Enable pressure coupling for equilibration (default: false)')
    parser.add_argument("-epi", "--eq_p_iso", type=lambda s: s.lower() in ['true','1','yes'], default=True,
                        help='Enable isotropic pressure for equilibration (default: true)', metavar="")
    parser.add_argument("-eax", "--eq_p_aniso_x", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along x for equilibration (default: false)', metavar="")
    parser.add_argument("-eay", "--eq_p_aniso_y", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along y for equilibration (default: false)', metavar="")
    parser.add_argument("-eaz", "--eq_p_aniso_z", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along z for equilibration (default: false)', metavar="")
    parser.add_argument("-egg", "--eq_genpsf", type=lambda s: s.lower() in ['true','1','yes'], default=True,
                        help='Generate PSF file during equilibration (default: true)', metavar="")
    parser.add_argument("-echk", "--eq_checkpoint", type=int, default=100000,
                        help='Intervals for saving checkpoint files during equilibration (default: 100000)', metavar="")

    # ---------------------------
    # Production simulation parameters
    # ---------------------------
    parser.add_argument("-pm", "--prod_mini", type=int, default=5000000,
                        help='Minimization steps for production (default: 5000000)', metavar="")
    parser.add_argument("-pn", "--prod_nstep", type=int, default=5000000,
                        help='Production simulation steps (default: 5000000)', metavar="")
    parser.add_argument("-pos", "--prod_nstout", type=int, default=5000,
                        help='Log output interval for production (default: 5000)', metavar="")
    parser.add_argument("-pod", "--prod_nstdcd", type=int, default=5000,
                        help='Trajectory output interval for production (default: 5000)', metavar="")
    parser.add_argument("-ps", "--prod_switching", action="store_true",
                        help='Enable production switching (default: false)')
    parser.add_argument("-pp", "--prod_press_couple", action="store_true",
                        help='Enable production pressure coupling (default: false)')
    parser.add_argument("-ppi", "--prod_p_iso", type=lambda s: s.lower() in ['true','1','yes'], default=True,
                        help='Enable isotropic pressure for production (default: true)', metavar="")
    parser.add_argument("-pax", "--prod_p_aniso_x", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along x for production (default: false)', metavar="")
    parser.add_argument("-pay", "--prod_p_aniso_y", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along y for production (default: false)', metavar="")
    parser.add_argument("-paz", "--prod_p_aniso_z", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Enable anisotropic pressure along z for production (default: false)', metavar="")
    parser.add_argument("-pg", "--prod_genpsf", type=lambda s: s.lower() in ['true','1','yes'], default=False,
                        help='Generate PSF file during production (default: false)', metavar="")
    parser.add_argument("-pchk", "--prod_checkpoint", type=int, default=100000,
                        help='Intervals for saving checkpoint files during production (default: 100000)', metavar="")

    # ---------------------------
    # Box and chain options
    # ---------------------------
    parser.add_argument("-b", "--box_type", type=str, choices=["cube", "slab"], default="cube",
                        help='Box type: "cube" or "slab" (default: "cube")', metavar="")
    parser.add_argument("-g", "--geometry", type=str, choices=["coil", "straight"], default="coil",
                        help='Initial chain geometry: "coil" or "straight" (default: "coil")', metavar="")
    parser.add_argument("-d", "--box_dim", type=str, default=None,
                        help='Custom box dimensions in nm as comma-separated values (e.g., "15.0,15.0,15.0"). Optional; if omitted, dimensions are computed based on box_type.',
                        metavar="")
    parser.add_argument("-ca", "--cube_add", type=float, default=0,
                        help='Extra nm for cube box dimensions (default: 0)', metavar="")
    parser.add_argument("-sa", "--slab_add", type=float, default=1000,
                        help='Extra nm added to Z dimension for slab box (default: 1000)', metavar="")
    parser.add_argument("-cd", "--chain_dist", type=float, default=0.5,
                        help='Chain separation in nm (default: 0.5)', metavar="")
    parser.add_argument("-cp", "--chain_pos", type=str, choices=["center", "random"], default="center",
                        help='Chain position (default: "center")', metavar="")

    # ---------------------------
    # Output file names
    # ---------------------------
    parser.add_argument("-tp", "--temp_pdb", type=str, default="temp.pdb",
                        help='Temporary PDB filename (default: "temp.pdb")', metavar="")
    parser.add_argument("-ip", "--initial_pdb", type=str, default="initial.pdb",
                        help='Simulation-ready initial PDB filename (default: "initial.pdb")', metavar="")
    parser.add_argument("-pf", "--psf_file", type=str, default="topol.psf",
                        help='PSF filename (default: "topol.psf")', metavar="")
    parser.add_argument("-ed", "--eq_dcd", type=str, default="eqm.dcd",
                        help='Equilibration trajectory filename (default: "eqm.dcd")', metavar="")
    parser.add_argument("-el", "--eq_log", type=str, default="eqm.log",
                        help='Equilibration log filename (default: "eqm.log")', metavar="")
    parser.add_argument("-epdb", "--eq_pdb", type=str, default="eqm.pdb",
                        help='Equilibration PDB filename (default: "eqm.pdb")', metavar="")
    parser.add_argument("-pd", "--prod_dcd", type=str, default="prod.dcd",
                        help='Production trajectory filename (default: "prod.dcd")', metavar="")
    parser.add_argument("-pl", "--prod_log", type=str, default="prod.log",
                        help='Production log filename (default: "prod.log")', metavar="")
    parser.add_argument("-ppdb", "--prod_pdb", type=str, default="prod.pdb",
                        help='Production PDB filename (default: "prod.pdb")', metavar="")
    parser.add_argument("-epar", "--eq_param", type=str, default="eqm.dat",
                        help='Equilibration parameter file (default: "eqm.dat")', metavar="")
    parser.add_argument("-ppar", "--prod_param", type=str, default="prod.dat",
                        help='Production parameter file (default: "prod.dat")', metavar="")

    # ---------------------------
    # Logging
    # ---------------------------
    parser.add_argument("-l", "--log", type=str, default="simulation.log",
                        help='Main simulation log file (default: "simulation.log")', metavar="")
    
    # ---------------------------
    # Output YAML file name
    # ---------------------------
    parser.add_argument("-o", "--output", type=str, default="config.yaml",
                        help='Output YAML file name (default: "config.yaml")', metavar="")

    args = parser.parse_args()

    # Process box_dim if provided: convert comma-separated string to list of floats.
    if args.box_dim:
        try:
            args.box_dim = [float(x.strip()) for x in args.box_dim.split(',')]
            if len(args.box_dim) != 3:
                raise ValueError
        except ValueError:
            parser.error("box_dim must be three comma-separated numbers, e.g., '15.0,15.0,15.0'")

    # Build configuration dictionary in sequential order.
    config = {
        "seq_file": args.seq_file,
        "codes_path": args.codes_path,
        "resources": args.resources,
        "forcefield": args.forcefield,
        "temp": args.temp,
        "press": args.press,
        "time_step": args.time_step,
        "friction": args.friction,
        "eq_mini": args.eq_mini,
        "eq_nstep": args.eq_nstep,
        "eq_nstout": args.eq_nstout,
        "eq_nstdcd": args.eq_nstdcd,
        "eq_switching": args.eq_switching,
        "eq_press_couple": args.eq_press_couple,
        "eq_p_iso": args.eq_p_iso,
        "eq_p_aniso_x": args.eq_p_aniso_x,
        "eq_p_aniso_y": args.eq_p_aniso_y,
        "eq_p_aniso_z": args.eq_p_aniso_z,
        "eq_genpsf": args.eq_genpsf,
        "eq_chk_int": args.eq_checkpoint,
        "prod_mini": args.prod_mini,
        "prod_nstep": args.prod_nstep,
        "prod_nstout": args.prod_nstout,
        "prod_nstdcd": args.prod_nstdcd,
        "prod_switching": args.prod_switching,
        "prod_press_couple": args.prod_press_couple,
        "prod_p_iso": args.prod_p_iso,
        "prod_p_aniso_x": args.prod_p_aniso_x,
        "prod_p_aniso_y": args.prod_p_aniso_y,
        "prod_p_aniso_z": args.prod_p_aniso_z,
        "prod_genpsf": args.prod_genpsf,
        "prod_chk_int": args.prod_checkpoint,
        "box_type": args.box_type,
        "geometry": args.geometry,
        "box_dim": args.box_dim,
        "cube_add": args.cube_add,
        "slab_add": args.slab_add,
        "chain_dist": args.chain_dist,
        "chain_pos": args.chain_pos,
        "temp_pdb": args.temp_pdb,
        "initial_pdb": args.initial_pdb,
        "psf_file": args.psf_file,
        "eq_dcd": args.eq_dcd,
        "eq_log": args.eq_log,
        "eq_pdb": args.eq_pdb,
        "prod_dcd": args.prod_dcd,
        "prod_log": args.prod_log,
        "prod_pdb": args.prod_pdb,
        "eq_param": args.eq_param,
        "prod_param": args.prod_param,
        "log": args.log
    }
    
    # Write the configuration dictionary to YAML file with keys in the defined order.

    if Path(args.output).exists():
        base, ext = os.path.splitext(args.output)
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        os.rename(args.output, f"{base}_{timestamp}{ext}")
        print(f"Existing configuration file backed up as {base}_{timestamp}{ext}")

    with open(args.output, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False, sort_keys=False)
    print(f"Configuration file '{args.output}' generated successfully.")

if __name__ == "__main__":
    main()
