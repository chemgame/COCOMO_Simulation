#!/usr/bin/env python
"""
Simulation Pipeline Driver

Usage:
  simulation.py --config <config.yaml>

All simulation parameters are provided via the YAML configuration file.
Defaults are assigned based on the provided sample config.yaml file.
"""

import os
import sys
import argparse
import logging
import yaml
import tempfile
from pathlib import Path
from datetime import datetime
import importlib

logger = logging.getLogger(__name__)

# --- Helper Functions ---
def error_exit(msg):
    logger.error(msg)
    sys.exit(1)

def check_file(fname, errmsg):
    if not Path(fname).exists():
        error_exit(errmsg)

def write_parameter_file(filename, content):
    with open(filename, "w") as f:
        f.write(content)
    logger.info("Wrote parameter file: %s", filename)

def read_sequences(seq_file, seq2len_mod):
    chain_lengths = []
    box_min_vals = []
    with open(seq_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                error_exit(f"Invalid line in sequences file: {line}")
            seq = parts[2]
            chain_lengths.append(len(seq))
            box_min_vals.append(seq2len_mod.calc_len(seq))
    return chain_lengths, box_min_vals

def compute_box_dimensions(box_type, cube_add, slab_add, max_box):
    if box_type == "cube":
        return max_box + cube_add, max_box + cube_add, max_box + cube_add
    else:
        return max_box, max_box, max_box + slab_add

def run_simulation(param_file, expected_dcd, cocomo_mod):
    try:
        cocomo_mod.run_sim(param_file)
        check_file(expected_dcd, f"Expected DCD file '{expected_dcd}' not found.")
    except Exception as err:
        error_exit(f"Simulation error with '{param_file}': {err}")

def process_trajectory(in_dcd, top_file, temp_pdb, out_pdb, chain_lengths, dcd2pdbFrame_mod, genpdb_cocomo_mod):
    dcd2pdbFrame_mod.extract_frame(in_dcd, top_file, temp_pdb, dumpframe=-1)
    check_file(temp_pdb, "Temporary PDB not created during trajectory processing.")
    genpdb_cocomo_mod.gen_pdb(temp_pdb, out_pdb, chain_lengths=chain_lengths, elastic_ch=" ")
    check_file(out_pdb, f"PDB file '{out_pdb}' not produced.")
    os.remove(temp_pdb)

def get_temp_pdb_file():
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tf:
        return tf.name

def load_config(config_file):
    try:
        with open(config_file) as cf:
            config = yaml.safe_load(cf)
    except Exception as e:
        error_exit(f"Failed to load configuration file '{config_file}': {e}")

    # Default configuration values
    defaults = {
        "seq_file": "sequences.dat",
        "codes_path": "/home/tb/samukher/CODE.BANK/bin/",
        "resources": "CUDA",
        "forcefield": "cocomo2",
        "temp": 300.0,
        "press": 1.0,
        "time_step": 0.02,
        "friction": 0.01,
        "eq_mini": 5000000,
        "eq_nstep": 500000,
        "eq_nstout": 5000,
        "eq_nstdcd": 5000,
        "eq_switching": False,
        "eq_press_couple": False,
        "eq_p_iso": True,
        "eq_p_aniso_x": False,
        "eq_p_aniso_y": False,
        "eq_p_aniso_z": False,
        "eq_genpsf": True,
        "eq_chk_int": 100000,
        "prod_mini": 5000000,
        "prod_nstep": 5000000,
        "prod_nstout": 5000,
        "prod_nstdcd": 5000,
        "prod_switching": False,
        "prod_press_couple": False,
        "prod_p_iso": True,
        "prod_p_aniso_x": False,
        "prod_p_aniso_y": False,
        "prod_p_aniso_z": False,
        "prod_genpsf": False,
        "prod_chk_int": 100000,
        "box_type": "cube",
        "geometry": "coil",
        "box_dim": None,
        "cube_add": 0,
        "slab_add": 1000,
        "chain_dist": 0.5,
        "chain_pos": "center",
        "temp_pdb": "temp.pdb",
        "initial_pdb": "initial.pdb",
        "psf_file": "topol.psf",
        "eq_dcd": "eqm.dcd",
        "eq_log": "eqm.log",
        "eq_pdb": "eqm.pdb",
        "prod_dcd": "prod.dcd",
        "prod_log": "prod.log",
        "prod_pdb": "prod.pdb",
        "eq_param": "eqm.dat",
        "prod_param": "prod.dat",
        "log": "simulation.log"
    }
    for key, val in defaults.items():
        config.setdefault(key, val)
    return config

# --- Argument Parsing ---
def parse_args():
    desc = """
    COCOMO Simulation Pipeline
    (A sample configuration file can be generated using the gen_config.py script)
    """
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="© Saumyak Mukherjee"
    )
    parser.add_argument("-c", "--config", type=str, default="config.yaml",
                        help="Path to YAML configuration file (default: config.yaml)", metavar="")
    return parser.parse_args()

# --- Main Function ---
def main():
    args = parse_args()
    config = load_config(args.config)

    # Setup logging: clear any existing handlers to avoid duplicate console output.
    root_logger = logging.getLogger()
    if root_logger.hasHandlers():
        root_logger.handlers.clear()
    log_file = config.get("log", "simulation.log")
    if Path(log_file).exists():
        base, ext = os.path.splitext(log_file)
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        os.rename(log_file, f"{base}_{timestamp}{ext}")
        print(f"Existing log file backed up as {base}_{timestamp}{ext}")
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s',
                        filename=log_file,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
    logging.getLogger().addHandler(console)

    # Ensure simulation modules are accessible.
    codes_dir = os.path.abspath(config["codes_path"])
    if codes_dir not in sys.path:
        sys.path.insert(0, codes_dir)
    logger.info("Using simulation modules from: %s", codes_dir)

    # Get simulation parameters from configuration.
    seq_file     = config["seq_file"]
    box_type     = config["box_type"]
    geometry     = config.get("geometry", "coil")
    chain_pos    = config.get("chain_pos", "center")
    temp         = config["temp"]
    resources    = config["resources"]
    press        = config["press"]
    time_step    = config["time_step"]
    friction     = config["friction"]

    eq_mini      = config["eq_mini"]
    eq_nstep     = config["eq_nstep"]
    eq_nstout    = config["eq_nstout"]
    eq_nstdcd    = config["eq_nstdcd"]
    eq_switch    = config["eq_switching"]
    eq_press     = config["eq_press_couple"]
    eq_p_iso     = config["eq_p_iso"]
    eq_p_aniso_x = config["eq_p_aniso_x"]
    eq_p_aniso_y = config["eq_p_aniso_y"]
    eq_p_aniso_z = config["eq_p_aniso_z"]
    eq_genpsf    = config["eq_genpsf"]
    eq_chk_int   = config["eq_chk_int"]

    prod_mini      = config["prod_mini"]
    prod_nstep     = config["prod_nstep"]
    prod_nstout    = config["prod_nstout"]
    prod_nstdcd    = config["prod_nstdcd"]
    prod_switch    = config["prod_switching"]
    prod_press     = config["prod_press_couple"]
    prod_p_iso     = config["prod_p_iso"]
    prod_p_aniso_x = config["prod_p_aniso_x"]
    prod_p_aniso_y = config["prod_p_aniso_y"]
    prod_p_aniso_z = config["prod_p_aniso_z"]
    prod_genpsf    = config["prod_genpsf"]
    prod_chk_int   = config["prod_chk_int"]

    cube_add    = config["cube_add"]
    slab_add    = config["slab_add"]
    chain_dist  = config["chain_dist"]

    # Box dimensions: use custom values if provided; otherwise compute based on sequences.
    if "box_dim" in config and config["box_dim"]:
        box_x, box_y, box_z = config["box_dim"]
        logger.info("Using custom box dimensions from config: X=%.2f, Y=%.2f, Z=%.2f", box_x, box_y, box_z)
    else:
        box_x = box_y = box_z = None  # Will be computed after processing sequences.
    
    # Dynamically import simulation modules.
    seq2len_mod = importlib.import_module("seq2len")
    gen_slab_mod = importlib.import_module("gen_slab")
    genpdb_cocomo_mod = importlib.import_module("genpdb_cocomo")
    dcd2pdbFrame_mod = importlib.import_module("dcd2pdbFrame")

    # Determine forcefield and import the corresponding module.
    forcefield = config.get("forcefield", "cocomo2")
    if forcefield == "cocomo":
        cocomo_mod = importlib.import_module("cocomo")
    elif forcefield == "cocomo2":
        cocomo_mod = importlib.import_module("cocomo2")
    else:
        error_exit(f"Error in forcefield selection: {forcefield} is not correct.")
    logger.info(f"Selected forcefield: {forcefield}")
    
    # Process sequences.
    chain_lengths, box_min_vals = read_sequences(seq_file, seq2len_mod)
    if not config.get("box_dim"):
        max_box_val = max(box_min_vals)
        box_x, box_y, box_z = compute_box_dimensions(box_type, cube_add, slab_add, max_box_val)
        logger.info("Computed box dimensions (nm): X=%.2f, Y=%.2f, Z=%.2f", box_x, box_y, box_z)
    
    # Generate initial simulation box.
    temp_box_pdb = get_temp_pdb_file()
    try:
        gen_slab_mod.gen_box(seq_file, (box_x, box_y, box_z), temp_box_pdb, chain_distance=chain_dist, geometry=geometry, position=chain_pos)
        check_file(temp_box_pdb, f"{temp_box_pdb} not created by gen_slab")
    except Exception as err:
        error_exit(f"Error in gen_slab: {err}")

    # Generate simulation-ready PDB.
    try:
        genpdb_cocomo_mod.gen_pdb(temp_box_pdb, config["initial_pdb"], chain_lengths=chain_lengths, elastic_ch=" ")
        check_file(config["initial_pdb"], f"{config['initial_pdb']} not created by genpdb_cocomo")
        os.remove(temp_box_pdb)
    except Exception as err:
        error_exit(f"Error in genpdb_cocomo: {err}")
    
    # Create equilibration parameter file.
    eq_content = "\n".join([
        f"resources    = {resources}",
        f"temp         = {temp}",
        f"press        = {press}",
        f"dt           = {time_step}",
        f"fric         = {friction}",
        f"mini_nstep   = {eq_mini}",
        f"nstep        = {eq_nstep}  # 10 ns",
        f"nstout       = {eq_nstout}",
        f"nstdcd       = {eq_nstdcd}",
        f"pdb_in       = {config['initial_pdb']}",
        f"psf_out      = {config['psf_file']}",
        f"dcd_out      = {config['eq_dcd']}",
        f"log_out      = {config['eq_log']}",
        f"SWITCHING    = {eq_switch}",
        f"press_couple = {eq_press}",
        f"p_iso        = {eq_p_iso}",
        f"p_aniso_x    = {eq_p_aniso_x}",
        f"p_aniso_y    = {eq_p_aniso_y}",
        f"p_aniso_z    = {eq_p_aniso_z}",
        f"genpsf       = {eq_genpsf}",
        f"chk_interval = {eq_chk_int}"
    ])
    write_parameter_file(config["eq_param"], eq_content)
    logger.info("Running equilibration simulation")
    run_simulation(config["eq_param"], config["eq_dcd"], cocomo_mod)
    
    # Process equilibration trajectory with a new temporary pdb file.
    temp_traj_eq = get_temp_pdb_file()
    process_trajectory(config["eq_dcd"], config["initial_pdb"], temp_traj_eq, config["eq_pdb"],
                         chain_lengths, dcd2pdbFrame_mod, genpdb_cocomo_mod)
    
    # Create production parameter file.
    prod_content = "\n".join([
        f"resources    = {resources}",
        f"temp         = {temp}",
        f"press        = {press}",
        f"dt           = {time_step}",
        f"fric         = {friction}",
        f"mini_nstep   = {prod_mini}",
        f"nstep        = {prod_nstep}  # 100 ns",
        f"nstout       = {prod_nstout}",
        f"nstdcd       = {prod_nstdcd}",
        f"pdb_in       = {config['eq_pdb']}",
        f"psf_out      = {config['psf_file']}",
        f"dcd_out      = {config['prod_dcd']}",
        f"log_out      = {config['prod_log']}",
        f"SWITCHING    = {prod_switch}",
        f"press_couple = {prod_press}",
        f"p_iso        = {prod_p_iso}",
        f"p_aniso_x    = {prod_p_aniso_x}",
        f"p_aniso_y    = {prod_p_aniso_y}",
        f"p_aniso_z    = {prod_p_aniso_z}",
        f"genpsf       = {prod_genpsf}",
        f"chk_interval = {prod_chk_int}"
    ])
    write_parameter_file(config["prod_param"], prod_content)
    logger.info("Running production simulation")
    run_simulation(config["prod_param"], config["prod_dcd"], cocomo_mod)
    
    # Process production trajectory with a new temporary pdb file.
    temp_traj_prod = get_temp_pdb_file()
    process_trajectory(config["prod_dcd"], config["initial_pdb"], temp_traj_prod, config["prod_pdb"],
                         chain_lengths, dcd2pdbFrame_mod, genpdb_cocomo_mod)
    
    logger.info("Simulation complete. Production PDB: %s", config["prod_pdb"])

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')
    main()
