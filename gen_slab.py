#!/usr/bin/env python
"""
Module: gen_slab.py
Description:
    Generate a PDB file with protein/RNA chains arranged in a simulation box,
    and recenter the entire system so that its center aligns with the geometric center 
    of the simulation box.
    
    Supports:
      -g/--geometry   : "straight" (linear rods) or "coil" (self-avoiding random walk)
      -p/--position   : "center" or "random"

Key Point:
    If a single line in seq_file has chain_count > 1, *all* those chains use the *same* chain ID.
    Only a new line in the sequence file triggers a different chain ID.

Usage:
    ./gen_slab.py -s sequences.dat -box 20 20 30 -o output.pdb -d 0.1 -g coil -p random
"""

import numpy as np
import MDAnalysis as mda
import argparse
import sys
import warnings
import logging

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# ---------------------
# Constants
# ---------------------
STEP_LENGTH = 0.38       # nm for coil geometry steps
MIN_CENTER_SEP = 0.7     # nm: chain centers must be at least this far apart
CENTER_MARGIN = 2.0      # nm: margin for "center" positioning in x-y plane
PDB_MAX_COORD = 9999.9   # Å: maximum safe coordinate for PDB format

# ---------------------
# Residue Mappings
# ---------------------
aa_map = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

rna_map = {
    'A': 'ADE', 'C': 'CYT', 'G': 'GUA', 'U': 'URA'
}

# ---------------------
# Self-Avoiding Random Walk
# ---------------------
def self_avoiding_random_walk(n_steps, step_length, min_dist=0.35, max_attempts=1000):
    """
    Generate a self-avoiding random walk in 3D.
    
    :param n_steps: number of residues
    :param step_length: distance (nm) between consecutive residues
    :param min_dist: minimal separation (nm) between new residue and existing ones
    :param max_attempts: maximum attempts to find a valid next step
    :return: (N,3) array of residue positions in nm
    """
    positions = [np.array([0.0, 0.0, 0.0])]
    attempts = 0
    while len(positions) < n_steps:
        if attempts >= max_attempts:
            raise ValueError("Failed to generate self-avoiding coil after maximum attempts.")
        direction = np.random.randn(3)
        direction /= np.linalg.norm(direction)
        new_pos = positions[-1] + step_length * direction
        if all(np.linalg.norm(new_pos - np.array(positions), axis=1) >= min_dist):
            positions.append(new_pos)
            attempts = 0
        else:
            attempts += 1
    return np.array(positions)

# ---------------------
# Main Box Generation
# ---------------------
def gen_box(seq_file, box_size, out_file, chain_distance, geometry, position):
    """
    Generate a simulation box PDB file based on a sequence file,
    recenter the complete system, and write the final output.

    :param seq_file: Path to the sequences file (each line: [PRO|RNA] chain_count sequence).
    :param box_size: (Lx, Ly, Lz) in nm
    :param out_file: Output PDB file name
    :param chain_distance: Distance between two chains (nm) if exactly two
    :param geometry: "straight" (linear rods) or "coil" (self-avoiding random walk)
    :param position: "center" or "random" positioning of chains in the box.
    """
    if any(dim <= 0 for dim in box_size):
        sys.exit("ERROR: box dimensions must be positive values.")

    # --- 1) Parse sequences file ---
    seq_types, chain_counts, sequences = [], [], []
    with open(seq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                sys.exit(f"ERROR: invalid line in {seq_file}: '{line}' (3 columns needed).")
            seq_type, ccount_str, seq_str = parts
            try:
                ccount = int(ccount_str)
            except ValueError:
                sys.exit(f"ERROR: chain_count must be int, got '{ccount_str}'.")
            if seq_type not in ("PRO", "RNA"):
                sys.exit(f"ERROR: sequence type must be 'PRO' or 'RNA', got '{seq_type}'.")
            seq_types.append(seq_type)
            chain_counts.append(ccount)
            sequences.append(seq_str)

    tot_chains = sum(chain_counts)
    tot_atoms = sum(len(seq) * cnt for seq, cnt in zip(sequences, chain_counts))
    Lx, Ly, Lz = box_size

    logger.info("Total chains: %d", tot_chains)
    logger.info("Total atoms : %d", tot_atoms)
    logger.info("Box (nm)    : %.3f x %.3f x %.3f", Lx, Ly, Lz)

    # --- 2) Assign Chain IDs Based on Sequence Lines ---
    alphabets = [chr(a) for a in range(ord('A'), ord('Z') + 1)] + \
                [chr(a) for a in range(ord('a'), ord('z') + 1)]
    chain_id_list = []
    for i, ccount in enumerate(chain_counts):
        chain_id_for_line = alphabets[i % len(alphabets)]
        chain_id_list.extend([chain_id_for_line] * ccount)
    
    if len(chain_id_list) != tot_chains:
        sys.exit("ERROR in chain ID assignment; mismatch in total chain count.")

    # --- 3) Determine chain centers ---
    if position == "center":
        if tot_chains == 2:
            centers = np.array([
                [-chain_distance / 2, 0.0, Lz / 2],
                [ chain_distance / 2, 0.0, Lz / 2]
            ])
        else:
            centers_2d = []
            while len(centers_2d) < tot_chains:
                x_rand = np.random.rand(1000) * (Lx - CENTER_MARGIN) - (Lx - CENTER_MARGIN) / 2
                y_rand = np.random.rand(1000) * (Ly - CENTER_MARGIN) - (Ly - CENTER_MARGIN) / 2
                new_pos = np.column_stack((x_rand, y_rand))
                
                if centers_2d:
                    dists = np.linalg.norm(
                        new_pos[:, None] - np.array(centers_2d)[None, :],
                        axis=-1
                    )
                    new_pos = new_pos[np.all(dists > MIN_CENTER_SEP, axis=1)]
                
                needed = tot_chains - len(centers_2d)
                centers_2d.extend(new_pos[:needed])
            
            centers_2d = np.array(centers_2d)
            z_vals = np.full(len(centers_2d), Lz / 2)
            centers = np.column_stack((centers_2d, z_vals))
    else:
        centers = []
        while len(centers) < tot_chains:
            new_pos = np.random.rand(1000, 3) * [Lx, Ly, Lz] - [Lx/2, Ly/2, Lz/2]
            
            if centers:
                dists = np.linalg.norm(
                    new_pos[:, None] - np.array(centers)[None, :],
                    axis=-1
                )
                new_pos = new_pos[np.all(dists > MIN_CENTER_SEP, axis=1)]
            
            needed = tot_chains - len(centers)
            centers.extend(new_pos[:needed])
        centers = np.array(centers)

    # --- 4) Build Coordinates for Each Chain ---
    ch_ids, atom_names, resnames, resids, pos = [], [], [], [], []
    total_chain_counter = 0  # global chain counter
    for seq_idx, (seq_type, ccount, seq_str) in enumerate(zip(seq_types, chain_counts, sequences)):
        mapping = aa_map if seq_type == "PRO" else rna_map
        
        for _ in range(ccount):
            seq_len = len(seq_str)
            center_x, center_y, center_z = centers[total_chain_counter]
            
            if geometry == "coil":
                coords = self_avoiding_random_walk(seq_len, STEP_LENGTH)
                coords -= coords.mean(axis=0)
                coords += np.array([center_x, center_y, center_z])
            else:  # "straight"
                z_vals = np.linspace(-seq_len/2 * STEP_LENGTH, seq_len/2 * STEP_LENGTH, seq_len)
                rod = np.column_stack((np.zeros(seq_len), np.zeros(seq_len), z_vals))
                rod -= rod.mean(axis=0)
                rod += np.array([center_x, center_y, center_z])
                coords = rod
            
            ch_id = chain_id_list[total_chain_counter]
            ch_ids.extend([ch_id] * seq_len)
            atom_names.extend(["CG"] * seq_len)
            resnames.extend([mapping.get(r, "UNK") for r in seq_str])
            resids.extend(range(1, seq_len + 1))
            pos.append(coords)
            
            total_chain_counter += 1

    # Combine coordinates and convert units from nm to Å
    positions_nm = np.vstack(pos)  # (total_atoms, 3) in nm
    positions_A = positions_nm * 10.0  # nm -> Å

    # --- 5) Scale if needed to fit PDB format ---
    max_abs_coord = np.max(np.abs(positions_A))
    if max_abs_coord > PDB_MAX_COORD:
        scale_factor = PDB_MAX_COORD / max_abs_coord
        logger.warning("System exceeds PDB limit (%.1f Å). Scaling by %.5f.",
                       max_abs_coord, scale_factor)
        positions_A *= scale_factor
        Lx *= scale_factor
        Ly *= scale_factor
        Lz *= scale_factor

    tot_atoms = positions_A.shape[0]
    u = mda.Universe.empty(n_atoms=tot_atoms, n_residues=tot_atoms, n_segments=tot_chains,
                           atom_resindex=np.arange(tot_atoms), trajectory=True)

    u.add_TopologyAttr("names", atom_names)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", resids)
    u.add_TopologyAttr("chainIDs", ch_ids)

    u.atoms.positions = positions_A
    # Set box dimensions in Å (here, Lx, Ly, Lz are in nm so multiply by 10)
    u.dimensions = np.array([Lx * 10, Ly * 10, Lz * 10, 90.0, 90.0, 90.0])

    # --- 6) Center the Entire System ---
    # Compute the geometric center of the simulation box.
    box_dim = u.dimensions[:3]
    box_center = box_dim / 2.0
    # Use center_of_mass if masses are available; otherwise, default to center_of_geometry.
    #if np.sum(u.atoms.masses) > 0:
    #    system_center = u.atoms.center_of_mass()
    #else:
    system_center = u.atoms.center_of_geometry()
    shift = box_center - system_center
    logger.info("Applying recentering shift: %s", shift)
    u.atoms.positions += shift

    # --- 7) Write the Final PDB File ---
    with mda.Writer(out_file, multiframe=False) as pdb_writer:
        pdb_writer.write(u.atoms)

    logger.info("PDB file generated: %s", out_file)

# ---------------------
# Argparse
# ---------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a PDB file with protein/RNA chains arranged in a 3D simulation box "
                    "and recenter the system to the geometric center of the box.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s", "--seq_file",
        type=str,
        required=True,
        metavar="",
        help="Path to the sequence file. Each line: [type PRO/RNA] [chain_count] [sequence]."
    )
    parser.add_argument(
        "-o", "--out_file",
        type=str,
        required=True,
        metavar="",
        help="Output PDB file name."
    )
    parser.add_argument(
        "-box", "--box_size",
        type=float,
        nargs=3,
        required=True,
        metavar="",
        help="Box dimensions in nm (Lx Ly Lz)."
    )
    parser.add_argument(
        "-d", "--chain_distance",
        type=float,
        default=0.1,
        metavar="",
        help="Distance (nm) between 2 chains if exactly two. Default: 0.1 nm."
    )
    parser.add_argument(
        "-g", "--geometry",
        choices=["straight", "coil"],
        default="straight",
        metavar="",
        help="Geometry for chain generation: 'straight' or 'coil'. Default: 'straight'."
    )
    parser.add_argument(
        "-p", "--position",
        choices=["center", "random"],
        default="center",
        metavar="",
        help="Positioning of chains in the box: 'center' or 'random'. Default: 'center'."
    )
    return parser.parse_args()

# ---------------------
# Main
# ---------------------
def main():
    logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')
    args = parse_args()
    logger.info("Running gen_box with geometry=%s, position=%s.", args.geometry, args.position)

    gen_box(
        seq_file=args.seq_file,
        box_size=tuple(args.box_size),
        out_file=args.out_file,
        chain_distance=args.chain_distance,
        geometry=args.geometry,
        position=args.position
    )

if __name__ == "__main__":
    main()
