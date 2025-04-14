#!/usr/bin/env python3
"""
Module for calculating the mass-weighted total radius of gyration (Rg) from an MD trajectory.

This module provides a function, calculate_rg, which computes the radius of gyration for a given trajectory
and topology file using MDAnalysis. It returns the average Rg, time values, and corresponding Rg values.
Additionally, the module can be executed directly from the command-line.

Usage as a script:
  python gyrate.py -f trajectory.dcd -s topology.pdb -b 0 -e 1000 -dt 1

Usage as a module:
  from gyrate import calculate_rg
  avg_rg, times, rg_values = calculate_rg(
      trajectory="trajectory.dcd",
      topology="topology.pdb",
      begin=0.0,
      end=1000,
      skip=1,
      output="gyrate.dat"
  )
"""

import argparse
import sys
import os
from datetime import datetime
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm


def backup_file_if_exists(fname):
    """Rename the output file with a timestamp if it already exists."""
    if os.path.exists(fname):
        base, ext = os.path.splitext(fname)
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        new_name = f"{base}_{timestamp}{ext}"
        os.rename(fname, new_name)


def calculate_rg(trajectory, topology, begin=0.0, end=None, skip=1, output="gyrate.dat"):
    """
    Calculate the mass-weighted radius of gyration (Rg) for an MD simulation.

    Parameters:
        trajectory (str): Path to the trajectory file (e.g., DCD).
        topology (str): Path to the topology file (e.g., PDB or PSF) containing mass information.
        begin (float, optional): Start processing at this time (ps). Default is 0.0.
        end (float or None, optional): Stop processing after this time (ps). If None, process until the end.
        skip (int, optional): Process every nth frame. Default is 1.
        output (str, optional): Path for writing the Rg vs time output file. Default is "gyrate.dat".

    Returns:
        tuple: A tuple containing:
            - avg_rg (float): Average radius of gyration over processed frames.
            - times (np.ndarray): Array of time values (ps) corresponding to processed frames.
            - rg_values (np.ndarray): Array of Rg values (Å) corresponding to processed frames.
    """
    try:
        u = mda.Universe(topology, trajectory)
    except Exception as err:
        sys.exit(f"Error loading topology and trajectory: {err}")

    times = []
    rg_values = []
    skip_counter = 0

    # Process each frame with a progress bar.
    for ts in tqdm(u.trajectory, desc="Processing frames"):
        t_ps = ts.time if hasattr(ts, "time") and ts.time is not None else ts.frame

        if t_ps < begin:
            continue
        if end is not None and t_ps > end:
            break
        if skip_counter % skip != 0:
            skip_counter += 1
            continue
        skip_counter += 1

        rg = u.atoms.radius_of_gyration()
        times.append(t_ps)
        rg_values.append(rg)

    if not rg_values:
        raise ValueError("No frames were processed. Check your time range and skip parameters.")

    times = np.array(times)
    rg_values = np.array(rg_values)
    avg_rg = np.mean(rg_values)
    print(f"Average radius of gyration: {avg_rg:.6f} Å")

    backup_file_if_exists(output)

    try:
        with open(output, "w") as f:
            f.write(f"# Average radius of gyration: {avg_rg:.6f} Å\n")
            f.write("# Time (ps)\tRg (Å)\n")
            for t, rg in zip(times, rg_values):
                f.write(f"{t:.3f}\t{rg:.6f}\n")
    except Exception as err:
        sys.exit(f"Error writing output file: {err}")

    return avg_rg, times, rg_values


def get_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate the mass-weighted total radius of gyration (Rg) using MDAnalysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--trajectory", required=True, metavar="",
                        help="Trajectory file (e.g., DCD)")
    parser.add_argument("-s", "--topology", required=True, metavar="",
                        help="Topology file (PDB or PSF) containing mass information.")
    parser.add_argument("-o", "--output", default="gyrate.dat", metavar="",
                        help="Output file for Rg vs time (default: gyrate.dat)")
    parser.add_argument("-b", "--begin", type=float, default=0.0, metavar="",
                        help="Begin time in ps")
    parser.add_argument("-e", "--end", type=float, default=None, metavar="",
                        help="End time in ps. If not specified, process until the end of trajectory.")
    parser.add_argument("-dt", "--skip", type=int, default=1, metavar="",
                        help="Process every nth frame (default: 1)")
    return parser.parse_args()


def main():
    """Main function to parse arguments and calculate the radius of gyration."""
    args = get_args()
    try:
        calculate_rg(
            trajectory=args.trajectory,
            topology=args.topology,
            begin=args.begin,
            end=args.end,
            skip=args.skip,
            output=args.output,
        )
    except Exception as err:
        sys.exit(err)


if __name__ == "__main__":
    main()
