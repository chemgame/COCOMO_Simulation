#!/usr/bin/env python3
"""
Calculate the mass-weighted total radius of gyration (Rg) of the system using MDAnalysis.
The output file (rg_vs_time.dat) contains two columns: time (ps) and Rg (Å). The file begins
with a commented line giving the average Rg over the processed frames. The average value is also
printed to the console.

Usage:
  python rg_calculator.py -f trajectory.dcd -s topology.pdb -b 0 -e 1000 -dt 1

Arguments:
  -f/--trajectory   Trajectory file (e.g., DCD)
  -s/--topology     Topology file (PDB or PSF) containing mass information.
  -b/--begin        Begin time in ps (default: 0)
  -e/--end          End time in ps (if not provided, process until the end)
  -dt/--skip        Process every nth frame (skip frames in between; default: 1)
"""

import argparse
import sys
import os
from datetime import datetime
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm

def get_args():
    parser = argparse.ArgumentParser(
        description="Calculate total mass-weighted radius of gyration (Rg) of the system using MDAnalysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--trajectory", required=True, metavar="",
                        help="Trajectory file (e.g., DCD)")
    parser.add_argument("-s", "--topology", required=True, metavar="",
                        help="Topology file (PDB or PSF)")
    parser.add_argument("-o", "--output", default="gyrate.dat", metavar="",
                        help="Output Rg vs time file (default: gyrate.dat)")
    parser.add_argument("-b", "--begin", type=float, default=0.0, metavar="",
                        help="Begin time (ps)")
    parser.add_argument("-e", "--end", type=float, default=None, metavar="",
                        help="End time (ps). If not given, process to the end of trajectory.")
    parser.add_argument("-dt", "--skip", type=int, default=1, metavar="",
                        help="Process every nth frame (skip frames in between)")
    return parser.parse_args()

def backup_file_if_exists(fname):
    if os.path.exists(fname):
        base, ext = os.path.splitext(fname)
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        new_name = f"{base}_{timestamp}{ext}"
        os.rename(fname, new_name)

def main():
    args = get_args()

    # Load the topology and trajectory with MDAnalysis.
    try:
        u = mda.Universe(args.topology, args.trajectory)
    except Exception as err:
        sys.exit(f"Error loading topology and trajectory: {err}")

    # Prepare containers for time and Rg values.
    times = []
    rg_values = []
    skip_counter = 0

    # Process each frame with a progress bar.
    for ts in tqdm(u.trajectory, desc="Processing frames"):
        # Use provided time if available; otherwise, default to frame index.
        t_ps = ts.time if hasattr(ts, "time") and ts.time is not None else ts.frame
        # Skip frames before the specified begin time.
        if t_ps < args.begin:
            continue
        # Stop if we've exceeded the end time.
        if args.end is not None and t_ps > args.end:
            break
        # Apply frame skipping.
        if skip_counter % args.skip != 0:
            skip_counter += 1
            continue
        skip_counter += 1

        # Calculate the mass-weighted radius of gyration for all atoms.
        rg = u.atoms.radius_of_gyration()
        times.append(t_ps)
        rg_values.append(rg)

    # Check if any frames were processed.
    if len(rg_values) == 0:
        sys.exit("No frames were processed. Check your time range and skip parameters.")

    times = np.array(times)
    rg_values = np.array(rg_values)
    avg_rg = np.mean(rg_values)

    # Print the average Rg to the console.
    print(f"Average radius of gyration: {avg_rg:.6f} Å")

    # Output file name.
    output_filename = args.output
    # Backup existing file if needed.
    backup_file_if_exists(output_filename)

    # Write the results to file with average Rg included in the header.
    try:
        with open(output_filename, "w") as f:
            f.write(f"# Average radius of gyration: {avg_rg:.6f} Å\n")
            f.write("# Time (ps)\tRg (Å)\n")
            for t, rg in zip(times, rg_values):
                f.write(f"{t:.3f}\t{rg:.6f}\n")
    except Exception as err:
        sys.exit(f"Error writing output file: {err}")

if __name__ == "__main__":
    main()
