#!/usr/bin/env python
"""
center.py

This module provides functions to center a molecular system within the simulation box
based on a user-defined MDAnalysis atom selection. The main functionality is available
through the `center_system` function, which can be imported into other scripts.

Usage as a script:
    $ ./center.py -s topology.psf -f trajectory.dcd -sel "protein" -o centered.dcd

Usage as a module:
    >>> from center import center_system
    >>> from MDAnalysis import Universe
    >>> # Create a Universe (with a static structure or trajectory)
    >>> u = Universe("topology.psf", "trajectory.dcd")
    >>> summary = center_system(u, "protein", begin=0, end=-1, skip=1, output_file="centered.dcd")
    >>> print(summary)
"""

import argparse
import logging
import warnings

import numpy as np
from MDAnalysis import Universe, Writer
from tqdm import tqdm

# Suppress warnings for a clean output
warnings.filterwarnings("ignore")

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Center a molecular system by shifting the selected atoms' center to the box center."
    )
    parser.add_argument("-f", "--trajectory", metavar="",
                        help="Trajectory file (e.g., dcd or xtc). Omit for a static structure")
    parser.add_argument("-s", "--topology", required=True, metavar="",
                        help="Topology file (e.g., psf or tpr)")
    parser.add_argument("-sel", "--atom-selection", required=True, metavar="",
                        help="MDAnalysis atom selection string")
    parser.add_argument("-o", "--output", required=True, metavar="",
                        help="Output file name. Format will match the input trajectory if provided")
    parser.add_argument("-b", "--begin", type=int, default=0, metavar="",
                        help="Index of the starting frame for processing")
    parser.add_argument("-e", "--end", type=int, default=-1, metavar="",
                        help="Index of the ending frame for processing (-1 for last frame)")
    parser.add_argument("-sk", "--skip", type=int, default=1, metavar="",
                        help="Skip interval between frames")
    return parser.parse_args()


def process_frame(u, sel):
    """
    Process the current frame.
    
    Computes the center of the simulation box and calculates the center of
    the selected atoms (using center-of-mass if masses are available; otherwise,
    the geometric center). Shifts all atom positions such that the selection's
    center is aligned with the box's geometric center.

    Parameters:
        u (Universe): The MDAnalysis Universe.
        sel (AtomGroup): The selected atoms.

    Returns:
        numpy.ndarray: The shift vector applied to the atoms.
    """
    ts = u.trajectory.ts
    # The dimensions array is [Lx, Ly, Lz, alpha, beta, gamma]; use the first three
    box_dim = ts.dimensions[:3]
    box_center = box_dim / 2.0

    # Calculate the center of the selection (center_of_mass if possible; otherwise, geometric)
    if np.sum(sel.masses) > 0:
        sel_center = sel.center_of_mass()
    else:
        sel_center = sel.center_of_geometry()

    # Compute the shift vector so that the selection is centered in the box.
    shift = box_center - sel_center
    # Apply the shift to all atoms (vectorized operation).
    u.atoms.positions += shift
    return shift


def center_system(universe, selection, begin=0, end=-1, skip=1, output_file=None):
    """
    Center the molecular system based on the specified atom selection.
    
    This function processes either a static structure or a trajectory.
    It shifts the entire system so that the center (or center of mass, if available)
    of the selected atoms aligns with the geometric center of the simulation box.
    If a trajectory is provided, it can process a specific range of frames.
    
    Parameters:
        universe (Universe): An MDAnalysis Universe instance.
        selection (str or AtomGroup): An MDAnalysis atom selection string or an already-selected AtomGroup.
        begin (int): Index of the starting frame to process (only applicable for trajectories).
        end (int): Index of the ending frame to process (-1 means the last frame; only for trajectories).
        skip (int): Frame skip interval (only for trajectories).
        output_file (str): Filename for the output. The output format will match the input.
                           If None, no output file is written.
    
    Returns:
        dict: A summary containing the number of processed frames and the average shift vector.
    """
    # Allow the selection to be either a string or an MDAnalysis AtomGroup.
    if isinstance(selection, str):
        sel = universe.select_atoms(selection)
    else:
        sel = selection

    logging.info("Number of atoms selected: %d", len(sel))

    # Determine frame processing indices if a trajectory is available.
    if hasattr(universe.trajectory, "n_frames") and universe.trajectory.n_frames > 1:
        total_frames = universe.trajectory.n_frames
        start = begin
        end = total_frames if end == -1 or end > total_frames else end
        frame_indices = range(start, end, skip)
        n_frames = len(frame_indices)
    else:
        frame_indices = [0]
        n_frames = 1

    # Set up the Writer if an output_file is given.
    writer = Writer(output_file, universe.atoms.n_atoms) if output_file is not None else None

    logging.info("Processing %d frame(s)...", n_frames)
    total_shift = np.zeros(3)
    count = 0

    # Process each frame with a progress bar.
    for frame in tqdm(frame_indices, desc="Processing frames"):
        universe.trajectory[frame]  # update Universe to the current frame
        shift = process_frame(universe, sel)
        total_shift += shift
        if writer is not None:
            writer.write(universe.atoms)
        count += 1

    if writer is not None:
        writer.close()

    summary = {
        "frames_processed": count,
        "average_shift": total_shift / count if count > 0 else total_shift
    }
    logging.info("Processed %d frame(s)", count)
    logging.info("Average shift applied: %s", summary["average_shift"])
    return summary


def main():
    """
    Command-line interface for centering the molecular system.
    """
    args = parse_args()

    # Create the MDAnalysis Universe.
    if args.trajectory:
        u = Universe(args.topology, args.trajectory)
    else:
        u = Universe(args.topology)
    
    # Center the system using the provided selection and frame parameters.
    center_system(u, args.atom_selection, begin=args.begin, end=args.end,
                  skip=args.skip, output_file=args.output)


if __name__ == "__main__":
    main()
