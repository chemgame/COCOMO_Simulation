#!/usr/bin/env python
"""
Module: dcd2pdbFrame.py
Description: Extract a given time frame from a .dcd trajectory as a PDB file.
Usage:
    As a standalone script:
        ./dcd2pdbFrame.py -i input.dcd -top topology.pdb -o output.pdb -dump 100.0
    Or import and call extract_frame(infile, intop, outfile, dumpframe=-1)
"""

import argparse
import os
import sys
import logging
import warnings
import MDAnalysis as mda

warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)

def extract_frame(infile, intop, outfile, dumpframe=-1):
    """Extract a frame from the DCD trajectory and write it as a PDB file."""
    try:
        u = mda.Universe(intop, infile)
    except Exception as e:
        logger.error("Unable to load trajectory and topology: %s", e)
        sys.exit(1)
    if dumpframe == -1:
        ts = len(u.trajectory) - 1
    else:
        dt = u.trajectory.dt
        ts = int(dumpframe / dt)
        if ts >= len(u.trajectory):
            logger.error("Time frame %s ps exceeds trajectory length.", dumpframe)
            sys.exit(1)
    try:
        u.trajectory[ts]
        with mda.Writer(outfile) as W:
            W.write(u)
    except Exception as e:
        logger.error("Error writing output file: %s", e)
        sys.exit(1)
    logger.info("Extracted frame written to %s", outfile)

def main():
    parser = argparse.ArgumentParser(
        description="Extract a time frame from a .dcd trajectory as a pdb file",
        epilog="Copyright reserved by Saumyak Mukherjee"
    )
    group = parser.add_argument_group("Required arguments")
    group.add_argument('-i', '--infile', type=str, required=True,
                       help='Input DCD file (str)', metavar='')
    group.add_argument('-top', '--intop', type=str, required=True,
                       help='Input topology file (str)', metavar='')
    group.add_argument('-o', '--outfile', type=str, required=True,
                       help='Output PDB file (str)', metavar='')
    group.add_argument('-dump', '--dumpframe', type=float, default=-1,
                       help='Time frame to extract in ps (float) [default: -1]', metavar='')
    args = parser.parse_args()
    extract_frame(args.infile, args.intop, args.outfile, args.dumpframe)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')
    main()
