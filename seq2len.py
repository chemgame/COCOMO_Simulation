#!/usr/bin/env python
"""
Module: seq2len.py
Description: Estimate minimum box length from a protein sequence.
Usage:
    As a standalone script:
        ./seq2len.py -s TENSTSAPAAKPKRAKASKKSTDHPPKYSDMI
    Or import and call calc_len(seq)
"""

import argparse
import logging
import sys

logger = logging.getLogger(__name__)

def calc_len(seq):
    """Calculate the minimum box length from the protein sequence using residue radii."""
    radii = {'A': 0.2845, 'R': 0.3567, 'N': 0.3150, 'D': 0.3114,
             'C': 0.3024, 'Q': 0.3311, 'E': 0.3279, 'G': 0.2617,
             'H': 0.3338, 'I': 0.3360, 'L': 0.3363, 'K': 0.3439,
             'M': 0.3381, 'F': 0.3556, 'P': 0.3187, 'S': 0.2927,
             'T': 0.3108, 'W': 0.3754, 'Y': 0.3611, 'V': 0.3205,
             'a': 0.4220, 'c': 0.4110, 'g': 0.4255, 'u': 0.4090}
    seq = seq.strip()
    for aa in seq:
        if aa not in radii:
            logger.error("Unknown amino acid: %s", aa)
            sys.exit(1)
    length = sum(radii[aa] for aa in seq)
    return 2 * length

def main():
    parser = argparse.ArgumentParser(
        description="Estimate minimum box length from a protein sequence",
        epilog="Copyright reserved by Saumyak Mukherjee"
    )
    parser.add_argument('-s', '--sequence', type=str, required=True,
                        help='Protein sequence (str)', metavar='')
    args = parser.parse_args()
    chainlength = calc_len(args.sequence)
    print(f"{chainlength:.4f}")

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')
    main()
