#!/usr/bin/env python
"""
Module: genpdb_cocomo.py
Description: Generate an input PDB file for COCOCMO simulation.
Usage:
    As a standalone script:
        echo "100\n200" | ./genpdb_cocomo.py -i input.pdb -o output.pdb -el "A"
    Or import and call gen_pdb(in_file, out_file, chain_lengths, elastic_ch)
"""

import mdtraj as md
import argparse
import os
import sys
import logging
import tqdm
from datetime import datetime

logger = logging.getLogger(__name__)

def gen_pdb(in_file, out_file, chain_lengths=None, elastic_ch=" "):
    """
    Generate a COCOCMO-compatible PDB file from an input PDB.
    
    Parameters:
      in_file (str): Input PDB file.
      out_file (str): Output PDB file.
      chain_lengths (list of int): List of chain lengths (number of residues).
                                   If None, the function will attempt to read from stdin or prompt interactively.
      elastic_ch (str): Chains to consider for elastic bonds (default: " ").
    """
    if not os.path.exists(in_file):
        logger.error("Input file %s does not exist!", in_file)
        sys.exit(1)
    
    # Backup output if exists
    if os.path.exists(out_file):
        base, ext = os.path.splitext(out_file)
        timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        new_name = f"{base}_{timestamp}{ext}"
        os.rename(out_file, new_name)
        logger.info("%s exists; backed up to %s", out_file, new_name)
    
    traj = md.load(in_file)
    top = traj.topology
    
    # Get unique chain IDs
    chains = []
    for chain in top.chains:
        if chain.chain_id not in chains:
            chains.append(chain.chain_id)
    
    # Get chain lengths from stdin if not provided
    if chain_lengths is None:
        if not sys.stdin.isatty():
            try:
                chain_lengths = [int(line.strip()) for line in sys.stdin if line.strip()]
            except Exception as e:
                logger.error("Error reading chain lengths from stdin: %s", e)
                sys.exit(1)
        else:
            chain_lengths = []
            for ch in chains:
                val = input(f'Number of residues in chain {ch}: ')
                chain_lengths.append(int(val))
    
    if len(chain_lengths) != len(chains):
        logger.error("Mismatch: %d chain lengths provided for %d chains.", len(chain_lengths), len(chains))
        sys.exit(1)
    nres = {ch: chain_lengths[i] for i, ch in enumerate(chains)}
    
    segidx = 0
    residx = 1
    amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                   "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                   "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    rna_bases = {"A": "ADE", "C": "CYT", "G": "GUA", "U": "URA"}
    
    # Modify each atom's properties in the topology
    for atom in tqdm.tqdm(top.atoms, desc=f'Writing {out_file}'):
        atom.name = 'CG'
        if atom.residue.chain.chain_id in elastic_ch:
            segtype = 'E'
            if atom.residue.name in rna_bases:
                atom.residue.name = rna_bases[atom.residue.name]
        else:
            if atom.residue.name in amino_acids:
                segtype = 'P'
            elif atom.residue.name in rna_bases:
                segtype = 'R'
                atom.residue.name = rna_bases[atom.residue.name]
            else:
                segtype = 'X'
        segment_id = f'{segtype}{segidx:03d}'
        atom.residue.segment_id = segment_id
        atom.residue.resSeq = residx
        residx += 1
        if residx > nres[atom.residue.chain.chain_id]:
            residx = 1
            segidx += 1
    traj.save(out_file)
    logger.info("PDB file %s written successfully!", out_file)

def main():
    parser = argparse.ArgumentParser(
        description="Generate input PDB for COCOCMO simulation",
        epilog="Copyright reserved by Saumyak Mukherjee"
    )
    grp_io = parser.add_argument_group("I/O")
    grp_io.add_argument('-i', '--infile', type=str, required=True,
                        help='Input PDB file (str)', metavar='')
    grp_io.add_argument('-o', '--outfile', type=str, default='cocomo_input.pdb',
                        help='Output PDB file (str) [default: cocomo_input.pdb]', metavar='')
    grp_net = parser.add_argument_group("Elastic network")
    grp_net.add_argument('-el', '--elastic', type=str, default=' ',
                        help='Chains for elastic bonds (str) [default: " "]', metavar='')
    args = parser.parse_args()
    gen_pdb(args.infile, args.outfile, chain_lengths=None, elastic_ch=args.elastic)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')
    main()
