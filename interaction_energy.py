#!/usr/bin/env python3

import warnings
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm
import matplotlib.pyplot as plt
from itertools import combinations

# Force‐field constants
kbond       = 4184.0      # kJ/mol/nm^2
l0_pro      = 0.38        # nm
l0_rna      = 0.5         # nm
theta0      = np.deg2rad(180.0)
kangle_pro  = 4.184       # kJ/mol/rad^2
kangle_rna  = 5.021       # kJ/mol/rad^2

kappa       = 1.0         # Debye length (nm)
# no elastic bonds in this version

ion_pi_pp   = 0.30        # cation-pi protein-protein epsilon
ion_pi_pr   = 0.20        # cation-pi protein-RNA epsilon
pi_pi_pp    = 0.10        # pi-pi protein-protein epsilon

# Coarse‐grained bead parameters (full dictionary from cocomo2.py)
ff_param = {
    'ALA': {'mass':71.079, 'charge':0.0, 'radius':0.2845, 'epsilon':0.295, 'azero':0.0002},
    'ARG': {'mass':157.197,'charge':1.0,'radius':0.3567,'epsilon':0.176,'azero':0.0},
    'ASN': {'mass':114.104,'charge':0.0,'radius':0.3150,'epsilon':0.176,'azero':0.0},
    'ASP': {'mass':114.080,'charge':-1.0,'radius':0.3114,'epsilon':0.176,'azero':0.0},
    'CYS': {'mass':103.139,'charge':0.0,'radius':0.3024,'epsilon':0.295,'azero':0.0002},
    'GLN': {'mass':128.131,'charge':0.0,'radius':0.3311,'epsilon':0.176,'azero':0.0},
    'GLU': {'mass':128.107,'charge':-1.0,'radius':0.3279,'epsilon':0.176,'azero':0.0},
    'GLY': {'mass':57.052, 'charge':0.0,'radius':0.2617,'epsilon':0.295,'azero':0.0002},
    'HIS': {'mass':137.142,'charge':0.0,'radius':0.3338,'epsilon':0.176,'azero':0.0},
    'ILE': {'mass':113.160,'charge':0.0,'radius':0.3360,'epsilon':0.295,'azero':0.0002},
    'LEU': {'mass':113.160,'charge':0.0,'radius':0.3363,'epsilon':0.295,'azero':0.0002},
    'LYS': {'mass':129.183,'charge':1.0,'radius':0.3439,'epsilon':0.176,'azero':0.0},
    'MET': {'mass':131.193,'charge':0.0,'radius':0.3381,'epsilon':0.295,'azero':0.0002},
    'PHE': {'mass':147.177,'charge':0.0,'radius':0.3556,'epsilon':0.295,'azero':0.0002},
    'PRO': {'mass':98.125, 'charge':0.0,'radius':0.3187,'epsilon':0.295,'azero':0.0002},
    'SER': {'mass':87.078, 'charge':0.0,'radius':0.2927,'epsilon':0.176,'azero':0.0},
    'THR': {'mass':101.105,'charge':0.0,'radius':0.3108,'epsilon':0.176,'azero':0.0},
    'TRP': {'mass':186.214,'charge':0.0,'radius':0.3754,'epsilon':0.295,'azero':0.0002},
    'TYR': {'mass':163.176,'charge':0.0,'radius':0.3611,'epsilon':0.295,'azero':0.0002},
    'VAL': {'mass':99.133, 'charge':0.0,'radius':0.3205,'epsilon':0.295,'azero':0.0002},
    'ADE': {'mass':315.697,'charge':-1.0,'radius':0.4220,'epsilon':0.41, 'azero':0.05},
    'CYT': {'mass':305.200,'charge':-1.0,'radius':0.4110,'epsilon':0.41, 'azero':0.05},
    'GUA': {'mass':345.200,'charge':-1.0,'radius':0.4255,'epsilon':0.41, 'azero':0.05},
    'URA': {'mass':305.162,'charge':-1.0,'radius':0.4090,'epsilon':0.41, 'azero':0.05},
}

def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', '--trajectory', required=True, metavar="",
                        help='Trajectory file (DCD)')
    parser.add_argument('-s', '--topology',   required=True, metavar="",
                        help='Topology file (PSF)')
    parser.add_argument('-c', '--pdbfile', required=True, metavar="",
                        help='Structure file (PDB)')
    parser.add_argument('-p', '--pairs',      required=True, metavar="",
            help='Comma-separated chain pairs, e.g. A:A,B:C, or P000:P001, etc.')
    parser.add_argument('-m', '--mode',       default='nonbonded', metavar="",
                        choices=['bonded', 'nonbonded', 'both'],
                        help='Energy mode to compute. Default = nonbonded')
    parser.add_argument('-o', '--output',     default='interaction_energy.dat', metavar="",
                        help='Output file name. Default = interaction_energy.dat')
    return parser.parse_args()

def get_chain_indices(pdb, cid):
    u_ = mda.Universe(pdb)
    if len(cid) == 1:
        ag = u_.select_atoms(f"chainID {cid}")
    else:
        ag = u_.select_atoms(f"segid {cid}")
    if len(ag) == 0:
        raise ValueError(f"No atoms found for segment or chain {cid}")
    return ag.indices

def prepare_bonded_terms(u):
    u.trajectory[0]
    bond_terms = []
    angle_terms = []
    # native bonds
    for b in u.bonds:
        i, j = b.atoms
        resn = i.resname
        r0   = l0_pro if len(resn) == 3 else l0_rna
        bond_terms.append((i.index, j.index, kbond, r0))
    # backbone angles
    for seg in {a.segid for a in u.atoms}:
        idxs = [a.index for a in u.select_atoms(f"segid {seg}")]
        for a, b, c in zip(idxs, idxs[1:], idxs[2:]):
            resn = u.atoms[b].resname
            kval = kangle_pro if len(resn) == 3 else kangle_rna
            angle_terms.append((a, b, kval, theta0))
    return bond_terms, angle_terms

def make_exclusions(bond_terms):
    excl = set()
    for i, j, *_ in bond_terms:
        excl.add((i, j)); excl.add((j, i))
    return excl

def bonded_energy(pos, idx1, idx2, bond_terms, angle_terms):
    eb = 0.0
    if np.array_equal(idx1, idx2):
        s = set(idx1)
        for i, j, k, r0 in bond_terms:
            if i in s and j in s:
                d = np.linalg.norm(pos[i] - pos[j])
                eb += 0.5 * k * (d - r0)**2
        for i, j, k, th0 in angle_terms:
            if {i, j, k} <= s:
                v1 = pos[i] - pos[j]; v2 = pos[k] - pos[j]
                c  = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                theta = np.arccos(np.clip(c, -1, 1))
                eb += 0.5 * k * (theta - th0)**2
    return eb

def nonbonded_energy(pos, idx1, idx2,
                     charges, radii, epsilons, azeros,
                     exclusions, box,
                     is_ion, is_aromatic, is_rna):
    P1, P2 = pos[idx1], pos[idx2]
    D = distance_array(P1, P2, box=box)
    if np.array_equal(idx1, idx2):
        np.fill_diagonal(D, np.inf)    
    mask = D < 3.0
    # exclude bonded and self
    for i, j in exclusions:
        if i in idx1 and j in idx2:
            a = idx1.tolist().index(i); b = idx2.tolist().index(j)
            mask[a, b] = False
    if np.array_equal(idx1, idx2):
        for ii in range(len(idx1)):
            mask[ii, ii] = False
    # Debye-Hückel mixing
    A1 = np.sqrt(azeros[idx1]); A2 = np.sqrt(azeros[idx2])
    Z1 = np.sign(charges[idx1]) * np.sqrt(0.75 * np.abs(charges[idx1]))
    Z2 = np.sign(charges[idx2]) * np.sqrt(0.75 * np.abs(charges[idx2]))
    Qmix = A1[:,None] * A2[None,:] + (Z1[:,None] + Z2[None,:])
    Udh   = Qmix * np.exp(-D / kappa) / D
    # 10-5 Lennard-Jones
    sigma_i = radii * 2 * 2**(-1/6)
    S       = 0.5 * (sigma_i[idx1][:,None] + sigma_i[idx2][None,:])
    E0      = np.sqrt(epsilons[idx1][:,None] * epsilons[idx2][None,:])
    Ulj     = 4 * E0 * ((S/D)**10 - (S/D)**5)
    # special overrides
    mask_ip = np.outer(is_ion[idx1],      is_aromatic[idx2])
    mask_pp = np.outer(is_aromatic[idx1], is_aromatic[idx2])
    mask_ir = np.outer(is_ion[idx1],      is_rna[idx2])
    Ulj[mask_ip] = 4 * ion_pi_pp * ((S[mask_ip]/D[mask_ip])**10 - (S[mask_ip]/D[mask_ip])**5)
    Ulj[mask_pp] = 4 * pi_pi_pp  * ((S[mask_pp]/D[mask_pp])**10 - (S[mask_pp]/D[mask_pp])**5)
    Ulj[mask_ir] = 4 * ion_pi_pr * ((S[mask_ir]/D[mask_ir])**10 - (S[mask_ir]/D[mask_ir])**5)
    return np.where(mask, Udh + Ulj, 0.0).sum()

def main():
    args = parse_args()

    warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message="DCDReader currently makes independent timesteps*"
    )

    # load Universe
    u = mda.Universe(args.topology, args.trajectory)
    if not hasattr(u.trajectory, '_auxs'):
        u.trajectory._auxs = {}

    # verify all residues present
    res_set = {a.resname for a in u.atoms}
    missing = res_set - set(ff_param.keys())
    if missing:
        raise KeyError(f"Missing ff_param entries for: {sorted(missing)}")

    # static bead arrays
    charges  = np.array([ff_param[a.resname]['charge']  for a in u.atoms])
    radii    = np.array([ff_param[a.resname]['radius']  for a in u.atoms])
    epsilons = np.array([ff_param[a.resname]['epsilon'] for a in u.atoms])
    azeros   = np.array([ff_param[a.resname]['azero']   for a in u.atoms])

    # bonded & angle terms
    bond_terms, angle_terms = prepare_bonded_terms(u)
    exclusions              = make_exclusions(bond_terms)

    # chain‐pair indices
    chain_pairs = []
    for token in args.pairs.split(','):
        c1, c2 = token.split(':')
        idx1 = get_chain_indices(args.pdbfile, c1)
        idx2 = get_chain_indices(args.pdbfile, c2)
        chain_pairs.append(((c1, c2), idx1, idx2))

    # group masks
    resnames    = np.array([a.resname for a in u.atoms])
    is_ion      = np.isin(resnames, ['ARG','LYS'])
    is_aromatic = np.isin(resnames, ['PHE','TRP','TYR'])
    is_rna      = np.isin(resnames, ['ADE','CYT','GUA','URA'])

    # PBC box [lx,ly,lz,alpha,beta,gamma] in nm/deg
    box = u.dimensions.copy()
    box[:3] /= 10.0

    # loop over frames
    n_frames = len(u.trajectory)
    all_results = []
    for ts in tqdm(range(n_frames)):
        u.trajectory[ts]
        pos = u.atoms.positions / 10.0
        frame_res = {}
        for (c1, c2), idx1, idx2 in chain_pairs:
            eb = bonded_energy(pos, idx1, idx2, bond_terms, angle_terms)
            en = nonbonded_energy(pos, idx1, idx2,
                                  charges, radii, epsilons, azeros,
                                  exclusions, box,
                                  is_ion, is_aromatic, is_rna)
            if args.mode == 'bonded':
                E = eb
            elif args.mode == 'nonbonded':
                E = en
            else:
                E = eb + en
            frame_res[(c1, c2)] = E
        all_results.append(frame_res)

    times = [u.trajectory[i].time for i in range(n_frames)]

    # write out .dat files
    for i, ((c1, c2), _, _) in enumerate(chain_pairs):
        Ets = np.array([fr[(c1, c2)] for fr in all_results])
        mu, sigma = Ets.mean(), Ets.std()/np.sqrt(len(Ets))
        fname = {args.output}
        with open(fname, 'w') as f:
            f.write(f"# mean={mu:.3f} kJ/mol   std={sigma:.3f} kJ/mol\n")
            f.write("# time(ps)    energy(kJ/mol)\n")
            for t, e in zip(times, Ets):
                f.write(f"{t:.3f}    {e:.6f}\n")

if __name__ == "__main__":
    main()
