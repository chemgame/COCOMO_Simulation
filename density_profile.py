#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
import MDAnalysis as mda
import warnings
from scipy.optimize import curve_fit

warnings.filterwarnings('ignore')

# Argument parsing
parser = argparse.ArgumentParser(description="Generate averaged density profiles along the z-axis for each protein chain")
parser.add_argument('-trj' , type=str  , required=True, metavar='', help="Path to the trajectory file (DCD)")
parser.add_argument('-top' , type=str  , required=True, metavar='', help="Path to the topology file (PDB)")
parser.add_argument('-ch'  , type=str  , nargs='+', required=True, metavar='', help="List of protein names")
parser.add_argument('-bins', type=int  , default=500, metavar='', help="Number of bins along the z-axis (Default = %(default)s)")
parser.add_argument('-b'   , type=int  , default=0, metavar='', help="Start frame for analysis (Default = %(default)s)")
parser.add_argument('-e'   , type=int  , default=None, metavar='', help="End frame for analysis (Default = last frame)")
parser.add_argument('-o'   , type=str  , metavar='', help="Output file name for the plot (Default = dens_ChainA-ChainB-...jpg)")

args = parser.parse_args()

# Load trajectory and topology
traj = md.load(args.trj, top=args.top)

# Determine the number of frames
total_frames = traj.n_frames

# Parameters
bins = args.bins
b = args.b
e = args.e if args.e is not None else total_frames  # Set e to the last frame if not provided
chains_input = args.ch

u = mda.Universe(args.top, args.trj)
dt = u.trajectory.dt

# Set up z-axis bins
z_box = traj.unitcell_lengths[0, 2]  # in nm
z_range = np.linspace(0, z_box, bins + 1)
dz = z_box / bins
z_centers = (z_range[:-1] + z_range[1:]) / 2

# Identify chains
chains = {chain: traj.topology.select(f'chainid == {chain_index}') 
          for chain_index, chain in enumerate(chains_input)}

# Initialize densities storage
densities = {chain: np.zeros(bins) for chain in chains.keys()}
total_density = np.zeros(bins)

# Masses dictionary (in mg)
masses = {
    'ALA': 71.079, 'ARG': 157.197, 'ASN': 114.104, 'ASP': 114.080, 'CYS': 103.139, 
    'GLN': 128.131, 'GLU': 128.107, 'GLY': 57.052, 'HIS': 137.142, 'ILE': 113.160, 
    'LEU': 113.160, 'LYS': 129.183, 'MET': 131.193, 'PHE': 147.177, 'PRO': 98.125, 
    'SER': 87.078, 'THR': 101.105, 'TRP': 186.214, 'TYR': 163.176, 'VAL': 99.133
}
amu_to_mg = 1.66053906660e-21
masses = {k: v * amu_to_mg for k, v in masses.items()}

# Calculate densities
for frame_idx, frame in enumerate(tqdm(traj[b:e], desc='Processing frames')):
    for chain, indices in chains.items():
        z_pos = frame.xyz[0, indices, 2]  # positions in nm
        masses_per_residue = np.array([masses[frame.topology.atom(index).residue.name] for index in indices])
        hist, _ = np.histogram(z_pos, bins=z_range, weights=masses_per_residue)
        densities[chain] += hist
        total_density += hist

# Normalize densities
vol_per_bin = dz * traj.unitcell_lengths[0, 0] * traj.unitcell_lengths[0, 1]  # volume in nm^3
vol_per_bin_ml = vol_per_bin * 1e-21  # convert to mL

for chain in densities:
    densities[chain] /= (e - b) * vol_per_bin_ml  # convert to mg/mL
total_density /= (e - b) * vol_per_bin_ml

# Fit total density to a double hyperbolic tangent function
def double_tanh_fit(z, rho0, z0, z1, width):
    return rho0 * (np.tanh((z - z0) / width) - np.tanh((z - z1) / width) + 1) / 2

popt, _ = curve_fit(double_tanh_fit, z_centers, total_density, 
                    p0=[np.max(total_density), z_box / 4, 3 * z_box / 4, z_box / 10])
rho0, z0, z1, width = popt
thickness = z1 - z0  # Thickness of the slab

# Plotting
fig, ax = plt.subplots(dpi=300, figsize=(8, 6))

for chain, density in densities.items():
    ax.plot(z_centers, density, label=f'{chain}')

# Plot total density and fit
ax.plot(z_centers, total_density, 'k--', label='Total Density')
#ax.plot(z_centers, double_tanh_fit(z_centers, *popt), 'k--', lw = 0.5, label=f'Fit (Thickness = {thickness:.2f} nm)')

# Labels and legend
ax.set_xlabel('Z (nm)', fontsize=15)
ax.set_ylabel('Density (mg/mL)', fontsize=15)
ax.legend(fontsize=12, frameon=False, loc='upper right')

# Save the figure
if not args.o:
    out = f'dens-prof_{'-'.join(chains_input)}.jpg'
else:
    out = args.o
fig.savefig(out, format="jpg", dpi=300, bbox_inches='tight')

# Output thickness
print(f"Estimated thickness of the protein slab: {thickness:.2f} nm")
