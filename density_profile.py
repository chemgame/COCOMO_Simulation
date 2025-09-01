#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
import MDAnalysis as mda
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore')

# Argument parsing
parser = argparse.ArgumentParser(description="Generate averaged density profiles along the z-axis for each protein chain")
parser.add_argument('-trj'  , type=str  , required=True, metavar='', help="Path to the trajectory file (DCD)")
parser.add_argument('-top'  , type=str  , required=True, metavar='', help="Path to the topology file (PDB)")
parser.add_argument('-ch'   , type=str  , nargs='+', required=True, metavar='', help="List of protein names")
parser.add_argument('-bins' , type=int  , default=500, metavar='', help="Number of bins along the z-axis (Default = %(default)s)")
parser.add_argument('-b'    , type=int  , default=0, metavar='', help="Start frame for analysis (Default = %(default)s)")
parser.add_argument('-e'    , type=int  , default=None, metavar='', help="End frame for analysis (Default = last frame)")
parser.add_argument('-o'    , type=str  , metavar='', help="Output file name for the plot (Default = dens_ChainA-ChainB-...jpg)")
parser.add_argument('-batch', type=int, default=50, metavar='', help="Number of frames to process at once (Default = %(default)s)")

args = parser.parse_args()

# Load trajectory and topology
try:
    traj = md.load(args.trj, top=args.top)
except:
    # Fallback for PSF files - convert via MDAnalysis
    u = mda.Universe(args.top, args.trj)
    # Convert to temporary PDB for MDTraj
    temp_pdb = 'temp_topology.pdb'
    u.select_atoms('all').write(temp_pdb)
    traj = md.load(args.trj, top=temp_pdb)
    import os
    os.remove(temp_pdb)

# Determine frame range
total_frames = traj.n_frames
b = args.b
e = args.e if args.e is not None else total_frames
chains_input = args.ch

# Set up z-axis bins
z_box = traj.unitcell_lengths[0, 2]  # in nm
bins = args.bins
z_edges = np.linspace(0, z_box, bins + 1)
dz = z_box / bins
z_centers = (z_edges[:-1] + z_edges[1:]) / 2

# Residue masses in mg
residue_masses = {
    'ALA': 71.079, 'ARG': 157.197, 'ASN': 114.104, 'ASP': 114.080, 'CYS': 103.139, 
    'GLN': 128.131, 'GLU': 128.107, 'GLY': 57.052, 'HIS': 137.142, 'ILE': 113.160, 
    'LEU': 113.160, 'LYS': 129.183, 'MET': 131.193, 'PHE': 147.177, 'PRO': 98.125, 
    'SER': 87.078, 'THR': 101.105, 'TRP': 186.214, 'TYR': 163.176, 'VAL': 99.133
}
amu_to_mg = 1.66053906660e-21

# Pre-compute all chain data
chain_indices = []
chain_masses = []
chain_starts = [0]
chain_names = []

for chain_index, chain_name in enumerate(chains_input):
    indices = traj.topology.select(f'chainid == {chain_index}')
    masses = np.array([residue_masses.get(traj.topology.atom(idx).residue.name, 0.0) * amu_to_mg 
                      for idx in indices])
    
    chain_indices.extend(indices)
    chain_masses.extend(masses)
    chain_starts.append(len(chain_indices))
    chain_names.append(chain_name)

# Convert to numpy arrays for vectorized operations
all_indices = np.array(chain_indices, dtype=int)
all_masses = np.array(chain_masses)
chain_starts = np.array(chain_starts)

# Initialize density arrays
densities = {chain: np.zeros(bins) for chain in chain_names}
total_density = np.zeros(bins)

# Pre-compute volume normalization
vol_per_bin = dz * traj.unitcell_lengths[0, 0] * traj.unitcell_lengths[0, 1]  # nm^3
vol_per_bin_ml = vol_per_bin * 1e-21  # convert to mL
normalization_factor = 1.0 / ((e - b) * vol_per_bin_ml)

# Ultra-fast vectorized calculation
batch_size = args.batch
n_frames_to_process = e - b

for batch_start in tqdm(range(b, e, batch_size), desc='Processing batches'):
    batch_end = min(batch_start + batch_size, e)
    actual_batch_size = batch_end - batch_start
    
    # Load batch of frames - shape: (batch_size, n_atoms, 3)
    batch_coords = traj.xyz[batch_start:batch_end]
    
    # Extract z-coordinates for all relevant atoms - shape: (batch_size, n_relevant_atoms)
    z_coords = batch_coords[:, all_indices, 2]
    
    # Vectorized histogram calculation across all frames and atoms
    # Flatten to process all frame-atom combinations at once
    z_flat = z_coords.flatten()  # shape: (batch_size * n_relevant_atoms,)
    masses_tiled = np.tile(all_masses, actual_batch_size)  # repeat masses for each frame
    
    # Single vectorized histogram call
    total_hist, _ = np.histogram(z_flat, bins=z_edges, weights=masses_tiled)
    total_density += total_hist
    
    # Calculate per-chain histograms using advanced indexing
    for i, chain_name in enumerate(chain_names):
        start_idx = chain_starts[i]
        end_idx = chain_starts[i + 1]
        
        # Extract z-coords and masses for this chain
        chain_z = z_coords[:, start_idx:end_idx].flatten()
        chain_masses_batch = np.tile(all_masses[start_idx:end_idx], actual_batch_size)
        
        # Vectorized histogram for this chain
        chain_hist, _ = np.histogram(chain_z, bins=z_edges, weights=chain_masses_batch)
        densities[chain_name] += chain_hist

# Apply normalization (vectorized)
total_density *= normalization_factor
for chain in densities:
    densities[chain] *= normalization_factor

# Fit total density to double hyperbolic tangent function
def double_tanh_fit(z, rho0, z0, z1, width):
    return rho0 * (np.tanh((z - z0) / width) - np.tanh((z - z1) / width) + 1) / 2

try:
    popt, _ = curve_fit(double_tanh_fit, z_centers, total_density, 
                        p0=[np.max(total_density), z_box / 4, 3 * z_box / 4, z_box / 10])
    rho0, z0, z1, width = popt
    thickness = z1 - z0
except:
    thickness = 0.0
    print("Warning: Fitting failed")

# Optimized plotting
fig, ax = plt.subplots(dpi=200, figsize=(8, 6))

# Plot all chains at once using vectorized operations
colors = plt.cm.Dark2(np.linspace(0, 1, len(chain_names)))
for i, (chain, density) in enumerate(densities.items()):
    ax.plot(z_centers, density, label=chain, lw=2, color=colors[i])

# Plot total density
ax.plot(z_centers, total_density, 'k--', label='Total', lw=1)

# Styling
ax.set_xlabel('Z (nm)', fontsize=20)
ax.set_ylabel('Density (mg/mL)', fontsize=20)
ax.legend(fontsize=15, frameon=True, loc='upper right', facecolor='lightgrey')
ax.tick_params(labelsize=15)
ax.grid(True, lw=0.5, ls='--')

# Save figure
if not args.o:
    out = f'dens-prof_{"-".join(chains_input)}.jpg'
else:
    out = args.o

fig.savefig(out, format="jpg", dpi=300, bbox_inches='tight')

# Output result
print(f"Estimated thickness of the protein slab: {thickness:.2f} nm")
