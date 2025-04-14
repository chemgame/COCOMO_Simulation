#!/usr/bin/env python
"""
Module: cocomo.py
Description: Run the COCOCMO MD simulation.
Usage:
    As a standalone script:
        ./cocomo.py -prm params.dat
    Or import and call run_sim(prm_file)
Note: DO NOT change any code related to force fields or simulation protocol.
"""

import sys
import os
import re
import glob
import logging
from datetime import datetime
import numpy as np
import argparse
import tqdm
import MDAnalysis as mda

from openmm.unit import *
from openmm import *
from openmm.app import *
import parmed as pmd

# Set up logger
logger = logging.getLogger(__name__)

# ----------------------------
# Custom Checkpoint Reporter
# ----------------------------
class CustomCheckpointReporter(CheckpointReporter):
    def __init__(self, file_template, reportInterval):
        super().__init__(file_template, reportInterval)
        self.file_template = file_template
        self._reportInterval = reportInterval

    def report(self, simulation, state):
        current_step = simulation.currentStep
        file_name = os.path.join('chkfiles', self.file_template.format(current_step))
        with open(file_name, 'wb') as f:
            f.write(simulation.context.createCheckpoint())

# ----------------------------
# Read parameters from file
# ----------------------------
def read_params(file_path):
    """Read simulation parameters from the given file."""
    params = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Remove comments and whitespace
                line = re.sub(r'#.*', '', line).strip()
                if '=' in line:
                    key, value = map(str.strip, line.split('=', 1))
                    if value.lower() in ['true', 'false']:
                        value = value.lower() == 'true'
                    else:
                        try:
                            num_val = float(value)
                            if num_val.is_integer():
                                value = int(num_val)
                            else:
                                value = num_val
                        except ValueError:
                            pass
                    params[key] = value
        globals().update(params)
        return params
    except Exception as e:
        logger.error("Error reading parameters: %s", e)
        raise

# ----------------------------
# Run the MD simulation
# ----------------------------
def run_sim(prm_file):
    """Run the COCOCMO MD simulation using parameters from the given file."""
    # Create checkpoint directory if needed
    if not os.path.exists('chkfiles'):
        os.makedirs('chkfiles')
    
    # Check for parameter file
    if not os.path.exists(prm_file):
        logger.error("Parameter file %s does not exist!", prm_file)
        sys.exit(1)
    params = read_params(prm_file)

    base_name = os.path.splitext(os.path.basename(prm_file))[0]

    # Check that the input PDB file exists
    if not os.path.exists(pdb_in):
        logger.error("PDB file %s does not exist!", pdb_in)
        sys.exit(1)
    
    # Back up files if already present
    if genpsf:
        if os.path.exists(psf_out):
            base, ext = os.path.splitext(psf_out)
            timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
            new_name = f"{base}_{timestamp}{ext}"
            os.rename(psf_out, new_name)
            logger.info("%s already exists; backed up to %s", psf_out, new_name)
    
    for fname in [dcd_out, log_out]:
        if os.path.exists(fname):
            base, ext = os.path.splitext(fname)
            timestamp = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
            new_name = f"{base}_{timestamp}{ext}"
            os.rename(fname, new_name)
            logger.info("%s already exists; backed up to %s", fname, new_name)

    # ----------------------------
    # FORCE FIELD PARAMETERS (DO NOT MODIFY)
    # ----------------------------
    kbond = 4184
    l0_pro = 0.38
    l0_rna = 0.5
    theta0 = 180
    kangle_pro = 4.184 
    kangle_rna = 5.021
    eps_polar = 0.40
    eps_nopol = 0.41
    eps_rna = 0.41
    cationpi_propro = 0.30
    cationpi_prorna = 0.20
    kappa = 1.00
    kelastic = 500
    min_elastic = 0.5
    max_elastic = 0.9

    ff_param = {
        'ALA': {'mass':  71.079, 'charge':  0.0, 'radius': 0.2845, 'epsilon': 0.41, 'azero': 0.00}, 
        'ARG': {'mass': 157.197, 'charge':  1.0, 'radius': 0.3567, 'epsilon': 0.40, 'azero': 0.05}, 
        'ASN': {'mass': 114.104, 'charge':  0.0, 'radius': 0.3150, 'epsilon': 0.40, 'azero': 0.05}, 
        'ASP': {'mass': 114.080, 'charge': -1.0, 'radius': 0.3114, 'epsilon': 0.40, 'azero': 0.05}, 
        'CYS': {'mass': 103.139, 'charge':  0.0, 'radius': 0.3024, 'epsilon': 0.40, 'azero': 0.05}, 
        'GLN': {'mass': 128.131, 'charge':  0.0, 'radius': 0.3311, 'epsilon': 0.40, 'azero': 0.05}, 
        'GLU': {'mass': 128.107, 'charge': -1.0, 'radius': 0.3279, 'epsilon': 0.40, 'azero': 0.05}, 
        'GLY': {'mass':  57.052, 'charge':  0.0, 'radius': 0.2617, 'epsilon': 0.41, 'azero': 0.00}, 
        'HIS': {'mass': 137.142, 'charge':  0.0, 'radius': 0.3338, 'epsilon': 0.40, 'azero': 0.05}, 
        'ILE': {'mass': 113.160, 'charge':  0.0, 'radius': 0.3360, 'epsilon': 0.41, 'azero': 0.00}, 
        'LEU': {'mass': 113.160, 'charge':  0.0, 'radius': 0.3363, 'epsilon': 0.41, 'azero': 0.00}, 
        'LYS': {'mass': 129.183, 'charge':  1.0, 'radius': 0.3439, 'epsilon': 0.40, 'azero': 0.05}, 
        'MET': {'mass': 131.193, 'charge':  0.0, 'radius': 0.3381, 'epsilon': 0.41, 'azero': 0.00}, 
        'PHE': {'mass': 147.177, 'charge':  0.0, 'radius': 0.3556, 'epsilon': 0.41, 'azero': 0.00}, 
        'PRO': {'mass':  98.125, 'charge':  0.0, 'radius': 0.3187, 'epsilon': 0.41, 'azero': 0.00}, 
        'SER': {'mass':  87.078, 'charge':  0.0, 'radius': 0.2927, 'epsilon': 0.40, 'azero': 0.05}, 
        'THR': {'mass': 101.105, 'charge':  0.0, 'radius': 0.3108, 'epsilon': 0.40, 'azero': 0.05}, 
        'TRP': {'mass': 186.214, 'charge':  0.0, 'radius': 0.3754, 'epsilon': 0.41, 'azero': 0.00}, 
        'TYR': {'mass': 163.176, 'charge':  0.0, 'radius': 0.3611, 'epsilon': 0.41, 'azero': 0.00}, 
        'VAL': {'mass':  99.133, 'charge':  0.0, 'radius': 0.3205, 'epsilon': 0.41, 'azero': 0.00}, 
        'ADE': {'mass': 315.697, 'charge': -1.0, 'radius': 0.4220, 'epsilon': 0.41, 'azero': 0.05}, 
        'CYT': {'mass': 305.200, 'charge': -1.0, 'radius': 0.4110, 'epsilon': 0.41, 'azero': 0.05}, 
        'GUA': {'mass': 345.200, 'charge': -1.0, 'radius': 0.4255, 'epsilon': 0.41, 'azero': 0.05}, 
        'URA': {'mass': 305.162, 'charge': -1.0, 'radius': 0.4090, 'epsilon': 0.41, 'azero': 0.05}
    }

    # ----------------------------
    # Build topology from input PDB
    # ----------------------------
    def read_pdb(filename):
        atoms = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.split()
                if line[0] == 'ATOM':
                    atoms.append([line[2], line[3], line[-2]])
        atoms = np.array(atoms, dtype='object')
        return atoms

    positions = PDBFile(pdb_in).positions
    atom_list = read_pdb(pdb_in)
    segnames = [atom_list[i,2] != atom_list[i+1,2] for i in range(-1, atom_list.shape[0]-1)]
    if (atom_list[-1,2] == atom_list[0,2]):
        segnames = [atom_list[0,2]]
    else:
        segnames = atom_list[segnames,2]

    top = topology.Topology()
    for seg in segnames:
        chain = top.addChain(seg)
        for atm in atom_list:
            if atm[2] == seg:
                residue = top.addResidue(atm[1], chain)
                top.addAtom(atm[0], element=element.carbon, residue=residue)
        atom_chain = list(chain.atoms())
        for i in range(len(atom_chain)-1):
            top.addBond(atom_chain[i], atom_chain[i+1])
    
    u = mda.Universe(pdb_in)
    boxx, boxy, boxz = u.dimensions[:3] / 10.0
    a = Quantity(np.zeros([3]), nanometers)
    a[0] = boxx * nanometers
    b = Quantity(np.zeros([3]), nanometers)
    b[1] = boxy * nanometers
    c = Quantity(np.zeros([3]), nanometers)
    c[2] = boxz * nanometers
    box = (a, b, c)
    top.setPeriodicBoxVectors(box)

    # ----------------------------
    # Add elastic bonds based on distance criteria
    # ----------------------------
    def elastic(top, positions, min_dist, max_dist):
        min_dist *= unit.nanometers
        max_dist *= unit.nanometers
        bond_distances = {}
        for chain in tqdm.tqdm(top.chains(), desc='Adding elastic bonds'):
            if not chain.id.startswith('E'):
                continue
            atoms = list(chain.atoms())
            for i in range(len(atoms)):
                for j in range(i + 3, len(atoms)):
                    dist = unit.norm(positions[atoms[i].index] - positions[atoms[j].index])
                    if min_dist < dist < max_dist:
                        top.addBond(atoms[i], atoms[j])
                        bond_distances[(atoms[i].index, atoms[j].index)] = dist.value_in_unit(unit.nanometers)
        return bond_distances

    bonddist_elastic = elastic(top, positions, min_elastic, max_elastic)
    logger.info("System topology:\n%s", top)

    # ----------------------------
    # Build the OpenMM system
    # ----------------------------
    system = openmm.System()
    for i in atom_list[:,1]:
        system.addParticle(ff_param[i]['mass']*unit.amu)
    system.setDefaultPeriodicBoxVectors(a, b, c)

    if SWITCHING == True:
        R_ON  = 2.9
        R_OFF = 3.1
    else:
        R_OFF = 3.0

    if genpsf:
        structure = pmd.openmm.load_topology(top, system, xyz=positions)
        structure.save(psf_out)

    # Bond force
    f_bond = openmm.HarmonicBondForce()
    for bond in top.bonds():
        if bond[0].residue.name in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
            f_bond.addBond(bond[0].index, bond[1].index, l0_pro*nanometer, kbond*kilojoules_per_mole/(nanometer**2))
        elif bond[0].residue.name in ['ADE','CYT','GUA','URA']:
            f_bond.addBond(bond[0].index, bond[1].index, l0_rna*nanometer, kbond*kilojoules_per_mole/(nanometer**2))
    system.addForce(f_bond)

    # Elastic network bond force
    f_elastic_bond = openmm.HarmonicBondForce()
    for (atom1, atom2), dist in bonddist_elastic.items():
        f_elastic_bond.addBond(atom1, atom2, dist*nanometer, kelastic*kilojoules_per_mole/(nanometer**2))
    system.addForce(f_elastic_bond)

    # Angle force
    f_angle = openmm.HarmonicAngleForce()
    for atoms in [[i for i in top.atoms() if i.residue.chain.id == seg] for seg in segnames]:
        if atoms[0].residue.name in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
            for i in range(len(atoms)-2):
                f_angle.addAngle(atoms[i].index, atoms[i+1].index, atoms[i+2].index, theta0*degrees, kangle_pro*kilojoule_per_mole/(radian**2))
        elif atoms[0].residue.name in ['ADE','CYT','GUA','URA']:
            for i in range(len(atoms)-2):
                f_angle.addAngle(atoms[i].index, atoms[i+1].index, atoms[i+2].index, theta0*degrees, kangle_rna*kilojoule_per_mole/(radian**2))
    system.addForce(f_angle)

    # Long-range electrostatic force
    k0 = kappa*nanometer
    if SWITCHING == True:
        equation1  = "select( step(r_on-r), longrange+sdel, switch ); "
        equation1 += "longrange = (A+Z)*exp(-r/K0)/r; "
        equation1 += "sdel = sk*((1/r_on)^3-1/(r_off)^3)^2 - (A+Z)/r_on*exp(-r_on/K0); "
        equation1 += "switch = sk*((1/r)^3-1/(r_off)^3)^2; "
        equation1 += "sk = -longrange_deriv_Ron/switch_deriv_Ron; "
        equation1 += "longrange_deriv_Ron = -1*(A+Z)*exp(-r_on/K0)/r_on*(1/K0+1/r_on); "
        equation1 += "switch_deriv_Ron = 6*(1/r_on^3-1/r_off^3)*1/r_on^4; "
        equation1 += "A=A1*A2; "
        equation1 += "Z=Z1+Z2 "
    else:
        equation1  = "(A+Z)/r*exp(-r/K0); "
        equation1 += "A=A1*A2; "
        equation1 += "Z=Z1+Z2"
    
    force1 = CustomNonbondedForce(equation1)
    force1.addGlobalParameter("K0", k0)
    force1.addPerParticleParameter("A")
    force1.addPerParticleParameter("Z")
    force1.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    force1.setCutoffDistance(R_OFF*nanometer)
    if SWITCHING == True:
        force1.addGlobalParameter("r_on", R_ON)
        force1.addGlobalParameter("r_off", R_OFF)
    for atom in top.atoms():
        charge_val = np.sqrt(0.75*np.abs(ff_param[atom.residue.name]['charge']))*np.sign(ff_param[atom.residue.name]['charge'])
        force1.addParticle([charge_val*nanometer*kilojoule/mole, ff_param[atom.residue.name]['azero']*(nanometer*kilojoule/mole)**(1/2)])
    force1.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
    system.addForce(force1)

    # Short-range LJ force
    if SWITCHING == True:
        equation2  = "select( step(r_on-r), shortrange+sdel, switch ); "
        equation2 += "shortrange = 4*epsilon*((sigma/r)^10-(sigma/r)^5); "
        equation2 += "sdel = sk*((1/r_on)^3-1/(r_off)^3)^2 - 4*epsilon*((sigma/r_on)^10-(sigma/r_on)^5); "
        equation2 += "switch = sk*((1/r)^3-1/(r_off)^3)^2; "
        equation2 += "sk = -shortrange_deriv_Ron/switch_deriv_Ron; "
        equation2 += "shortrange_deriv_Ron = 4*epsilon*(-10*(sigma/r_on)^11*1/r_on+5*(sigma/r_on)^5*1/r_on ); "
        equation2 += "switch_deriv_Ron = 6*(1/r_on^3-1/r_off^3)*1/r_on^4; "
        equation2 += "sigma=0.5*(sigma1+sigma2); "
        equation2 += "epsilon=sqrt(epsilon1*epsilon2)"
    else:
        equation2  = "4*epsilon*((sigma/r)^10-(sigma/r)^5); "
        equation2 += "sigma=0.5*(sigma1+sigma2); "
        equation2 += "epsilon=sqrt(epsilon1*epsilon2)"
    
    force2 = CustomNonbondedForce(equation2)
    force2.addPerParticleParameter("sigma")
    force2.addPerParticleParameter("epsilon")
    force2.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    force2.setCutoffDistance(R_OFF*nanometer)
    if SWITCHING == True:
        force2.addGlobalParameter("r_on", R_ON)
        force2.addGlobalParameter("r_off", R_OFF)
    for atom in top.atoms():
        force2.addParticle([ff_param[atom.residue.name]['radius']*2*2**(-1/6)*nanometer, ff_param[atom.residue.name]['epsilon']*kilojoule/mole])
    force2.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
    system.addForce(force2)

    # Protein-protein cation-pi force
    if any(atom.residue.name in ['ARG', 'LYS'] for atom in top.atoms()) and any(atom.residue.name in ['PHE','TRP','TYR'] for atom in top.atoms()):
        force3 = CustomNonbondedForce(equation2)
        force3.addPerParticleParameter("sigma")
        force3.addPerParticleParameter("epsilon")
        force3.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        force3.setCutoffDistance(R_OFF*nanometer)
        if SWITCHING == True:
            force3.addGlobalParameter("r_on", R_ON)
            force3.addGlobalParameter("r_off", R_OFF)
        for atom in top.atoms():
            force3.addParticle([ff_param[atom.residue.name]['radius']*2*2**(-1/6)*nanometer, cationpi_propro*kilojoule/mole])
        force3.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
        arg_lys  = [atom.index for atom in top.atoms() if atom.residue.name in ['ARG', 'LYS']]
        aromatic = [atom.index for atom in top.atoms() if atom.residue.name in ['PHE','TRP','TYR']]
        force3.addInteractionGroup(arg_lys, aromatic)
        system.addForce(force3)

    # Protein-RNA cation-pi force
    if any(atom.residue.name in ['ARG', 'LYS'] for atom in top.atoms()) and any(atom.residue.name in ['ADE','CYT','GUA','URA'] for atom in top.atoms()):
        force4 = CustomNonbondedForce(equation2)
        force4.addPerParticleParameter("sigma")
        force4.addPerParticleParameter("epsilon")
        force4.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        force4.setCutoffDistance(R_OFF*nanometer)
        if SWITCHING == True:
            force4.addGlobalParameter("r_on", R_ON)
            force4.addGlobalParameter("r_off", R_OFF)
        for atom in top.atoms():
            force4.addParticle([ff_param[atom.residue.name]['radius']*2*2**(-1/6)*nanometer, cationpi_prorna*kilojoule/mole])
        force4.createExclusionsFromBonds([(i[0].index, i[1].index) for i in top.bonds()], 1)
        arg_lys  = [atom.index for atom in top.atoms() if atom.residue.name in ['ARG', 'LYS']]
        nucleic = [atom.index for atom in top.atoms() if atom.residue.name in ['ADE','CYT','GUA','URA']]
        force4.addInteractionGroup(arg_lys, nucleic)
        system.addForce(force4)

    # Remove center-of-mass motion
    system.addForce(openmm.CMMotionRemover())

    # ----------------------------
    # Simulation setup
    # ----------------------------
    integrator = LangevinIntegrator(temp*kelvin, fric/picosecond, dt*picoseconds)
    if press_couple:
        if p_iso:
            logger.info("Isotropic pressure coupling enabled")
            system.addForce(MonteCarloBarostat(press*bar, temp*kelvin))
        elif p_aniso_x:
            logger.info("Anisotropic pressure coupling enabled along x-axis")
            system.addForce(MonteCarloAnisotropicBarostat((press, press, press)*bar, temp*kelvin, True, False, False))
        elif p_aniso_y:
            logger.info("Anisotropic pressure coupling enabled along y-axis")
            system.addForce(MonteCarloAnisotropicBarostat((press, press, press)*bar, temp*kelvin, False, True, False))
        elif p_aniso_z:
            logger.info("Anisotropic pressure coupling enabled along z-axis")
            system.addForce(MonteCarloAnisotropicBarostat((press, press, press)*bar, temp*kelvin, False, False, True))
        else:
            logger.error("Error in pressure coupling parameter!")
            sys.exit(1)

    platform = Platform.getPlatformByName(resources)
    if resources == 'CUDA':
        prop = dict(CudaPrecision='mixed')
        simulation = Simulation(top, system, integrator, platform, prop)
    else:
        simulation = Simulation(top, system, integrator, platform)

    chk_freq = params.get("chk_interval")
    chk_template = f"{base_name}_{{}}.chk"
    checkpoint_loaded = False

    try:
        chk_files = glob.glob(os.path.join('chkfiles', f"{base_name}_*.chk"))
        #chk_files = glob.glob(os.path.join('chkfiles', f"*.chk"))
        if chk_files:
            latest_chk = max(chk_files, key=lambda x: int(x.split('_')[1].split('.')[0]))
            logger.info("Loading checkpoint from %s", latest_chk)
            with open(latest_chk, 'rb') as f:
                simulation.context.loadCheckpoint(f.read())
            checkpoint_loaded = True
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        logger.info("No checkpoint found; starting new simulation.")
        simulation.context.setPositions(positions)

    if not checkpoint_loaded:
        logger.info("Initial energy: %s", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        logger.info("Starting energy minimization (%s steps)...", mini_nstep)
        simulation.minimizeEnergy(tolerance=100.0*kilojoules_per_mole/nanometer, maxIterations=mini_nstep)
        logger.info("Post-minimization energy: %s", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        simulation.context.setVelocitiesToTemperature(temp*kelvin, 870516298)

    simulation.reporters.append(DCDReporter(dcd_out, nstdcd))
    simulation.reporters.append(StateDataReporter(log_out, nstout, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, volume=True,
                                                    density=True, progress=True, remainingTime=True, speed=True,
                                                    totalSteps=nstep, separator='\t'))
    simulation.reporters.append(CustomCheckpointReporter(chk_template, chk_freq))

    logger.info("Running simulation for %s steps...", nstep)
    remaining_steps = nstep - simulation.currentStep
    simulation.step(remaining_steps)
    logger.info("Final potential energy: %s", simulation.context.getState(getEnergy=True).getPotentialEnergy())

    #final_step = simulation.currentStep
    #final_chk = chk_template.format(final_step)
    #with open(final_chk, 'wb') as f:
    #    f.write(simulation.context.createCheckpoint())
    #logger.info("Final checkpoint saved to %s", final_chk)

# ----------------------------
# Main CLI entry point
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Run COCOCMO MD simulation",
        epilog="Copyright reserved by Saumyak Mukherjee"
    )
    parser.add_argument('-prm', '--parameters', type=str, default='params.dat',
                        help='Input parameter file (str) [default: params.dat]', metavar='')
    args = parser.parse_args()
    try:
        run_sim(args.parameters)
    except Exception as e:
        logger.error("Simulation failed: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')
    main()
