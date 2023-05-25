"""
    Simulation_utils
    --------------------------------------------------------------------------------

    Utils for running CALVADOS simulations.

    --------------------------------------------------------------------------------
"""


import os
import shutil
from datetime import datetime as dt
import pandas as pd
import numpy as np

import mdtraj as md
from simtk import openmm, unit
from simtk.openmm import app, XmlSerializer

from utils import log as logger
from residues import residues
from conditions import conditions


#························································································#
#································ C A L V A D O S ·······································#
#························································································#

# CALVADOS 2 model is default
AH_cutoff = 2 # nm
DH_cutoff = 4 # nm
residues['AH_lambda'] = residues[f'CALVADOS2']


#························································································#
#······························ S I M U L A T I O N ·····································#
#························································································#

def openmm_simulate(dir: str, steps: int, top_path: str=None, sequence: str=None, boxlength: float=200, eqsteps: int=1000,cond: str='default', platform=None, stride: int=3000, verbose=False, log=True, savechk=True) -> None:
    """
    
    Takes a sequence (single-chain) or a topology (n-chain) as well as simulation specifications,
    runs a a CALVADOS coarse-grained simulation.

    Results of simulations are saved in the specified directory.

    Simulation can be run under different conditions, with options found in `conditions`

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `dir`: `str`
            Directory for simulation files

        `steps`: `int`
            Number of steps to run the simulation for (10 fs steps)

        `top_path`: `str`
            Path to the topology to use as input for simulation;
            Overwrites input to `sequence`

        `sequence`: `str`
            A sequence to submit for single-chain simulation;
            Ignored if `top_path` is provided

        `boxlength`: `float`
            The side length of the simulation (cubic) box [nm];
            Only used for generating an arbitrary topology for if `sequence` is provided.

        `eqsteps`: `int`
            Number of steps to subtract from trajectory as 'equilibration steps'

        `cond`: `str`
            The standard conditions to run the simulation with; see `conditions` for choices.

        `platform`: `str`
            Platform name to use for `openmm.Platform.getPlatformByName()`; 
            Specifies use of CPU or GPU (`CUDA`)

        `stride`: `int`
            The frame sampling frequency; 
            Number of steps between snapshots

        `verbose`: `bool`
            Whether to print log messages to stdout

        `log`: `bool`
            Whether to write log messages to `<dir>/simulate.log`

        `savechk`: `bool`
            Whether to save a checkpoint after the simulation ends

    """
    
    # Setting up logging
    log = logger(write=log, print=verbose, file=f'{dir}/simulate.log')
    log.message(f"[{dt.now()}] SIMULATION '{dir}' ")
    log.message(f"[{dt.now()}] Sequence: {sequence}")

    # Getting conditions and residue data
    log.message(f"[{dt.now()}] Preparing simulation with '{cond}' conditions")
    condition = conditions.loc[cond]

    # Calculating histidine charge based on Henderson-Hasselbalch equation
    H_pKa = 6
    residues.loc['H','q'] = 1. / (1 + 10**(condition.pH - H_pKa))

    # Generating arbitrary topology from sequence if topology not provided
    if top_path:
        shutil.copyfile(top_path, f'{dir}/top.pdb')
    elif sequence:
        top_path = f'{dir}/top.pdb'
        generate_save_topology(sequence, boxlength, top_path)
    else:
        # Throwing error if neither topology or sequence is provided
        raise TypeError("Either `sequence` or `top_path` must be provided to simulate_openmm_multiple()!")

    # Loading topology
    top = app.pdbfile.PDBFile(top_path)
    unitcell_vector = top.getTopology().getUnitCellDimensions()

    # Retrieving residue parameters for all chains in topology
    seqs = extract_sequences(md.load(top_path).topology)

    # Initiating OpenMM system
    system = openmm.System()

    # Defining simulation box from topology
    x = unit.Quantity(np.zeros([3]), unit.nanometers)
    x[0] = unitcell_vector.x * unit.nanometers
    y = unit.Quantity(np.zeros([3]), unit.nanometers)
    y[1] = unitcell_vector.y * unit.nanometers
    z = unit.Quantity(np.zeros([3]), unit.nanometers)
    z[2] = unitcell_vector.z * unit.nanometers
    system.setDefaultPeriodicBoxVectors(x, y, z)
    
    # Initialising energy objects
    r_0=0.38
    k=8033
    hb = openmm_harmonic_bond(r_0, k)
    ah = openmm_ashbaugh_hatch(epsilon_factor=condition.eps_factor)
    dh = openmm_debye_huckel(T=condition.temp, c=condition.ionic)
    
    # Adding residues to system and setting particle-specific parameters
    for _, res in seqs.iterrows():
        system.addParticle(res.MW*unit.amu)
        ah.addParticle([res.AH_sigma*unit.nanometer, res.AH_lambda*unit.dimensionless])
        dh.addParticle([res.q])

    # Setting residue interactions for bonded residues
    for i in seqs.index[:-1]:
        if seqs.loc[i, 'chain'] == seqs.loc[i+1, 'chain']:
            hb.addBond(i, i+1, r_0*unit.nanometer, k*unit.kilojoules_per_mole/(unit.nanometer**2))
            ah.addExclusion(i, i+1)
            dh.addExclusion(i, i+1)
        
    # Adding energy objects
    for energy_term in [hb, ah, dh]:
        system.addForce(energy_term)

    # Serialising system
    serialized_system = XmlSerializer.serialize(system)
    with open(f'{dir}/system.xml','w') as file:
        file.write(serialized_system)
    
    # Setting integrator and platform
    friction = 0.01
    stepsize = 0.010 # 10 fs timestep
    integrator = openmm.LangevinIntegrator(condition.temp*unit.kelvin, friction/unit.picosecond, stepsize*unit.picosecond) 
    platform = openmm.Platform.getPlatformByName(platform)

    # Initiating simulation object
    simulation = app.simulation.Simulation(top.topology, system, integrator, platform) #, dict(CudaPrecision='mixed')) 

    # Checking for checkpoint to start from
    check_point = f'{dir}/restart.chk'
    if os.path.isfile(check_point):
        log.message(f"[{dt.now()}] Reading from check point file '{check_point}'")
        simulation.loadCheckpoint(check_point)
        simulation.reporters.append(app.dcdreporter.DCDReporter(f'{dir}/traj.dcd', stride, append=True))
        eqsteps = 0
    
    # Else start from scratch
    else:
        log.message(f"[{dt.now()}] Starting from scratch")
        simulation.context.setPositions(top.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(app.dcdreporter.DCDReporter(f'{dir}/pretraj.dcd', stride))

    # Setting up OpenMM log
    simulation.reporters.append(app.statedatareporter.StateDataReporter(
        f'{dir}/traj.log',
        int(stride),
        potentialEnergy=True,
        temperature=True,
        step=True,
        speed=True,
        elapsedTime=True,
        separator='\t'))

    # Running simulation
    log.message(f"[{dt.now()}] Running simulation of {steps * stepsize / 1000} ns ({steps} steps)")
    simulation.step(steps)

    # Saving final checkpoint
    if savechk:
        log.message(f"[{dt.now()}] Saving check point in '{check_point}'")
        simulation.saveCheckpoint(check_point)

    # Generating trajectory without equilibration
    log.message(f"[{dt.now()}] Saving formatted trajectory in '{dir}/traj.dcd'")
    save_dcd(traj_path=f'{dir}/pretraj.dcd', top_path=f'{dir}/top.pdb', file_path=f'{dir}/traj.dcd', eqsteps=eqsteps)
    os.remove(f'{dir}/pretraj.dcd')


#························································································#
#·································· T O P O O G Y ·······································#
#························································································#

def extract_sequences(top: md.Topology, res: pd.DataFrame=residues.copy()) -> pd.DataFrame:
    """
    
    Takes a MDTraj topology, generates a DataFrame containing CALVADOS parameters for each residue in each chain of the topology.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `top`: `md.Topology`
            A MDTraj topology
    
        `res`: `pandas.DataFrame`
            A `residues` template DataFrame

    Returns
    -------

        `seqs`: `pandas.DataFrame`
            A DataFrame containing residue-specific parameters for each residue in topology;
            See Schema below

    --------------------------------------------------------------------------------

    Schema of DataFrame
    -------------------
    Each row represents a residue in the topology

        `[index]`: `int`
            An integer representing the order of the residue in the topology
        
        `chain`: `int`
            The chain ID of the residue

        `res`: `int`
            The residue ID of the residue in the chain
        
        `aa`: `str`
            The one letter code of the residue amino acid

        `MW`, `AH_lambda`, `AH_sigma`, `q`
            See `residues`

    """

    # Preparing dataframe for storing sequences
    seqs = pd.DataFrame(columns=['chain', 'res', 'aa', 'MW', 'AH_lambda', 'AH_sigma', 'q'])

    # Looping over chains
    n = 0
    res = res.set_index('one')
    for chain in top.chains:

        # Looping over residues
        for i, aa in enumerate(chain.residues):
            seqs.loc[n] = {'chain': chain.index, 'res': aa.index, 'aa': aa.name} | res.loc[aa.name].to_dict()
            n += 1

        # Modifying terminal residues of chain
        nter = seqs[seqs.chain == chain.index].iloc[0].name
        cter = seqs[seqs.chain == chain.index].iloc[-1].name
        seqs.loc[nter, 'q'] += 1
        seqs.loc[cter, 'q'] -= 1
        seqs.loc[nter, 'MW'] += 2
        seqs.loc[cter, 'MW'] += 16

    return seqs

#························································································#
def generate_save_topology(seq, boxlength: float, file_path: str) -> None:
    """
    
    Takes a sequence, generates a MDTraj topology of the sequence as a chain and saves it to a .pfb-file.

    The topology will be of the chain as a straight line with X=0, Y=0, and spanning the Z-dimension where it is centered.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str|list`
            An amino acid sequence

        `boxlength`: `float`
            The length of the simulation box sides; Used to center the sequence.
            Restraints are placed to prevent the chain from surpassing half the `boxlength` in length.

        `file_path`: `str`
            The path to save the file as (include `.pdb` suffix).

    """
    
    # Generating initial residue coordinates; straight chain centered on the z-axis
    N = len(seq)
    rise_per_residue = 0.38
    assert N * rise_per_residue < boxlength/2, f"Box dimensions is too small compared to sequence! ({boxlength} nm)"
    coordinates = [[0, 0, boxlength/2 + rise_per_residue*(i - N/2)] for i in range(N)]

    # Generating topology
    top = md.Topology()
    chain = top.add_chain()
    for aa in seq:
        res = top.add_residue(aa, chain)
        top.add_atom(aa, element=md.element.carbon, residue=res)
    for i in range(chain.n_atoms-1):
        top.add_bond(chain.atom(i),chain.atom(i+1))

    # Saving topology
    traj = md.Trajectory(xyz=np.array(coordinates).reshape(N, 3), topology=top, time=0, unitcell_lengths=[boxlength]*3, unitcell_angles=[90]*3)
    traj.save_pdb(file_path)
    fix_pdb_conects(file_path)

#························································································#
def merge_topologies(trajs: list, boxlength: float=None) -> md.Trajectory:
    """

    Takes a list of topologies in the format of single-frame MDTraj trajectories,
    merges the topologies into one.

    The unitcell will match that of the first topology, if no new box length is specified.

    Make sure input topologies doesn't overlap; Consider translating structures by
    modifying the md.Trajectory.xyz attribute of the topology trajectory.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `tops`: `list[md.Trajectory]`
            MDTraj single-frame trajectories of the topologies to combine;
            (I.e. load using `md.load("xxx.pdb")`)

        `boxlength`: `float`
            An argument to overwride the box size

    Returns
    -------

        `merged_top`: `md.Trajectory`
            A single-frame trajectory containing the merged topologies

    """

    # Merging topology objects
    merged_traj = trajs[0]
    for traj in trajs[1:]:
        merged_traj = merged_traj.stack(traj)

    # Setting new box dimensions
    if boxlength:
        merged_traj.unitcell_lengths = [[boxlength, boxlength, boxlength]]

    # Centering coordinates
    merged_traj.center_coordinates()
    merged_traj.xyz += merged_traj.unitcell_lengths[0,0]/2
    
    return merged_traj

#························································································#
def fix_pdb_conects(pdb_path: str) -> None:
    """
    
    Takes the path to a PDB-file,
    overwrites the file such that the CONECT lines matches the specified chains in the model.

    PDB-files created by MDTraj.Trajectory.save("xxx.pdb") appears to not create correct CONECT statements,
    which disturbs the function of topology-based methods.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `pdb_path`: `str`
            The path to the .pdb file to proces

    """

    # Loading file for reading
    assert '.pdb' in pdb_path
    with open(pdb_path, 'r') as file:
        lines = file.readlines()

        # Assembling chains
        chains = []
        seq = []
        for line in lines:
            if "ATOM" in line:
                seq += [line.split()[1]]
            elif "TER" in line:
                seq += [line.split()[1]]
                chains.append(seq)
                seq = []
        chains = [sorted([int(atom) for atom in chain]) for chain in chains]

        # Assembling new CONECTs
        conects = []
        for chain in chains:
            for i, atom in enumerate(chain):
                if i == 0:
                    conect = [atom, atom+1]
                elif i == len(chain) - 1:
                    conect = [atom, atom-1]
                else:
                    conect = [atom, atom-1, atom+1]
                conects.append(conect)

        # Formatting new CONECT lines
        conect_lines = []
        for conect in conects:
            conect_line = 'CONECT' + ''.join([' '*(5-len(str(i))) + str(i) for i in conect])
            conect_lines.append(conect_line + '\n')

    # Overwriting  file
    with open(pdb_path, 'w') as file:
        while 'CONECT' not in lines[0]:
            file.write(lines.pop(0))
        file.writelines(conect_lines)
        file.write("END")


#························································································#
#···································· O P E N M M ·······································#
#························································································#

def openmm_harmonic_bond(r_0=0.38, k=8033) -> openmm.HarmonicBondForce:
    """
    
    Sets up a harmonic bond restraint energy term,
    returns the corresponding OpenMM HarmonicBondForce object (which needs particle-specific interactions set).

    The object has periodic boundary conditions enabled.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `r_0`: `float`
            The resting distance in the model [nm]

        `k`: `float`
            The energy penalty term in the model [kJ/(mol·nm^2)]

    Returns
    -------

        `hb`: `openmm.HarmonicBondForce`
            An OpenMM HarmonicBondForce object

    """

    # Initiating model object
    hb = openmm.HarmonicBondForce()

    # Enable periodic boundary conditions
    hb.setUsesPeriodicBoundaryConditions(True)

    return hb

#························································································#
def openmm_ashbaugh_hatch(epsilon_factor: float, r_cutoff=AH_cutoff) -> openmm.CustomNonbondedForce:
    """
    
    Sets up a Ashbaugh-Hatch energy term,
    returns the corresponding OpenMM CustomNonbondedForce object (which needs particle-specific parameters/interactions set).

    The object has periodic boundary conditions enabled.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `epsilon_factor`: `float`
            For calculating Lennard-Jones epsilon parameter (see `ah_parameters`)

        `r_cutoff`: `float`
            The r_cutoff distance for the interaction [nm]


    Returns
    -------

        `ah`: `openmm.CustomNonbondedForce`
            An OpenMM CustomNonbondedForce object of the energy term

    """

    # Initiating model object
    energy_expression = 'select(step(r-2^(1/6)*s),4*epsilon*l*((s/r)^12-(s/r)^6-shift),4*epsilon*((s/r)^12-(s/r)^6-l*shift)+epsilon*(1-l))'
    parameter_expression = 's=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/ah_cutoff)^12-(0.5*(s1+s2)/ah_cutoff)^6'
    ah = openmm.CustomNonbondedForce(energy_expression+';'+parameter_expression)

    # Calculating Lennard-Jones epsilon parameter
    epsilon = ah_parameters(epsilon_factor)

    # Adding global parameters
    ah.addGlobalParameter('epsilon', epsilon*unit.kilojoules_per_mole)
    ah.addGlobalParameter('ah_cutoff', r_cutoff*unit.nanometer)

    # Specifying position specific parameters sigma and lambda
    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')
    
    # Set periodic boundary conditions with specified r_cutoff
    ah.setForceGroup(0)
    ah.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setCutoffDistance(r_cutoff*unit.nanometer)

    return ah

#························································································#
def openmm_debye_huckel(T: float, c: float, r_cutoff=DH_cutoff) -> openmm.CustomNonbondedForce:
    """
    
    Sets up a Debye-Hückel energy term,
    returns the corresponding OpenMM CustomNonbondedForce object (which needs particle-specific parameters/interactions set).

    The object has periodic boundary conditions enabled.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `T`: `float`
            The absolute temperature [°K]

        `c`: `float`
            The ionic strength [mM]

        `r_cutoff`: `float`
            The r_cutoff distance for the interaction [nm]


    Returns
    -------

        `dh`: `openmm.CustomNonbondedForce`
            An OpenMM CustomNonbondedForce object of the energy term

    """

    # Initiating model object
    energy_expression = 'q*yukawa_epsilon*(exp(-yukawa_kappa*r)/r-shift)'
    parameter_expression = 'q=q1*q2; shift=exp(-yukawa_kappa*dh_cutoff)/dh_cutoff'
    dh = openmm.CustomNonbondedForce(energy_expression+';'+parameter_expression)

    # Calculating the inverse Debye length kappa and
    # the Debye-Hückel potential coefficients epsilon
    # (Also known as Yukawa kappa/epsilon)
    yukawa_kappa, yukawa_epsilon = dh_parameters(T, c)

    # Adding global parameters
    dh.addGlobalParameter('yukawa_kappa', yukawa_kappa/unit.nanometer)
    dh.addGlobalParameter('yukawa_epsilon', yukawa_epsilon*unit.nanometer*unit.kilojoules_per_mole)
    dh.addGlobalParameter('dh_cutoff', r_cutoff*unit.nanometer)

    # Specifying position specific parameter q
    dh.addPerParticleParameter('q')
    
    # Set periodic boundary conditions with specified r_cutoff
    dh.setForceGroup(1)
    dh.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    dh.setCutoffDistance(r_cutoff*unit.nanometer)

    return dh


#························································································#
#······························ P A R A M E T E R S ·····································#
#························································································#

def ah_parameters(epsilon_factor: float) -> float:
    """
    
    Calculates the Lennard-Jones potential epsilon parameter.
    Used to calculate the Ashbaugh-Hatch potential.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `epsilon_factor`: `float`
            TODO The scaling factor for the ??? epsilon value []

    Returns
    -------

        `lj_epsilon`: `float`
            The Lennard-Jones epsilon parameter [kJ/mol]

    """
    # TODO Where does this value 4.184 come from?
    lj_epsilon = 4.184 * epsilon_factor

    return lj_epsilon

#························································································#
def dh_parameters(T: float, c: float) -> tuple:
    """
    
    Calculates the Yukawa epsilon and kappa parameters.
    Used to calculate the Debye-Hückel potential.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `T`: `float`
            Absolute temperature [°K]

        `c`: `float`
            Ionic strength of the solution [M]/[mol/L]

    Returns
    -------

        `kappa`: `float`
            The inverse Debye-Hückel length [1/nm]; 
            used for the scaling of the exponential term in the Debye-Hückel equation

        `epsilon`: `float`
            The coefficient in the Debye-Hückel equation [nm·kJ/mol]

    """

    # Setting constants
    e = 1.6021766e-19       # C             | Elementary charge
    R = 8.3145              # J/(mol·°K)    | Ideal gas constant
    N_A = 6.022e+23         # 1/mol         | Avogadro's constant
    eps_0 = 8.854188e-12    # F/m           | Vacuum permittivity
    pi = np.pi
    eps_r = 5321*(T**-1) + 233.76 - 0.9297*(T) + 0.1417*1e-2*(T**2) - 0.8292*1e-6*(T**3) # [unitless] | Emperical scalar

    # Calculating Bjerrum length
    B = N_A * (e**2) / (4*pi*eps_0*eps_r*R*T) # m
    B *= 1e9 # nm (from m)

    # Calculating inverse Debye-Hückel length
    c *= 1e-24 # mol/nm^3 (from mol/dm^3)
    yukawa_kappa = float(np.sqrt(8*pi*B*N_A*c)) # nm 

    # Calculating the coefficient of the Debye–Hückel equation
    yukawa_epsilon = N_A * (e**2) / (4*np.pi*eps_0*eps_r) # m·J/mol
    yukawa_epsilon *= 10**6 # nm·kJ/mol

    return yukawa_kappa, yukawa_epsilon


#························································································#
#·························· P O S T P R O C E S S I N G ·································#
#························································································#

def save_dcd(traj_path: str, top_path: str, file_path: str, eqsteps: int=0) -> None:
    """
    
    Generates coordinate and trajectory in convenient formats.
    Recenters trajectory frames and removes initial equilibritation steps.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj_path`: `str`
            Path to raw trajectory .dcd file

        `top_path`: `str`
            Path to topology .pdb file

        `file_path`: `str`
            Path to save new processed trajectory .dcd file to

    """

    # Loading trajectory
    traj = md.load(traj_path, top=top_path)

    # Applying periodic boundary conditions to the molecules in each frame of the trajectory
    traj = traj.image_molecules(anchor_molecules=[set(traj.top.chain(0).atoms)], make_whole=True)

    # Centering in box
    traj.center_coordinates()
    traj.xyz += traj.unitcell_lengths[0,0]/2

    # Filtering out equilibration from final trajectory
    tocut = eqsteps
    traj[int(tocut):].save_dcd(file_path)
