"""
    Simulation_utils
    --------------------------------------------------------------------------------

    Utils for running CALVADOS simulations.

    --------------------------------------------------------------------------------
"""


import os
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
def openmm_simulate(sequence: str, boxlength: float, dir: str, steps: int, eqsteps: int=1000,cond: str='default', platform=None, stride: int=3000, verbose=False, log=True, savechk=True) -> None:
    """
    
    Takes a sequence and simulation specifications,
    runs a single-chain CALVADOS coarse-grained simulation.

    Results of simulations are saved in the specified directory.

    Simulation can be run under different conditions, with options found in `conditions`

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `sequence`: `str`
            A sequence to submit for simulation

        `boxlength`: `float`
            The side length of the simulation (cubic) box [nm]

        `dir`: `str`
            Directory for simulation files

        `steps`: `int`
            Number of steps to run the simulation for (10 fs steps)

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
    
    log = logger(write=log, print=verbose, file=f'{dir}/simulate.log')
    log.message(f"[{dt.now()}] SIMULATION '{dir}' ")
    log.message(f"[{dt.now()}] Sequence: {sequence}")

    # Getting conditions and residue data
    log.message(f"[{dt.now()}] Preparing simulation with '{cond}' conditions")
    condition = conditions.loc[cond]

    # Formating terminal residues as special residue types
    sequence, residues = format_terminal_res(sequence)

    # Calculating histidine charge based on Henderson-Hasselbalch equation
    H_pKa = 6
    residues.loc['H','q'] = 1. / (1 + 10**(condition.pH - H_pKa))


    # Initiating OpenMM system
    system = openmm.System()

    # Defining simulation box 
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = boxlength * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = boxlength * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = boxlength * unit.nanometers
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Generating and loading topology
    top_path = f'{dir}/top.pdb'
    generate_save_topology(sequence, boxlength, top_path)
    top = app.pdbfile.PDBFile(top_path)

    # Adding molecular weights
    for mw in map(lambda aa: residues.MW[aa], sequence):
        system.addParticle(mw*unit.amu)
    
    # Initialising energy objects
    r_0=0.38
    k=8033
    hb = openmm_harmonic_bond(r_0, k)
    ah = openmm_ashbaugh_hatch(epsilon_factor=condition.eps_factor)
    dh = openmm_debye_huckel(T=condition.temp, c=conditions.ionic)
    
    # Setting particle-specific parameters
    for aa in sequence:
        aa = residues.loc[aa]
        ah.addParticle([aa.AH_sigma*unit.nanometer, aa.AH_lambda*unit.dimensionless])
        dh.addParticle([aa.q])

    # Setting residue interaction pairs
    N = len(sequence)
    for i in range(N-1):
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
    
    # Setting integrator and CPU/GPU
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
        simulation.reporters.append(app.dcdreporter.DCDReporter(f'{dir}/pretraj.dcd', stride, append=True))
        eqsteps = 0
    
    # Else start from scratch
    else:
        log.message(f"[{dt.now()}] Starting from scratch")
        simulation.context.setPositions(top.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(app.dcdreporter.DCDReporter(f'{dir}/pretraj.dcd', stride))

    # Setting up log
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
#····························· P R E P A R A T I O N ····································#
#························································································#

def format_terminal_res(seq, res: pd.DataFrame=residues.copy()):
    """
    
    Takes a sequence and a `residues` DataFrame, modifies the sequence with special terminal residue types 'X' and 'Z'
    for the N- and C-terminal respectively.
    Returns the modified sequence and the modified `residues` DataFrame.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str|list`
            An amino acid sequence

        `res`: `pandas.DataFrame`
            A `residues` DataFrame

    Returns
    -------

        `seq`: `str`
            The modified sequence with terminal 'X'/'Z' residues

        `res`: `pandas.DataFrame`
            A modified `residues` DataFrame with 'X'/'Z' residue types

    """

    # Gettning standard residue data
    res = res.set_index('one')

    # Adding new residue types, and using original terminal residues as templates
    res.loc['X'] = res.loc[seq[0]].copy()
    res.loc['Z'] = res.loc[seq[-1]].copy()
    res.loc['X','MW'] += 2
    res.loc['Z','MW'] += 16
    res.loc['X','q'] += 1
    res.loc['Z','q'] -= 1

    # Modfiying sequence
    seq = list(seq)
    seq[0] = 'X'
    seq[-1] = 'Z'
    seq = ''.join(seq)

    return seq, res


#························································································#
def generate_save_topology(seq, boxlength: float, file_path: str) -> None:
    """
    
    Takes a sequence, generates a MDTraj topology of the sequence as a chain and saves it to a .pfb-file.

    The topology will be of the chain as a straight line on the X-Y plane, centered in the Z-dimension.
    
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
    traj = md.Trajectory(xyz=np.array(coordinates).reshape(N, 3), topology=top, time=0, unitcell_lengths=[boxlength]*3, unitcell_angles=[90,90,90])
    traj.save_pdb(file_path)


#························································································#
#·································· M O D E L S ·········································#
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

def save_dcd(traj_path: str, top_path: str, file_path: str, eqsteps: int=1000) -> None:
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
