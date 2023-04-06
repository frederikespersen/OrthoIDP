"""
    Evolve
    --------------------------------------------------------------------------------

    Script for running Francesco Pesce's IDP evolution algorithm.
    Takes shell arguments, orchestrates preparations, and initiates simulation.

    --------------------------------------------------------------------------------

    Run `python evolve.py -h` for documentation on arguments.

"""


import argparse
import os
import sys
import random
import pandas as pd
import shutil
import numpy as np
import mdtraj as md
from time import time


#························································································#

# Setting up argument parser
# Monte Carlo parameters
parser = argparse.ArgumentParser(prog="Evolve", description="Runs an IDP shuffling-evolution algorithm for target observable value")
parser.add_argument('-r', '--restart',
                    type=int,
                    required=False,
                    default=None,
                    help="which generation of a previous run to restart from (default: 'None') ['None': start new evolution; '-1' from the last]")
parser.add_argument('-d', '--dir',
                    type=str,
                    required=False,
                    default='evolution',
                    help="the name of the directory to put simulation results in within results/")
parser.add_argument('-f', '--fasta',
                    type=str,
                    required=True,
                    help="path to FASTA file with sequence to simulate")
parser.add_argument('-m', '--measure',
                    type=str,
                    required=False,
                    default='Rg',
                    help="what structural obervable to use as evolutionary constraint. Options: ['Rg': Radius of gyration; 'nu': Scaling exponent]")
parser.add_argument('-t', '--target_value',
                    type=float,
                    required=True,
                    help="the target value of the structural observable indicated in --measure")
parser.add_argument('-L', '--L_at_half_acceptance',
                    type=float,
                    required=True,
                    help="???")
parser.add_argument('-a', '--simulated_annealing',
                    action='store_true',
                    required=False,
                    help="Whether to simulate annealing by progressively lowering the Monte Carlo control parameter")
parser.add_argument('-s', '--source',
                    type=str,
                    required=False,
                    default=None,
                    help="path source code for simulate utils (default: Same dir as evolve.py)")
# Simulation parameters
parser.add_argument('-b', '--boxlength',
                    type=float,
                    required=False,
                    default=200,
                    help="boxlength [nm] for simulation cube (default_ 200 nm)")
parser.add_argument('-c', '--conditions',
                    type=str,
                    required=False,
                    default='default',
                    help="physical conditions to simulate under (default: 'default')")
parser.add_argument('-n', '--steps',
                    type=int,
                    required=False,
                    default=100000000,
                    help="the number of timesteps [10 fs] to run the simulation for (default: 100,000,000 / 1 µs)")
parser.add_argument('-p', '--platform',
                    type=str,
                    required=False,
                    default='CUDA',
                    help="computional platform to use (default: CUDA)")

# Parsing arguments
args = parser.parse_args()

restart = args.restart
dir = args.dir
measure = args.measure
L_at_half_acceptance = args.L_at_half_acceptance
target = args.target_value
source_path = args.source
cond = args.conditions


#························································································#
#································ P R E P A R A T I O N ·································#
#························································································#

# Timing
t0 = time()

# Setting source path as the dir of evolve.py (this file) if nothing else is specified
if source_path is None:
    source_path = '/'.join(__file__.split('/')[:-1])
sys.path.append(source_path)

# Importing source code
import evolve_utils
from utils import log as logger
from utils import read_fasta
from simulate_utils import openmm_simulate

#························································································#

# Extracting fasta file
seq = read_fasta(args.fasta, just_seq=True)

#························································································#

# Setting results directory
dir = f'results/{dir}'

# Creating directory for evolution
os.makedirs(dir, exist_ok=True)
os.chdir(dir)

# Setting logging
log = logger(write=True, print=False, file='evolution.log', timestamp=True)

#························································································#

# Standard settings
store_filename = 'evolution.pkl'
max_pool_size = 10

#························································································#

# Setting measure for observable
compute_observable = evolve_utils.structural_measures[measure]

#························································································#

log.message("")

# Starting from scratch (Generation 0)
if restart is None:
    log.message(f"STARTING FROM GENERATION 0")
    log.message(f"Input sequence: {seq}")
    g = 0

    # Calculating initial Monte Carlo control parameter
    c = -1 * args.L_at_half_acceptance / np.log(0.5)

    # Running intial simulation
    if not os.path.isdir(f'g{g}'):
        log.message(f"Simulating input sequence")
        os.makedirs(f'g{g}')
        openmm_simulate(sequence=seq, dir=f'g{g}', boxlength=args.boxlength, cond=args.conditions, steps=args.steps, platform=args.platform)
    else:
        log.message(f"Input simulation found")

    # Calculating observable
    traj = md.load_dcd(f'g{g}/traj.dcd', f'g{g}/top.pdb')
    obs = compute_observable(seq, traj)
    log.message(f"Calculated observable '{measure}' to {obs:.4f}")
    
    # Initialising DataFrame for storing results
    store = pd.DataFrame(columns=['sequence','observable','simulate','mc','c'])
    store.loc[g] = {'sequence': seq,
                    'observable': obs,
                    'simulate': True,
                    'mc': True,
                    'c': c}

    # Pickling initial generation
    store.to_pickle(store_filename)

    # Calculating energy of each frame for reweighting
    E = evolve_utils.compute_total_energy(seq, traj)

    # Saving trajectory in temporary Series pool
    pool = pd.DataFrame({'traj':traj}, index = [g])

# Restarting
else:

    # Unpickling previous evolution
    store = pd.read_pickle(store_filename)

    # Restarting from latest generation
    if restart == -1:
        restart = store.index[-1]

    log.message(f"RESTARTING FROM GENERATION {restart}")
    
    # Removing generations after restart if exists
    if store[store.index > restart].simulate.sum() > 0:
        log.message(f"Removing simulation data from generations {', '.join([g for g in store.index if g > restart])}")
        for g in store.index:
            if g > restart:
                if store.simulate[g]:
                    shutil.rmtree('g'+str(g))
                store = store.drop(g)

    # Loading trajectories and saving them in temporary DataFrame pool
    log.message("Loading previous trajectories for reweighting")
    pool = pd.DataFrame(columns=['traj'])
    for g in store[store.simulate].index[-max_pool_size:]:
        traj = md.load_dcd(f'g{g}/traj.dcd', f'g{g}/top.pdb')
        pool.loc[g] = {'traj': traj}

    # Using latest control parameter value
    c = store.c.iloc[-1]

log.message(f"Targeting a '{measure}' of {target}")
log.message(f"Starting evolution with Monte Carlo control parameter of {c:.6f}")

#························································································#

# Calculating energies for all sequence in the pool using all trajectories
# Creating dataframe with rows of used sequences and columns of used trajectories

# Looping over trajectories
for g in pool.index:
    traj = pool.traj[g]

    # Looping over sequences
    Etots = []
    for seq in store.sequence[pool.index]:

        # Computing energy for each combination
        Etots.append(evolve_utils.compute_total_energy(seq, traj, cond))

    # Saving results in DataFrame
    pool[f'{g}'] = Etots

#························································································#

# Starting after restart, else after generation 0
if restart is None:
    start = 0
else:
    start = restart
start += 1

# Timing
t = time() - t0
log.message(f"Preparations finished in {t/60:.2f} minutes")


#························································································#
#·································· E V O L U T I O N ···································#
#·································· A L G O R I T H M ···································#
#························································································#

# Looping over generations
for g in range(start, 100000):

    # Timing
    t0 = time()

    #···························· I N I T I A L I S I N G ·······························#

    # Generating new sequence and initialising new entry 
    log.message("")
    log.message(f"INITIALIZING GENERATION {g}")
    last_seq = store[store.mc].sequence.iloc[-1]
    seq = evolve_utils.swap_sequence(last_seq)
    log.message(f"Generated sequence: {seq}")

    # Checking whether the sequence has been previously generated
    identical_seqs = store[seq == store.sequence]

    #···················· C O M P U T I N G   O B S E R V A B L E ·······················#

    # Checking whether the sequence has been previously simulated
    if sum(identical_seqs.simulate) > 0:
        
        # Copying old entry
        log.message(f"Copying previously simulated identical entry in generation {identical_seqs.index[0]}")
        previous_entry = identical_seqs[identical_seqs.simulate].iloc[0]
        store.loc[g] = {'sequence': seq,
                        'observable': previous_entry.observable,
                        'simulate': False,
                        'mc': False,
                        'c': c}
    
    # Else, attempting to calculate observable by reweighting / simulation
    else:
        log.message("No previous simulation found for sequence")
        log.message("Evaluating whether reweighting is possible")
        weights, eff_frames = evolve_utils.compute_MBAR(seq, pool, cond)
        
        # Check whether there is sufficient effective frames for reweighting
        if eff_frames >= 20000:
            log.message(f"Reweighting (Effective frames: {eff_frames:.0f})")

            # Calculating observable of the sequence using previous trajectories
            obs_pool = [compute_observable(seq, traj) for traj in pool.traj]
            obs = np.average(obs_pool, weights=weights)

            # Filling out entry
            store.loc[g] = {'sequence': seq,
                            'observable': obs,
                            'simulate': False,
                            'mc': False,
                            'c': c}

        # Else simulate as last resort
        else:
            log.message(f"Not enough frames for reweighting (Effective frames: {eff_frames})")
            log.message("Simulating sequence")
            os.makedirs(f'g{g}')
            openmm_simulate(sequence=seq, dir=f'g{g}', boxlength=args.boxlength, cond=args.conditions, steps=args.steps, platform=args.platform)

            # Calculating observable
            traj = md.load_dcd(f'g{g}/traj.dcd', f'g{g}/top.pdb')
            obs = compute_observable(seq, traj)
            
            # Filling out entry
            store.loc[g] = {'sequence': seq,
                            'observable': obs,
                            'simulate': True,
                            'mc': False,
                            'c': c}

            # Updating pool with new simulation
            if len(pool) == max_pool_size:
                earlist_g = pool.index[0]
                pool.drop(index=earlist_g, columns=str(earlist_g), inplace=True)
            pool.loc[g] = [traj] + [evolve_utils.compute_total_energy(seq, t, cond) for t in pool.traj]
            pool[f'{g}'] = [evolve_utils.compute_total_energy(s, traj, cond) for s in store.sequence[pool.index]]

    log.message(f"Calculated observable '{measure}' to {store.observable[g]:.4f}")

    #··················· M O N T E   C A R L O   E V O L U T I O N ······················#

    # Calculating cost function
    last_g = store[store.mc].index[-1]
    L = abs(store.observable[g] - target) - abs(store.observable[last_g] - target)

    # Calculating probability of acceptance
    p = min(1, np.exp(-L/c))

    # Monte Carlo sampling
    if random.random() <= p:

        # Accepting move
        store.mc[g] = True
        log.message(f"Accepted evolved generation with a probability of {p:.4%}")
    
    else:
        # Rejecting move
        store.mc[g] = False
        log.message(f"Rejected evolved generation with a probability of {p:.4%}")

    log.message(f"MC acceptance rate is currently at {store.mc.sum()/len(store):.2%}%")

    #································ P I C K L I N G ···································#

    # Pickling current evolution
    store.to_pickle(store_filename)

    # Timing
    t = time() - t0
    if t < 3600:
        log.message(f"Generation finished in {t/60:.2f} minutes")
    else:
        log.message(f"Generation finished in {t/3600:.2f} hours")
        
    #······························· A N N E A L I N G ··································#

    # Checking for annealing
    if args.simulated_annealing:
        
        # Counting how many times the current control parameter value was used
        nc = sum(c == store.c)

        # Annealing if value count is large compared to sequence length
        if nc >= len(seq)*2:
            c *= 0.99
            log.message("")
            log.message("SIMULATING ANNEALING")
            log.message(f"Lowering control parameter by 1% to {c}")
    
