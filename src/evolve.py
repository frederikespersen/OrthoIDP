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
import pandas as pd
import shutil
import numpy as np
import mdtraj as md

import evolve_utils
import analyse_utils
from data_utils import read_fasta
from simulate_utils import simulate


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="Evolve", description="Runs an IDP shuffling-evolution algorithm for target observable value")
parser.add_argument('-f', '--fasta',
                    type=str,
                    required=True,
                    help="path to FASTA file with sequence to simulate")
parser.add_argument('-d', '--dir',
                    type=str,
                    required=False,
                    help="directory to store results in (default uses name of fasta file)")
parser.add_argument('-r', '--restart',
                    type=int,
                    required=False,
                    help="which step of a previous run to restart from (default: None - start new evolution)")
parser.add_argument("--L_at_half_acceptance", type=float)

# Parsing arguments
args = parser.parse_args()
fasta_path = args.fasta
dir = args.dir
restart = args.restart

#························································································#

# Extracting fasta file
init_sequence = read_fasta(fasta_path, just_seq=True)

#························································································#

# Using fasta file name as evolution directory if nothing else is specified
if dir is None:
    dir = fasta_path.split('/')[-1].replace('.fasta', '')

# Creating directory for evolution
os.makedirs(dir, exist_ok=True)

#························································································#

# Standard settings
store_filename = 'evolution.pkl'
boxlength = 200
steps = 51000000
cond = 'default'
compute_observable = analyse_utils.compute_rg
max_generations = 50000

#························································································#

# Starting from scratch (Generation 0)
if restart is None:
    evolve_utils.log(f'INITIALIZING INPUT SEQUENCE AND SIMULATE', header=True)
    g = 0
    
    # Calculating initial Monte Carlo control parameter
    c = -1 * args.L_at_half_acceptance / np.log(0.5)

    # Initialising DataFrame for storing results
    store = pd.DataFrame(columns=['sequence','observable','simulate','mc','c','distances','mask'])
    store.iloc[g] = dict(fasta=init_sequence, obs=None, simulate=True, mc=True, c=c, distances=None, mask=None)

    # Running intial simulation
    if os.path.isdir('g0'):
        evolve_utils.log('Start sequence simulation available, skipping simulation.')
    else:
        simulate(sequence=init_sequence, boxlength=boxlength, dir="g0", steps=steps, eqsteps=1000, cond=cond, vmodel=3, platform='CODA')

    # Calculating mean observable
    init_traj = md.load_dcd('g0/traj.dcd', 'g0/top.pdb')
    store.obs[g] = compute_observable(init_sequence, init_traj).mean()

    # Calculating distances
    distances, mask = analyse_utils.compute_distances('g0')
    store.distances[g] = distances
    store.mask[g] = mask

    # Pickling initial generation
    store.to_pickle(store_filename)

# Restarting
else: 
    evolve_utils.log(f'RESTARTING FROM GENERATION {restart}', header=True)

    # Loading previous evolution
    store = pd.read_pickle(store_filename)

    # Removing generations after restart
    for g in store.index:
        if g > restart:
            if store.simulate[g] == True:
                shutil.rmtree('g'+str(g))
            store = store.drop(g)

#························································································#

# Starting after restart, else after generation 0
if restart is None:
    start = 0 + 1
else:
    start = restart + 1

# Looping over generations
for g in range(start, max_generations):
    evolve_utils.log(f"\nINITIALIZE GENERATION {g}")

    # Initialising pool
    
