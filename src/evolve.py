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
from process_data import read_fasta
from simulate_utils import simulate
from analyse_results import compute_rg


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
compute_observable = compute_rg

#························································································#

# Starting from scratch
if restart is None:
    g = 0
    evolve_utils.log(f'INITIALIZING INPUT SEQUENCE AND SIMULATE', header=True)

    # Calculating initial Monte Carlo control parameter
    c = -1 * args.L_at_half_acceptance / np.log(0.5)

    # Initialising DataFrame for storing results
    store = pd.DataFrame(columns=['sequence','observable','simulate','mc','mc_cp'])
    store.iloc[0] = dict(fasta=init_sequence, obs=None, simulate=True, mc=True, mc_cp=c)

    # Running intial simulation
    if os.path.isdir('g0'):
        evolve_utils.log('Start sequence simulation available, skipping simulation.')
    else:
        simulate(sequence=init_sequence, boxlength=boxlength, dir="g0", steps=steps, eqsteps=1000, cond=cond, vmodel=3, platform='CODA')

    # Calculating mean observable
    init_traj = md.load_dcd('g0/traj.dcd', 'g0/top.pdb')
    store.obs[0] = compute_observable(init_sequence, init_traj).mean()

    # Calculating total energy for reweighting
    
    pool_d = [d]
    pool_mask = [mask]
    pool_ndx = [0]
    emat = calcEtot(residues,d,parameters.loc[ID],store.fasta[0], mask)
    len_pool_i = len(emat)
    store.to_pickle('./evolution.pkl')

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

    # Start from last control parameter value
    c = store.mc_cp.iloc[-1]


#························································································#

