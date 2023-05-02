"""
    Evolve
    --------------------------------------------------------------------------------

    Script for running a sequence evolution algorithm.
    Takes shell arguments, orchestrates preparations, and initiates simulation.

    --------------------------------------------------------------------------------

    Run `python evolve.py -h` for documentation on arguments.

"""


import argparse
import os
import sys
import pandas as pd
import numpy as np
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
                    help="path to FASTA file with sequence to evolve")
parser.add_argument('-m', '--measure',
                    type=str,
                    required=False,
                    default='kappa',
                    help="what structural obervable to use as evolutionary constraint. Options: ['kappa']")
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
parser.add_argument('-g', '--max_gen',
                    type=int,
                    required=False,
                    default=100000,
                    help="The maximum amount of generations to generate (default: 100000 generations)")
parser.add_argument('-i', '--identical_generations',
                    type=bool,
                    required=False,
                    default=True,
                    help="Whether identical generations (sequences) are allowed or should be skipped (default: True); Consider ")
parser.add_argument('-s', '--source',
                    type=str,
                    required=False,
                    default=None,
                    help="path source code for utils (default: Same dir as evolve.py)")
parser.add_argument('-p', '--save_period',
                    type=int,
                    required=False,
                    default=10,
                    help="The saving period, i.e. how many generations to generate before saving")

# Parsing arguments
args = parser.parse_args()
restart = args.restart
dir = args.dir
measure = args.measure
L_at_half_acceptance = args.L_at_half_acceptance
target = args.target_value
source_path = args.source



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
store_filename = 'evolution.csv'

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

    # Calculating observable
    obs = compute_observable(seq)
    log.message(f"Calculated observable '{measure}' to {obs:.4f}")
    
    # Initialising DataFrame for storing results
    store = pd.DataFrame(columns=['sequence','observable','mc','c'])
    store.index.name = 'generation'
    store.loc[g] = {'sequence': seq,
                    'observable': obs,
                    'mc': True,
                    'c': c}

    # Saving initial generation
    store.to_csv(store_filename)

# Restarting
else:

    # Loading previous evolution
    store = pd.read_csv(store_filename)

    # Restarting from latest generation
    if restart == -1:
        restart = store.index[-1]

    log.message(f"RESTARTING FROM GENERATION {restart}")
    
    # Removing generations after restart if exists
    log.message(f"Removing simulation data from generations {', '.join([g for g in store.index if g > restart])}")
    store = store.drop([g for g in store.index if g > restart])

    # Using latest control parameter value
    c = store.c.iloc[-1]

log.message(f"Targeting a '{measure}' of {target}")
log.message(f"Starting evolution with Monte Carlo control parameter of {c:.6f}")

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
n = 0
for g in range(start, args.max_gen):

    # Timing
    t0 = time()

    #···························· I N I T I A L I S I N G ·······························#

    # Generating new sequence and initialising new entry 
    log.message("")
    log.message(f"INITIALIZING GENERATION {g}")
    last_seq = store[store.mc].sequence.iloc[-1]
    seq = evolve_utils.swap_sequence(last_seq)

    # Optionally checking whether the sequence has been previously generated:
    if args.identical_generations:
        while seq in store.sequence:
            seq = evolve_utils.swap_sequence(last_seq)
    
    log.message(f"Generated sequence: {seq}")

    #···················· C O M P U T I N G   O B S E R V A B L E ·······················#

    # Checking whether the sequence has been previously generated
    if seq in store.sequence:
        
        # Copying old entry
        log.message(f"Copying previous identical entry in generation {store[store.sequence  == seq].index[0]}")
        former_g = store[store.sequence  == seq].iloc[0]
        store.loc[g] = {'sequence': seq,
                        'observable': former_g.observable,
                        'mc': False,
                        'c': c}
    
    # Else calculating observable anew
    else:
        
        # Calculating observable
        obs = compute_observable(seq)
            
        # Filling out entry
        store.loc[g] = {'sequence': seq,
                        'observable': obs,
                        'mc': False,
                        'c': c}

    log.message(f"Calculated observable '{measure}' to {store.observable[g]:.4f}")

    #··················· M O N T E   C A R L O   E V O L U T I O N ······················#

    # Determining whether new generation is accepted
    last_g = store[store.mc].index[-1]
    accept, p = evolve_utils.mc_move(store.observable.loc[g], store.observable.loc[last_g], target, store.c.loc[g])

    if accept:

        # Accepting move
        store.mc[g] = True
        log.message(f"Accepted evolved generation with a probability of {p:.4%}")
    
    else:
        # Rejecting move
        store.mc[g] = False
        log.message(f"Rejected evolved generation with a probability of {p:.4%}")

    log.message(f"MC acceptance rate is currently at {store.mc.sum()/len(store):.2%}%")

    #································ P I C K L I N G ···································#

    # Saving current evolution
    n += 1
    if n >= args.save_period:
        log.message(f"Saving results to '{store_filename}' after {n} generations")
        n = 0
        store.to_csv(store_filename)

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
    