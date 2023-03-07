"""
    Simulation
    --------------------------------------------------------------------------------

    [Template]
    Script for running a simulation.
    Takes shell arguments, orchestrates preparations, and initiates simulation.
    Template variables are marked like `<variable>`.

    --------------------------------------------------------------------------------

    Run `python simulate.py -h` for documentation on arguments.

"""


import argparse
import os
import sys

sys.path.append("<src path>")
from simulate_utils import simulate
from conditions import conditions
from process_data import read_fasta


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="Simulate", description="Runs a single-chain CALVADOS simulation of an input sequence")
parser.add_argument('-f', '--fasta',
                    type=str,
                    required=True,
                    help="path to FASTA file with sequence to simulate")
parser.add_argument('-p', '--platform',
                    type=str,
                    required=True,
                    help="computional platform to use (i.e. CUDA)")
parser.add_argument('-b', '--boxlength',
                    type=float,
                    required=False,
                    default=200,
                    help="boxlength [nm] for simulation cube (default 200 nm)")
parser.add_argument('-c', '--conditions',
                    type=str,
                    required=False,
                    default='default',
                    choices=conditions['name'],
                    help="physical conditions to simulate under (default: 'default')")
parser.add_argument('-d', '--dir',
                    type=str,
                    required=False,
                    help="directory to store results in (default uses name of fasta file)")
parser.add_argument('-s', '--steps',
                    type=int,
                    required=False,
                    default=202000000,
                    help="the number of timesteps [10 fs] to run the simulation for (default: 202000000)")

# Parsing arguments
args = parser.parse_args()
fasta_path = args.fasta
platform = args.platform
boxlength = args.boxlength
cond = args.conditions
dir = args.dir
steps = args.steps

#························································································#

# Extracting fasta file
sequence = read_fasta(fasta_path, just_seq=True)

#························································································#

# Using fasta file name as simulation directory if nothing else is specified
if dir is None:
    dir = fasta_path.split('/')[-1].replace('.fasta', '')

# Creating directory for simulation
os.makedirs(dir, exist_ok=True)

#························································································#

# Starting simulation
simulate(sequence, boxlength=boxlength, dir=dir, steps=steps, eqsteps=1000, cond=cond, platform=platform)

#························································································#