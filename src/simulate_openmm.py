"""
    Simulate_OpenMM
    --------------------------------------------------------------------------------

    Script for running an OpenMM simulation.
    Takes shell arguments, orchestrates preparations, and initiates simulation.

    --------------------------------------------------------------------------------

    Run `python openmm_simulate.py -h` for documentation on arguments.

"""


import argparse
import os
import sys


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="Simulate", description="Runs a single-chain CALVADOS simulation of an input sequence")
parser.add_argument('-f', '--fasta',
                    type=str,
                    required=True,
                    help="path to FASTA file with sequence to simulate")
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
parser.add_argument('-d', '--dir',
                    type=str,
                    required=False,
                    help="directory to store results in (default uses name of fasta file)")
parser.add_argument('-n', '--steps',
                    type=int,
                    required=False,
                    default=202000000,
                    help="the number of timesteps [10 fs] to run the simulation for (default: 202000000)")
parser.add_argument('-p', '--platform',
                    type=str,
                    required=False,
                    default='CUDA',
                    help="computional platform to use (default: CUDA)")
parser.add_argument('-s', '--source',
                    type=str,
                    required=False,
                    help="path source code for simulate utils (default: Same dir as simulate.py)")

# Parsing arguments
args = parser.parse_args()
fasta_path = args.fasta
source_path = args.source
platform = args.platform
boxlength = args.boxlength
cond = args.conditions
dir = args.dir
steps = args.steps

#························································································#

# Setting source path as the dir of simulate.py (this file) if nothing else is specified
if source_path is None:
    source_path = '/'.join(__file__.split('/')[:-1])
sys.path.append(source_path)

# Importing source code
from simulate_utils import openmm_simulate as simulate
from utils import read_fasta

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
simulate(sequence, boxlength=boxlength, dir=dir, steps=steps, eqsteps=1000, cond=cond, platform=platform, verbose=False, log=True)

#························································································#