"""
    Simulation
    --------------------------------------------------------------------------------

    [Template]
    Script for running a simulation.
    Takes shell arguments, orchestrates preparations, and initiates simulation.
    May be experiment specific.

    --------------------------------------------------------------------------------
"""


import argparse
import os
from simulation_utils import simulate


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser()
parser.add_argument("--fasta", type=str)

# Parsing arguments
args = parser.parse_args()
fasta_path = args.fasta

#························································································#

# Extracting fasta file
with open(fasta_path, 'r') as fasta:
    seq_lines = fasta.readlines()[1:].strip()
    sequence = ''.join(seq_lines)

#························································································#

# Using fasta file name as simulation directory
dir = fasta_path.split('/')[-1].replace('.fasta', '')

# Creating directory for simulation
os.makedirs(dir, exist_ok=True)

#························································································#

# Starting simulation
simulate(sequence, boxlength=200, dir=dir, steps=202000000, eqsteps=1000, cond='default', platform='CUDA')

#························································································#