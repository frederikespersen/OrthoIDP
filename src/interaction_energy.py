"""
    Interaction energy
    --------------------------------------------------------------------------------

    Script for calculating the interaction energy between two trajectory selections
    Takes shell arguments, calculates interaction energy for each frame,
    and saves results to file.
    Calculates using CALVADOS force field.

    Includes options to compute collective variables like center of mass distance
    and minimal intermolecular (between selections) distance

    --------------------------------------------------------------------------------

    Run `python interaction_energy.py -h` for documentation on arguments.

"""


import argparse
import os
import sys
import mdtraj as md
import pandas as pd
import numpy as np


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="Interaction energy", description="Calculates interaction energy for each frame in trajectory")

# Input arguments
parser.add_argument('-t', '--traj',
                    type=str,
                    required=True,
                    help="The path to the trajectory DCD file")
parser.add_argument('-p', '--top',
                    type=str,
                    required=True,
                    help="The path to the topology PDB file")
parser.add_argument('-a', '--sel1',
                    type=str,
                    required=True,
                    help="The first selection of the trajectory (Lik 'chainid 0')")
parser.add_argument('-b', '--sel2',
                    type=str,
                    required=True,
                    help="The second selection of the trajectory (Lik 'chainid 0')")
parser.add_argument('-c', '--conditions',
                    type=str,
                    required=True,
                    help="Which conditions to calculate energy for [Default: 'defualt']")
parser.add_argument('-o', '--output',
                    type=str,
                    required=True,
                    help="The path to the file to output results to (.csv, overwritten)")

# Energy arguments
parser.add_argument('-a', '--ashbaugh_hatch',
                    type=bool,
                    required=False,
                    default=True,
                    help="Whether to include Ashbaugh-Hatch potential in energy calculation [Default: True]")
parser.add_argument('-d', '--debye_huckel',
                    type=bool,
                    required=False,
                    default=True,
                    help="Whether to include Debye-Hückel potential in energy calculation [Default: True]")
parser.add_argument('-h', '--harmonic',
                    type=bool,
                    required=False,
                    default=False,
                    help="Whether to include Harmonic potential in energy calculation [Default: False]")

# Collective variables
parser.add_argument('-m', '--com',
                    action='store_true',
                    required=False,
                    help="Whether to calculate center of mass distance as collective variable")
parser.add_argument('-i', '--minimum_inter',
                    action='store_true',
                    required=False,
                    help="Whether to calculate the minimum interresidue distance as collective variable")

# Other arguments
parser.add_argument('-s', '--source',
                    type=str,
                    required=False,
                    help="Path source code for utils (default: Same dir as interaction_energy.py)")

# Parsing arguments
args = parser.parse_args()

#························································································#

# Setting source path as the dir of simulate.py (this file) if nothing else is specified
if args.source is None:
    source_path = '/'.join(__file__.split('/')[:-1])
else:
    source_path = args.source
sys.path.append(source_path)

# Importing source code
import analyse_utils
from conditions import conditions

#························································································#

# Creating directory for output
dir = '/'.join(args.output.split('/')[:-1])
os.makedirs(dir, exist_ok=True)

#························································································#

# Loading trajectory
traj = md.load_dcd(args.traj, args.top)

#························································································#

# Calculating residue distances between selections
pairs = traj.top.select_pairs(args.sel1, args.sel2)
distances = md.compute_distances(traj, pairs)

#························································································#

# Computing energy
E = analyse_utils.compute_energy(traj, args.conditions, distances=distances,
                                 ah=args.ashbaugh_hatch,
                                 dh=args.debye_huckel,
                                 hb=args.harmonic,
                                 pairs_ij=pairs)

# Assembling to DataFrame
potentials = []
if args.ashbaugh_hatch:
    potentials.append("Ashbaugh-Hatch [kJ/mol]")
if args.debye_huckel:
    potentials.append("Debye-Hückel [kJ/mol]")
if args.harmonic:
    potentials.append("Harmonic bond [kJ/mol]")
df = pd.DataFrame(E.sum(axis=2), columns=potentials)

#························································································#

# Computing minimum interresidue distance
if args.minimum_inter:
    df['Minimum interresidue distance [nm]'] = distances.min(axis=1)

#························································································#

# Computing center of mass for each frame (Optional)
if args.com:
    com1 = analyse_utils.compute_com(traj.atom_slice(np.unique(pairs[:,0])))
    com2 = analyse_utils.compute_com(traj.atom_slice(np.unique(pairs[:,1])))
    com_diff = (((com1-com2)**2).sum(axis=1))**0.5
    df['Center of mass distance [nm]'] = com_diff

#························································································#

# Saving results
df.to_csv(args.output)