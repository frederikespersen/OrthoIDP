"""
    R0_scan
    --------------------------------------------------------------------------------

    Script for computing the Chi squared value for a range of r0-values for a trajectory

    --------------------------------------------------------------------------------
"""


# Imports
import argparse
import os
import mdtraj as md
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from datetime import datetime as dt
from csv import writer

import sys
sys.path.append("../../src")
from utils import log as logger


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="R0 Scan", description="Runs a scan of an R0 range for a trajectory")
parser.add_argument('-d', '--dir',
                    type=str,
                    required=True,
                    help="Directory containing trajectory (traj.dcd) and topology (top.pdb)")
parser.add_argument('-s', '--start',
                    type=float,
                    required=True,
                    help="The minimum r0 value to check")
parser.add_argument('-e', '--end',
                    type=float,
                    required=True,
                    help="The maximum r0 value to check")
parser.add_argument('-i', '--increment',
                    type=float,
                    required=True,
                    help="The increment in the range of r0 values to check")
parser.add_argument('-o', '--output',
                    type=str,
                    required=True,
                    help="The .csv file to append the results to")
args = parser.parse_args()

#························································································#

# Loading trajectory
traj = md.load_dcd(f'{args.dir}traj.dcd', f'{args.dir}top.pdb')
id = args.dir.split('/')[-2]

#························································································#

# Defining chi2 computation
def r0_chi2(r0, traj):
    pairs = traj.top.select_pairs('all','all')
    ij = pairs[:,1] - pairs[:,0]

    # Calculating interresidue cartesian distances across frames
    d = md.compute_distances(traj, pairs)
    d_mean = d.mean(axis=0)
    d_mean_err = d.std(axis=0, ddof=0)/np.sqrt(d.shape[0])

    # Defining model
    model = lambda x, v: r0 * np.power(x,v)
    p0 = [0.6]
    
    # Fitting model
    ij_cutoff=10
    mask = ij>ij_cutoff
    popt, pcov = curve_fit(model, ij[mask], d_mean[mask], p0=p0, sigma=d_mean_err[mask])
    v = popt[0]

    # Calculating chi2
    y = d_mean[mask]
    y_fit = model(ij[mask], v)
    chi2 = (((y_fit-y)/y)**2).sum()

    return chi2

#························································································#

# Defining scanning range
scan = np.arange(args.start, args.end + args.increment, args.increment)

#························································································#

# Looping over R0
results = [id]
for r0 in scan:
    results.append(r0_chi2(r0, traj))

# Saving results
with open(args.output, 'a') as file:
    writer(file).writerow(results)

#························································································#