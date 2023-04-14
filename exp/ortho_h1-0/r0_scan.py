"""
    R0_scan
    --------------------------------------------------------------------------------

    Script for scanning for the best common R0 value for computing the scaling exponent
    across all orthologs.

    --------------------------------------------------------------------------------
"""


# Imports
import os
import mdtraj as md
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from datetime import datetime as dt

import sys
sys.path.append("../../src")
from utils import log as logger


#························································································#

log = logger(write=True, print=False, file=f'results/r0_scan.log')

#························································································#

# Finding ortholog simulations and loading trajectories
orthologs = [id for id in os.listdir('results') if os.path.isdir(f'results/{id}')]

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

r0_min = 0.5
r0_max = 0.6
scan_len = 101
scan = np.linspace(r0_min, r0_max, scan_len)

#························································································#

results = np.empty((len(orthologs), scan_len))

# Looping over trajectories
for i, id in enumerate(orthologs):
    traj = md.load_dcd(f'results/{id}/traj.dcd', f'results/{id}/top.pdb')
    log.message(f"[{dt.now()}] FITTING SCAN FOR TRAJEJCTORY OF {id} [{i+1}/{len(orthologs)}]")

    # Looping over R0
    for j, r0 in enumerate(scan):
        results[i, j] = r0_chi2(r0, traj)

    # Saving results
    np.savetxt("results/r0_scan.csv", results, delimiter=',')
    log.message(f"[{dt.now()}] Saved to results/r0_scan.csv")

#························································································#