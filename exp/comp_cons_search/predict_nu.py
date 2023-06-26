"""
    Predict nu
    --------------------------------------------------------------------------------

    Parallel prediction of nu for ortholog sequences using SVR model.

    --------------------------------------------------------------------------------

    Run `python predict_nu.py -h` for documentation on arguments.

"""


import argparse
import os
import sys
import numpy as np
import pandas as pd
import itertools
from localcider.sequenceParameters import SequenceParameters
import joblib

sys.path.append('../../src')
from residues import residues


#························································································#

# Setting up argument parser
parser = argparse.ArgumentParser(prog="Predict nu", description="Predict nu using SVR model")

# Input arguments
parser.add_argument('-n', '--split',
                    type=int,
                    required=True,
                    help="The amount of subdivisions to split the input DataFrame into")
parser.add_argument('-x', '--split_idx',
                    type=int,
                    required=True,
                    help="The index of subdivision to treat in this script")
parser.add_argument('-i', '--input',
                    type=str,
                    required=True,
                    help="The .csv file containing the sequences")
parser.add_argument('-h', '--header',
                    type=str,
                    required=True,
                    help="The header for the column in the .csv input file that contains sequences")
parser.add_argument('-o', '--output',
                    type=str,
                    required=True,
                    help="The output file for results")

# Parsing arguments
args = parser.parse_args()

#························································································#

# Creating directory for output
dir = '/'.join(args.output.split('/')[:-1])
if dir == '': 
    dir = '.'
else:
    os.makedirs(dir, exist_ok=True)

#························································································#

# Getting residue data
_residues = residues.set_index('one')

# Loading SVR model
model = joblib.load('svr_model.joblib')

# Calculating features
def svr_features(seq):

    # Getting residue data
    res = _residues.loc[list(seq)]
    l = res.AH_lambda.to_numpy()
    q = res.q.to_numpy()

    # Fixing charges, for both FCR/SCD (q) and kappa (_seq)
    q[0] += 1
    q[-1] -= 1
    aa_q = ['A', 'K', 'K' 'D', 'D'] # q = 0 => A;   q = 1,2 => K;   q = -1,-2 => D
    _seq = aa_q[q[0]] + seq[1:-1] + aa_q[q[-1]]

    # Calculating position pairs for sequence decorator measures
    N = len(seq)
    ij_pairs = np.array(list(itertools.combinations(range(len(seq)),2)))
    ij_dist = np.diff(ij_pairs).flatten().astype(float)

    # Calculating features
    scd = np.sum(q[ij_pairs].prod(axis=1)*np.sqrt(ij_dist))/N
    shd = (l[ij_pairs].sum(axis=1)/np.abs(ij_dist)).sum()/N
    kappa = SequenceParameters(_seq).get_kappa()
    fcr = abs(q).mean()
    mean_lambda = l.mean()
    features = np.array([scd, shd, kappa, fcr, mean_lambda])

    return features

# Predicting
def svr_nu(features):
    nu = model.predict(features)
    return nu

#························································································#

# Loading data
data = pd.read_csv(args.input)

# Splitting into subdivision
data = np.array_split(data,args.split)[args.split_idx]

#························································································#

# Calculating features
features = np.stack(data[args.header].apply(svr_features))

# Predicting nu
nu = svr_nu(features)

#························································································#

# Assembling results
results = pd.concat([data,
                     pd.DataFrame(features, columns=['scd', 'shd', 'kappa', 'fcr', 'mean_lambda'], index=data.index),
                     pd.DataFrame(nu, columns=['nu'], index=data.index)],
                     axis=1)

# Appending results
results.to_csv('data/idr_orthologs_processed.csv', mode='a', header=False)
