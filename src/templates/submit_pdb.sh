#!/bin/bash

###################################################################################
# SUBMIT SINGLE PDB Template
# --------------------------
# Template for a submitting a single topology in a pdb file for simulation
# Fill out fields marked by '<field>', and place the file in the corresponding ~/exp/< >/ folder.
###################################################################################

#SBATCH --job-name=<jobname>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t <walltime>
#SBATCH -o results/<jobname>.out
#SBATCH -e results/<jobname>.err


# Getting job input file
input_file="data/<filename>.pdb"

# Displaying job info
echo "[`date`] STARTED Job ID: $SLURM_JOB_ID | Input: $input_file"

# <server> env settings
source <source>
conda <env>

# Submitting simulation
python ../../src/simulate_openmm_top.py -t $input_file <options>

echo "[`date`] FINISHED Job ID: $SLURM_JOB_ID | Input: $input_file"
