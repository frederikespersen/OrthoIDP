#!/bin/bash

###################################################################################
# SUBMIT SINGLE FASTA Template
# ----------------------------
# Template for a submitting a single sequence in a fasta file for simulation
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
input_file="data/<filename>.fasta"

# Displaying job info
echo "[`date`] STARTED Job ID: $SLURM_JOB_ID | Input: $input_file"

# <server> env settings
source <source>
conda <env>

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file <options>

echo "[`date`] FINISHED Job ID: $SLURM_JOB_ID | Input: $input_file"
