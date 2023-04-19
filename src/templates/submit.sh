#!/bin/bash

###################################################################################
# SUBMIT Template
# ---------------
# Template for a simulation submission script; Fill out fields marked by '<field>'.
# Place the file in the corresponding ~/exp/< >/ folder.
###################################################################################

#SBATCH --job-name=<jobname>
#SBATCH --partition=<partition>
#SBATCH --array=<range>%<simultaneous>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t <walltime>
#SBATCH -o results/out
#SBATCH -e results/err


# Getting job input file
input_files=($(ls data/*.fasta))
input_file=${input_files[$SLURM_ARRAY_TASK_ID]}

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source <source>
conda <env>

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file <options>

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
