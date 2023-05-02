#!/bin/bash
#SBATCH --job-name=cons_seqevo
#SBATCH --partition=sbinlab_ib2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 168:00:00            # 1 week
#SBATCH --mem=16G
#SBATCH -o results/cons.out
#SBATCH -e results/cons.err

# Getting job input file
input_file="data/H1-0_AVG.fasta"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/evolve.py --dir cons --fasta $input_file --measure kappa --target 0.15 --L_at_half_acceptance 0.001 --simulated_annealing

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
