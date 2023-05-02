#!/bin/bash
#SBATCH --job-name=cons_seqevo
#SBATCH --partition=sbinlab_ib2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00            # 1/2 day
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
python ../../src/evolve.py --restart -1 --dir cons --fasta $input_file --measure kappa --target 0.1889910274293869 --L_at_half_acceptance 0.01 --simulated_annealing

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
