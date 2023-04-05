#!/bin/bash
#SBATCH --job-name=h1-0_max_evo
#SBATCH --partition=sbinlab_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 168:00:00            # 1 week
#SBATCH --mem=16G
#SBATCH -o results/max.out
#SBATCH -e results/max.err

# Getting job input file
input_file="data/H1-0_AVG.fasta"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/evolve.py -d max -f $input_file --measure "Rg" --target_value 10 --L_at_half_acceptance 0.01 --simulated_annealing

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
