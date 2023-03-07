#!/bin/bash

## Setting SBATCH parameters
#SBATCH --job-name=idpevo_h1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 168:00:00            # 1 week
#SBATCH -o out
#SBATCH -e err
#SBATCH --partition=sbinlab_gpu
#SBATCH --gres=gpu:1

## Pointing to source code for Python dependencies
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

## Running simulation script
python evolve.py --target_Obs 0 --restart -1 --fasta /home/fknudsen/OrthoIDP/exp/idpevo_h1-0/data/H1-0_AVG.fasta --L_at_half_acceptance 0.01 --MC_move swap_full --simulated_annealing