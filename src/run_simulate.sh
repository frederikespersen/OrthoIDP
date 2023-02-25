#!/bin/bash

##
##   Run_simulation
##   -------
##   [Template]
##   Takes a .fasta file, runs a simulation using the sequence and puts results in directory with the same name as the .fasta file.
##


## Setting SBATCH parameters [Note: specify partition for cluster]
#SBATCH --job-name=sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 500:00:00
#SBATCH -o out
#SBATCH -e err
#SBATCH --partition=sbinlab_gpu
#SBATCH --gres=gpu:1


## Pointing to source code (Using Francesco's)
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm


## Running simulation script
python sim.py --fasta $1