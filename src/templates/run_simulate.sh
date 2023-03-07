#!/bin/bash

##   Run_simulate
##   -------
##   [Template]
##   Takes a .fasta file, runs a simulation using the sequence and puts results in directory with the same name as the .fasta file.
##   Template variables are marked like `<variable>`.


## Setting SBATCH parameters
#SBATCH --job-name=sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t <walltime-hours>:00:00
#SBATCH -o out
#SBATCH -e err
#SBATCH --partition=<partition>
#SBATCH --gres=gpu:1


## Pointing to source code for Python dependencies
source <source>
conda activate <env>


## Running simulation script
python <simulate_submit> --fasta $1