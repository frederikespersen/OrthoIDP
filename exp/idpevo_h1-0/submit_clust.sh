#!/bin/bash
#SBATCH --job-name=avg_clust
#SBATCH --partition=sbinlab_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 168:00:00
#SBATCH -o results/avg_clust.out
#SBATCH -e results/avg_clust.err


# Get the current condition
input_file="data/H1-0_AVG_CLUST.fasta"
output_dir="avg_clust"

# Displaying job info
echo "[`date`] STARTED Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file -d $output_dir

echo "[`date`] FINISHED Job ID: $SLURM_JOB_ID | Input: $input_file"