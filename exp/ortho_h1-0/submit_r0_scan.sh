#!/bin/bash
#SBATCH --job-name=r0_scan
#SBATCH --partition=sbinlab_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -t 168:00:00
#SBATCH -o results/r0_scan.out
#SBATCH -e results/r0_scan.err

# Displaying job info
echo "[`date`] STARTED Job ID: $SLURM_JOB_ID"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting scan
python r0_scan.py

echo "[`date`] FINISHED Job ID: $SLURM_JOB_ID"
