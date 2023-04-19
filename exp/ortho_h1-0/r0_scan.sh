#!/bin/bash
#SBATCH --job-name=r0_scan
#SBATCH --partition=sbinlab_ib2
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --array=0-188
#SBATCH --mem=10G
#SBATCH -t 168:00:00
#SBATCH -o results/r0_scan.out
#SBATCH -e results/r0_scan.err


# Getting job input directories
input_dirs=($(ls -d $5))
input_dir=${input_dirs[$SLURM_ARRAY_TASK_ID]}

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_dir"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting r0 scan
python r0_scan.py -d $input_dir -s $1 -e $2 -i $3 -o $4

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_dir"