#!/bin/bash
#SBATCH --job-name=e1a_linkerss
#SBATCH --partition=sbinlab_ib2
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --array=0-188
#SBATCH -t 168:00:00
#SBATCH -o results/out
#SBATCH -e results/err


# Getting job input file
input_files=($(ls data/*.fasta))
input_file=${input_files[$SLURM_ARRAY_TASK_ID]}

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm.py -f $input_file -p CPU

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
