#!/bin/bash
#SBATCH --job-name=ortho_h1-0
#SBATCH --ntasks=1
#SBATCH --partition=sbinlab_gpu
#SBATCH --array=0-190%10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
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
python ../../src/simulate.py -f $input_file

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file"
