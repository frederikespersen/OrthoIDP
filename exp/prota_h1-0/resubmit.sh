#!/bin/bash
#SBATCH --job-name=tc_prota_h1-0
#SBATCH --partition=sbinlab_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 168:00:00
#SBATCH -o results/resubmit.out
#SBATCH -e results/resubmit.err

# Getting job input file
input_top="H1-0_PROTA_WT"
input_cond="default"
input_file="data/$input_top.pdb"
output_dir="two_chain/$input_cond/$input_top"

# Displaying job info
echo "[`date`] STARTED Job ID: $SLURM_JOB_ID | Input: $input_file"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_top.py -t $input_file -c $input_cond -d $output_dir -n 1000000000

echo "[`date`] FINISHED Job ID: $SLURM_JOB_ID | Input: $input_file"