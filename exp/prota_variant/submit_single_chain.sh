#!/bin/bash
#SBATCH --job-name=sc_prota_variant
#SBATCH --partition=sbinlab_gpu
#SBATCH --array=0-19%20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o results/single_chain.out
#SBATCH -e results/single_chain.err


# Setting input options
input_seqs=(H1-0_VAR_k0.60 H1-0_VAR_k0.21 H1-0_VAR_k0.85 H1-0_VAR_k0.40 H1-0_VAR_k0.07 H1-0_VAR_k0.70 H1-0_VAR_k0.50 H1-0_VAR_k0.30 H1-0_VAR_k0.92 H1-0_VAR_k0.75 H1-0_VAR_k0.55 H1-0_VAR_k0.14 H1-0_VAR_k0.65 H1-0_VAR_k0.80 H1-0_VAR_k0.46 H1-0_VAR_k0.89 H1-0_VAR_k0.36 H1-0_VAR_k0.18 H1-0_VAR_k0.27 H1-0_VAR_k0.11)
input_cond="ionic_240"

# Get the current sequence
input_seq=${input_seqs[$SLURM_ARRAY_TASK_ID]}
input_file="data/$input_seq.fasta"
output_dir="single_chain/$input_cond/$input_seq"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file -c $input_cond -d $output_dir -n 5000000

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"