#!/bin/bash
#SBATCH --job-name=sc_comb_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-49%50
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00                 # 2 days
#SBATCH -o results/single_chain.out
#SBATCH -e results/single_chain.err


# Setting input options
input_seqs=(H1-0_VAR_Rg3.43_k0.62 H1-0_VAR_Rg3.63_k0.54 H1-0_VAR_Rg3.74_k0.42 H1-0_VAR_Rg3.83_k0.39 H1-0_VAR_Rg3.89_k0.34 H1-0_VAR_Rg3.96_k0.37 H1-0_VAR_Rg4.03_k0.20 H1-0_VAR_Rg4.18_k0.27 H1-0_VAR_Rg4.26_k0.23 H1-0_VAR_Rg3.46_k0.58 H1-0_VAR_Rg3.65_k0.47 H1-0_VAR_Rg3.76_k0.38 H1-0_VAR_Rg3.85_k0.30 H1-0_VAR_Rg3.92_k0.25 H1-0_VAR_Rg3.97_k0.25 H1-0_VAR_Rg4.05_k0.30 H1-0_VAR_Rg4.19_k0.21 H1-0_VAR_Rg4.27_k0.33 H1-0_VAR_Rg3.51_k0.50 H1-0_VAR_Rg3.68_k0.41 H1-0_VAR_Rg3.79_k0.33 H1-0_VAR_Rg3.87_k0.26 H1-0_VAR_Rg3.92_k0.33 H1-0_VAR_Rg3.98_k0.28 H1-0_VAR_Rg4.07_k0.25 H1-0_VAR_Rg4.20_k0.17H1-0_VAR_Rg3.57_k0.45 H1-0_VAR_Rg3.72_k0.36 H1-0_VAR_Rg3.81_k0.27 H1-0_VAR_Rg3.87_k0.42 H1-0_VAR_Rg3.94_k0.22 H1-0_VAR_Rg3.99_k0.22 H1-0_VAR_Rg4.09_k0.19 H1-0_VAR_Rg4.23_k0.22H1-0_VAR_Rg3.60_k0.48 H1-0_VAR_Rg3.72_k0.49 H1-0_VAR_Rg3.82_k0.45 H1-0_VAR_Rg3.88_k0.37 H1-0_VAR_Rg3.94_k0.28 H1-0_VAR_Rg3.99_k0.34 H1-0_VAR_Rg4.15_k0.18 H1-0_VAR_Rg4.24_k0.18H1-0_VAR_Rg3.62_k0.42 H1-0_VAR_Rg3.72_k0.54 H1-0_VAR_Rg3.83_k0.34 H1-0_VAR_Rg3.89_k0.29 H1-0_VAR_Rg3.95_k0.32 H1-0_VAR_Rg4.02_k0.26 H1-0_VAR_Rg4.15_k0.21 H1-0_VAR_Rg4.24_k0.26)
input_cond="default"

# Get the current sequence
input_seq=${input_seqs[$SLURM_ARRAY_TASK_ID]}
input_file="data/$input_seq.fasta"
output_dir="single_chain/$input_cond/$input_seq"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# ROBUST env settings
source /home/fknudsen/.bashrc
conda activate orthoidp

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file -c $input_cond -d $output_dir -n 10000000 -b 100 -p CPU

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"