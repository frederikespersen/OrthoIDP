#!/bin/bash
#SBATCH --job-name=sc_rg_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-19%20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00                 # 2 days
#SBATCH -o results/single_chain.out
#SBATCH -e results/single_chain.err


# Setting input options
input_seqs=(H1-0_VAR_Rg3.44  H1-0_VAR_Rg3.61  H1-0_VAR_Rg3.75  H1-0_VAR_Rg3.86  H1-0_VAR_Rg3.96  H1-0_VAR_Rg4.07  H1-0_VAR_Rg4.22 H1-0_VAR_Rg3.50  H1-0_VAR_Rg3.66  H1-0_VAR_Rg3.79  H1-0_VAR_Rg3.90  H1-0_VAR_Rg3.99  H1-0_VAR_Rg4.15  H1-0_VAR_Rg4.25 H1-0_VAR_Rg3.56  H1-0_VAR_Rg3.72  H1-0_VAR_Rg3.83  H1-0_VAR_Rg3.93  H1-0_VAR_Rg4.03  H1-0_VAR_Rg4.19)
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
python ../../src/simulate_openmm_fasta.py -f $input_file -c $input_cond -d $output_dir -n 10000000 -b 100

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"