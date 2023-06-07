#!/bin/bash
#SBATCH --job-name=tc_k2_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-19%20
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00                # 2 days
#SBATCH -o results/two_chain.out
#SBATCH -e results/two_chain.err


# Setting input options
input_tops=(H1-0_VAR_k0.07_PROTA_WT_25nm  H1-0_VAR_k0.27_PROTA_WT_25nm  H1-0_VAR_k0.40_PROTA_WT_25nm  H1-0_VAR_k0.55_PROTA_WT_25nm  H1-0_VAR_k0.70_PROTA_WT_25nm H1-0_VAR_k0.11_PROTA_WT_25nm  H1-0_VAR_k0.30_PROTA_WT_25nm  H1-0_VAR_k0.46_PROTA_WT_25nm  H1-0_VAR_k0.60_PROTA_WT_25nm  H1-0_VAR_k0.75_PROTA_WT_25nm H1-0_VAR_k0.14_PROTA_WT_25nm  H1-0_VAR_k0.36_PROTA_WT_25nm  H1-0_VAR_k0.50_PROTA_WT_25nm  H1-0_VAR_k0.65_PROTA_WT_25nm)
input_cond="ionic_240"

# Get the current sequence
input_top=${input_tops[$SLURM_ARRAY_TASK_ID]}
input_file="data/$input_top.pdb"
output_dir="two_chain/$input_cond/$input_top"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# ROBUST env settings
source /home/fknudsen/.bashrc
conda activate orthoidp

# Submitting simulation
python ../../src/simulate_openmm_top.py -t $input_file -c $input_cond -d $output_dir -n 200000000

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"