#!/bin/bash
#SBATCH --job-name=ec_rg_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-19%20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o results/energy_calc.out
#SBATCH -e results/energy_calc.err


# Setting input options
input_tops=(H1-0_VAR_k0.60_PROTA_WT_25nm H1-0_VAR_k0.21_PROTA_WT_25nm H1-0_VAR_k0.85_PROTA_WT_25nm H1-0_VAR_k0.40_PROTA_WT_25nm H1-0_VAR_k0.07_PROTA_WT_25nm H1-0_VAR_k0.70_PROTA_WT_25nm H1-0_VAR_k0.50_PROTA_WT_25nm H1-0_VAR_k0.30_PROTA_WT_25nm H1-0_VAR_k0.92_PROTA_WT_25nm H1-0_VAR_k0.75_PROTA_WT_25nm H1-0_VAR_k0.55_PROTA_WT_25nm H1-0_VAR_k0.14_PROTA_WT_25nm H1-0_VAR_k0.65_PROTA_WT_25nm H1-0_VAR_k0.80_PROTA_WT_25nm H1-0_VAR_k0.46_PROTA_WT_25nm H1-0_VAR_k0.89_PROTA_WT_25nm H1-0_VAR_k0.36_PROTA_WT_25nm H1-0_VAR_k0.18_PROTA_WT_25nm H1-0_VAR_k0.27_PROTA_WT_25nm H1-0_VAR_k0.11_PROTA_WT_25nm)
input_cond="ionic_240"

# Getting arguments
input_top=${input_tops[$SLURM_ARRAY_TASK_ID]}
input_dir="results/two_chain/$input_cond/$input_top"
input_traj="$input_dir/traj.dcd"
input_top="$input_dir/top.pdb"
output_file="$input_dir/interaction_energy.csv"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_traj; $input_top"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/interaction_energy.py -t $input_traj -p $input_top -x "chainid 0" -y "chainid 1" -c $input_cond -o $output_file --com --minimum_inter

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_traj; $input_top"