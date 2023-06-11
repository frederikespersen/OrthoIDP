#!/bin/bash
#SBATCH --job-name=ec_comb_prota_variant
#SBATCH --partition=qgpu
#SBATCH --array=0-49%50
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o results/energy_calc.out
#SBATCH -e results/energy_calc.err


# Setting input options
input_tops=(H1-0_VAR_Rg3.43_k0.62_PROTA_WT_25nm H1-0_VAR_Rg3.63_k0.54_PROTA_WT_25nm H1-0_VAR_Rg3.74_k0.42_PROTA_WT_25nm H1-0_VAR_Rg3.83_k0.39_PROTA_WT_25nm H1-0_VAR_Rg3.89_k0.34_PROTA_WT_25nm H1-0_VAR_Rg3.96_k0.37_PROTA_WT_25nm H1-0_VAR_Rg4.03_k0.20_PROTA_WT_25nm H1-0_VAR_Rg4.18_k0.27_PROTA_WT_25nm H1-0_VAR_Rg4.26_k0.23_PROTA_WT_25nm H1-0_VAR_Rg3.46_k0.58_PROTA_WT_25nm H1-0_VAR_Rg3.65_k0.47_PROTA_WT_25nm H1-0_VAR_Rg3.76_k0.38_PROTA_WT_25nm H1-0_VAR_Rg3.85_k0.30_PROTA_WT_25nm H1-0_VAR_Rg3.92_k0.25_PROTA_WT_25nm H1-0_VAR_Rg3.97_k0.25_PROTA_WT_25nm H1-0_VAR_Rg4.05_k0.30_PROTA_WT_25nm H1-0_VAR_Rg4.19_k0.21_PROTA_WT_25nm H1-0_VAR_Rg4.27_k0.33_PROTA_WT_25nm H1-0_VAR_Rg3.51_k0.50_PROTA_WT_25nm H1-0_VAR_Rg3.68_k0.41_PROTA_WT_25nm H1-0_VAR_Rg3.79_k0.33_PROTA_WT_25nm H1-0_VAR_Rg3.87_k0.26_PROTA_WT_25nm H1-0_VAR_Rg3.92_k0.33_PROTA_WT_25nm H1-0_VAR_Rg3.98_k0.28_PROTA_WT_25nm H1-0_VAR_Rg4.07_k0.25_PROTA_WT_25nm H1-0_VAR_Rg4.20_k0.17_PROTA_WT_25nm H1-0_VAR_Rg3.57_k0.45_PROTA_WT_25nm H1-0_VAR_Rg3.72_k0.36_PROTA_WT_25nm H1-0_VAR_Rg3.81_k0.27_PROTA_WT_25nm H1-0_VAR_Rg3.87_k0.42_PROTA_WT_25nm H1-0_VAR_Rg3.94_k0.22_PROTA_WT_25nm H1-0_VAR_Rg3.99_k0.22_PROTA_WT_25nm H1-0_VAR_Rg4.09_k0.19_PROTA_WT_25nm H1-0_VAR_Rg4.23_k0.22_PROTA_WT_25nm H1-0_VAR_Rg3.60_k0.48_PROTA_WT_25nm H1-0_VAR_Rg3.72_k0.49_PROTA_WT_25nm H1-0_VAR_Rg3.82_k0.45_PROTA_WT_25nm H1-0_VAR_Rg3.88_k0.37_PROTA_WT_25nm H1-0_VAR_Rg3.94_k0.28_PROTA_WT_25nm H1-0_VAR_Rg3.99_k0.34_PROTA_WT_25nm H1-0_VAR_Rg4.15_k0.18_PROTA_WT_25nm H1-0_VAR_Rg4.24_k0.18_PROTA_WT_25nm H1-0_VAR_Rg3.62_k0.42_PROTA_WT_25nm H1-0_VAR_Rg3.72_k0.54_PROTA_WT_25nm H1-0_VAR_Rg3.83_k0.34_PROTA_WT_25nm H1-0_VAR_Rg3.89_k0.29_PROTA_WT_25nm H1-0_VAR_Rg3.95_k0.32_PROTA_WT_25nm H1-0_VAR_Rg4.02_k0.26_PROTA_WT_25nm H1-0_VAR_Rg4.15_k0.21_PROTA_WT_25nm H1-0_VAR_Rg4.24_k0.26_PROTA_WT_25nm)
input_cond="ionic_330"

# Getting arguments
input_top=${input_tops[$SLURM_ARRAY_TASK_ID]}
input_dir="results/two_chain/$input_cond/$input_top"
input_traj="$input_dir/traj.dcd"
input_top="$input_dir/top.pdb"
output_file="$input_dir/interaction_energy.csv"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_traj; $input_top"

# ROBUST env settings
source /home/fknudsen/.bashrc
conda activate orthoidp

# Submitting simulation
python ../../src/interaction_energy.py -t $input_traj -p $input_top -x "chainid 0" -y "chainid 1" -c $input_cond -o $output_file --com --minimum_inter

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_traj; $input_top"