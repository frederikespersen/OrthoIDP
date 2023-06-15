#!/bin/bash
#SBATCH --job-name=ec_para_prota_variant
#SBATCH --partition=qgpu
#SBATCH --array=0-9%10
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o results/energy_calc.out
#SBATCH -e results/energy_calc.err


# Setting input options
input_tops=(H1-0_PROTA_WT_25nm H1-1_PROTA_WT_25nm H1-3_PROTA_WT_25nm H1-5_PROTA_WT_25nm H1-7_PROTA_WT_25nm H1-10_PROTA_WT_25nm H1-2_PROTA_WT_25nm H1-4_PROTA_WT_25nm H1-6_PROTA_WT_25nm H1-8_PROTA_WT_25nm)
input_cond="ionic_290"

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