#!/bin/bash
#SBATCH --job-name=tc_k2_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-19%20
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00               # 2 days
#SBATCH -o results/two_chain.out
#SBATCH -e results/two_chain.err


# Setting input options
input_tops=(H1-0_VAR_Rg3.44_PROTA_WT_25nm H1-0_VAR_Rg3.61_PROTA_WT_25nm H1-0_VAR_Rg3.75_PROTA_WT_25nm H1-0_VAR_Rg3.86_PROTA_WT_25nm H1-0_VAR_Rg3.96_PROTA_WT_25nm H1-0_VAR_Rg4.07_PROTA_WT_25nm H1-0_VAR_Rg4.22 H1-0_VAR_Rg3.50_PROTA_WT_25nm H1-0_VAR_Rg3.66_PROTA_WT_25nm H1-0_VAR_Rg3.79_PROTA_WT_25nm H1-0_VAR_Rg3.90_PROTA_WT_25nm H1-0_VAR_Rg3.99_PROTA_WT_25nm H1-0_VAR_Rg4.15_PROTA_WT_25nm H1-0_VAR_Rg4.25 H1-0_VAR_Rg3.56_PROTA_WT_25nm H1-0_VAR_Rg3.72_PROTA_WT_25nm H1-0_VAR_Rg3.83_PROTA_WT_25nm H1-0_VAR_Rg3.93_PROTA_WT_25nm H1-0_VAR_Rg4.03_PROTA_WT_25nm H1-0_VAR_Rg4.19)
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