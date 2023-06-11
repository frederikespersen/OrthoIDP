#!/bin/bash
#SBATCH --job-name=tc_para_prota
#SBATCH --partition=qgpu
#SBATCH --array=0-9%10
#SBATCH --gres=gpu:1
#SBATCH -t 72:00:00              # 3 days
#SBATCH -o results/two_chain.out
#SBATCH -e results/two_chain.err


# Setting input options
input_tops=(H1-0_PROTA_WT_25nm H1-1_PROTA_WT_25nmH1-3_PROTA_WT_25nmH1-5_PROTA_WT_25nmH1-7_PROTA_WT_25nm H1-10_PROTA_WT_25nmH1-2_PROTA_WT_25nmH1-4_PROTA_WT_25nmH1-6_PROTA_WT_25nmH1-8_PROTA_WT_25nm)
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
python ../../src/simulate_openmm_top.py -t $input_file -c $input_cond -d $output_dir -n 1000000000

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"