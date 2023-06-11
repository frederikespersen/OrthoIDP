#!/bin/bash
#SBATCH --job-name=ec_prota_h1-0
#SBATCH --partition=qgpu
#SBATCH --array=0-8%9
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH -t 168:00:00
#SBATCH -o results/energy_calc.out
#SBATCH -e results/energy_calc.err


# Setting input options
input_conds=("default" "Borgia_in_silico" "ionic_165" "ionic_180" "ionic_205" "ionic_240" "ionic_290" "ionic_330" "ionic_340")

# Getting arguments
input_cond=${input_conds[$SLURM_ARRAY_TASK_ID]}
input_dir="results/two_chain_25nm/$input_cond/H1-0_PROTA_WT_25nm"
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