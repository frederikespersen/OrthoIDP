#!/bin/bash
#SBATCH --job-name=tc_prota_h1-0
#SBATCH --partition=sbinlab_ib2
#SBATCH --array=0-17%18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -t 168:00:00
#SBATCH -o results/three_chain.out
#SBATCH -e results/three_chain.err


# Setting input options
input_tops=("2xH1-0_PROTA_WT" "H1-0_2xPROTA_WT")
input_conds=("default" "Borgia_in_silico" "ionic_165" "ionic_180" "ionic_205" "ionic_240" "ionic_290" "ionic_330" "ionic_340")

# Calculate the total number of tasks
num_tops=${#input_tops[@]}
num_conds=${#input_conds[@]}

# Calculate the current topology and condition indices
top_idx=$(($SLURM_ARRAY_TASK_ID % $num_tops))
cond_idx=$((($SLURM_ARRAY_TASK_ID - $top_idx) / $num_tops))

# Get the current topology and condition
input_top=${input_tops[$top_idx]}
input_cond=${input_conds[$cond_idx]}
input_file="data/$input_top.pdb"
output_dir="three_chain/$input_cond/$input_top"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_top.py -t $input_file -c $input_cond -d $output_dir  -p CPU

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"