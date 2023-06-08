#!/bin/bash
#SBATCH --job-name=tc_prota_h1-0
#SBATCH --partition=qgpu
#SBATCH --array=0-8%9
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o results/two_chain_25nm.out
#SBATCH -e results/two_chain_25nm.err


# Setting input options
input_top="H1-0_PROTA_WT_25nm"
input_conds=("default" "Borgia_in_silico" "ionic_165" "ionic_180" "ionic_205" "ionic_240" "ionic_290" "ionic_330" "ionic_340")

# Get the current condition
input_cond=${input_conds[$SLURM_ARRAY_TASK_ID]}
input_file="data/$input_top.pdb"
output_dir="two_chain_25nm/$input_cond/$input_top"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# ROBUST env settings
source /home/fknudsen/.bashrc
conda activate orthoidp

# Submitting simulation
python ../../src/simulate_openmm_top.py -t $input_file -c $input_cond -d $output_dir -n 1000000000

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"