#!/bin/bash
#SBATCH --job-name=is_prota_h1-0
#SBATCH --partition=sbinlab_gpu
#SBATCH --array=0-13%14
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o results/ionic_scan.out
#SBATCH -e results/ionic_scan.err


# Setting input options
input_seqs=("H1-0_WT" "PROTA_WT")
input_conds=("ionic_165" "ionic_180" "ionic_205" "ionic_240" "ionic_290" "ionic_330" "ionic_340")

# Calculate the total number of tasks
num_seqs=${#input_seqs[@]}
num_conds=${#input_conds[@]}

# Calculate the current sequence and condition indices
seq_idx=$(($SLURM_ARRAY_TASK_ID % $num_seqs))
cond_idx=$((($SLURM_ARRAY_TASK_ID - $seq_idx) / $num_seqs))

# Get the current sequence and condition
input_seq=${input_seqs[$seq_idx]}
input_cond=${input_conds[$cond_idx]}
input_file="data/$input_seq.fasta"
output_dir="ionic_scan/$input_cond/$input_seq"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file -c $input_cond -d $output_dir

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"