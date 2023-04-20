#!/bin/bash

###################################################################################
# SUBMIT COMBINATIONS Template
# --------------------------------
# Template for a submitting combinations of simulation parameters, like sequences and conditions.
# Fill out fields marked by '<field>', and place the file in the corresponding ~/exp/< >/ folder.
###################################################################################

#SBATCH --job-name=<jobname>
#SBATCH --partition=<partition>
#SBATCH --array=<no. combinations>%<simultaneous>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t <walltime>
#SBATCH -o results/<jobname>.out
#SBATCH -e results/<jobname>.err


# Setting input options for sequences and conditions
input_seqs=("H1-0_WT" "PROTA_WT")
input_conds=("ionic_165" "ionic_180" "ionic_205" "ionic_240" "ionic_290" "ionic_330" "ionic_340")

# Calculate the number of sequences and conditions
num_seqs=${#input_seqs[@]}
num_conds=${#input_conds[@]}

# Calculate the current sequence and condition indices
seq_idx=$(($SLURM_ARRAY_TASK_ID % $num_seqs))
cond_idx=$((($SLURM_ARRAY_TASK_ID - $seq_idx) / $num_seqs))

# Get the current sequence and condition
input_seq=${input_seqs[$seq_idx]}
input_cond=${input_conds[$cond_idx]}
input_file="data/$input_seq.fasta"
output_dir="<jobname>/$input_cond/$input_seq"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/simulate_openmm_fasta.py -f $input_file -c $input_cond -d $output_dir

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_cond"