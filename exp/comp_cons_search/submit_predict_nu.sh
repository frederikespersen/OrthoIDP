#!/bin/bash
#SBATCH --job-name=nu_predict
#SBATCH --partition=sbinlab_ib2
#SBATCH --array=0-9999%100
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -o results/predict_nu.out
#SBATCH -e results/predict_nu.err


# Getting arguments
split = 10000
split_idx = $SLURM_ARRAY_TASK_ID
input_file="data/idr_orthologs.csv"
output_file="data/idr_orthologs_processed.csv"

# Creating output file
if
    [ $split_idx==0 ]
then
    echo ,idr,idr_full_len,idr_ortholog,ortholog_full_len,ortholog_seq,scd,shd,kappa,fcr,mean_lambda,nu > $output_file
fi

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $split; $split_idx"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python predict_nu.py -n $split -x $split_idx -i $input_file -o $output_file -h "ortholog_seq"

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $split; $split_idx"