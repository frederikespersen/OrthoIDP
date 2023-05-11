#!/bin/bash
#SBATCH --job-name=r_seqevo
#SBATCH --partition=sbinlab_ib2
#SBATCH --array=0-10%11
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00            # 1/2 day
#SBATCH -o results/range.out
#SBATCH -e results/range.err

# Getting job input options
input_seqs=(H1-0_AVG_MAX_KAPPA H1-0_AVG_MAX_KAPPA H1-0_AVG_MAX_KAPPA H1-0_AVG_MAX_KAPPA H1-0_AVG_MAX_KAPPA H1-0_AVG_MIN_KAPPA H1-0_AVG_MIN_KAPPA H1-0_AVG_MIN_KAPPA H1-0_AVG_MIN_KAPPA H1-0_AVG_MIN_KAPPA H1-0_AVG_MIN_KAPPA)
input_targets=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

# Getting kappa target
input_target=${input_targets[$SLURM_ARRAY_TASK_ID]}
input_seq=${input_seqs[$SLURM_ARRAY_TASK_ID]}
input_file="data/$input_seq.fasta"

# Displaying job info
echo "[`date`] STARTED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_target"

# DeiC env settings
source /groups/sbinlab/fpesce/.bashrc
conda activate openmm

# Submitting simulation
python ../../src/evolve.py --dir range/$input_target --fasta $input_file --measure kappa --target $input_target --L_at_half_acceptance 0.01 --simulated_annealing -g 50000

echo "[`date`] FINISHED Job Array ID: $SLURM_ARRAY_TASK_ID | Job ID: $SLURM_JOB_ID | Input: $input_file; $input_target"
