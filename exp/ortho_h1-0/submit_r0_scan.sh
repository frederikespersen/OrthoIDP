#!/bin/bash
#SBATCH --job-name=r0_scan
#SBATCH --partition=sbinlab_ib2
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -t 1:00:00
#SBATCH -o results/r0_scan.out
#SBATCH -e results/r0_scan.err

# Setting scan parameters
start=5
end=6
increment=0.01

# Creating new .csv file for appending results with initial header line
output="results/r0_scan.csv"
scan=`seq -s "," $start $increment $end`
header="id,$scan"
echo $header > $output

# Submitting scans
sbatch r0_scan.sh $start $end $increment $output "../ortho_h1-0/results/*/"