# Source code

## Description
This directory contains utility scripts for both data orchestration, analysis, and simulations.

## Contents
This directory contains:
- **`sim`**: Scripts for simulations
  - `run_simulate.sh`: **Template** shell script for submitting a simulation with Slurm (*.gitignored*)
  - `simulate.py`: **Temlate** python script for submitting a simulation
  - `simulate_utils.py`: Utilitary functions for `simulate.py`
  - `residues.py`: Residue-level data and parameters as a Python object
  - `conditions.py`: Standard condition sets as a Python object
- `process_data.py`: Script for generating input data
- `analyse_data.py`: Script for analysing input data
- `analyse_results.py`: Script for analysing output results
