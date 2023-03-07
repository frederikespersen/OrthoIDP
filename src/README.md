# Source code

## Description
This directory contains utility scripts for both data orchestration, analysis, and simulations.

## Contents
This directory contains:
- `simulate_utils.py`: Utilitary functions for `simulate.py`
- `simulate.py`: Python script for submitting a simulation
- `evolve_utils.py`: Utilitary functions for `evolve.py`
- `evolve.py`: Python script for submitting an evolution
- `residues.py`: Residue-level data and parameters as a Python object
- `conditions.py`: Standard condition sets foorsimulatin as a Python object
- `process_data.py`: Script for generating input data
- `analyse_data.py`: Script for analysing input data
- `analyse_results.py`: Script for analysing output results
- **`templates/`**: Templates to fill out for submitting simulations
    - `run_simulate.sh`: **Template** Shell script for submitting a simulation with Slurm