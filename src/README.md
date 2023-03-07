# Source code

## Description
This directory contains utility scripts for both data orchestration, analysis, and simulations.

## Contents
This directory contains:
- `data_utils.py`: Script for generating input data
- `analyse_utils.py`: Script for analysis
- `simulate.py`: Script for submitting a simulation
- `simulate_utils.py`: Utilitary functions for `simulate.py`
- `evolve.py`: Script for submitting an evolution
- `evolve_utils.py`: Utilitary functions for `evolve.py`
- `residues.py`: Residue-level data and parameters as a Python object
- `conditions.py`: Standard condition sets foorsimulatin as a Python object
- **`templates/`**: Templates to fill out for submitting simulations
    - `run_simulate.sh`: **Template** Shell script for submitting a simulation with Slurm
