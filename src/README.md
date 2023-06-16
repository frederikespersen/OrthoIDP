# Source code

## Description
This directory contains utility scripts for both data orchestration, analysis, and simulations.

## Contents
This directory contains:
- `utils.py`: Script with general utils
- `data_utils.py`: Script for generating input data
- `analyse_utils.py`: Script for analysis
- `simulate_utils.py`: Script for simulations
- `simulate_openmm_fasta.py`: Submission script for submitting an single-chain OpenMM CALVADOS simulation using a single `.fasta`  sequence as input
- `simulate_openmm_top.py`: Submission script for submitting an single-chain OpenMM CALVADOS simulation using a `.pdb` topology as input
- `evolve.py`: Submission script for submitting an evolution
- `evolve_utils.py`: Utilitary functions for `evolve.py`
- `residues.py`: Residue-level data and parameters as a Python object
- `conditions.py`: Standard condition sets foorsimulatin as a Python object
- `pretraj_fixer.py`: Script for postprocessing prematurely terminated simulations
- **`templates/`**: Templates to fill out for submitting simulations
    - `submit_*.sh`: **Template** Shell script for submitting a simulation with Slurm
