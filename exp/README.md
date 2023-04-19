# Experiments

## Description
This directory contains data from separate experiments. The experiments all involve simulating IDRs using CALVADOS.

Most of the input processes of this directory is orchestrated either on remote servers or in `~/Data.ipynb` and `~/Analyse.ipynb`.


## Contents
This directory contains:
- **`initial/`**: Initial test run of 3 variants of 2 human histones (H2B, H1.0) to verify histones as a model system
- **`ortho_h1-0/`**: Single-chain simulations of orthologs of histone H1-0
- **`para_h1`**: Single-chain simulations of human histone H1 paralogs
- **`idp-evo_h1-0/`**: Using a $R_g$ evolution algorithm on an average histone H1.0-sequence
- **`e1a_linkers`**: Recreating the work of González-Foulet [Nature, 2022] using CALVADOS2
- **`prota_h1-0`**: Two-chain simulations of ProTα and histone H1.0

Each experiment is structures as:
- **`<experiment>/**`**: Files from a specific experiment
    - **`data/`**: Input data for simulations and analysis
    - **`results/`**: Output data from simulations (*.gitignored*)
        - **`<variant_id>/`**: Variant specific simulation (*.gitignored*)
    - `<experiment>.ipynb`: Notebook for orchestrating data and analysis
    - `<experiment>.json`: Metadata on experiment (See below)
    - `submit*.sh`: SLURM submission script(s)

## Experiment metadata
Metadata on the IDR sequences used for the experiment is stored in .json files named like the experiment (i.e. `~/exp/initial/initial.json`).

The file has the following schema:
- `templates`: The intrinsically disordered proteins to use for the analysis
    - `<protein_id>`: A custom identifier to be use for the protein and its variants
        - `uniprot_id`: The UniProt ID to use for the protein
        - `name`: The name of the protein in UniProt
        - `description`: The description of the protein in UniProt
        - `species`: The species from which the protein originates
        - `sequence`: The full protein sequence
- `data`: Metadata on the input data
    - `<sequence_id>`: Like `<protein_id>`, but sometimes with a suffix if several variants exists
        - `template`: The template `<protein_id>` for the sequence
        - `sequence`: The sequence in the .fasta file with the corresponding `<sequence_id>`
        - `region`: The general location of the extracted IDR, like NTD, CTD, or INT (Internal)
        - `location`: The precise location of the IDR in the full length protein
        - `variant`: Which variant of the `<protein_id>` that is used
