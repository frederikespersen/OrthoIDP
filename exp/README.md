# Results

## Description
This directory contains data from separate experiments.
Most of the input processes of this directory is orchestrated either on remote servers or in ``~/Data.ipynb`` and ``~/Analyse.ipynb``.


## Contents
This directory contains:
- **``initial``**: Initial test run of 3 variants of 2 histones to verify histones as a model system
- **``ortho``**: [*To be completed*]


Each experiment is structures as:
- **``<experiment>**``**: Files from a specific experiment
    - **``data``**: Input data for simulations and analysis
    - **``results``**: Output data from simulations and analysis
        - **``<variant_id>``**: Variant specific simulation
    - ``<experiment>.json``: Metadata on experiment (See below)


## Experiment metadata
Metadata on the IDR sequences used for the experiment is stored in .json files named like the experiment (i.e. ``~/exp/initial/initial.json``).

The file has the following schema:
- ``accession_numbers``: Contains the accession numbers to be used, typically RefSeq / MANE sequences
    - ``<protein_id>``: A custom identifier to be use for the protein and its variants
- ``data``: Metadata on the input data
    - ``<sequence_id>``: Like ``<protein_id>``, but sometimes with a suffix if several variants exists
        - ``template``: The template ``<protein_id>`` for the sequence
        - ``name``: The GenBank description of the protein
        - ``species``: The species from which the protein originates
        - ``sequence``: The sequence in the .fasta file with the corresponding ``<sequence_id>``
        - ``region``: The general location of the extracted IDR, like NTD, CTD, or INT (Internal)
        - ``location``: The precise location of the IDR in the full length protein
        - ``variant``: Which variant off the ``<protein_id>`` that is used

**To first generate data for the experiment the file must have entries in ``accession_numbers``.**