# Data

## Description
This directory contains input data for the project.
Most of the input processes of this directory is orchestrated in ``~/Data.ipynb`` using ``~/src/process_data.py``.

## Contents
This directory contains:
- **``seqs``**: Sequence files
    - **``raw``**: Genbank format files of raw protein sequences from accession numbers
    - **``idr``**: FASTA format files of longest IDR region of proteins in ``raw``
    - **``var``**: FASTA format files (one-line-sequence) of variants of sequences in ``ìdr``
    - ``seqs/seqs.json``: Metadata on files in ``seqs/raw``and ``seqs/idr``
    - ``seqs/vars.json``: Metadata on files in ``seqs/var``
