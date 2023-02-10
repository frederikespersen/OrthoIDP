# Data
This directory contains:
- **`histones/`**: FASTA sequences of terminal histone IDRs
- `Data.ipynb`: Notebook for orchestrating data download and processing of sequences
- `accession_numbers.json`: JSON containing accession numbers for proteins (see format below)

## Format of accession_numbers.json
The JSON has the following structure:
- `[Group]`: The directory for the data
  - `[Protein name]`: The name for the data file (followed by `terminal`)
    - `acc_num`: The accession number of the chosen RefSeq for the protein
    - `terminal`: Whether to use the first (`NTD`) or last (`CTD`) disordered region in the protein
