"""
    Analyse_data
    --------------------------------------------------------------------------------

    Utils for analysing sequences prior to simulation.

    --------------------------------------------------------------------------------
"""


import json
import pandas as pd
from localcider.sequenceParameters import SequenceParameters


#························································································#
#································· G E N E R A L ········································#
#························································································#

def load_fasta_seq(fasta_path: str) -> tuple[str]:
    """

    Takes a FASTA file path, returns the sequence, id, and description.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        fasta_path: str
            Path to .fasta file

    Returns
    -------

        `seq`: `str`
            The sequence in the FASTA file

        `id`: `str`
            The id in the FASTA file

        `desc`: `str`
            The description in the FASTA file

    """

    # Loading file
    with open(fasta_path, 'r') as file:
        lines = file.readlines()

    # Removing whitespace
    lines = [line.strip() for line in lines]

    # Assembling data (Whether standard FASTA format or one-line-sequence FASTA format)
    seq = ''.join(lines[1:])
    id = lines[0].split(' ')[0][1:]
    desc = ' '.join(lines[0].split(' ')[1:])

    return seq, id, desc


#························································································#
def load_metadata(metadata_path: str) -> pd.DataFrame:
    """
    
    Takes the path to a meta data .json-file, returns the meta data loaded into a DataFrame.

    Row granularity is variants in 'data' field of .json file.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `metadata_path`: `str`
            The path to the meta deta .json file to load in

    Returns
    -------

        `metadata`: `pandas.DataFrame`
            A DataFrame containing the meta data fields

    """

    # Loading meta data file
    with open(metadata_path, 'r') as file:
        metadata = json.load(file)
    
    # Parsing data and templates fields
    data = pd.DataFrame(metadata['data']).transpose()
    templates = pd.DataFrame(metadata['templates']).transpose()

    # Joining templates on to data (unique fields only)
    unique_fields = list(set(templates) - set(data))
    metadata = data.join(templates[unique_fields])

    return metadata


#························································································#
#····························· A M I N O   A C I D S ····································#
#························································································#

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
"""

A list of the one-letter codes of the naturally occuring amino acids.

"""


#························································································#
def amino_acid_content(seqs: str|list|pd.Series) -> pd.DataFrame:
    """

    Takes one or more sequences, returns a DataFrame of the frequencies of amino acids.

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seqs`: `str | list | pandas.Series`
            Sequence(s) to calculate frequencies for

    Returns
    -------

        `freqs`: `pandas.DataFrame`
            A DataFrame with the amino acid frequency for each amino acid in the sequence(s).

    """

    # Formatting sequence(s) as a pd.Series
    if type(seqs) == str:
        seqs = pd.Series([seqs])
    elif type(seqs) == list:
        seqs = pd.Series(seqs)
    
    # Initiating DataFrame
    freqs = pd.DataFrame(index=seqs.index)

    # Calculating frequencies
    seqs = seqs.str.upper()
    for aa in amino_acids:
        freqs[aa] = seqs.apply(lambda seq: len(list(filter(lambda s: s == aa, seq))) / len(seq))

    return freqs


#························································································#
#··································· C I D E R ··········································#
#························································································#

def cider_parameters(seqs: str|list|pd.Series) -> pd.DataFrame:
    """

    Takes one or more sequences, returns a DataFrame of CIDER parameters.

    --------------------------------------------------------------------------------

    More on CIDER from PappuLab:
    - [CIDER](http://pappulab.wustl.edu/CIDER/about/)
    - [localCIDER](http://pappulab.github.io/localCIDER/)

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seqs`: `str | list | pandas.Series`
            Sequence(s) to calculate parameters for

    Returns
    -------

        `params`: `pandas.DataFrame`
            A DataFrame with select CIDER parameters calcualted for the sequence(s)

    """

    # Formatting sequence(s) as a pd.Series
    if type(seqs) == str:
        seqs = pd.Series([seqs])
    elif type(seqs) == list:
        seqs = pd.Series(seqs)
    
    # Mapping sequence(s) to a SequenceParameteres object
    seqs = seqs.map(SequenceParameters)

    # Initiating DataFrame
    params = pd.DataFrame(index=seqs.index)

    # Calculating parameters
    params['kappa'] = seqs.apply(lambda SeqOb: SeqOb.get_kappa())
    params['FCR'] = seqs.apply(lambda SeqOb: SeqOb.get_FCR())
    params['NCPR'] = seqs.apply(lambda SeqOb: SeqOb.get_NCPR())
    params['Hydrophobicity'] = seqs.apply(lambda SeqOb: SeqOb.get_mean_hydropathy())
    params['Frac. dis. prom.'] = seqs.apply(lambda SeqOb: SeqOb.get_fraction_disorder_promoting())

    return params

