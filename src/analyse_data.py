"""
    --------------------------------------------------------------------------------
    ``ANALYSE_DATA``
    --------------------------------------------------------------------------------

    Utils for analysing sequences prior to simulation.

    --------------------------------------------------------------------------------
"""


import pandas as pd
from localcider.sequenceParameters import SequenceParameters


#························································································#
#································· G E N E R A L ········································#
#························································································#

def load_fasta_seq(fasta_path: str) -> tuple[str]:
    """
    Takes a FASTA file path, returns the id, description, and encoded sequence.

    --------------------------------------------------------------------------------

    :param fasta_path: Path to .fasta file
    :return: The sequence, id, and description in the FASTA file as a tuple

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
#··································· C I D E R ··········································#
#························································································#

def cider_parameters(seq: str, name=0) -> pd.DataFrame:
    """
    Takes a sequence, returns a DataFrame of its CIDER parameters.

    --------------------------------------------------------------------------------

    More on CIDER from PappuLab:
    - [CIDER](http://pappulab.wustl.edu/CIDER/about/)
    - [localCIDER](http://pappulab.github.io/localCIDER/)

    --------------------------------------------------------------------------------

    :param seq: Sequence to calculate parameters for
    :param name: Index of the single row in the DataFrame
    :return: A single-row DataFrame with select CIDER parameters

    """
    
    # Creating a SequenceParameters object from sequence
    SeqOb = SequenceParameters(seq)

    # Initiating DataFrame
    params = pd.DataFrame(index=[name])

    # Calculating parameters
    params['kappa'] = SeqOb.get_kappa()
    params['FCR'] = SeqOb.get_FCR()
    params['NCPR'] = SeqOb.get_NCPR()
    params['Hydrophobicity'] = SeqOb.get_mean_hydropathy()
    params['Frac. dis. prom.'] = SeqOb.get_fraction_disorder_promoting()

    return params

