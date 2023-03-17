"""
    Analyse_utils
    --------------------------------------------------------------------------------

    Utils for analysis.

    --------------------------------------------------------------------------------
"""


import json
import pandas as pd
from localcider.sequenceParameters import SequenceParameters
import mdtraj as md
import numpy as np

from residues import residues
import data_utils


#························································································#
#································· G E N E R A L ········································#
#························································································#

def load_fasta_seq(fasta_path: str) -> tuple:
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
amino_acid_types = {
    'A': 'Hydrophobic',
    'C': 'Polar',
    'D': 'Negative',
    'E': 'Negative',
    'F': 'Hydrophobic',
    'G': 'Special',
    'H': 'Polar',
    'I': 'Hydrophobic',
    'K': 'Positive',
    'L': 'Hydrophobic',
    'M': 'Hydrophobic',
    'N': 'Polar',
    'P': 'Special',
    'Q': 'Polar',
    'R': 'Positive',
    'S': 'Polar',
    'T': 'Polar',
    'V': 'Hydrophobic',
    'W': 'Hydrophobic',
    'Y': 'Polar'
}
"""

A dict of the general types (Hydrophobic, Polar, Positive, Negative, Special) of amino acids.

--------------------------------------------------------------------------------

Schema
------

    `<AA>`: One-letter code for amino acid

        `<type>`: The general type of the amino acid


"""


#························································································#
def amino_acid_content(seqs) -> pd.DataFrame:
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
def average_sequence(seqs) -> str:
    """
    
    Takes a list or Series of sequences, returns an average sequence.
    
    Averaging occurs by first generating the average amino acid counts and then 
    assembling a random sequence from that.
    The average amino acid counts are found by multiplying the average amino acid
    frequency of each amino acid by the average sequence length.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seqs`: `list|pandas.Series`
            Sequences to average

    Returns
    -------

        `avg`: `str`
            The generated average series
    
    """

    # Setting datatype to Series
    seqs = pd.Series(seqs)

    # Finding average amino acid frequencies
    freqs = amino_acid_content(seqs).mean()

    # Finding average length
    N = seqs.str.len().mean()

    # Finding average amino acid counts
    counts = freqs * N

    # Assembling random sequence using counts
    avg = ''
    for aa, c in counts.items():
        avg += aa * round(c)
    avg = data_utils.shuffle_seq(avg, seed=1)

    return avg


#························································································#
#··································· C I D E R ··········································#
#························································································#

def cider_parameters(seqs) -> pd.DataFrame:
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


#························································································#
#······························ T R A J E C T O R Y ·····································#
#························································································#

def compute_rg(seq, traj: md.Trajectory) -> np.ndarray:
    """
    
    Takes a trajectory and the original sequence for a simulation,
    returns the radius of gyration for each frame.

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seq`: `str|list``
            The sequence for the trajectory

        `traj`: `md.Trajectory``
            An OpenMM Trajectory object

    Returns
    -------

        `rg`: `np.ndarray`
            An array of radius of gyration value for each frame

    """

    # Recreating sequence and residue-type data
    seq, res = data_utils.format_terminal_res(seq)

    # Getting molecular weights for sequence
    mass = np.array([res.loc[aa, 'MW'] for aa in seq])

    # Calculate and return a Rg Series
    rg = md.compute_rg(traj, mass)

    return rg

#························································································#
def compute_distances(traj: md.Trajectory) -> tuple:
    """
    
    Takes a trajectory for a simulation,
    returns the distance between all residue pairs for each frame and
    a mask array for non-bonded pairs.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory``
            An OpenMM Trajectory object

    Returns
    -------

        `distances`: `np.ndarray`
            An array of distances between all residue; Dimensions: (`frame`, `pair`)

        `mask`: `np.array`
            A mask array for whether a pair is bonded

    """

    # Selecting all pairs of residues
    pairs = traj.top.select_pairs('all','all')

    # Finding pairs that corresponds to bonds (difference of ID number in [-1,0,1]) (Not used for AH- and DH-potentials)
    mask = np.abs(pairs[:,0]-pairs[:,1]) > 1
    pairs = pairs[mask]

    # Computing pairs between distances over trajectory
    distances = md.compute_distances(traj, pairs).astype(np.float32)

    # Setting threshold for distances for energy calculations
    distance_threshold = 4.
    distances[distances > distance_threshold] = np.inf

    return distances, mask

#························································································#
