"""
    Data_utils
    --------------------------------------------------------------------------------

    Utils for downloading and preprocessing sequences prior to simulation.

    --------------------------------------------------------------------------------
"""


import requests
import io
import random
import pandas as pd
from Bio import SeqIO
from Bio import SwissProt
from Bio.SwissProt import Record

from residues import residues
import analyse_utils


#························································································#
#······························· D O W N L O A D I N G ··································#
#························································································#

def get_protein_xml_record(uniprot_id: str) -> Record:
    """

     Takes a protein UniProt ID,
     returns a SeqRecord object of the corresponding UniProt .xml file.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `uniprot_id`: `str`
            A UniProt ID
    
    Returns
    -------

        `record`: `Bio.SeqRecord`
            A UniProt .xml file as a ecord object retrived with the `uniprot_id`

    """

    # Fetching UniProt .xml record
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml"
    response = requests.get(url)
    content = io.StringIO(response.content.decode("utf-8"))
    record = SeqIO.read(content, 'uniprot-xml')

    return record


#························································································#
def get_protein_txt_record(uniprot_id: str) -> Record:
    """

     Takes a protein UniProt ID,
     returns a Bio.SwissProt.Record object of the corresponding UniProt .txt (SwissProt) file.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `uniprot_id`: `str`
            A UniProt ID
    
    Returns
    -------

        `record`: `Bio.SwissProt.Record`
            A UniProt .txt (SwissProt) file as a ecord object retrived with the `uniprot_id`

    """

    # Fetching UniProt .txt (SwissProt) record
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    response = requests.get(url)
    content = io.StringIO(response.content.decode("utf-8"))
    record = list(SwissProt.parse(content))[0]

    return record


#························································································#
def get_protein_metadata(uniprot_id: str) -> tuple:
    """
    
    Takes a protein UniProt ID,
    Returns the UniProt data on the name, description, species of origin, and full-length sequence of the protein.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `uniprot_id`: `str`
            A UniProt ID

    Returns
    -------

        `name`: `str`
            The UniProt name of the protein

        `description`: `str`
            The UniProt description of the protein

        `organism`: `str`
            The species of origin (organism) of the protein

        `sequence`: `str``
            The protein sequence

    """
    
    # Fetching UniProt XML record
    record = get_protein_xml_record(uniprot_id)

    # Extracting metadata
    name = record.name
    description = record.description
    organism = record.annotations["organism"]
    sequence = str(record.seq)

    return name, description, organism, sequence


#························································································#
def get_protein_idr(uniprot_id: str, i_idr: int=0, length_order=False, region_threshold: int=0) -> str:
    """

    Takes an intrinsically disordered protein (IDP) UniProt ID,
    extracts a specified intrinsically disordered region (IDR) of the protein.
    Returns the sequence, location, and region (NTD, INT, CTD) of the IDR.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `uniprot_id`: `str`
            A UniProt ID, corresponding to a protein with an IDR
    
        `i_idr`: `int`
            The index of the disordered region in the protein
            [0 = first/NTD; -1 =last/CTD]
            (See `length order` below)
    
        `length_order`: `bool`
            Whether to sort the disordered regions by descending length before choosing with `i_idr`

        `region_threshold`: `int`
            The difference in residue positions between N/C-terminal and region end to classify NTD/CTD;
            Default 4 (I.e. within 0-4 residues of NTD/CTD residue)

    Returns
    -------
    
        `seq`: `str`
            The extracted IDR sequence

        `location`: `str`
            The positionwise location of the IDR int the sequence [i:j] (1-indexed)

        `region`: `str`
            The general location of the IDR, either N-terminal (NTD), internal (INT), or C-terminal (CTD)

    """

    # Fetching UniProt XML record
    record = get_protein_txt_record(uniprot_id)

    # Finding disordered regions
    idrs = []
    for feature in record.features:
        if 'Disordered' in feature.qualifiers.values():
            start = feature.location.start
            end = feature.location.end
            idrs.append((start, end))
    if not idrs:
        raise ValueError(f"No disordered regions were found for UniProt ID {uniprot_id}!")

    # Choosing IDR
    if length_order:
        idrs.sort(key=lambda p: p[1] - p[0], reverse=True)
    idr = idrs[i_idr]

    # Parsing sequence, location, and description
    seq = str(record.sequence[idr[0]:idr[1]])
    location = f'{idr[0]}:{idr[1]}'

    # Describing whether region is terminal
    region = 'INT'
    if idr[0] <= region_threshold:
        region = "NTD"
    elif len(record.sequence) - idr[1] <= region_threshold:
        region = "CTD"

    return seq, location, region


#························································································#
#····························· P R E P R O C E S S I N G·································#
#························································································#

variant_types = {
    "wt": {
        "name": "Wild type",
        "function": lambda seq, seed: seq    
    },
    "rand": {
        "name": "Randomly shuffled",
        "function": lambda seq, seed: shuffle_seq(seq, seed)},
    "clust": {
        "name": "Terminally clustered charges",
        "function": lambda seq, seed: cluster_seq(seq, ['K', 'R'], ['D', 'E'], seed)}
}
"""

A dictionary containing descriptions and mapping functions for generating variants.

--------------------------------------------------------------------------------

Schema
------

`<variant_id>`: Short name of variant

        `name`: Description of variant

        `function`: Lambda function for generating variant from sequence and seed

"""


#························································································#
def generate_variant(seq: str, variant: str, seed=None) -> str:
    """

    Takes a sequence, 
    generates and returns a variant of the sequence.
    
    --------------------------------------------------------------------------------

    Possible variants
    -----------------
    (See `variant_types`):
    - `wt`: Wild type
    - `rand`: Randomly shuffled
    - `clust`: Positive charges clustered in C-terminal end, negative charges in N-terminal end

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str`
            A protein sequence
    
        `variant`: `str`
            The variant to generate; 
            See `variant_types` for options
    
        `seed`: `int`
            Seed for random events
    
    Returns
    -------
        
        `seq`: `str`
            The variant sequence

    """

    # Getting variant transformation
    try:
        map_func = variant_types[variant]["function"]
    except KeyError:
        raise ValueError(f"Variant type '{variant}' is not valid, see documentation!")

    # Applying transformation
    seq = map_func(seq, seed)

    return seq
    

#························································································#
def shuffle_seq(seq: str, seed=None) -> str:
    """

    Takes a sequence, shuffles it randomly.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str`
            Sequence to be shuffled
    
        `seed`: `int`
            Seed for random event
    
    Returns
    -------
        
        `seq`: `str`
            The shuffled sequence

    """
    
    # Prepping
    seq = list(seq)

    # Shuffling
    random.seed(seed)
    random.shuffle(seq)

    # Joining
    seq = ''.join(seq)

    return seq


#························································································#
def cluster_seq(seq: str, ngroup: list, cgroup: list, seed=None, mc_threshold=1.) -> str:
    """

    Takes a sequence, clusters two groups of distinct amino acids in oppposite 
    ends of a sequence by switiching residue positions in the sequence.
    
    Optional parameter for whether or not to do Monte Carlo criteria for switching.

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seq`: `str`
            Sequence to be clustered
    
        `ngroup`: `list[str]`
            List of amino acids (one-letter codes) to cluster in the N-terminal end of the sequence
    
        `cgroup`: `list[str]`
            List of amino acids (one-letter codes) to cluster in the C-terminal end of the sequence
    
        `seed`: `int`
            Seed for random events
    
        `mc_threshold`: `float`
            Minimum random value in [0:1] for switching to occur
            (1 = 100% chance of switching to cluster, i.e. no randomness)
    
    Returns
    -------
    
        `seq`: `str`
            The clustered sequence

    """

    # Prepping
    random.seed(seed)
    seq = list(seq)
    assert (mc_threshold <= 1) and (mc_threshold > 0), "Monte Carlo criteria must be 0 < x <= 1!"

    # Finding positions of N-terminal group residues
    i_ngroup = [i for i, aa in enumerate(seq) if aa in ngroup]
    
    # Loop over sequence for N-terminal group
    for i, aa in enumerate(seq):

        # If all downstream group residues are sorted, stop
        if len(i_ngroup) == 0:
            break

        # If N-terminal group residue, pass and remove group residue from indexed list (i.e. 'downsteam residues')
        if aa in ngroup:
            i_ngroup.remove(i)

        # If not, switch with downstream N-terminal group residue
        else:
            
            # Monte Carlo criteria
            if random.random() < mc_threshold:

                # Switch residue with random downstream group residue
                j = random.choice(i_ngroup)
                seq[i], seq[j] = seq[j], seq[i]
                i_ngroup.remove(j)
        
    # Reversing sequence and repeating for C-terminal group residues
    seq.reverse()
    i_cgroup = [i for i, aa in enumerate(seq) if aa in cgroup]

    # Reverse sequence and loop over it for C-terminal group
    for i, aa in enumerate(seq):

        # If all downstream group residues are sorted, stop
        if len(i_cgroup) == 0:
            break

        # If group residue, pass and remove group residue from indexed list (i.e. 'downsteam residues')
        if aa in cgroup:
            i_cgroup.remove(i)

        # If not, switch with downstream group residue
        else:
            
            # Monte Carlo criteria
            if random.random() < mc_threshold:

                # Switch residue with random downstream group residue
                j = random.choice(i_cgroup)
                seq[i], seq[j] = seq[j], seq[i]
                i_cgroup.remove(j)

    # Re-reversing sequence
    seq.reverse()
    seq = ''.join(seq)

    return seq

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
    freqs = analyse_utils.amino_acid_content(seqs).mean()

    # Finding average length
    N = seqs.str.len().mean()

    # Finding average amino acid counts
    counts = freqs * N

    # Assembling random sequence using counts
    avg = ''
    for aa, c in counts.items():
        avg += aa * round(c)
    avg = shuffle_seq(avg, seed=1)

    return avg