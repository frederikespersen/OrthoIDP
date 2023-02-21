"""
    Process_data
    --------------------------------------------------------------------------------

    Utils for downloading and preprocessing sequences prior to simulation.

    --------------------------------------------------------------------------------
"""


import os
import random
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = 'tgw325@alumni.ku.dk'


#························································································#
#······································ R A W ···········································#
#························································································#

def get_protein_gp(acc_num: str) -> SeqRecord:
    """

    Takes a protein accession number, like RefSeq or UniProt ID,
    returns the corresponding GenPept file as a SeqRecord object.

    --------------------------------------------------------------------------------
 
    Parameters
    ----------
    
        `acc_num`: `str`
            A protein identifier, like a RefSeq accession number or an UniProt ID


    Returns
    -------

        `record`: `Bio.SeqRecord`
            The GenPept file as a SeqRecord object

    """

    # Fetching file as GenPept record
    handle = Entrez.efetch(db='protein',id=acc_num, rettype='gp', retmode='text')
    record = SeqIO.read(handle, 'genbank')

    return record


#························································································#
#····························· P R E P R O C E S S I N G·································#
#························································································#

def extract_idr(gp, i_idr: int=0, length_order=False) -> tuple[str]:
    """

    Takes a protein GenPept file as a SeqRecord object,
    extracts a specified IDR region of the protein.
    Returns the sequence, region (NTD, INT, CTD), and location of the IDR.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `gp`: `Bio.SeqRecord`
            A GenPept file as a SeqRecord object
    
        `i_idr`: `int`
            The index of the disordered region in the protein
            [0 = first/NTD; -1 =last/CTD]
            (See `length order` below)
    
        `length_order`: `bool`
            Whether to sort the disordered regions by descending length before choosing with `i_idr`

    Returns
    -------
    
        `seq`: `str`
            The extracted IDR sequence

        `region`: `str`
            The general location of the IDR, either N-terminal (NTD), internal (INT), or C-terminal (CTD)

        `location`: `str`
            The positionwise location of the IDR int he format [i:j] (0-indexed)

    """

    # Loading GenPept record
    record = gp

    # Finding disordered regions
    idrs = []
    for feature in record.features:
            if feature.type == 'Region':
                if 'Disordered' in feature.qualifiers['region_name'][0]:
                    idrs.append(feature)
    
    # Choosing IDR
    if length_order:
        idrs.sort(key=lambda f: len(f), reverse=True)
    idr = idrs[i_idr]
    seq = idr.extract(record.seq)
    location = idr.location

    # Checking whether domain is terminal
    region = 'INT'
    if 0 in location:
        region = "NTD"
    elif len(record.seq)-1 in location:
        region = "CTD"

    return str(seq), region, str(location)


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

`<variant_id>`

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

        `seq`: 
            `str`Sequence to be shuffled
    
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
def cluster_seq(seq: str, ngroup: list[str], cgroup: list[str], seed=None, mc_threshold=1.) -> str:
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

