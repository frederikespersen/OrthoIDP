"""
    ``PROCESS_DATA``
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

def get_protein_gp(acc_num: str, filename: str, dir='data/seqs/raw', verbose=False) -> str:
    """
    Takes a RefSeq protein accession number, downloads it in genbank format.

    --------------------------------------------------------------------------------

    :param ``acc_num``: Accession number of the RefSeq protein
    :param ``filename``: Filename without .gp suffix
    :param ``dir``: Directory to save file in
    :param ``verbose``: Whether to print file actions
    :return: Path to file

    """

    # Fetching file as genbank record
    handle = Entrez.efetch(db='protein',id=acc_num, rettype='gp', retmode='text')
    record = SeqIO.read(handle, 'genbank')

    # Saving as genbank file
    filepath = "/".join([dir, filename + '.gp'])
    os.makedirs(dir, exist_ok=True)
    with open(filepath, 'w') as file:
        SeqIO.write(record, file, 'genbank')
    if verbose:
        print(f"The RefSeq {acc_num} has been downloaded as '{filepath}'")

    return filepath


#························································································#
#····························· P R E P R O C E S S I N G·································#
#························································································#

def extract_idr_fasta(gp_path: str, i_idr: int=0, length_order=False, fasta_dir='data/seqs/idr', fasta_id='', fasta_desc='', verbose=False) -> tuple[str]:
    """
    Takes a genbank protein file path, extracts a specified IDR region of the protein in FASTA format.

    --------------------------------------------------------------------------------

    :param ``gp_path``: Path to the .gp file
    :param ``i_idr``: The index of the disordered region in the protein [0 = first/NTD; -1 =last/CTD] (See ``length order``though)
    :param ``length_order``: Whether to sort the disordered regions by descending length before choosing with ``i_idr``
    :param ``fasta_dir``: Directory to save file in
    :param ``fasta_id``: The ID to use for the fasta file
    :param ``fasta_desc``: The description to use for the fasta file
    :param ``verbose``: Whether to print file actions
    :return: Path to file and location of IDR

    """

    # Loading genbank file
    with open(gp_path, 'r') as file:
        records = list(SeqIO.parse(file, 'genbank'))

        # Getting record
        assert len(records) != 0, f"No record contained in {gp_path}!"
        assert len(records) == 1, f"More than one record contained in {gp_path}!"
        record = records[0]

    # Finding disordered regions
    idrs = []
    for feature in record.features:
            if feature.type == 'Region':
                if 'Disordered' in feature.qualifiers['region_name'][0]:
                    idrs.append(feature)
    
    # Choosing IDR
    if length_order:
        idrs.sort(key=lambda f: len(f), reverse=True)
    idr_feature = idrs[i_idr]
    idr_seq = idr_feature.extract(record.seq)
    idr_loc = idr_feature.location

    # Checking whether domain is terminal
    idr_pos = 'IDR'
    if 0 in idr_loc:
        idr_pos = "NTD"
    elif len(record.seq)-1 in idr_loc:
        idr_pos = "CTD"
    
    # Creating SeqIO object
    name = gp_path.split('/')[-1][:-3] + '_' + idr_pos
    if fasta_id == '':
        fasta_id = name
    assert ' ' not in fasta_id, f"Fasta identifier must not contain whitespace! {fasta_id}"
    record = SeqRecord(idr_seq, id=str(fasta_id), description=str(fasta_desc))

    # Saving as FASTA file
    filepath = '/'.join([fasta_dir, name + '.fasta'])
    os.makedirs(fasta_dir, exist_ok=True)
    with open(filepath, 'w') as file:
        SeqIO.write(record, file, 'fasta')
    if verbose:
        print(f"The IDR in {idr_loc} of {gp_path} has been extracted to '{filepath}'")

    return filepath, str(idr_loc)


#························································································#
#······································ F I N A L ·······································#
#························································································#

variant_types = {
    "wt": {
        "name": "Wild-type",
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

Schema:

``<variant_id>``

        ``name``: Description of variant

        ``function``: Lambda function for generating variant from sequence and seed

"""


#························································································#
def generate_variant_fasta(fasta_path: str, variant: str, filename: str, dir='data/seqs/var', seed=None) -> str:
    """
    Takes a FASTA IDR file path, 
    generates a variant of the sequence and saves it in FASTA format with a 1-line sequence.
    
    --------------------------------------------------------------------------------

    Possible variants (See ``variant_types``):
    - ``wt``: Wild type
    - ``rand``: Randomly shuffled (See ``shuffle_seq()``)
    - ``clust``: Positive charges clustered in C-terminal end, negative charges in N-terminal end

    --------------------------------------------------------------------------------

    :param fasta_path: Path to the .fasta file
    :param variant: The variant to generate, choose from: ``wt``, ``rand``, or ``clust``
    :param filename: Filename without .fasta suffix
    :param dir: Directory to save file in
    :param seed: Seed for random events
    :return: Path to file

    """

    # Loading FASTA file
    with open(fasta_path, 'r') as file:
        records = list(SeqIO.parse(file, 'fasta'))

        # Getting record
        assert len(records) != 0, f"No record contained in {fasta_path}!"
        assert len(records) == 1, f"More than one record contained in {fasta_path}!"
        record = records[0]

    # Generating variant
    try:
        map_func = variant_types[variant]["function"]
    except KeyError:
        raise ValueError(f"Variant type '{variant}' is not valid, see documentation!")
    seq = map_func(str(record.seq), seed)

    # Saving as (one-line-sequence) FASTA file
    filepath = '/'.join([dir, filename + '.fasta'])
    os.makedirs(dir, exist_ok=True)
    with open(filepath, 'w') as file:
        file.write('>' + filename + '\n')
        file.write(seq)

    return filepath
    

#························································································#
def shuffle_seq(seq: str, seed=None) -> str:
    """
    Takes a sequence, shuffles it randomly.

    --------------------------------------------------------------------------------

    :param seq: Sequence to be shuffled
    :param seed: Seed for random event
    :return: The shuffled sequence

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

    :param seq: Sequence to be clustered
    :param ngroup: List of amino acids to cluster in the N-terminal end of the sequence
    :param cgroup: List of amino acids to cluster in the C-terminal end of the sequence
    :param seed: Seed for random events
    :param mc_threshold: Minimum random value in [0:1] for switching to occur (1 = 100% chance of switching to cluster, i.e. no randomness)
    :return: The clustered sequence

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

