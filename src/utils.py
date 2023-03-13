"""
    Utils
    --------------------------------------------------------------------------------

    General utils.

    --------------------------------------------------------------------------------
"""


#························································································#
def log(message: str, logfile='log.log', header=False) -> None:
    """
    
    Takes a message, writes it to a log file.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `message`: `str`
            A string to output to a log file

        `logfile`: `str`
            The path to the log file to output to

        `header`: `bool``
            Whether the logged message is a header, which will add a line beneath it
    
    """

    # Writing to log file
    with open(logfile, 'a+') as file:
        file.write(message + '\n')

        # Creating line for header
        if header:
            file.write('\n' + '-'*len(message) + '\n')


#························································································#
def read_fasta(path: str, just_seq: bool=False) -> dict|str|list:
    """
    
    Takes a path to a FASTA file, returns the sequence(s) contained herein.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `path`: `str`
            A path to a readable FASTA file

        `just_seq`: `bool`
            Whether to return a dict of ID(s) and sequence(s) (default, `False`) or just the sequence(s) (`True`)

    Returns
    -------

        `seqs`: `dict|str|list`
            The sequences of the FASTA file, either in a dict with ID(s) as key(s) or as s string (one sequence) / list (several sequences)
        
    """

    # Reading file
    with open(path, 'r') as fasta:
        lines = [line.strip() for line in fasta.readlines()]
    
    # Looping over file lines
    seqs = {}
    id = ''
    for line in lines:

        # Finding IDs
        if line[0] == '>':
            # Setting ID for subsequent sequence
            id = line[1:].split(' ')[0]
            seqs[id] = ''

        # Finding sequences
        else:
            # Appending to last ID
            seqs[id] += line
    
    # Formatting results:
    if just_seq:

        # As list
        seqs = list(seqs.values())
        if len(seqs) == 1:
            # As string
            seqs = seqs[0]

    return seqs