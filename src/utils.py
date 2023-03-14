"""
    Utils
    --------------------------------------------------------------------------------

    General utils.

    --------------------------------------------------------------------------------
"""


#························································································#
class log():
    """
    A simple class for logging.
    """

    def __init__(self, write: bool, print: bool, file: str='log.log') -> None:
        """

        Takes arguments for whether to write outputs to file and/or print to stdout, as well as a filepath to write log to.

        --------------------------------------------------------------------------------

        Parameters
        ----------

            `write`: `bool`
                Wether to write log messages to a file (`file`)

            `print`: `bool`
                Wether to print log messages to stdout

            `file`: `str`
                Which file to write log messages to (if `write = True`)

        Methods
        -------

            `message`
                Method for writing messages to log

        """

        # Sets attributes
        self.write = write
        self.print = print
        self.logfile = file
    
    def message(self, message: str, header: bool=False) -> None:
        """

        Takes a message, logs it according to object settings.

        --------------------------------------------------------------------------------

        Parameters
        ----------

            `message`: `str`
                Message to log.

            `header`: `bool`
                Whether to treat the message like a header (new paragraph, underlined)

        """

        # Printing
        if self.print:
            print(message)
        
        # Writing
        if self.write:

            # Writing to log file
            with open(self.logfile, 'a+') as file:

                if header:
                    print('\n')

                file.write(message + '\n')

                if header:
                    file.write('-'*len(message) + '\n')


#························································································#
def read_fasta(path: str, just_seq: bool=False):
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