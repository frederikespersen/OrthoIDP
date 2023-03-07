"""
    Evolve_utils
    --------------------------------------------------------------------------------

    Utils for running Francesco Pesce's IDP evolution algorithm.

    --------------------------------------------------------------------------------
"""


from residues import residues


#························································································#
#································· G E N E R A L  ·······································#
#························································································#

def log(message: str, logfile='log.txt', header=False) -> None:
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
            file.write('-'*len(message) + '\n')
