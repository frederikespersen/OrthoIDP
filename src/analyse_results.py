"""
    Analyse_results
    --------------------------------------------------------------------------------

    Utils for analysing simulation results.

    --------------------------------------------------------------------------------
"""


import pandas as pd
import mdtraj as md
import numpy as np

from sim.residues import residues

#························································································#
#······························ T R A J E C T O R Y ·····································#
#························································································#

def calculate_rg(seq: str|list, traj: md.Trajectory) -> np.ndarray:
    """
    
    Takes a trajectory and the original sequence for the simulation,
    returns a the radius of gyration for each frame.

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
    seq, res = format_terminal_res(seq, residues.copy())

    # Getting molecular weights for sequence
    mass = np.array([res.loc[aa, 'MW'] for aa in seq])

    # Creating mass-frame matrix (I.e. sequence molecular weights for each frame)
    framemass = np.tile(mass, (traj.n_frames, 1))

    # Calculate and return a Rg Series
    rg = md.compute_rg(traj, framemass)

    return rg


#························································································#
# FIXME Temporary: Import error when loading from simulate_utils because of improper simtk installation
def format_terminal_res(seq: str|list, res: pd.DataFrame):
    """
    
    Takes a sequence and a `residues` DataFrame, modifies the sequence with special terminal residue types 'X' and 'Z'
    for the N- and C-terminal respectively.
    Returns the modified sequence and the modified `residues` DataFrame.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str|list`
            An amino acid sequence

        `res`: `pandas.DataFrame`
            A `residues` DataFrame

    Returns
    -------

        `seq`: `str`
            The modified sequence with terminal 'X'/'Z' residues

        `res`: `pandas.DataFrame`
            A modified `residues` DataFrame with 'X'/'Z' residue types

    """

    # Gettning standard residue data
    res.set_index('one', inplace=True)

    # Adding new residue types, and using original terminal residues as templates
    res.loc['X'] = res.loc[seq[0]].copy()
    res.loc['Z'] = res.loc[seq[-1]].copy()
    res.loc['X','MW'] += 2
    res.loc['Z','MW'] += 16
    res.loc['X','q'] += 1
    res.loc['Z','q'] -= 1

    # Modfiying sequence
    seq = list(seq)
    seq[0] = 'X'
    seq[-1] = 'Z'
    seq = ''.join(seq)

    return seq, res