"""
    Analyse_results
    --------------------------------------------------------------------------------

    Utils for analysing simulation results.

    --------------------------------------------------------------------------------
"""


import pandas as pd
import mdtraj as md
import numpy as np

import simulate_utils


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
    seq, res = simulate_utils.format_terminal_res(seq)

    # Getting molecular weights for sequence
    mass = np.array([res.loc[aa, 'MW'] for aa in seq])

    # Creating mass-frame matrix (I.e. sequence molecular weights for each frame)
    framemass = np.tile(mass, (traj.n_frames, 1))

    # Calculate and return a Rg Series
    rg = md.compute_rg(traj, framemass)

    return rg
