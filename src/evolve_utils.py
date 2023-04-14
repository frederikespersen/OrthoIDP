"""
    Evolve_utils
    --------------------------------------------------------------------------------

    Utils for running Francesco Pesce's IDP evolution algorithm.

    --------------------------------------------------------------------------------
"""


import random
import mdtraj as md
import numpy as np
import pandas as pd
import pymbar

import analyse_utils
import simulate_utils
from simulate_utils import format_terminal_res
from conditions import conditions


#························································································#
#···································· G E N E R A L ·····································#
#························································································#

def load_g_trajectory(g: int) -> md.Trajectory:
    """
    
    Takes the number of a simulated generation,
    returns the MDTRaj trajectory of the simulation.
    
    --------------------------------------------------------------------------------

    Parameters
    ----------

        `g`: `int`
            The number of the generation to load a trajectory for;
            Must be simulated generation

    Returns
    ----------

        `traj`: `md.Trajectory`
            The simulation trajectory
    
    """

    # Defining simulation directory
    dir = f'g{g}'

    # Checking whether exists

    # Loading trajectory
    try:
        traj = md.load_dcd(dir+'traj.dcd', dir+'top.pdb')
    except OSError:
        raise ValueError(f"Generation {g} does not have simulation data!")

    return traj

#························································································#
#································ O B S E R V A B L E S ·································#
#························································································#

structural_measures = {
    'Rg': lambda seq, traj: analyse_utils.compute_rg(seq, traj).mean(),
    'nu': lambda seq, traj: analyse_utils.compute_scaling_exponent(traj)
}
"""

A dictionary containing structural measures for use in Metropolis criterion.
All measures are set to take a sequence and trajectory as arguments.

--------------------------------------------------------------------------------

Schema
------

`<measure_id>`: Function to evaluate measure from sequence and OpenMM trajectory.

"""


#························································································#
#································ M O N T E   C A R L O ·································#
#························································································#

def swap_sequence(seq: str) -> str:
    """

    Takes a sequence,
    swaps to positions randomly

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str`
            The sequence to swap

    Returns
    ----------

        `swap`: `float`
            The swapped sequence
    
    """

    # Converting to list
    swap = list(seq)

    # Choosing two random positions
    pos1 = random.randint(0,len(seq)-1)
    pos2 = random.randint(0,len(seq)-1)

    # Swapping sequences
    swap[pos1], swap[pos2] = seq[pos2], seq[pos1]

    # Converting back to string
    swap = ''.join(swap)

    return swap


#························································································#
#································ R E W E I G H T I N G ·································#
#························································································#

def compute_MBAR(seq, pool, cond='default'):
    """
    
    Takes a new sequence, previous trajectories and the energy of their frames,
    returns MBAR weights and the estimated number of effective frames.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `seq`: `str|list`
            A new sequence to be reweighted against previous trajectories

        `pool`: `pd.DataFrame`
            A DataFrame containing one column 'traj' and fields representing total energy per frame;
            Row indices represent the generation of the sequence used for calculating energy, and
            column names represent the generation of the trajectory used for calculating energy.

        `cond`: `str`
            The standard conditions that the simulations were run with; see `conditions` for choices.

    Returns
    ----------

        `weights`: `np.ndarray[float]`
            The weight of each of the trajectories for reweighting

        `eff_frames`: `float`
            The effective number of independent frames

    
    """
    # Getting energies from previous simulations
    Etots = pool.drop(columns='traj').to_numpy()
    Etots = np.array([np.concatenate(row) for row in Etots])

    # Calculating the energy of the sequence using previous trajectories
    Etot = []
    for traj in pool.traj:
        Etot.append(analyse_utils.compute_energy(seq, traj, cond).sum(axis=0))
    Etot = np.array(Etot).flatten()
    Etots = np.vstack([Etots, Etot])

    # Normalising energies by RT
    R = 8.3145 * 1e-3 # kJ/mol·K
    T = conditions.loc[cond].temp # K
    Etots /= (R * T)

    # Calculating MBAR weights
    N_frame = len(pool.iloc[-1,-1])
    N_frames = np.full((len(pool)), N_frame)
    mbar = pymbar.MBAR(Etots, N_frames)
    weights = mbar.W_nk[...,-1]

    # Calculating the effective number of frames (Using Kullback-Leibler divergence)
    norm_weights = weights / weights.sum()
    norm_weights = norm_weights[np.where(norm_weights > 1e-50)]
    rel_entropy = np.sum(norm_weights * np.log(norm_weights * weights.size))
    rel_eff_frames = np.exp(-rel_entropy)
    eff_frames = rel_eff_frames * N_frames.sum()

    # Formatting weights with generation indices
    weights = pd.Series(weights, index=pool.index)

    return weights, eff_frames
