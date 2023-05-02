"""
    Evolve_utils
    --------------------------------------------------------------------------------

    Utils for running Francesco Pesce's IDP evolution algorithm.

    --------------------------------------------------------------------------------
"""


import random
import numpy as np
from localcider.sequenceParameters import SequenceParameters


#························································································#
#································ O B S E R V A B L E S ·································#
#························································································#

structural_measures = {
    'kappa': lambda seq: SequenceParameters(seq).get_kappa()
}
"""

A dictionary containing measures for use in Metropolis criterion.
All measures are set to take a sequence as arguments.

--------------------------------------------------------------------------------

Schema
------

`<measure_id>`: Function to evaluate measure from sequence.

"""


#························································································#
#································ M O N T E   C A R L O ·································#
#························································································#

def swap_sequence(seq: str) -> str:
    """

    Takes a sequence,
    swaps two positions randomly

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
def mc_move(obs_new: float, obs_old: float, target: float, c: float) -> bool:
    """

    Takes the observable of two states, a target, and a control parameter to scale the cost function by,
    returns whether to accept the new state by the Metropolis criterion as well as the acceptance probability.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `obs_new`: `float`
            The observable of the new state

        `obs_old`: `float`
            The observable of the old state

        `target`: `float`
            The target observable given Metropolis criterion
        
         `c`: `float`
            A parameter to scale the cost function by; Influences acceptance rate

    Returns
    ----------

        `accept`: `bool`
            Whether to accept the new state

        `p`: `float`
            The acceptance probability
    
    """
    
    # Calculating cost function
    L = abs(obs_new - target) - abs(obs_old - target)

    # Calculating probability of acceptance by Metropolis criterion
    p = min(1, np.exp(-L/c))

    # Monte Carlo sampling
    if random.random() <= p:

        # Accepting move
        accept = True
    
    else:
        # Rejecting move
        accept = False

    return accept, p
