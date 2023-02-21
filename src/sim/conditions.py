"""
    Conditions
    --------------------------------------------------------------------------------

    Standard condition settings as a Python object.

    --------------------------------------------------------------------------------
"""


import pandas as pd


#························································································#

conditions = pd.DataFrame(
[
    {
        "name":"default",
        "eps_factor":0.2,
        "temp":293,
        "pH":7.4,
        "ionic":0.2
    }
]
)
"""

A DataFrame containing standard condition setups.

--------------------------------------------------------------------------------

Setups
------
    `default`
        Default conditions.

--------------------------------------------------------------------------------

Fields
------

    `name`: `str`
        Setup name

    `eps_factor`: `float`
        TODO The solvent ??? permittivity [?]
    
    `temp`: `float`
        The absolute temperature [°K]
    
    `pH`: `float`
        The solvent pH []
    
    `ionic`: `float`
        The ionic strength of the solution [mM]

"""