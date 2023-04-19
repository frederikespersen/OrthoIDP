"""
    Conditions
    --------------------------------------------------------------------------------

    Standard condition settings as a Python object.

    --------------------------------------------------------------------------------
"""


import pandas as pd
import argparse


#························································································#

conditions = pd.DataFrame(
[
    {
        "name":"default",
        "eps_factor":0.2,
        "temp":298,
        "pH":7.0,
        "ionic":0.15
    },
    {
        "name":"Borgia_in_silico",
        "eps_factor":0.2,
        "temp":300,
        "pH":6.0,
        "ionic":0.165
    }
]
).set_index('name')
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
        The ionic strength of the solution [M]

"""

#························································································#

if __name__ == '__main__':
    
    # Setting up option for producing .csv file
    parser = argparse.ArgumentParser(prog="Conditions", description="Generates a .csv of the condition options")
    parser.add_argument('--csv',
                        action='store_true',
                        required=False,
                        help="whether to generate a .csv file of the conditions")
    
    # Producing .csv
    if parser.parse_args().csv:
        conditions.to_csv("conditions.csv")

    # Printing
    print(conditions)