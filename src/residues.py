"""
    Residues
    --------------------------------------------------------------------------------

    Residue data as a Python object.

    --------------------------------------------------------------------------------
"""


import pandas as pd


#························································································#

residues = pd.DataFrame(
[
    {
        "one": "R",
        "three": "ARG",
        "MW": 156.19,
        "lambdas": 0.7307624767517166,
        "sigmas": 0.6559999999999999,
        "q": 1,
        "CALVADOS1": 0.7249915947715212,
        "CALVADOS2": 0.7307624767517166
    },
    {
        "one": "D",
        "three": "ASP",
        "MW": 115.09,
        "lambdas": 0.0416040480605567,
        "sigmas": 0.5579999999999999,
        "q": -1,
        "CALVADOS1": 0.0291821237763497,
        "CALVADOS2": 0.0416040480605567
    },
    {
        "one": "N",
        "three": "ASN",
        "MW": 114.1,
        "lambdas": 0.4255859009787713,
        "sigmas": 0.568,
        "q": 0,
        "CALVADOS1": 0.4383272997027284,
        "CALVADOS2": 0.4255859009787713
    },
    {
        "one": "E",
        "three": "GLU",
        "MW": 129.11,
        "lambdas": 0.0006935460962935,
        "sigmas": 0.5920000000000001,
        "q": -1,
        "CALVADOS1": 0.0061002816086497,
        "CALVADOS2": 0.0006935460962935
    },
    {
        "one": "K",
        "three": "LYS",
        "MW": 128.17,
        "lambdas": 0.1790211738990582,
        "sigmas": 0.636,
        "q": 1,
        "CALVADOS1": 0.0586171731586979,
        "CALVADOS2": 0.1790211738990582
    },
    {
        "one": "H",
        "three": "HIS",
        "MW": 137.14,
        "lambdas": 0.4663667290557992,
        "sigmas": 0.608,
        "q": 0,
        "CALVADOS1": 0.4651948082346978,
        "CALVADOS2": 0.4663667290557992
    },
    {
        "one": "Q",
        "three": "GLN",
        "MW": 128.13,
        "lambdas": 0.3934318551056041,
        "sigmas": 0.602,
        "q": 0,
        "CALVADOS1": 0.3268188050525212,
        "CALVADOS2": 0.3934318551056041
    },
    {
        "one": "S",
        "three": "SER",
        "MW": 87.08,
        "lambdas": 0.4625416811611541,
        "sigmas": 0.518,
        "q": 0,
        "CALVADOS1": 0.4648570130065605,
        "CALVADOS2": 0.4625416811611541
    },
    {
        "one": "C",
        "three": "CYS",
        "MW": 103.14,
        "lambdas": 0.5615435099141777,
        "sigmas": 0.5479999999999999,
        "q": 0,
        "CALVADOS1": 0.610362354303913,
        "CALVADOS2": 0.5615435099141777
    },
    {
        "one": "G",
        "three": "GLY",
        "MW": 57.05,
        "lambdas": 0.7058843733666401,
        "sigmas": 0.45,
        "q": 0,
        "CALVADOS1": 0.7012713677972457,
        "CALVADOS2": 0.7058843733666401
    },
    {
        "one": "T",
        "three": "THR",
        "MW": 101.11,
        "lambdas": 0.3713162976273964,
        "sigmas": 0.562,
        "q": 0,
        "CALVADOS1": 0.5379777613307019,
        "CALVADOS2": 0.3713162976273964
    },
    {
        "one": "A",
        "three": "ALA",
        "MW": 71.07,
        "lambdas": 0.2743297969040348,
        "sigmas": 0.504,
        "q": 0,
        "CALVADOS1": 0.0011162643859539,
        "CALVADOS2": 0.2743297969040348
    },
    {
        "one": "M",
        "three": "MET",
        "MW": 131.2,
        "lambdas": 0.5308481134337497,
        "sigmas": 0.618,
        "q": 0,
        "CALVADOS1": 0.7458993420826714,
        "CALVADOS2": 0.5308481134337497
    },
    {
        "one": "Y",
        "three": "TYR",
        "MW": 163.18,
        "lambdas": 0.9774611449343455,
        "sigmas": 0.6459999999999999,
        "q": 0,
        "CALVADOS1": 0.9950108229594324,
        "CALVADOS2": 0.9774611449343455
    },
    {
        "one": "V",
        "three": "VAL",
        "MW": 99.13,
        "lambdas": 0.2083769608174481,
        "sigmas": 0.5860000000000001,
        "q": 0,
        "CALVADOS1": 0.4185006852559869,
        "CALVADOS2": 0.2083769608174481
    },
    {
        "one": "W",
        "three": "TRP",
        "MW": 186.22,
        "lambdas": 0.9893764740371644,
        "sigmas": 0.6779999999999999,
        "q": 0,
        "CALVADOS1": 0.9844235478393932,
        "CALVADOS2": 0.9893764740371644
    },
    {
        "one": "L",
        "three": "LEU",
        "MW": 113.16,
        "lambdas": 0.6440005007782226,
        "sigmas": 0.618,
        "q": 0,
        "CALVADOS1": 0.5563020305733198,
        "CALVADOS2": 0.6440005007782226
    },
    {
        "one": "I",
        "three": "ILE",
        "MW": 113.16,
        "lambdas": 0.5423623610671892,
        "sigmas": 0.618,
        "q": 0,
        "CALVADOS1": 0.6075268330845265,
        "CALVADOS2": 0.5423623610671892
    },
    {
        "one": "P",
        "three": "PRO",
        "MW": 97.12,
        "lambdas": 0.3593126576364644,
        "sigmas": 0.5559999999999999,
        "q": 0,
        "CALVADOS1": 0.3729641853599348,
        "CALVADOS2": 0.3593126576364644
    },
    {
        "one": "F",
        "three": "PHE",
        "MW": 147.18,
        "AH_lambda": 0.8672358982062975,
        "AH_sigma": 0.636,
        "q": 0,
        "CALVADOS1": 0.9216959832175944,
        "CALVADOS2": 0.8672358982062975
    }
]
)
"""

A DataFrame containing residue data for CALVADOS simulation.

--------------------------------------------------------------------------------

Fields
------

    `one`: `str`
        Residue names (one-letter codes)

    `three`: `str`
        Residue names (three-letter codes)

    `MW`: `float`
        Molecular weight of residue [g/mol]

    `AH_lambda`: `float`
        Lambda parameter in Ashbaugh-Hatch potential interactions quantifying residue hydrophobicity [];
        Default set to CALVADOS2 parameters

    `AH_sigma`: `float`
        Sigma parameter in Ashbaugh-Hatch potential interactions quantifying residue size [nm]
    
    `q`: `int`
        Residue unit charge [e]

    `CALVADOS1`: `float`
        CALVADOS "stickyness" parameter for Ashbaugh-Hatch potential lambda value (Version 1) []

    `CALVADOS1`: `float`
        CALVADOS "stickyness" parameter for Ashbaugh-Hatch potential lambda value (Version 2) []

"""

#························································································#

if __name__ == '__main__':
    print(residues)