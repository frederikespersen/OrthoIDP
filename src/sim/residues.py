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
        "one":"R",
        "three":"ARG",
        "MW":156.19,
        "AH_lambda":0.7249915948,
        "AH_sigma":0.656,
        "q":1,
        "HPS":0.0,
        "HPSUrry":0.478824,
        "AVG":0.1582837898,
        "M1":0.7249915948,
        "M2":0.8139478491,
        "M3":0.7233336788
    },
    {
        "one":"D",
        "three":"ASP",
        "MW":115.09,
        "AH_lambda":0.0291821238,
        "AH_sigma":0.558,
        "q":-1,
        "HPS":0.378378,
        "HPSUrry":0.214119,
        "AVG":0.2074558136,
        "M1":0.0291821238,
        "M2":0.0731916654,
        "M3":0.0017064814
    },
    {
        "one":"N",
        "three":"ASN",
        "MW":114.1,
        "AH_lambda":0.4383272997,
        "AH_sigma":0.568,
        "q":0,
        "HPS":0.432432,
        "HPSUrry":0.508236,
        "AVG":0.2648094726,
        "M1":0.4383272997,
        "M2":0.0581789695,
        "M3":0.1596122196
    },
    {
        "one":"E",
        "three":"GLU",
        "MW":129.11,
        "AH_lambda":0.0061002816,
        "AH_sigma":0.592,
        "q":-1,
        "HPS":0.459459,
        "HPSUrry":-0.08,
        "AVG":0.2217395834,
        "M1":0.0061002816,
        "M2":0.0042107191,
        "M3":0.0224500142
    },
    {
        "one":"K",
        "three":"LYS",
        "MW":128.17,
        "AH_lambda":0.0586171732,
        "AH_sigma":0.636,
        "q":1,
        "HPS":0.513514,
        "HPSUrry":0.302354,
        "AVG":0.1851426522,
        "M1":0.0586171732,
        "M2":0.1804654043,
        "M3":0.0948078101
    },
    {
        "one":"H",
        "three":"HIS",
        "MW":137.14,
        "AH_lambda":0.4651948082,
        "AH_sigma":0.608,
        "q":0,
        "HPS":0.513514,
        "HPSUrry":0.684707,
        "AVG":0.4091773174,
        "M1":0.4651948082,
        "M2":0.5112064894,
        "M3":0.4869673179
    },
    {
        "one":"Q",
        "three":"GLN",
        "MW":128.13,
        "AH_lambda":0.3268188051,
        "AH_sigma":0.602,
        "q":0,
        "HPS":0.513514,
        "HPSUrry":0.478824,
        "AVG":0.2546409647,
        "M1":0.3268188051,
        "M2":0.4346511118,
        "M3":0.4678384204
    },
    {
        "one":"S",
        "three":"SER",
        "MW":87.08,
        "AH_lambda":0.464857013,
        "AH_sigma":0.518,
        "q":0,
        "HPS":0.594595,
        "HPSUrry":0.508236,
        "AVG":0.3724420532,
        "M1":0.464857013,
        "M2":0.4771931962,
        "M3":0.4872252949
    },
    {
        "one":"C",
        "three":"CYS",
        "MW":103.14,
        "AH_lambda":0.6103623543,
        "AH_sigma":0.548,
        "q":0,
        "HPS":0.594595,
        "HPSUrry":0.56706,
        "AVG":0.7691022569,
        "M1":0.6103623543,
        "M2":0.8476010023,
        "M3":0.3998242232
    },
    {
        "one":"G",
        "three":"GLY",
        "MW":57.05,
        "AH_lambda":0.7012713678,
        "AH_sigma":0.45,
        "q":0,
        "HPS":0.648649,
        "HPSUrry":0.49353,
        "AVG":0.4573157523,
        "M1":0.7012713678,
        "M2":0.7898753832,
        "M3":0.7841272883
    },
    {
        "one":"T",
        "three":"THR",
        "MW":101.11,
        "AH_lambda":0.5379777613,
        "AH_sigma":0.562,
        "q":0,
        "HPS":0.675676,
        "HPSUrry":0.508236,
        "AVG":0.4254603148,
        "M1":0.5379777613,
        "M2":0.2311010548,
        "M3":0.2737417379
    },
    {
        "one":"A",
        "three":"ALA",
        "MW":71.07,
        "AH_lambda":0.0011162644,
        "AH_sigma":0.504,
        "q":0,
        "HPS":0.72973,
        "HPSUrry":0.522942,
        "AVG":0.5331954368,
        "M1":0.0011162644,
        "M2":0.0054961934,
        "M3":0.0030752979
    },
    {
        "one":"M",
        "three":"MET",
        "MW":131.2,
        "AH_lambda":0.7458993421,
        "AH_sigma":0.618,
        "q":0,
        "HPS":0.837838,
        "HPSUrry":0.596471,
        "AVG":0.7255656169,
        "M1":0.7458993421,
        "M2":0.7584982013,
        "M3":0.9928888134
    },
    {
        "one":"Y",
        "three":"TYR",
        "MW":163.18,
        "AH_lambda":0.995010823,
        "AH_sigma":0.646,
        "q":0,
        "HPS":0.864865,
        "HPSUrry":0.817059,
        "AVG":0.6030785548,
        "M1":0.995010823,
        "M2":0.99682662,
        "M3":0.9844421087
    },
    {
        "one":"V",
        "three":"VAL",
        "MW":99.13,
        "AH_lambda":0.4185006853,
        "AH_sigma":0.586,
        "q":0,
        "HPS":0.891892,
        "HPSUrry":0.584707,
        "AVG":0.7797027423,
        "M1":0.4185006853,
        "M2":0.4049719035,
        "M3":0.4277711496
    },
    {
        "one":"W",
        "three":"TRP",
        "MW":186.22,
        "AH_lambda":0.9844235478,
        "AH_sigma":0.678,
        "q":0,
        "HPS":0.945946,
        "HPSUrry":0.92,
        "AVG":0.76512858,
        "M1":0.9844235478,
        "M2":0.9108322527,
        "M3":0.7527633945
    },
    {
        "one":"L",
        "three":"LEU",
        "MW":113.16,
        "AH_lambda":0.5563020306,
        "AH_sigma":0.618,
        "q":0,
        "HPS":0.972973,
        "HPSUrry":0.640589,
        "AVG":0.8027772272,
        "M1":0.5563020306,
        "M2":0.5599323165,
        "M3":0.3351707808
    },
    {
        "one":"I",
        "three":"ILE",
        "MW":113.16,
        "AH_lambda":0.6075268331,
        "AH_sigma":0.618,
        "q":0,
        "HPS":0.972973,
        "HPSUrry":0.625883,
        "AVG":0.8626948823,
        "M1":0.6075268331,
        "M2":0.4009209849,
        "M3":0.6867378047
    },
    {
        "one":"P",
        "three":"PRO",
        "MW":97.12,
        "AH_lambda":0.3729641854,
        "AH_sigma":0.556,
        "q":0,
        "HPS":1.0,
        "HPSUrry":0.678824,
        "AVG":0.4372146457,
        "M1":0.3729641854,
        "M2":0.3724021223,
        "M3":0.4706966794
    },
    {
        "one":"F",
        "three":"PHE",
        "MW":147.18,
        "AH_lambda":0.9216959832,
        "AH_sigma":0.636,
        "q":0,
        "HPS":1.0,
        "HPSUrry":0.74353,
        "AVG":0.8206231056,
        "M1":0.9216959832,
        "M2":0.9028267586,
        "M3":0.8709039993
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
        Switched out by CALVADOS parameters M1, M2, or M3

    `AH_sigma`: `float`
        Sigma parameter in Ashbaugh-Hatch potential interactions quantifying residue size [nm]
    
    `q`: `int`
        Residue charge [e]

    `HPS`: `float`
        Lambda parameter in HPS model []

    `HPSUrry`: `float`
        Lambda parameter in HPS-Urry model []

    `AVG`: `float`
        TODO ??? [?]

    `M1`: `float`
        CALVADOS "stickyness" parameter for Ashbaugh-Hatch potential lambda value (Version 1) []

    `M2`: `float`
        CALVADOS "stickyness" parameter for Ashbaugh-Hatch potential lambda value (Version 2) []

    `M3`: `float`
        CALVADOS "stickyness" parameter for Ashbaugh-Hatch potential lambda value (Version 3) []

"""
