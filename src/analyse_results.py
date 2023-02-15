"""
    Analyse_results
    --------------------------------------------------------------------------------

    Utils for analysing simulation results.

    --------------------------------------------------------------------------------
"""


import pandas as pd
import mdtraj as md


#························································································#
#··························· S I M U L A T I O N   S P E C S ····························#
#························································································#

def simulation_time(trajlog_path: str) -> float:
    """
    
    Takes a trajectory log file path, returns the duration of the simulation.
    
    Assumes log file is tab-delimited and uses the ``Elapsed Time (s)``-column.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        ``trajlog_path``: ``str``
            Path to traj.log file

    Returns
    -------
    
        ``time``: ``float``
            The time of the simulation in hours

    """

    # Loading log-file
    log = pd.read_csv(trajlog_path, delimiter='\t')

    # Finding time of last log entry
    time = log['Elapsed Time (s)'].iloc[-1]

    # Returning time in hours
    return time / 3600

