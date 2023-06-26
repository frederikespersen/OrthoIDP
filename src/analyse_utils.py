"""
    Analyse_utils
    --------------------------------------------------------------------------------

    Utils for analysis.

    --------------------------------------------------------------------------------
"""


import json
import itertools
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from localcider.sequenceParameters import SequenceParameters
import mdtraj as md
import numpy as np

import simulate_utils
from conditions import conditions

plt.rcParams["font.family"] = "STIXGeneral"


#························································································#
#································· G E N E R A L ········································#
#························································································#

def load_fasta_seq(fasta_path: str) -> tuple:
    """

    Takes a FASTA file path, returns the sequence, id, and description.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        fasta_path: str
            Path to .fasta file

    Returns
    -------

        `seq`: `str`
            The sequence in the FASTA file

        `id`: `str`
            The id in the FASTA file

        `desc`: `str`
            The description in the FASTA file

    """

    # Loading file
    with open(fasta_path, 'r') as file:
        lines = file.readlines()

    # Removing whitespace
    lines = [line.strip() for line in lines]

    # Assembling data (Whether standard FASTA format or one-line-sequence FASTA format)
    seq = ''.join(lines[1:])
    id = lines[0].split(' ')[0][1:]
    desc = ' '.join(lines[0].split(' ')[1:])

    return seq, id, desc


#························································································#
def load_metadata(metadata_path: str, join=True) -> pd.DataFrame:
    """
    
    Takes the path to a meta data .json-file, returns the meta data loaded into a DataFrame.

    Row granularity is variants in 'data' field of .json file.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `metadata_path`: `str`
            The path to the meta deta .json file to load in

        `join`: `bool`
            Whether to join template metadata to data metadata;
            Only viable if N:1 cardinality for data:template

    Returns
    -------

        `metadata`: `pandas.DataFrame`
            A DataFrame containing the meta data fields
            (If `join=True`, a tuple of two dataframes (data; templates) are returned instead)

    """

    # Loading meta data file
    with open(metadata_path, 'r') as file:
        metadata = json.load(file)
    
    # Parsing data and templates fields
    data = pd.DataFrame(metadata['data']).transpose()
    templates = pd.DataFrame(metadata['templates']).transpose()

    # Joining templates on to data (unique fields only)
    if join:
        unique_fields = list(set(templates) - set(data))
        metadata = data.join(templates[unique_fields], on='template')
    else:
        metadata = (data, templates)

    return metadata


#························································································#
#····························· A M I N O   A C I D S ····································#
#························································································#

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
"""

A list of the one-letter codes of the naturally occuring amino acids.

"""


#························································································#
amino_acid_types = pd.Series({
    'A': 'Hydrophobic',
    'C': 'Polar',
    'D': 'Negative',
    'E': 'Negative',
    'F': 'Hydrophobic',
    'G': 'Special',
    'H': 'Polar',
    'I': 'Hydrophobic',
    'K': 'Positive',
    'L': 'Hydrophobic',
    'M': 'Hydrophobic',
    'N': 'Polar',
    'P': 'Special',
    'Q': 'Polar',
    'R': 'Positive',
    'S': 'Polar',
    'T': 'Polar',
    'V': 'Hydrophobic',
    'W': 'Hydrophobic',
    'Y': 'Polar'
})
"""

A pd.Series of the general types (Hydrophobic, Polar, Positive, Negative, Special) of amino acids.

--------------------------------------------------------------------------------

Schema
------

    `<AA>`: One-letter code for amino acid

        `<type>`: The general type of the amino acid


"""


#························································································#
def amino_acid_content(seqs) -> pd.DataFrame:
    """

    Takes one or more sequences, returns a DataFrame of the frequencies of amino acids.

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seqs`: `str | list | pandas.Series`
            Sequence(s) to calculate frequencies for

    Returns
    -------

        `freqs`: `pandas.DataFrame`
            A DataFrame with the amino acid frequency for each amino acid in the sequence(s).

    """

    # Formatting sequence(s) as a pd.Series
    if type(seqs) == str:
        seqs = pd.Series([seqs])
    elif type(seqs) == list:
        seqs = pd.Series(seqs)
    
    # Initiating DataFrame
    freqs = pd.DataFrame(index=seqs.index)

    # Calculating frequencies
    seqs = seqs.str.upper()
    for aa in amino_acids:
        freqs[aa] = seqs.apply(lambda seq: len(list(filter(lambda s: s == aa, seq))) / len(seq))

    return freqs


#························································································#
#··································· C I D E R ··········································#
#························································································#

def cider_parameters(seqs) -> pd.DataFrame:
    """

    Takes one or more sequences, returns a DataFrame of CIDER parameters.

    --------------------------------------------------------------------------------

    More on CIDER from PappuLab:
    - [CIDER](http://pappulab.wustl.edu/CIDER/about/)
    - [localCIDER](http://pappulab.github.io/localCIDER/)

    --------------------------------------------------------------------------------

    Parameters
    ----------
    
        `seqs`: `str | list | pandas.Series`
            Sequence(s) to calculate parameters for (`str`, `list` interpreted single sequence)

    Returns
    -------

        `params`: `pandas.DataFrame`
            A DataFrame with select CIDER parameters calcualted for the sequence(s)

    """

    # Formatting sequence(s) as a pd.Series
    if type(seqs) == str:
        seqs = pd.Series([seqs])
    elif type(seqs) == list:
        seqs = pd.Series([''.join(seqs)])
    
    # Mapping sequence(s) to a SequenceParameteres object
    seqs = seqs.map(SequenceParameters)

    # Initiating DataFrame
    params = pd.DataFrame(index=seqs.index)

    # Calculating parameters
    params['kappa'] = seqs.apply(lambda SeqOb: SeqOb.get_kappa())
    params['FCR'] = seqs.apply(lambda SeqOb: SeqOb.get_FCR())
    params['NCPR'] = seqs.apply(lambda SeqOb: SeqOb.get_NCPR())
    params['Hydrophobicity'] = seqs.apply(lambda SeqOb: SeqOb.get_mean_hydropathy())
    params['Frac. dis. prom.'] = seqs.apply(lambda SeqOb: SeqOb.get_fraction_disorder_promoting())

    return params


#························································································#
#······························ T R A J E C T O R Y ·····································#
#························································································#

def log_duration(log_path: str) -> float:
    """
    
    Takes the path to a simulation log file,
    returns the wall (real) time duration of the simulation.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `log_path`: `str`
            The path to the traj.log file of the simulation

    Returns
    ----------

        `duration`: `float`
            The duration of the simulation in wall time [h]
    
    """

    # Reading log file
    log = pd.read_csv(log_path, delimiter='\t')
    
    # Finding the time stamp of the last step
    duration = log["Elapsed Time (s)"].iloc[-1]

    # Converting to hours
    duration /= 3600

    return duration


#························································································#

def compute_energy(traj: md.Trajectory, cond='default', ah=True, dh=True, hb=False, pairs_ij=None, distances=None) -> np.ndarray:
    """
    
    Takes a trajectory,
    returns the energy of each residue pair of frame calculated using the CALVADOS model.
    Includes options to choose which potentials to calculate.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

        `cond`: `str`
            The standard conditions to run the simulation with; see `conditions` for choices.
        
        `ah`: `bool`
            Whether to calculate Ashbaugh-Hatch potential (AH)

        `db`: `bool`
            Whether to calculate Debye-Hückel potential (DH)

        `hb`: `bool`
            Whether to calculate harmonic bond potential (HB)

        `pairs_ij`: `np.darray`
            An array specifying the residue IDs to calculate the energy for (Created using `md.Topology.select_pairs()`);
            Default is all residue pairs.
            
        `distances`: `np.darray`
            Precalculated distances from md.distances();
            Default is to calculate from scratch
            
    Returns
    ----------

        `E`: `np.ndarray[float]; [n_potentials: n_frames: n_pairs]`
            The energy of each pair in each frame of the trajectory evaluated using the CALVADOS model;
            Order of potentials is always: AH, DH, HB

    """

    # Deriving sequence used for simulation
    residues = simulate_utils.extract_sequences(traj.topology)

    # Defining pairs
    if pairs_ij is None:
        pairs_ij = traj.top.select_pairs('all','all')

    # Calculating pairwise distances
    if distances is None:
        distances = md.compute_distances(traj, pairs_ij)


    # Defining energy calculation functions
    def calc_lj(r, s, e):
        return 4*e*((s/r)**12-(s/r)**6)
    def calc_ah(r, s, l, e):
        return np.where(r<=s*np.power(2,1/6), calc_lj(r, s, e)+e*(1-l), l*calc_lj(r, s, e))
    def calc_dh(r, q, ye, yk):
        return q*ye*(np.exp(-yk*r)/r)
    def calc_hb(r, k, r_0):
        return 1/2 * k * (r - r_0)**2
    

    # Preparing global parameters
    cond = conditions.loc[cond]
    r_0=0.38
    k=8033
    e = simulate_utils.ah_parameters(cond.eps_factor)
    yukawa_kappa, yukawa_epsilon = simulate_utils.dh_parameters(cond.temp, cond.ionic)

    # Preparing residue pair specific parameters
    s = (residues.AH_sigma.loc[pairs_ij[:, 0]].to_numpy() + residues.AH_sigma.loc[pairs_ij[:, 1]].to_numpy())/2
    l = (residues.AH_lambda.loc[pairs_ij[:, 0]].to_numpy() + residues.AH_lambda.loc[pairs_ij[:, 1]].to_numpy())/2
    q = (residues.q.loc[pairs_ij[:, 0]].to_numpy() * residues.q.loc[pairs_ij[:, 1]].to_numpy())
    chain = residues.chain.to_numpy()
    bonded = ((chain[pairs_ij[:,1]] == chain[pairs_ij[:,0]]) & (pairs_ij[:,1] - pairs_ij[:,0] == 1))

    # Formatting parameters to trajectory
    r = distances
    s = s[np.newaxis, :]
    l = l[np.newaxis, :]
    q = q[np.newaxis, :]
    bonded = bonded[np.newaxis, :]

    # Calculating energies
    energies = []
    if ah:
        ah = (calc_ah(r, s, l, e) - calc_ah(simulate_utils.AH_cutoff, s, l, e))
        ah = np.where(~bonded, ah, 0)
        energies.append(ah)
    if dh:
        dh = (calc_dh(r, q, yukawa_epsilon, yukawa_kappa) - calc_dh(simulate_utils.DH_cutoff, q, yukawa_epsilon, yukawa_kappa))
        dh = np.where(~bonded, dh, 0)
        energies.append(dh)
    if hb:
        hb = calc_hb(r, k, r_0)
        hb = np.where(bonded, hb, 0)
        energies.append(hb)
    E = np.stack(energies)

    return E

#························································································#

def compute_com(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the center of mass (COM) of each frame in the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `com`: `np.ndarray[float]`
            The center of mass for each frame in the trajectory

    """

    # Deriving sequence used for simulation
    residues = simulate_utils.extract_sequences(traj.topology)

    # Computing center of mass for each frame
    masses = residues.MW.values
    com = np.sum(traj.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()

    return com

#························································································#

def compute_gyration_tensor(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the gyration tensor of each frame in the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `q`: `np.ndarray[float]`
            The gyration tensor for each frame in the trajectory

    """

    # Deriving sequence used for simulation
    residues = simulate_utils.extract_sequences(traj.topology)

    # Computing center of mass for each frame
    com = compute_com(traj)

    # Computing gyration tensor for each frame
    masses = residues.MW.values
    si = traj.xyz - com[:,np.newaxis,:]
    q = np.einsum('jim,jin->jmn', si*masses.reshape(1,-1,1),si)/masses.sum()

    return q

#························································································#

def compute_end_to_end(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the end-to-end distance of each frame in the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `Re`: `np.ndarray[float]`
            The end-to-end distance for each frame in the trajectory

    """

    # Computing end-to-end distance for each frame
    Re = md.compute_distances(traj, [[0, traj.n_atoms-1]])

    return Re

#························································································#

def compute_rg(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the radius of gyration for each frame in the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `Rg`: `np.ndarray[float]`
            The radius of gyration for each frame in the trajectory

    """

    # Computing gyration tensor for each frame
    q = compute_gyration_tensor(traj)

    # Computing radius of gyration for each frame:
    tr_q = np.trace(q, axis1=1, axis2=2)
    Rg = np.sqrt(tr_q)

    return Rg

#························································································#

def compute_asphericity(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the asphericity (as defined in Aronovitz & Nelson 1986) averaged over the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `Delta`: `np.ndarray[float]`
            The asphericity for each frame in the trajectory

    """

    # Computing gyration tensor for each frame
    q = compute_gyration_tensor(traj)

    # Computing traceless gyration tensor for each frame
    tr_q = np.trace(q, axis1=1, axis2=2)
    tr_q_mean = tr_q/3
    q_hat = q - tr_q_mean.reshape(-1,1,1)*np.identity(3).reshape(-1,3,3)

    # Computing asphericity for each frame
    tr_q_hat_sq = np.trace(q_hat**2, axis1=1, axis2=2)
    Delta = 3/2*tr_q_hat_sq.mean()/(tr_q.mean()**2)

    return Delta

#························································································#

def compute_prolateness(traj: md.Trajectory):
    """
    
    Takes a trajectory,
    computes the prolateness (as defined in Aronovitz & Nelson 1986) averaged over the trajectory.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object

    Returns
    ----------

        `S`: `np.ndarray[float]`
            The prolateness for each frame in the trajectory

    """

    # Computing gyration tensor for each frame
    q = compute_gyration_tensor(traj)

    # Computing traceless gyration tensor for each frame
    tr_q = np.trace(q, axis1=1, axis2=2)
    tr_q_mean = tr_q/3
    q_hat = q - tr_q_mean.reshape(-1,1,1)*np.identity(3).reshape(-1,3,3)

    # Computing prolateness for each frame
    S = 27*np.linalg.det(q_hat).mean()/(tr_q.mean()**3)

    return S

#························································································#

def compute_scaling_exponent(traj: md.Trajectory, r0_fix: float=0.518, ij_cutoff=5, plot=False) -> tuple:
    """
    
    Takes a trajectory for a simulation,
    returns the fitted scaling exponent,
    which describes the power-law relationsship between residue pair sequence distance and cartesian distance.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            An OpenMM Trajectory object
        
        `r0_fix`: `float`
            A fixed r0 parameter for fitting scaling exponent; If set to `None`, the parameter will be fitted instead.

        `ij_cutoff`: `int`
            The minimum interresidue sequence distance |i j| for data to have to use for fitting model.

        `plot`: `bool`
            Whether to plot the fit.

    Returns
    -------

        `v`: `float`
            The fitted scaling exponent []

        `v_err`: `float`
            The error on the fitted scaling exponent []

        `r0`: `float``
            The (optionally) fitted r0 [nm]

        `r0_err`: `float`
            The error on the fitted r0 [nm]

    """

    # Finding pairs and calculating interresidue sequence distances
    pairs = traj.top.select_pairs('all','all')
    ij = pairs[:,1] - pairs[:,0]

    # Calculating interresidue cartesian distances as quadratic mean
    d = md.compute_distances(traj, pairs)
    d_mean = np.sqrt(np.square(d).mean(axis=0))
    d_mean_err = d.std(axis=0, ddof=0)/np.sqrt(d.shape[0])

    # Defining model (dependent on fixed r0 or not)
    if r0_fix:
        model = lambda x, v: r0_fix * np.power(x,v)
        p0 = [0.5]
    else:
        model = lambda x, v, r0: r0 * np.power(x,v)
        p0 = [0.5, 0.68]
    
    # Fitting model
    mask = ij>ij_cutoff
    popt, pcov = curve_fit(model, ij[mask], d_mean[mask], p0=p0, sigma=d_mean_err[mask])

    # Extracting fitted value(s)
    if r0_fix:
        v, r0 = popt[0], r0_fix
        v_err, r0_err = np.sqrt(np.diag(pcov))[0], None
    else:
        v, r0 = popt
        v_err, r0_err = np.sqrt(np.diag(pcov))
    
    # Plotting
    if plot:
        plt.scatter(ij, d_mean, alpha=0.01, color='darkred')
        plt.plot(np.unique(ij), model(np.unique(ij), *popt), c='black')
        plt.xlabel("$| i - j |$")
        plt.ylabel("$\sqrt{\ \overline{{r_{i,j}}^2}\ }$")
    
    return v, v_err, r0, r0_err

#························································································#

def compact_frame(traj) -> md.Trajectory:
    """
    
    Takes a trajectory,
    returns the one frame of the trajectory with the lowest Rg (Most compact).

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj`: `md.Trajectory`
            A MDTraj trajectory

    Returns
    -------

        `compact_frame`: `md.Trajectory`
            An OpenMM Trajectory object containing one frame, corresponding to the most compact

    """

    # Loading sequence
    Rg = compute_rg(traj)

    # Selecting most compact frame
    compact_frame = traj[Rg == Rg.min()]

    return compact_frame


#························································································#
#···························· B I O C H E M I S T R Y ···································#
#························································································#

def compute_Kd(energy, com_diff, T, bins, plot=True) -> float:
    """
    
    Takes arrays of total interaction energy (in kJ/mol) and center of mass distances (in nm),
    computes and returns $K_d$.

    Takes option to plot a scatter of data alongside binning used for computation.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `energy`: `np.darray[float]`
            A one-dimensional array of the total interaction energy between two species;
            In kJ/mol units
        
        `com_diff`: `np.darray[float]`
            A one-dimensional array of the center-of-mass distance between two species;
            In nm units

        `T`: `float`
            The absolute temperature;
            In K units

        `bins`: `int`
            The amount of bins to split the data into

        `plot`: `bool`
            Whether to plot a scatter of the input data and the computed bins

    Returns
    -------

        `Kd`: `float`
            The dissociation constant Kd;
            in nM units

    """
    
    # Binning center of mass difference
    bin_edges = np.linspace(0, np.max(com_diff), bins + 1)
    bin_indices = np.digitize(com_diff, bin_edges)

    # Binning energy
    mean_energy = np.zeros(bins)
    err_energy = np.zeros(bins)
    for i in range(bins):
        bin_mask = (bin_indices == i + 1)
        mean_energy[i] = np.mean(energy[bin_mask])
        err_energy[i] = np.std(energy[bin_mask])/np.sqrt(bin_mask.sum())

    # Defining constants
    R = 8.314462618e-3 # kJ/(mol·K)
    N_A = 6.022e+23 # 1/mol
    pi = np.pi

    # Filtering NaN-value bins off
    r = (((bin_edges[:-1] + bin_edges[1:])/2)*1e-8)[~np.isnan(mean_energy)] # dm
    E = mean_energy[~np.isnan(mean_energy)] # kJ/mol
    E_err = err_energy[~np.isnan(mean_energy)] # kJ/mol

    # Calculating effective potential mean force
    pmf = E + 2*R*T*np.log(r) # kJ/mol

    # Adjusting for offset in energy when no contact is made (i.e. longest CoM distance, if sampled properly)
    offset = pmf[-1] # kJ/mol
    pmf -= offset # kJ/mol

    # Calculating Kd
    Kd = 1/(4*pi*N_A*simpson((np.exp(-(pmf)/((R*T)**2)) * (r**2)), r)) # M
    Kd *= 1e9 # nM

    if plot:
        plt.plot(r*1e8, pmf, linestyle='--', color='black', label='Mean')
        plt.fill_between(r*1e8, pmf-E_err, pmf+E_err, alpha=0.2, color='darkred', label='± Error of mean')
        xlim = plt.xlim(left=0)
        plt.hlines(0, *xlim, linestyles='-', color='grey', linewidth=1, alpha=0.5)
        plt.scatter(com_diff, energy - energy[com_diff.idxmax()], alpha=0.05, s=1, color='darkred')
        plt.xlabel('Center of mass distance\n[nm]', fontsize=12)
        plt.ylabel('Interaction energy\n[kJ/mol]', fontsize=12)
        plt.show()
    
    return Kd