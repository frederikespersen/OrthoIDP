#························································································#
#························································································#

# Quick-fix script to reformat any prematurely terminated simulations

#························································································#
#························································································#


import os
import mdtraj as md
import argparse


#························································································#


def save_dcd(traj_path: str, top_path: str, file_path: str, eqsteps: int=0) -> None:
    """
    
    Generates coordinate and trajectory in convenient formats.
    Recenters trajectory frames and removes initial equilibritation steps.

    --------------------------------------------------------------------------------

    Parameters
    ----------

        `traj_path`: `str`
            Path to raw trajectory .dcd file

        `top_path`: `str`
            Path to topology .pdb file

        `file_path`: `str`
            Path to save new processed trajectory .dcd file to

    """

    # Loading trajectory
    traj = md.load(traj_path, top=top_path)

    # Applying periodic boundary conditions and centering around chain 0 in each frame of the trajectory
    traj = traj.image_molecules(anchor_molecules=[set(traj.top.chain(0).atoms)], make_whole=True)

    # Filtering out equilibration from final trajectory
    tocut = eqsteps
    traj[int(tocut):].save_dcd(file_path)


#························································································#


# Setting up argument parser
parser = argparse.ArgumentParser(prog="Pretraj-fixer", description="Quick-fix script to reformat any prematurely terminated simulations")

# Input arguments
parser.add_argument('-p', '--path',
                    type=str,
                    required=False,
                    default='.',
                    help="The path to search for pretraj.dcd files in")
parser.add_argument('-e', '--eqsteps',
                    type=int,
                    required=False,
                    default=1000,
                    help="The amount of equilibration steps to remove")
parser.add_argument('-t', '--test',
                    type=bool,
                    required=False,
                    default=False,
                    help="Whether this is to test which files are located")

# Parsing arguments
args = parser.parse_args()
path = args.path
eqsteps = args.eqsteps
test = args.test


#························································································#


# Finding pretraj directories
pretraj_dirs = []
for root, dirs, files in os.walk(path):
    if 'pretraj.dcd' in files:
        pretraj_dirs.append(root)

# Applying reformatting
for dir_path in pretraj_dirs:
    print("Fixing:", dir_path)
    if not test:
        try:
            save_dcd(traj_path=f'{dir_path}/pretraj.dcd', top_path=f'{dir_path}/top.pdb', file_path=f'{dir_path}/traj.dcd', eqsteps=eqsteps)
        except:
            print("FAILED TO FIX:", dir_path)