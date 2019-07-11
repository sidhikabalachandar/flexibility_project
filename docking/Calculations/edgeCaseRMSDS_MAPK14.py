"""
The purpose of this code is to collect the poses from crystal structures, glide docking, and combind docking
and aggregate into a single .mae file.
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 aggregate_ligand_poses.py
"""

import os
from docking.docking_class import Docking_Set
import time

def rmsd():

    docking_config = []
    rmsd_set_info = []
    ligand = "2BAK"
    name = '3GCP_to_{}'.format(ligand)
    folderDocking = '/home/users/sidhikab/MAPK14/{}'.format(name)
    grid_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids/{}/{}.zip'.format(ligand, ligand)
    prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/3GCP_lig/3GCP_lig.mae'

    ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/3GCP_lig.mae'

    docking_config.append({'folder': folderDocking,
                           'name': name,
                           'grid_file': grid_file,
                           'prepped_ligand_file': prepped_ligand_file,
                           'glide_settings': {}})
    rmsd_set_info.append({'folder' : folderDocking,
                          'name' : name,
                          'ligand_file' : ligand_file})
    run_config = {'run_folder': '/home/users/sidhikab/MAPK14/run',
                  'group_size': 5,
                  'partition': 'owners',
                  'dry_run': False}


    dock_set = Docking_Set()
    dock_set.run_rmsd_set(rmsd_set_info, run_config)
    rmsdError = True
    for i in range (1, 15):
        time.sleep(60)
        print("RMSD minute ", i)
        if checkRMSD(ligand):
            print('RMSD complete')
            rmsdError = False
            break
    if rmsdError:
        print("RMSD did not complete in 15 minutes - possible error")

def checkRMSD(ligand):
    filename = '/home/users/sidhikab/MAPK14/3GCP_to_{}/3GCP_to_{}_rmsd.csv'.format(ligand, ligand)
    return os.path.isfile(filename)
if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    rmsd()