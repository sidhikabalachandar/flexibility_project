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

def docking(grid):
    """
    :param data_path: (string) location of folder data
    :param save_location: (string) where to save the new .mae files, directory path
    :param docking_version: (string) what docking version to use
    :param scores_version: (string) what score version to sue
    :return:
    """

    ligands = os.listdir(grid)
    docking_config = []
    rmsd_set_info = []

    for ligand in ligands:
        for struc in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/MAPK14/{}'.format(name)
            grid_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids/{}/{}.zip'.format(struc, struc)
            prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(ligand, ligand)
            folderRmsd = 'home/users/sidhikab/MAPK14/{}'.format(ligand)
            ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/{}_lig.mae'

            docking_config.append({'folder': folderDocking,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {}})
    run_config = {'run_folder': '/home/users/sidhikab/MAPK14/run',
                  'group_size': 5,
                  'partition': 'rondror',
                  'dry_run': False}

    dock_set = Docking_Set()
    dock_set.run_docking_set(docking_config, run_config)
    # wait until docking is done to perform rmsd calculations
    for i in range(1, 15):
        if (all(dock_set.check_docking_set_done(docking_config))):
            print("Docking Completed")
            dock_set.run_rmsd_set(docking_config, run_config)
    print("Docking did not complete in 15 minutes - possible error")

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)