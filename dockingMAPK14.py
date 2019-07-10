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
import pickle
import time

def docking(grid):

    ligands = os.listdir(grid)
    docking_config = []
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/MAPK14/{}'.format(name)
            grid_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids/{}/{}.zip'.format(struc, struc)
            prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(ligand, ligand)

            ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/{}_lig.mae'.format(ligand)

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
    dock_set.run_docking_set(docking_config, run_config)
    error = True

    for i in range(1, 15):
        time.sleep(60)
        print("Docking Minute", i)

        if all(dock_set.check_docking_set_done(docking_config)):
            error = False
            print("Docking Completed")
            dock_set.run_rmsd_set(rmsd_set_info, run_config)
            pickleError = True

            for i in range(1, 15):
                time.sleep(60)
                print("RMSD Minute", i)

                if rmsd_complete(grid):
                    pickleError = False
                    outfile = open('/home/users/sidhikab/MAPK14_rmsds', 'wb')
                    pickle.dump(dock_set.get_docking_results(rmsd_set_info), outfile)
                    outfile.close()
                    print("RMSD completed")
                    break

            if pickleError:
                print("RMSD did not complete in 15 minutes - possible error")
            break
    if error:
         print("Docking did not complete in 15 minutes - possible error")

def rmsd_complete(grid):
    ligands = os.listdir(grid)
    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            rmsd_file = '/home/users/sidhikab/MAPK14/{}/{}_rmsd.csv'.format(name, name)
            if not os.path.isfile(rmsd_file):
                return False
    return True

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)