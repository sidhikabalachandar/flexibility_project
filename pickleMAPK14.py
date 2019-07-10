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
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/MAPK14/{}'.format(name)
            ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/{}_lig.mae'.format(ligand)

            rmsd_set_info.append({'folder' : folderDocking,
                                  'name' : name,
                                  'ligand_file' : ligand_file})


    dock_set = Docking_Set()

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