"""
The purpose of this code is to check if the rmsd calculations have been completed
It can be run on sherlock using
$ ~/miniconda/bin/python3 check_rmsd.py
"""

import os
from docking.docking_class import Docking_Set
import time
import pickle


MAX_NUM_LIGANDS = 25
GROUP_SIZE = 10


'''
Get docking configuration and rmsd set information for all ligands for a given protein
'''
def get_docking_info(protein):
    ligand_folder = "/scratch/PI/rondror/combind/bpp_data/" + protein + "/docking/grids"
    ligands = os.listdir(ligand_folder)
    if len(ligands) > MAX_NUM_LIGANDS:
        ligands = ligands[:MAX_NUM_LIGANDS]
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)

            ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/structures/ligands/{}_lig.mae'.format(protein,
                                                                                                         ligand)

            rmsd_set_info.append({'folder': folderDocking,
                                  'name': name,
                                  'ligand_file': ligand_file})
    return rmsd_set_info


'''
Check if rmsd calculation is finished
'''
def rmsd_complete(protein, grid):
    ligands = os.listdir(grid)
    done = True
    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            rmsd_file = '/home/users/sidhikab/{}/{}/{}_rmsd.csv'.format(protein, name, name)
            if not os.path.isfile(rmsd_file):
                print(rmsd_file)
                done = False
    return done


folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)
proteins = [proteins[0]]
dock_set = Docking_Set()

error = True
for i in range(0, 15):
    print("RMSD Minute", i * 5)
    all_done = True
    for protein in proteins:
        if protein[0] != '.':
            grid = "/scratch/PI/rondror/combind/bpp_data/{}/docking/grids".format(protein)
            if not rmsd_complete(protein, grid):
                print(protein, "not complete")
                all_done = False
    if all_done:
        error = False
        print("RMSD Completed")
        break

    time.sleep(300)
if error:
    print("RMSD did not complete in 70 minutes - possible error")