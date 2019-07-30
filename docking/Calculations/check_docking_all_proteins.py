"""
The purpose of this code is to check which proteins have finished docking
It can be run on sherlock using
$ ~/miniconda/bin/python3 check_docking_all_proteins.py
"""

import os
from docking.docking_class import Docking_Set
import time

MAX_NUM_LIGANDS = 25

'''
Get docking configuration and rmsd set information for all ligands for a given protein
'''
def get_docking_info(folder, protein):
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = os.listdir(ligand_folder)
    if len(ligands) > MAX_NUM_LIGANDS:
        ligands = ligands[:MAX_NUM_LIGANDS]

    print(len(ligands))
    docking_config = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
            grid_file = '{}/{}/docking/grids/{}/{}.zip'.format(folder, protein, struc, struc)
            prepped_ligand_file = '{}/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(folder, protein, ligand, ligand)

            docking_config.append({'folder': folderDocking,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {}})
    return docking_config
folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)
proteins = proteins[:4]
counter = 0

error = True
for i in range(0, 15):
    print("RMSD Minute", i * 5)
    all_done = True
    for protein in proteins:
        if protein[0] != '.':
            docking_config = get_docking_info(folder, protein)
            dock_set = Docking_Set()
            if all(dock_set.check_docking_set_done(docking_config)):
                counter += 1
                print(protein, " completed")
            else:
                print(protein, " not completed")
                all_done = False

    if all_done:
        break

    time.sleep(300)
print(counter, " out of ", len(proteins), " completed: ", counter * 100 / len(proteins), "% complete")