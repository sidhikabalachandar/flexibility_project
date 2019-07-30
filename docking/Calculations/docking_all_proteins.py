"""
The purpose of this code is to perform the dockings of a group of proteins
The docking is done between all ligands and all protein structures
If the protein has more than 25 ligands, only the first 25 ligands are considered
It can be run on sherlock using
$ ~/miniconda/bin/python3 docking_all_proteins.py
"""

import os
from docking.docking_class import Docking_Set
import time


MAX_NUM_LIGANDS = 25
GROUP_SIZE = 5


'''
Get docking configuration and rmsd set information for all ligands for a given protein
'''
def get_docking_info(folder, protein):
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = os.listdir(ligand_folder)
    if len(ligands) > MAX_NUM_LIGANDS:
        ligands = ligands[:MAX_NUM_LIGANDS]
    docking_config = []

    # for struc in ligands:
    struc = '1GHW'
    #     for ligand in ligands:
    ligand = '2JHO'
    name = '{}_to_{}'.format(ligand, struc)
    folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
    grid_file = '/scratch/PI/rondror/combind/bpp_data/{}/docking/grids/{}/{}.zip'.format(protein, struc, struc)
    prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(
        protein, ligand, ligand)

    docking_config.append({'folder': folderDocking,
                           'name': name,
                           'grid_file': grid_file,
                           'prepped_ligand_file': prepped_ligand_file,
                           'glide_settings': {}})
    return docking_config


folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)
proteins = [proteins[0]]

for protein in proteins:
    if protein[0] != '.':
        docking_config = get_docking_info(folder, protein)
        run_config = {'run_folder': '/home/users/sidhikab/all_protein_docking/{}/run'.format(protein),
                      'group_size': GROUP_SIZE,
                      'partition': 'owners',
                      'dry_run': False}

        dock_set = Docking_Set()
        dock_set.run_docking_set(docking_config, run_config)

error = True
for i in range(0, 15):
    time.sleep(300)
    print("Docking Minute", i * 5)
    all_done = True
    for protein in proteins:
        if protein[0] != '.':
            docking_config = get_docking_info(folder, protein)
            dock_set = Docking_Set()
            if not all(dock_set.check_docking_set_done(docking_config)):
                all_done = False
    if all_done:
        error = False
        print("Docking Completed")
        break
if error:
    print("Docking did not complete in 70 minutes - possible error")