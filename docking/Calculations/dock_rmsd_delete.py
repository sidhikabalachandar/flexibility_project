"""
The purpose of this code is to first  perform the dockings of a group of proteins
Then the rmsds are calculated
Finally the pose viewer files are deleted
It can be run on sherlock using
$ ~/miniconda/bin/python3 dock_rmsd_delete.py
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
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
            grid_file = '/scratch/PI/rondror/combind/bpp_data/{}/docking/grids/{}/{}.zip'.format(protein, struc, struc)
            prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(
                protein, ligand, ligand)
            ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/structures/ligands/{}_lig.mae'.format(protein, ligand)

            docking_config.append({'folder': folderDocking,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {},
                                   'ligand_file': ligand_file})
            rmsd_set_info.append({'folder': folderDocking,
                                  'name': name,
                                  'ligand_file': ligand_file})
    return (docking_config, rmsd_set_info)


'''
Check if rmsd calculation is finished
'''
def rmsd_complete(grid):
    ligands = os.listdir(grid)
    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            rmsd_file = '/home/users/sidhikab/MAPK14/{}/{}_rmsd.csv'.format(name, name)
            if not os.path.isfile(rmsd_file):
                return False
    return True


folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)
proteins = [proteins[4]]
dock_set = Docking_Set()

for protein in proteins:
    if protein[0] != '.':
        (docking_config, rmsd_set_info) = get_docking_info(folder, protein)
        run_config = {'run_folder': '/home/users/sidhikab/all_protein_docking/{}/run'.format(protein),
                      'group_size': GROUP_SIZE,
                      'partition': 'owners',
                      'dry_run': False}

        dock_set.run_docking_set(docking_config, run_config)

error = True
for i in range(0, 15):
    time.sleep(300)
    print("RMSD Minute", i * 5)
    all_done = True
    for protein in proteins:
        if protein[0] != '.':
            (docking_config, rmsd_set_info) = get_docking_info(folder, protein)
            grid = "/scratch/PI/rondror/combind/bpp_data/{}/docking/grids".format(protein)
            if not rmsd_complete(grid):
                all_done = False
    if all_done:
        error = False
        print("RMSD Completed")
        break
if error:
    print("RMSD did not complete in 70 minutes - possible error")