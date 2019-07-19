"""
The purpose of this code is to collect the RMSDs for the dockings of all proteins
The docking is done between all ligands and all protein structures
If the protein has more than 50 ligands, only the first 50 ligands are considered
It can be run on sherlock using
$ ~/miniconda/bin/python3 docking_all_proteins.py
"""

import os
from docking.docking_class import Docking_Set
import time


'''
Get docking configuration and rmsd set information for all ligands for a given protein
'''
def get_docking_info(protein):
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = os.listdir(ligand_folder)
    if len(ligands) > 50:
        ligands = ligands[:50]
    docking_config = []
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
            grid_file = '/scratch/PI/rondror/combind/bpp_data/{}/docking/grids/{}/{}.zip'.format(protein, struc, struc)
            prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(
                protein, ligand, ligand)

            ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/structures/ligands/{}_lig.mae'.format(protein,
                                                                                                         ligand)

            docking_config.append({'folder': folderDocking,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {}})
            rmsd_set_info.append({'folder': folderDocking,
                                  'name': name,
                                  'ligand_file': ligand_file})
    return (docking_config, rmsd_set_info)


folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)

for protein in proteins:
    if protein[0] != '.':
        (docking_config, rmsd_set_info) = get_docking_info(protein)
        run_config = {'run_folder': '/home/users/sidhikab/all_protein_docking/{}/run'.format(protein),
                      'group_size': 20,
                      'partition': 'owners',
                      'dry_run': False}

        dock_set = Docking_Set()
        if not all(dock_set.check_docking_set_done(docking_config)):
            dock_set.run_docking_set(docking_config, run_config)

error = True
for i in range(0, 15):
    time.sleep(300)
    print("Docking Minute", i * 5)
    all_done = True
    for protein in  proteins:
        if protein[0] != '.':
            (docking_config, rmsd_set_info) = get_docking_info(protein)
            if not all(dock_set.check_docking_set_done(docking_config)):
                all_done = False
    if all_done:
        error = False
        print("Docking Completed")
        break
if error:
    print("Docking did not complete in 70 minutes - possible error")