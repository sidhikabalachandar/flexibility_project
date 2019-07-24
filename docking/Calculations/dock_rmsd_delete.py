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
combind_root = "/scratch/PI/rondror/combind/bpp_data"


def get_docking_info(folder, protein):
    '''
    Get docking configuration and rmsd set information for all ligands for a given protein
    '''
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:MAX_NUM_LIGANDS] #sorted
    docking_config = []
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            output_folder = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
            root = '/scratch/PI/rondror/combind/bpp_data/'
            grid_file = folder + '{}/docking/grids/{}/{}.zip'.format(protein, struc, struc)
            prepped_ligand_file = folder + '{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(protein, ligand, ligand)
            ligand_file = folder + '{}/structures/ligands/{}_lig.mae'.format(protein, ligand)

            docking_config.append({'folder': output_folder,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {},
                                   'ligand_file': ligand_file})
    return docking_config


proteins = os.listdir(combind_root)
proteins = proteins[-1]

dock_set = Docking_Set()

for protein in proteins:
    if protein[0] != '.':
        (docking_config, rmsd_set_info) = get_docking_info(combind_root, protein)
        run_config = {'run_folder': '/home/users/sidhikab/all_protein_docking/{}/run'.format(protein),
                      'group_size': GROUP_SIZE,
                      'partition': 'owners',
                      'dry_run': False}

        dock_set.run_docking_rmsd_delete(docking_config, run_config)