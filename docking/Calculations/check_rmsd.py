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
def get_docking_info(ligands, protein):
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


folder = "/scratch/PI/rondror/combind/bpp_data/"
proteins = os.listdir(folder)
dock_set = Docking_Set()

all_done = True
for protein in proteins:
    if protein[0] != '.':
        print(protein)
        ligand_folder = "/scratch/PI/rondror/combind/bpp_data/" + protein + "/docking/grids"
        ligands = os.listdir(ligand_folder)
        if len(ligands) > MAX_NUM_LIGANDS:
            ligands = ligands[:MAX_NUM_LIGANDS]

        rmsd_set_info = get_docking_info(ligands, protein)
        bool_done_list = dock_set.check_rmsd_set_done(rmsd_set_info)
        name_incomplete = []
        counter = 0
        not_done_counter = 0
        for struc in ligands:
            for ligand in ligands:
                if not bool_done_list[counter]:
                    name = '{}_to_{}'.format(ligand, struc)
                    name_incomplete.append(name)
                    not_done_counter += 1
                counter += 1

        print(not_done_counter, "/", counter, "incomplete")
        print(name_incomplete)
