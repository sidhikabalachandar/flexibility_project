"""
The purpose of this code is to collect the poses from crystal structures, glide docking, and combind docking
and aggregate into a single .mae file.
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ~/miniconda/bin/python3 rmsd_aggregation.py
"""

import os
from docking.docking_class import Docking_Set
import pickle


'''
Get docking configuration and rmsd set information for all ligands for a given protein
'''
def get_docking_info(folder, protein, max_num_ligands):
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = os.listdir(ligand_folder)
    if len(ligands) > max_num_ligands:
        ligands = ligands[:max_num_ligands]
    rmsd_set_info = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            folderDocking = '/home/users/sidhikab/all_protein_docking/{}/{}'.format(protein, name)
            ligand_file = '/scratch/PI/rondror/combind/bpp_data/{}/structures/ligands/{}_lig.mae'.format(protein, ligand)
            rmsd_set_info.append({'folder': folderDocking,
                                  'name': name,
                                  'ligand_file': ligand_file})
    return rmsd_set_info


if __name__ == '__main__':
    folder = "/scratch/PI/rondror/combind/bpp_data"
    proteins = os.listdir('/home/users/sidhikab/all_protein_docking/')
    rmsds = {}
    dock_set = Docking_Set()
    max_num_ligands = 25

    for protein in proteins:
        if protein[0] != '.':
            rmsd_set_info = get_docking_info(folder, protein, max_num_ligands)
            docking_results = dock_set.get_docking_results(rmsd_set_info)
            struc_dict = {}
            for name, rmsd in docking_results.items():
                #name is ligand_to_struc
                ls = name.split('_to_')
                if ls[1] not in struc_dict:
                    struc_dict[ls[1]] = {}
                struc_dict[ls[1]][ls[0]] = rmsd
            rmsds[protein] = struc_dict

    outfile = open('/home/users/sidhikab/all_protein_docking/rmsds', 'wb')
    pickle.dump(rmsds, outfile)
    outfile.close()
    print("RMSD completed")