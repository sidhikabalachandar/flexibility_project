"""
The purpose of this code is to first  perform the dockings of a group of proteins
Then the rmsds are calculated
Finally the pose viewer files are deleted
It can be run on sherlock using
$ ~/miniconda/bin/python3 dock_rmsd_delete.py
"""

import os
import sys
from docking.docking_class import Docking_Set
import pickle


'''
Get docking configuration and rmsd set information for all ligands for a given protein
:param folder: path to root folder for protein data
:param protein: name of protein
:param max_ligands: maximum number of ligands to collect data for
:param output_folder_root: path to save location root
:return: a list of dictionaries containing docking configuration information for every pair of ligands to protein structures for a given protein
'''
def get_docking_info(folder, protein, max_ligands, output_folder_root):
    '''
    Get docking configuration and rmsd set information for all ligands for a given protein
    '''
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands] #sorted
    # ligands = ['3IW6_to_2BAL', '3O8P_to_2BAL', '3ROC_to_2BAL', '4DLI_to_2BAL', '4F9W_to_2BAL']
    docking_config = []

    # for name in ligands:
    for struc in ligands:
        for ligand in ligands:
            if struc != ligand:
        # ligand, struc = name.split('_to_')
                name = '{}_to_{}'.format(ligand, struc)
                output_folder = output_folder_root + '/{}/{}'.format(protein, name)
                grid_file = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/grids/{}/{}.zip'.format(name, name)
                prepped_ligand_file = folder + '/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(protein, ligand, ligand)
                ligand_file = folder + '/{}/structures/ligands/{}_lig.mae'.format(protein, ligand)
                docking_config.append({'folder': output_folder,
                                       'name': name,
                                       'grid_file': grid_file,
                                       'prepped_ligand_file': prepped_ligand_file,
                                       'glide_settings': {},
                                       'ligand_file': ligand_file})
    return docking_config


'''
Get the list of all proteins
:param combind_root: path to the combind root folder
:return: list of protein name strings
'''
def get_proteins(combind_root):
    '''
    Get the list of all proteins
    :param combind_root: path to the combind root folder
    :return: list of protein name strings
    '''
    proteins = sorted(os.listdir(combind_root))
    proteins = [p for p in proteins if p[0] != '.']
    print(proteins)
    return proteins


if __name__ == '__main__':
    max_ligands = 25
    combind_root = '/scratch/PI/rondror/combind/bpp_data'
    output_folder = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/mut_rmsds'
    result_folder = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/mut_rmsds'
    proteins = ['MAPK14']
    dock_set = Docking_Set()

    task = sys.argv[1]
    if task == 'run_dock':
        for protein in proteins:
            docking_config = get_docking_info(combind_root, protein, max_ligands, output_folder)
            run_config = {'run_folder': output_folder+'/{}/run'.format(protein),
                          'group_size': 15,
                          'partition': 'owners',
                          'dry_run': False}
            print(protein)
            dock_set.run_docking_rmsd_delete(docking_config, run_config, incomplete_only=True)

    if task == 'check':
        for protein in proteins:
            docking_config = get_docking_info(combind_root, protein, max_ligands, output_folder)
            done = dock_set.check_rmsd_set_done(docking_config)
            missing = [item[0]['name'] for item in zip(docking_config, done) if not item[1]]
            print('{}: Missing {}/{}'.format(protein, len(missing), len(docking_config)))
            print(missing)

    if task == 'results':
        rmsds = {}
        for protein in proteins:
            print(protein)
            docking_config = get_docking_info(combind_root, protein, max_ligands, output_folder)
            docking_results = dock_set.get_docking_results(docking_config)
            struc_dict = {}
            for name, rmsd in docking_results.items():
                # name is ligand_to_struc
                ls = name.split('_to_')
                if ls[1] not in struc_dict:
                    struc_dict[ls[1]] = {}
                struc_dict[ls[1]][ls[0]] = rmsd
            rmsds[protein] = struc_dict
        with open(result_folder + '/rmsds.pkl', 'wb') as outfile:
            pickle.dump(rmsds, outfile)
