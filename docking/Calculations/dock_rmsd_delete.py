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

def get_docking_info(folder, protein, max_ligands, output_folder_root):
    '''
    Get docking configuration and rmsd set information for all ligands for a given protein
    '''
    ligand_folder = folder + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands] #sorted
    docking_config = []

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            output_folder = output_folder_root + '/{}/{}'.format(protein, name)
            grid_file = folder + '/{}/docking/grids/{}/{}.zip'.format(protein, struc, struc)
            prepped_ligand_file = folder + '/{}/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(protein, ligand, ligand)
            ligand_file = folder + '/{}/structures/ligands/{}_lig.mae'.format(protein, ligand)
            docking_config.append({'folder': output_folder,
                                   'name': name,
                                   'grid_file': grid_file,
                                   'prepped_ligand_file': prepped_ligand_file,
                                   'glide_settings': {},
                                   'ligand_file': ligand_file})
    return docking_config

if __name__ == '__main__':
    max_ligands = 25
    combind_root = '/scratch/PI/rondror/combind/bpp_data'
    output_folder = '/home/users/lxpowers/projects/combind/all_docking'

    proteins = sorted(os.listdir(combind_root))
    proteins = [proteins[-1]]
    print(proteins)
    dock_set = Docking_Set()

    for protein in proteins:
        if protein[0] != '.':
            if sys.argv[1] == 'run_dock':
                docking_config = get_docking_info(combind_root, protein, max_ligands, output_folder)
                run_config = {'run_folder': output_folder+'/{}/run'.format(protein),
                              'group_size': 5,
                              'partition': 'owners',
                              'dry_run': True}

                dock_set.run_docking_rmsd_delete(docking_config, run_config)
            else:
                #check progress
                done = dock_set.check_rmsd_set_done(docking_config)
                missing = [item[0]['name'] for item in zip(docking_config, done) if not item[1]]
                print('{}: Missing {}/{}'.format(protein, len(missing), len(docking_config)))
                print(missing)
