"""
The purpose of this code is to delete all current pose viewer files to increase the amount of available space in the home directory
It can be run on sherlock using
$ ~/miniconda/bin/python3 delete_pose_viewer.py
"""

import os


MAX_NUM_LIGANDS = 25
GROUP_SIZE = 10

folder = "/scratch/PI/rondror/combind/bpp_data/"
proteins = os.listdir(folder)
proteins = proteins[:4]

for protein in proteins:
    if protein[0] != '.':
        ligand_folder = folder + protein + "/docking/grids"
        ligands = os.listdir(ligand_folder)
        if len(ligands) > MAX_NUM_LIGANDS:
            ligands = ligands[:MAX_NUM_LIGANDS]

        for struc in ligands:
            for ligand in ligands:
                name = '{}_to_{}'.format(ligand, struc)
                file_name = '/home/users/sidhikab/all_protein_docking/{}/{}/{}_pv.maegz'.format(protein, name, name)
                if os.path.exists(file_name):
                    os.remove(file_name)

print("All pose viewer files removed")