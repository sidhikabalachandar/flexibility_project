"""
The purpose of this code is to check if the rmsd calculations have been completed
It can be run on sherlock using
$ ~/miniconda/bin/python3 check_rmsd.py
"""

import os
from dock_rmsd_delete import get_docking_info
from docking.docking_class import Docking_Set

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
            docking_config = get_docking_info(combind_root, protein, max_ligands, output_folder)
            done = dock_set.check_rmsd_set_done(docking_config)
            missing = [item[0]['name'] for item in zip(docking_config, done) if not item[1]]
            print('{}: Missing {}/{}'.format(protein, len(missing),len(docking_config)))
            print(missing)

