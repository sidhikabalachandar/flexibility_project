"""
The purpose of this code is to collect the poses from crystal structures, glide docking, and combind docking
and aggregate into a single .mae file.
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 aggregate_ligand_poses.py
"""

import os
import pickle

def docking(grid):

    ligands = os.listdir(grid)

    outfile = open('/home/users/sidhikab/MAPK14_ligand_names', 'wb')
    pickle.dump(ligands, outfile)
    outfile.close()

    all_pairs = []
    for ligand in ligands:
        for struc in ligands:
            pair = '{}_to_{}'.format(ligand, struc)
            all_pairs.append(pair)

    outfile = open('/home/users/sidhikab/MAPK14_ligand_docking_names', 'wb')
    pickle.dump(all_pairs, outfile)
    outfile.close()

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)