"""
The purpose of this code is to collect the names of every pair of ligands and structures for MAPK14 and pickle the output
Format is ligand_to_struc
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ~/miniconda/bin/python3 ligand_pickler.py
"""

import os
import pickle


'''
collects the names of every pair of ligands and structures for MAPK14 and pickles the output
:param grid: path to the grid files of MAPK14
'''
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