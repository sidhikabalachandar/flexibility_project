"""
The purpose of this code is to list the names of every pair of ligands and structures for MAPK14
Format is ligand_to_struc
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ~/miniconda/bin/python3 ligand_pickler.py
"""

import os


'''
lists the names of every pair of ligands and structures for MAPK14
:param grid: path to the grid files of MAPK14
'''
def names(grid):
    ligands = os.listdir(grid)
    counter = 0

    for struc in ligands:
        for ligand in ligands:
            name = '{}_to_{}'.format(ligand, struc)
            print(counter, name)
            counter += 1

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    names(grid)