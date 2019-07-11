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