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
from docking.docking_class import Docking

def docking(grid):

    ligands = os.listdir(grid)
    dock = Docking("/home/users/sidhikab/MAPK14", "")
    for ligand in ligands:
        min_gscores = []
        for struc in ligands:
            pair = '{}_to_{}'.format(ligand, struc)
            gscores, emodels, rmsds = dock.get_gscores_emodels(pair)
            min_gscores.append((min(gscores), gscores.index(min(gscores)), struc))

        best = min(min_gscores)
        print("For ligand ", ligand, " best pose is pose ", best[1], " for protein corresponding to ligand ", best[2])
        print("The corresponding gscore is ", best[0])

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)