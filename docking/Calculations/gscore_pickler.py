"""
The purpose of this code is to collect the gscores for every pair of ligands and structures for MAPK14 and pickle the output
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ~/miniconda/bin/python3 gscore_pickler.py
"""

import os
from docking.docking_class import Docking
import pickle



'''
collects the gscores for every pair of ligands and structures for MAPK14 and pickles the output
:param grid: path to the grid files of MAPK14
'''
def docking(grid):
    ligands = os.listdir(grid)
    dock = Docking("/home/users/sidhikab/MAPK14", "")
    all_gscores = []
    for ligand in ligands:
        min_gscores = []
        for struc in ligands:
            pair = '{}_to_{}'.format(ligand, struc)
            gscores, emodels, rmsds = dock.get_gscores_emodels(pair)
            min_gscores.append(gscores[0])

        all_gscores.append(min_gscores)


    outfile = open('/home/users/sidhikab/MAPK14_gscores', 'wb')
    pickle.dump(all_gscores, outfile)
    outfile.close()

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)