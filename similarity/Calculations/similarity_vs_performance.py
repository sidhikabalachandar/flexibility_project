"""
This file obtains the tanimoto coefficient for every pair of structures for MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 similarity_vs_performance.py
"""

import os
from logger_class import logger
from schrodinger.structure import StructureReader
from schrodinger.application.canvas.fingerprint import CanvasFingerprintGenerator
from schrodinger.application.canvas.similarity import CanvasFingerprintSimilarity
import pickle


'''
Obtains the tanimoto coefficient for every pair of structures for MAPK14 and pickles them
:param grid: path to the grid files of the MAPK14 structures
:return:
'''
def docking(grid):

    ligands = os.listdir(grid)
    similarities = []
    lgr = logger()
    cfg = CanvasFingerprintGenerator(lgr, 'Radial')
    cfs = CanvasFingerprintSimilarity(lgr)
    cfs.setMetric('Tanimoto')

    for struc_ligand in ligands:
        for ligand in ligands:
            prepped_ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(ligand, ligand)
            prepped_struc_ligand_file = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(struc_ligand, struc_ligand)
            lig = cfg.generate(next(StructureReader(prepped_ligand_file)))
            struc_lig = cfg.generate(next(StructureReader(prepped_struc_ligand_file)))
            similarities.append(cfs.calculateSimilarity(lig, struc_lig))

    outfile = open('/home/users/sidhikab/MAPK14_similarities', 'wb')
    pickle.dump(similarities, outfile)
    outfile.close()

if __name__ == '__main__':
    grid = "/scratch/PI/rondror/combind/bpp_data/MAPK14/docking/grids"
    docking(grid)