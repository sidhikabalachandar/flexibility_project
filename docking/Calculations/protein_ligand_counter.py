"""
The purpose of this code is to count the number of ligands for each protein in the bpp_data folder
It can be run using
$ ~/miniconda/bin/python3 protein_ligand_counter.py
"""

import os

folder = "/scratch/PI/rondror/combind/bpp_data"
proteins = os.listdir(folder)

for protein in proteins:
    if protein[0] != '.':
        ligand_folder = folder + "/" + protein + "/docking/grids"
        ligands = os.listdir(ligand_folder)
        print(protein, " has ",len(ligands), " ligands")