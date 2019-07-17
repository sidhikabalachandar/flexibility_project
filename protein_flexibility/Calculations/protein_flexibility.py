"""
The purpose of this code is to collect the poses from crystal structures, glide docking, and combind docking
and aggregate into a single .mae file.
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 aggregate_ligand_poses.py
"""

import schrodinger.structutils.rmsd as r
from schrodinger.structure import StructureReader, StructureWriter
import os

def get_all_res(grid):
    ligands = os.listdir(grid)
    rmsds = []
    for struc in ligands:
        arr = []
        for ligand in ligands:
            #get struc's protein structure and atom list
            struc_struc_file = "{}/{}/{}_out.mae".format(grid, struc, struc)
            struc_structures = list(StructureReader(struc_struc_file))
            struc_struc = struc_structures[0]
            struc_molecules = list(struc_struc.molecule)
            struc_prot_atoms = []

            for mol in struc_molecules:
                if len(mol.residue) == 1:
                    struc_lig_atoms = mol.getAtomList()
                    struc_lig = struc_struc.extract(struc_lig_atoms)
                else:
                    struc_prot_atoms += mol.getAtomList()

            struc_prot = struc_struc.extract(struc_prot_atoms)

            #get lig's protein structure and atom list
            lig_struc_file = "{}/{}/{}_out.mae".format(grid, lig, lig)
            lig_structures = list(StructureReader(lig_struc_file))
            lig_struc = lig_structures[0]
            lig_molecules = list(lig_struc.molecule)
            lig_prot_atoms = []

            for mol in lig_molecules:
                if len(mol.residue) == 1:
                    lig_lig_atoms = mol.getAtomList()
                    lig_lig = lig_struc.extract(lig_lig_atoms)
                else:
                    lig_prot_atoms += mol.getAtomList()

            lig_prot = lig_struc.extract(lig_prot_atoms)

            #calculate RMSD
            #PROBLEM: atom lists must be the same length
            #this is not the case for many pairs
            #will this be solved if we calculate the rmsd of all residues within 5 A of the ligand?
            arr.append(r.calculate_in_place_rmsd(struc_prot, struc_prot_atoms, lig_prot, lig_prot_atoms))

        rmsds.append(arr)


if __name__ == '__main__':
    grid = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files'
    get_all_res(grid)