'''
This protocol can be used to find the amino acid sequence for each pair of structures
The result is stored in a pickled 2D array
The 2D array will be used for pairwise alignment

Store outputs in Data/Alignments
Store 1 alignment pickled file per protein

how to run this file:
ml load chemistry
ml load schrodinger
$SCHRODINGER/run python3 sequence_finder.py
'''

from schrodinger.structure import StructureReader
import os
import pickle


'''
Get the list of all proteins
:param combind_root: path to the combind root folder
:return: list of protein name strings
'''
def get_proteins(combind_root):
    proteins = sorted(os.listdir(combind_root))
    proteins = [p for p in proteins if p[0] != '.']
    print(proteins)
    return proteins


'''
Get the list of all ligands
:param protein: name of the protein
:param max_ligands: maximum number of ligands to analyze for each protein
:param combind_root: path to the combind root folder
:return: list of ligand name strings
'''
def get_ligands(protein, max_ligands, combind_root):
    ligand_folder = combind_root + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
    return ligands


'''
Get the amino acid sequence
:param file: .mae file for the structure
:return: the amino acid string for all amino acids in chain A
'''
def get_sequence_from_str(file):
    s = list(StructureReader(file))[0]
    str = ''
    for m in list(s.molecule):
        if len(m.residue) != 1:
            for r in list(m.residue):
                atom = list(r.atom)[0]
                if (atom.chain == "A"):  # to fix for MAPK14: 3GCU
                    str += atom.pdbcode
    return str

if __name__ == '__main__':
    max_ligands = 25
    combind_root = '/scratch/PI/rondror/combind/bpp_data'
    result_folder = '/home/users/lxpowers/projects/combind/flexibility/flexibility_project/protein_flexibility/Data'
    proteins = get_proteins(combind_root)

    for protein in proteins:
        print(protein)
        structure_folder = '{}/{}/structures/aligned_files'.format(combind_root, protein)
        pdb_ids = get_ligands(protein, max_ligands, combind_root)
        seqs = {}
        for pdb_id in pdb_ids:
            pdb_id = pdb_id[:4]
            ending = '/{}/{}_out.mae'.format(pdb_id, pdb_id)
            seqs[pdb_id] = get_sequence_from_str(structure_folder + ending)

        with open('{}/sequences/{}_sequences.pkl'.format(result_folder, protein), 'wb') as f:
            pickle.dump(seqs, f)
