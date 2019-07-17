'''
This protocol can be used to find the amino acid sequence for each pair of structures
The result is stored in a pickled 2D array
The 2D array will be used for pairwise alignment
'''
#how to run this file:
#ml load chemistry
#ml load schrodinger
#$SCHRODINGER/run python3 sequence_finder.py

from schrodinger.structure import StructureReader
import os
import pickle


folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files/'
ligands = os.listdir(folder)
ligands.remove("4DLI")

strs = []
for ligand in ligands:
	ending = '{}/{}_out.mae'.format(ligand, ligand)
	s = list(StructureReader(folder + ending))[0]
	str = ''
	for mol in list(s.molecule):
		if len(mol.residue) != 1:
			for res in list(mol.residue):
				str += list(res.atom)[0].pdbcode

	strs.append(str)

outfile = open('/home/users/sidhikab/MAPK14_amino_acid_strs', 'wb')
pickle.dump(strs, outfile)
outfile.close()