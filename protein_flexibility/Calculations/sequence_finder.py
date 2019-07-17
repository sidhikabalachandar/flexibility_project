from schrodinger.structure import StructureReader

#how to run this file:
#ml load chemistry
#ml load schrodinger
#$SCHRODINGER/run python3 compute_rmsd_schrodinger.py
#~/miniconda/bin/python3


'''
This protocol can be used to analyze each pair of structures
This will work better than defining the residues using all structures and ligands
because it will be more specific to one pair
'''

folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files/'
struc_1 = list(StructureReader(folder + '3HUC/3HUC_out.mae'))[0]
struc_2 = list(StructureReader(folder + '3GCS/3GCS_out.mae'))[0]

str_struc_1 = ''
for mol in list(struc_1.molecule):
	if len(mol.residue) != 1:
		for res in list(mol.residue):
			str_struc_1 += list(res.atom)[0].pdbcode

str_struc_2 = ''
for mol in list(struc_2.molecule):
	if len(mol.residue) != 1:
		for res in list(mol.residue):
			str_struc_2 += list(res.atom)[0].pdbcode

f = open('seq.txt', "w+")
f.write(str_struc_1)
f.write("\n")
f.write(str_struc_2)