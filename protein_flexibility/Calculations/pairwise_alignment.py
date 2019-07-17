'''
This protocol can be used to find the pairwise alignemnt between the amino acid strings of each pair of proteins
'''
#how to run this file:
# ~/miniconda/bin/python3 pairwise_alignment.py

from Bio import pairwise2
import os
import pickle


infile = open('MAPK14_amino_acid_strs','rb')
strs = pickle.load(infile)
infile.close()

folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files/'
ligands = os.listdir(folder)
ligands.remove("4DLI")

paired_strs = []
for i in range(len(ligands)):
	str_s1 = strs[i]
	arr = []

	for j in range(len(ligands)):
		str_s2 = strs[j]
		alignments = pairwise2.align.globalxx(str_s1, str_s2)
		arr.append((alignments[0][0], alignments[0][1]))

	paired_strs.append(arr)

outfile = open('/home/users/sidhikab/MAPK14_pairwise_alignment', 'wb')
pickle.dump(paired_strs, outfile)
outfile.close()