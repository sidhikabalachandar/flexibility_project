from Bio import pairwise2

#how to run this file:
# ~/miniconda/bin/python3


'''
This protocol can be used to analyze each pair of structures
This will work better than defining the residues using all structures and ligands
because it will be more specific to one pair
'''

with open('seq.txt', "r") as f:
	lines = f.readlines()
	str_struc_1 = lines[0]
	str_struc_2 = lines[1]

alignments = pairwise2.align.globalxx(str_struc_1, str_struc_2)
paired_str_struc_1 = alignments[0][0]
paired_str_struc_2 = alignments[0][1]

f = open('paired_seq.txt', "w+")
f.write(paired_str_struc_1)
#f.write("\n")
f.write(paired_str_struc_2)