'''
This protocol can be used to find the rmsd between every pair of protein structures
Only the residues within 4 angstroms of either structures' ligands are considered
'''
# how to run this file:
# ml load chemistry
# ml load schrodinger
# $SCHRODINGER/run python3 rmsd_calculator.py

from schrodinger.structure import StructureReader
import schrodinger.structutils.measure as measure
import schrodinger.structutils.rmsd as rmsd
import schrodinger.structutils.analyze as analyze
import os
import pickle


'''
This function gets the pdbcode, chain, resnum, and inscode of every residue in the protein structure
It ignores any residues associated with the ligand
'''
def get_all_res(s):
	r_list = []
	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				r_list.append((list(r.atom)[0].pdbcode, list(r.atom)[0].chain, list(r.atom)[0].resnum, list(r.atom)[0].inscode))
	return r_list


'''
Maps unique residue identifiers to list index in alignment string

alignment_string: (string) output from alignment program, contains one letter codes and dashes
	example: 'TE--S--T-'
r_list: list of unique identifiers of each residue in order of sequence
	number of residues in r_list must be equal to number of residues in alignment_string
'''
def map_residues_to_align_index(alignment_string, r_list):
	r_to_i_map = {}
	counter = 0
	for i in range(len(alignment_string)):
		if counter >= len(r_list):
			break
		if alignment_string[i] == r_list[counter][0]:
			r_to_i_map[r_list[counter]] = i
			counter += 1;
	return r_to_i_map


'''
This function inverses an input map
The keys become values and vice versa
'''
def inv_map(m):
	return {v: k for k, v in m.items()}


'''
This function gets the unique identifier for all residues within 4 angstroms of the ligand
'''
def get_res_near_ligand(s):
	cutoff = 4
	l = analyze.find_ligands(s)[0]
	#get atom indexes for the ligand
	close_a_indices = measure.get_atoms_close_to_subset(s, l.atom_indexes, cutoff)
	close_r_set = set([(s.atom[i].pdbcode, s.atom[i].chain, s.atom[i].resnum, s.atom[i].inscode) for i in close_a_indices])
	return close_r_set


'''
This function gets the atom list corresponding to a given list of unique residue identifiers from a given protein structure
'''
def get_atom_list(s, final_r_list):
	complete_a_list = []

	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				if (list(r.atom)[0].pdbcode, list(r.atom)[0].chain, list(r.atom)[0].resnum, list(r.atom)[0].inscode) in final_r_list:
					complete_a_list.append(r.getAtomList())

	return complete_a_list


folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files/'
ligands = os.listdir(folder)
ligands.remove("4DLI")

infile = open('MAPK14_pairwise_alignment','rb')
paired_strs = pickle.load(infile)
infile.close()

rmsds = []

for i in range(len(ligands)):
	ending_1 = '{}/{}_out.mae'.format(ligands[i], ligands[i])
	s1 = list(StructureReader(folder + ending_1))[0]

	arr = []

	for j in range(len(ligands)):
		#to monitor progress
		print(i, ligands[i], j, ligands[j])

		ending_2 = '{}/{}_out.mae'.format(ligands[j], ligands[j])
		s2 = list(StructureReader(folder + ending_2))[0]

		(paired_str_s1, paired_str_s2) = paired_strs[i][j]

		r_list_s1 = get_all_res(s1)
		r_list_s2 = get_all_res(s2)

		r_to_i_map_s1 = map_residues_to_align_index(paired_str_s1, r_list_s1)
		r_to_i_map_s2 = map_residues_to_align_index(paired_str_s2, r_list_s2)
		i_to_r_map_s1 = inv_map(r_to_i_map_s1)
		i_to_r_map_s2 = inv_map(r_to_i_map_s2)

		valid_r_s1 = get_res_near_ligand(s1)
		valid_r_s2 = get_res_near_ligand(s2)

		final_r_list_s1 = []
		final_r_list_s2 = []

		for r in valid_r_s1:
			s1index = r_to_i_map_s1[r]
			if paired_str_s1[s1index] == paired_str_s2[s1index]:
				if r not in final_r_list_s1:
					final_r_list_s1.append(r)
				if i_to_r_map_s2[s1index] not in final_r_list_s2:
					final_r_list_s2.append(i_to_r_map_s2[s1index])

		for r in valid_r_s2:
			s2index = r_to_i_map_s2[r]
			if paired_str_s2[s2index] == paired_str_s1[s2index]:
				if r not in final_r_list_s2:
					final_r_list_s2.append(r)
				if i_to_r_map_s1[s2index] not in final_r_list_s1:
					final_r_list_s1.append(i_to_r_map_s1[s2index])

		a_list_s1 = get_atom_list(s1, final_r_list_s1)
		a_list_s2 = get_atom_list(s2, final_r_list_s2)

		final_a_list_s1 = []
		final_a_list_s2 = []

		for k in range(len(a_list_s1)):
			if len(a_list_s1[k]) == len(a_list_s2[k]):
				final_a_list_s1 += a_list_s1[k]
				final_a_list_s2 += a_list_s2[k]

		arr.append(rmsd.superimpose(s1, final_a_list_s1, s2, final_a_list_s2))

	rmsds.append(arr)

print(rmsds)
outfile = open('/home/users/sidhikab/MAPK14_pairwise_struc_rmsds', 'wb')
pickle.dump(rmsds, outfile)
outfile.close()
