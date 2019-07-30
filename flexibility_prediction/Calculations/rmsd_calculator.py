'''
This protocol can be used to find the rmsd between the residues in the binding pocket of every pair of structures ofa protein
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
import csv


#Loop over each protein
#Loop over each starting ligand and structure
#Loop over each target ligand and structure
#Get the alignment string between these two structures
#Get all the residues in each structure
#Create a map of residues to alignment string index
#Create a map of alignment string index to residues
#Get all residues near the ligand
#Reduce to only the ligands that are paired between each structure


'''
This function gets the pdbcode, chain, resnum, and inscode of every residue in the protein structure
It ignores any residues associated with the ligand
'''
def get_all_res(s):
	r_list = []
	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				if list(r.atom)[0].chain == "A":
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
			counter += 1
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
def get_atoms(s, final_r_list):
	asl_list = []
	complete_a_list = []

	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				if (list(r.atom)[0].pdbcode, list(r.atom)[0].chain, list(r.atom)[0].resnum, list(r.atom)[0].inscode) in final_r_list:
					asl_list.append(r.getAsl())
					complete_a_list.append(r.getAtomList())
	return (asl_list, complete_a_list)


if __name__ == '__main__':
	folder = "/scratch/PI/rondror/combind/bpp_data/"
	#proteins = os.listdir(folder)
	proteins = ['MAPK14']
	protein_dict = {}
	infile = open('MAPK14_ASL_to_resinfo_dict', 'rb')
	ASL_to_feature = pickle.load(infile)
	infile.close()


	with open('res_flexibility.csv', 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow(['protein', 'start ligand', 'target ligand', 'rmsd', 'bfactor', 'normalized bfactor', 'normal variate bfactor', 'res name', 'secondary structure'])

		for protein in proteins:
			if protein[0] != '.':
				ligand_file = folder + protein + "/structures/aligned_files/"
				ligands = os.listdir(ligand_file)
				ligands.remove("4DLI")
				infile = open('MAPK14_pairwise_alignment','rb')
				paired_strs = pickle.load(infile)
				infile.close()
				rmsds = []

				for i, start in enumerate(ligands):
					ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, start, start)
					s1 = list(StructureReader(folder + ending_1))[0]

					arr = []

					for j, target in enumerate(ligands):

						if i != j:

							ending_2 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, target, target)
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

							(asl_list_s1, a_list_s1) = get_atoms(s1, final_r_list_s1)
							(asl_list_s2, a_list_s2) = get_atoms(s2, final_r_list_s2)

							for k in range(len(a_list_s1)):
								if len(a_list_s1[k]) == len(a_list_s2[k]):
									rmsd_val = rmsd.calculate_in_place_rmsd(s1, a_list_s1[k], s2, a_list_s2[k])
									feature = ASL_to_feature[protein][start][asl_list_s1[k]]
									writer.writerow([protein, start, target, rmsd_val, feature[0], feature[1], feature[2], feature[3], feature[4]])

