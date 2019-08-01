'''
This protocol can be used to find the rmsd between the residues in the binding pocket of every pair of structures ofa protein
Only the residues within 4 angstroms of either structures' ligands are considered
'''
# how to run this file:
# ml load chemistry
# ml load schrodinger
# $SCHRODINGER/run python3 rmsd_calculator.py

import schrodinger.structure as structure
import schrodinger.structutils.measure as measure
import schrodinger.structutils.rmsd as rmsd
import schrodinger.structutils.analyze as analyze
import os
import pickle
import csv
import sys
import time


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
			#print(r_list[counter], i)
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
def get_res_near_ligand(s, r_to_i_map):
	cutoff = 4
	ligand_alist = []
	for m in list(s.molecule):
		for r in list(m.residue):
			if list(r.atom)[0].chain == "L":
				ligand_alist += r.getAtomList()
	if ligand_alist == []:
		return 0
	#get atom indexes for the ligand
	close_a_indices = measure.get_atoms_close_to_subset(s, ligand_alist, cutoff)
	close_r_set = set({})
	for i in close_a_indices:
		if (s.atom[i].pdbcode, s.atom[i].chain, s.atom[i].resnum, s.atom[i].inscode) in r_to_i_map:
			close_r_set.add((s.atom[i].pdbcode, s.atom[i].chain, s.atom[i].resnum, s.atom[i].inscode))
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


def get_ligands(protein, max_ligands, combind_root):
    ligand_folder = combind_root + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
    return ligands


def get_proteins(combind_root):
	'''
	Get the list of all proteins
	:param combind_root: path to the combind root folder
	:return: list of protein name strings
	'''
	proteins = sorted(os.listdir(combind_root))
	proteins = [p for p in proteins if p[0] != '.']
	print(proteins)
	return proteins


def compute_protein_rmsds(protein, rmsd_file, ASL_to_feature, combind_root):
	max_ligands = 25


	with open(rmsd_file, 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow(['protein', 'start ligand', 'target ligand', 'rmsd', 'bfactor', 'normalized bfactor', 'normal variate bfactor', 'res name', 'secondary structure'])

		ligands = get_ligands(protein, max_ligands, combind_root)
		infile = open('../../protein_flexibility/Data/alignments/{}_alignment.pkl'.format(protein),'rb')
		paired_strs = pickle.load(infile)
		infile.close()

		for start in ligands:

			if start not in ASL_to_feature[protein]:
				print(protein, start, "not in dict")
				#continue
			ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, start, start)
			s1 = list(structure.StructureReader(combind_root + ending_1))[0]

			for target in ligands:

				if start != target:
					ending_2 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, target, target)
					s2 = list(structure.StructureReader(combind_root + ending_2))[0]

					if start < target:
						(paired_str_s1, paired_str_s2) = paired_strs[start][target]
					else:
						(paired_str_s2, paired_str_s1) = paired_strs[target][start]

					r_list_s1 = get_all_res(s1)
					r_list_s2 = get_all_res(s2)

					r_to_i_map_s1 = map_residues_to_align_index(paired_str_s1, r_list_s1)
					r_to_i_map_s2 = map_residues_to_align_index(paired_str_s2, r_list_s2)
					i_to_r_map_s1 = inv_map(r_to_i_map_s1)
					i_to_r_map_s2 = inv_map(r_to_i_map_s2)

					valid_r_s1 = get_res_near_ligand(s1, r_to_i_map_s1)
					valid_r_s2 = get_res_near_ligand(s2, r_to_i_map_s2)

					if valid_r_s1 == set({}):
						print(protein, start, "no residues close to the ligand")
						#continue

					if valid_r_s1 == 0:
						print(protein, target, "pose viewer file has no ligand")
						#continue

					if valid_r_s2 == set({}):
						print(protein, start, "no residues close to the ligand")
						#continue

					if valid_r_s2 == 0:
						print(protein, target, "pose viewer file has no ligand")
						#continue

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


if __name__ == '__main__':
	task = sys.argv[1]
	combind_root = '/scratch/PI/rondror/combind/bpp_data/'
	result_folder = '/home/users/sidhikab/flexibility_project/flexibility_prediction/Data'
	save_folder = result_folder+'/rmsds/'
	partition = 'rondror'
	infile = open('../Data/ASL_to_resinfo_dict', 'rb')
	ASL_to_feature = pickle.load(infile)
	infile.close()

	if task == 'all':
		proteins = get_proteins(combind_root)
		#submit jobs for each protein
		cmd = 'sbatch -p {} -t 1:00:00 -o {}_rmsd.out --wrap="$SCHRODINGER/run python3 rmsd_calculator.py  protein {}"'
		for prot_name in proteins:
			print(prot_name)
			if not os.path.exists(save_folder + '{}_rmsds.csv'.format(prot_name)):
				os.system(cmd.format(partition, save_folder + '/' + prot_name, prot_name))
				time.sleep(0.5)
			else:
				print("Exists")

	if task == 'protein':
		protein = sys.argv[2]
		rmsd_file = save_folder + '{}_rmsds.csv'.format(protein)
		print(protein, rmsd_file)
		compute_protein_rmsds(protein, rmsd_file, ASL_to_feature, combind_root)
