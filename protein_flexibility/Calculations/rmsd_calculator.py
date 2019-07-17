from schrodinger.structure import StructureReader
import schrodinger.structutils.measure as measure
import schrodinger.structutils.rmsd as rmsd
import schrodinger.structutils.analyze as analyze

#how to run this file:
#ml load chemistry
#ml load schrodinger
#$SCHRODINGER/run python3 compute_rmsd_schrodinger.py
#~/miniconda/bin/python3


def get_all_res(s):
	r_list = []
	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				r_list.append((list(r.atom)[0].pdbcode, list(r.atom)[0].chain, list(r.atom)[0].resnum, list(r.atom)[0].inscode))
	return r_list


def map_residues_to_align_index(alignment_string, r_list):
	'''
	Maps unique residue identifiers to list index in alignment string

	alignment_string: (string) output from alignment program, contains one letter codes and dashes
		example: 'TE--S--T-'
	residue_list: (list of anything) unique identifiers of each residue in order of sequence
		number of residues in residue_list must be equal to number of residues in alignment_string
	'''
	r_to_i_map = {}
	counter = 0
	for i in range(len(alignment_string)):
		if counter >= len(r_list):
			break
		if alignment_string[i] == r_list[counter][0]:
			r_to_i_map[r_list[counter]] = i
			counter += 1;
	return r_to_i_map


def inv_map(m):
	return {v: k for k, v in m.items()}

'''
This protocol can be used to analyze each pair of structures
This will work better than defining the residues using all structures and ligands
because it will be more specific to one pair
'''


def get_res_near_ligand(s):
	cutoff = 4
	l = analyze.find_ligands(s)[0]
	#get atom indexes for the land
	close_a_indices = measure.get_atoms_close_to_subset(s, l.atom_indexes, cutoff)
	close_r_set = set([(s.atom[i].pdbcode, s.atom[i].chain, s.atom[i].resnum, s.atom[i].inscode) for i in close_a_indices])
	return close_r_set


def get_atom_list(s, final_r_list):
	complete_a_list = []

	for m in list(s.molecule):
		if len(m.residue) != 1:
			for r in list(m.residue):
				if (list(r.atom)[0].pdbcode, list(r.atom)[0].chain, list(r.atom)[0].resnum, list(r.atom)[0].inscode) in final_r_list:
					complete_a_list.append(r.getAtomList())

	return complete_a_list


with open('paired_seq.txt', "r") as f:
	lines = f.readlines()
	paired_str_s1 = lines[0].strip()
	paired_str_s2 = lines[1].strip()


folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/aligned_files/'
s1 = list(StructureReader(folder+'3HUC/3HUC_out.mae'))[0]
s2 = list(StructureReader(folder+'3GCS/3GCS_out.mae'))[0]
'''
This procedure works for this pair of structures (1KV1 and 3GCS). For most of the others it fails. 
The reason is the residue numbering is not consistent. In some of the structures 
there are multiple residue 169 or 170s. This is because the sequences of some of the 
proteins are slightly different and during the structural alignment that was done on these 
structures, the sequence was changed automatically. We should talk to Joe about how to handle this 
because he was involved in the initial alignment. 
'''

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
		final_r_list_s1.append(r)
		final_r_list_s2.append(i_to_r_map_s2[s1index])

first_s1 = final_r_list_s1
first_s2 = final_r_list_s2

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

for i in range(len(a_list_s1)):
	if len(a_list_s1[i]) == len(a_list_s2[i]):
		final_a_list_s1 += a_list_s1[i]
		final_a_list_s2 += a_list_s2[i]

print('rmsd:', rmsd.superimpose(s1, final_a_list_s1, s2, final_a_list_s2))

# structures = list(StructureReader('1KV1_out.mae'))
# struc = structures[0]
# molecules = list(struc.molecule)  #http://content.schrodinger.com/Docs/r2015-4/python_api/api/schrodinger.structure.Structure-class.html
# lig = molecules[3]
# lig_struc = struc.extract(lig.getAtomList())
# protein_atom_list = molecules[0].getAtomList() + molecules[1].getAtomList() + molecules[2].getAtomList()
# prot_struc = lig_struc = struc.extract(protein_atom_list)

# valid_residue_id_list = sorted([a[1] for a in valid_res])
# valid_residue_id_list_str = [str(i) for i in valid_residue_id_list]
#
# #transform residue numbers to a pymol selection to check in our pymol session
# #this will create selection over residues near ligand
# pymol_res_sele = '+'.join(valid_residue_id_list_str)
# pymol_sele = 'sele (all_MAPK14_pv.1KV1 or all_MAPK14_pv.3GCS) and resi ' + pymol_res_sele
# print(pymol_sele)
#
# schrodinger_res_sele = ' '.join(valid_residue_id_list_str)
# selection_asl = 'NOT atom.element H and chain. A and res. '+schrodinger_res_sele
# print(selection_asl)
# conf_rmsd = rmsd.ConformerRmsd(struc_1, struc_2, asl_expr=selection_asl) # in place, heavy atom RMSD calc.
# conf_rmsd.renumber_structures = True
# st = conf_rmsd.prepareStructure()
# print(st)
# print('rmsd:', conf_rmsd.calculate())
