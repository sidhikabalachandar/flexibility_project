from schrodinger.structure import StructureReader
import schrodinger.structutils.measure as measure
import schrodinger.structutils.rmsd as rmsd
import schrodinger.structutils.analyze as analyze
#how to run this file:
#ml load chemistry
#ml load schrodinger
#$SCHRODINGER/run python3 compute_rmsd_schrodinger.py
#$SCHRODINGER/run ~/miniconda/bin/python3 compute_rmsd_schrodinger.py


'''
This protocol can be used to analyze each pair of structures
This will work better than defining the residues using all structures and ligands
because it will be more specific to one pair
'''

def get_atoms_near_ligand(struc):
	cutoff = 4
	lig = analyze.find_ligands(struc)[0]
	#get atom indexes for the ligand 
	close_atom_indices = measure.get_atoms_close_to_subset(struc, lig.atom_indexes, cutoff)
	close_residue_set = set([(struc.atom[i].pdbcode, struc.atom[i].resnum, struc.atom[i].chain) for i in close_atom_indices])
	return close_residue_set

def struc_has_residue(struc, query):
	try: 
		struc.findResidue(query)
		return True
	except: 
		return False

folder = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/processed_files/'
struc_1 = list(StructureReader(folder+'3HUC/3HUC_out.mae'))[0]
struc_2 = list(StructureReader(folder+'3GCS/3GCS_out.mae'))[0]
'''
This procedure works for this pair of structures (1KV1 and 3GCS). For most of the others it fails. 
The reason is the residue numbering is not consistent. In some of the structures 
there are multiple residue 169 or 170s. This is because the sequences of some of the 
proteins are slightly different and during the structural alignment that was done on these 
structures, the sequence was changed automatically. We should talk to Joe about how to handle this 
because he was involved in the initial alignment. 
'''

close_residue_set1 = get_atoms_near_ligand(struc_1)
close_residue_set2 = get_atoms_near_ligand(struc_2)
combined_residues = close_residue_set1.union(close_residue_set2)

#the selection for rmsd must have all residues present in both structures
#check for each residue id in residue id list if it is present in both structures, otherwise remove it
valid_res = []
for res in combined_residues:
	query = res[2]+':'+str(res[1])
	if struc_has_residue(struc_1, query) and struc_has_residue(struc_2, query):
		valid_res.append(res)

print(valid_res)
valid_residue_id_list = sorted([a[1] for a in valid_res])
valid_residue_id_list_str = [str(i) for i in valid_residue_id_list]

#transform residue numbers to a pymol selection to check in our pymol session 
#this will create selection over residues near ligand
pymol_res_sele = '+'.join(valid_residue_id_list_str)
pymol_sele = 'sele (all_MAPK14_pv.1KV1 or all_MAPK14_pv.3GCS) and resi ' + pymol_res_sele
print(pymol_sele)

schrodinger_res_sele = ' '.join(valid_residue_id_list_str)
selection_asl = 'NOT atom.element H and chain. A and res. '+schrodinger_res_sele
print(selection_asl)
conf_rmsd = rmsd.ConformerRmsd(struc_1, struc_2, asl_expr=selection_asl) # in place, heavy atom RMSD calc.
conf_rmsd.renumber_structures = True
st = conf_rmsd.prepareStructure()
print(st)
print('rmsd:', conf_rmsd.calculate())

