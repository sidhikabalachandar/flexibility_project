"""
The purpose of this code is to collect the b factor, residue identity, and secondary structure
of each of the residues in the ligand binding pocket of each structure of each protein
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 residue_feature_vector_creator.py
"""

import os
from schrodinger.structure import StructureReader, StructureWriter
import pickle
import statistics
import schrodinger.structutils.analyze as analyze
import sys
import time


#Loop over each protein
#Loop over each structure for each protein
#Loop over each residue near the ligand for each structure
#create a dictionary structured {protein : {ligand :  {ASL : bfactor, normalized bfactor, pdbcode, secondary structure} } }


'''
This function gets the average bfactor of a structure
'''
def bfactor_stats(s):
    bfactors = []
    for m in list(s.molecule):
        if len(m.residue) != 1:
            for r in list(m.residue):
                bfactors.append(r.temperature_factor)
    return (statistics.mean(bfactors), statistics.stdev(bfactors))

'''
This function gets all of the residues bfactors, name, and secondary structure
'''
def get_all_res(s, num_neighbors):
    (avg, sdev) = bfactor_stats(s)
    if sdev == 0:
        return None
    r_dict = {}
    for j, m in enumerate(list(s.molecule)):
        #print('Molecule', j)
        if len(m.residue) != 1:
            residues = list(m.residue)
            for i in range(len(residues)):
                if residues[i].secondary_structure == -1:
                    continue
                name = residues[i].pdbres
                num = residues[i].resnum
                bfactor = residues[i].temperature_factor
                normalized_bfactor = (residues[i].temperature_factor - avg) / sdev
                prevBfactor = (residues[(i - 1) % len(residues)].temperature_factor - avg) / sdev
                nextBfactor = (residues[(i + 1) % len(residues)].temperature_factor - avg) / sdev
                prev2Bfactor = (residues[(i - 2) % len(residues)].temperature_factor - avg) / sdev
                next2Bfactor = (residues[(i + 2) % len(residues)].temperature_factor - avg) / sdev
                mol_weight = sum(map(lambda x:x.atomic_weight, list(residues[i].atom)))
                sasa = analyze.calculate_sasa(s, residues[i].atom)
                secondary_structure =  residues[i].secondary_structure

                r_dict[residues[i].getAsl()] = (name, num, bfactor, normalized_bfactor, prevBfactor, prev2Bfactor, nextBfactor, next2Bfactor, mol_weight, sasa, secondary_structure)
    return r_dict


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

def get_ligands(protein, max_ligands, combind_root):
    ligand_folder = combind_root + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
    return ligands


def create_feature_vector(protein, pickle_file, combind_root):
    max_ligands = 25
    num_neighbors = 2
    ligands = get_ligands(protein, max_ligands, combind_root)
    ligand_dict = {}
    for ligand in ligands:
        print(ligand)
        ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, ligand, ligand)
        struc = list(StructureReader(combind_root + '/' + ending_1))[0]
        residues = get_all_res(struc, num_neighbors)
        if residues != None:
            ligand_dict[ligand] = residues

    outfile = open(pickle_file, 'wb')
    pickle.dump(ligand_dict, outfile)
    outfile.close()


if __name__ == '__main__':
    task = sys.argv[1]
    combind_root = '/scratch/PI/rondror/combind/bpp_data/'
    result_folder = '/home/users/sidhikab/flexibility_project/flexibility_prediction/Data'
    save_folder = result_folder + '/feature_vectors/'
    partition = 'rondror'

    if task == 'all':
        proteins = get_proteins(combind_root)
        #submit jobs for each protein
        cmd = 'sbatch -p {} -t 1:00:00 -o {}_rmsd.out --wrap="$SCHRODINGER/run python3 residue_feature_vector_creator.py  protein {}"'
        for prot_name in proteins:
            print(prot_name)
            if not os.path.exists(save_folder + prot_name):
                os.system(cmd.format(partition, save_folder + '/' + prot_name, prot_name))
                time.sleep(0.5)
            else:
                print("Exists")

    if task == 'protein':
        protein = sys.argv[2]
        pickle_file = save_folder + protein
        print(protein, pickle_file)
        create_feature_vector(protein, pickle_file, combind_root)