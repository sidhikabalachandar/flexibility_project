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
from schrodinger.protein import rotamers
import schrodinger.structutils.rmsd as rmsd
from schrodinger.structutils.interactions import steric_clash
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
def get_all_res(s, rot_s):
    cutoff = 2
    (avg, sdev) = bfactor_stats(s)
    if sdev == 0:
        return None
    r_dict = {}
    # for i in range(len(list(s.molecule))):
    i = 0
    m = list(s.molecule)[i]
    rot_m = list(rot_s.molecule)[i]
    residues = list(m.residue)
        # for j in range(len(list(m.residue))):
    j = 0
    r = list(m.residue)[j]
    rot_r = list(rot_m.residue)[j]
    if r.secondary_structure == -1:
        continue
    rotamer_lib = rotamers.Rotamers(rot_s, list(rot_r.atom)[0])
    a_ls = r.getAtomList()
    r_rmsd_ls = []
    rmsd_ls = []
    counter = 0
    for rotamer in list(rotamer_lib.rotamers):
        rotamer.apply()
        rot_a_ls = rot_r.getAtomList()
        no_r_a_ls = [x for x in rot_s.getAtomIndices() if x not in rot_a_ls]
        if steric_clash.clash_volume(rot_s, rot_a_ls, rot_s, no_r_a_ls) > cutoff:
            counter += 1
            r_rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))
        rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))
    print('Num rots without conflict = ', counter)

    name = r.pdbres
    num = r.resnum
    bfactor = r.temperature_factor
    normalized_bfactor = (r.temperature_factor - avg) / sdev
    prevBfactor = (residues[(i - 1) % len(residues)].temperature_factor - avg) / sdev
    nextBfactor = (residues[(i + 1) % len(residues)].temperature_factor - avg) / sdev
    prev2Bfactor = (residues[(i - 2) % len(residues)].temperature_factor - avg) / sdev
    next2Bfactor = (residues[(i + 2) % len(residues)].temperature_factor - avg) / sdev
    mol_weight = sum(map(lambda x: x.atomic_weight, list(r.atom)))
    num_rots = len(rotamer_lib.rotamers)
    avg_rot_rmsd = statistics.mean(rmsd_ls)
    num_r_rots = len(r_rmsd_ls)
    avg_r_rot_rmsd = statistics.mean(r_rmsd_ls)
    sasa = analyze.calculate_sasa(s, r.atom)
    secondary_structure = r.secondary_structure

    r_dict[r.getAsl()] = (name, num, bfactor, normalized_bfactor, prevBfactor, prev2Bfactor,
                          nextBfactor, next2Bfactor, mol_weight, num_rots, avg_rot_rmsd, num_r_rots,
                          avg_r_rot_rmsd, sasa, secondary_structure)
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
    ligands = get_ligands(protein, max_ligands, combind_root)
    ligand_dict = {}
    # for ligand in ligands:
    ligand = '1KV1'
    print(ligand)
    ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, ligand, ligand)
    s = list(StructureReader(combind_root + '/' + ending_1))[0]
    rot_s = list(StructureReader(combind_root + '/' + ending_1))[0]
    residues = get_all_res(s, rot_s)
    if residues != None:
        ligand_dict[ligand] = residues

    # outfile = open(pickle_file, 'wb')
    # pickle.dump(ligand_dict, outfile)
    # outfile.close()


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