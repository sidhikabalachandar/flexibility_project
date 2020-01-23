"""
The purpose of this code is to collect the features of each residues of each structure of each protein
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


'''
This function gets the mean and standard deviation of all of the bfactors of a structure
:param s: the protein structure 
:return: the mean and the standard deviation of the list of bfactors associated with the protein structure
'''
def bfactor_stats(s):
    bfactors = []
    for m in list(s.molecule):
        if len(m.residue) != 1:
            for r in list(m.residue):
                bfactors.append(r.temperature_factor)
    return (statistics.mean(bfactors), statistics.stdev(bfactors))


'''
This function gets all of the residues features
:param s: the protein structure 
:param rot_s: the protein structure that will be mutated to different rotamers
:return: the dictionary of each residue's ASL to that residue's features
'''
def get_all_res(s, rot_s):
    cutoff = 50
    (avg, sdev) = bfactor_stats(s)

    if sdev == 0:
        return None

    r_dict = {}
    exception_counter = 0

    for i in range(len(list(s.molecule))):
        m = list(s.molecule)[i]
        rot_m = list(rot_s.molecule)[i]
        residues = list(m.residue)

        for j in range(len(list(m.residue))):
            r = list(m.residue)[j]
            rot_r = list(rot_m.residue)[j]

            if r.secondary_structure == -1:
                continue

            name = r.pdbres
            num = r.resnum
            bfactor = r.temperature_factor
            normalized_bfactor = normalizedBFactor(residues, i, avg, sdev)
            prevBfactor = normalizedBFactor(residues, i - 1, avg, sdev)
            nextBfactor = normalizedBFactor(residues, i + 1, avg, sdev)
            prev2Bfactor = normalizedBFactor(residues, i - 2, avg, sdev)
            next2Bfactor = normalizedBFactor(residues, i + 2, avg, sdev)
            mol_weight = molecularWeight(r)
            (num_rots, avg_rot_rmsd, num_r_rots, avg_r_rot_rmsd) = rotamers(rot_s, rot_r, s, r, cutoff)
            sasa = analyze.calculate_sasa(s, r.atom)
            secondary_structure = r.secondary_structure
            #need to get the xyz position of residue
            position = r.getAlphaCarbon().xyz
            r_dict[r.getAsl()] = (name, num, bfactor, normalized_bfactor, prevBfactor, prev2Bfactor,
                                  nextBfactor, next2Bfactor, mol_weight, num_rots, avg_rot_rmsd, num_r_rots,
                                  avg_r_rot_rmsd, sasa, secondary_structure, position[0], position[1], position[2])

    print("Num exceptions =", exception_counter)
    return r_dict


'''
This function finds the normalized bfactor for a particular residue
:param residues: a list of all residues in the protein structure
:param index: the index of the particular residue in question
:param avg: the average bfactor over all residues in the protein structure
:param sdev: the standard deviation calculated over all residues in the protein structure
:return: the normalized bfactor value
'''
def normalizedBFactor(residues, index, avg, sdev):
    return (residues[index % len(residues)].temperature_factor - avg) / sdev

'''
This function finds the molecular weight of a particular residue
:param residue: a residue object corresponding to the residue in question
:return: the molecular weight of the residue
'''
def molecularWeight(residue):
    return sum(map(lambda x: x.atomic_weight, list(residue.atom)))


'''
This function finds the rotamer information of a particular residue
:param rot_s: a structure object corresponding to the residue in question that will be mutated to the different rotamers
:param rot_r: a structure object corresponding to the residue in question that will not be mutated to the different rotamers
:param s: a structure object corresponding to the protein structure that will be mutated to the different rotamers
:param r: a structure object corresponding to the protein structure that will not be mutated to the different rotamers
:param cutoff: the cutoff for what is considered an acceptable rotamer option
               if the clash of the rotamer and protein structure is less than cutoff then it is considered a viable rotamer option
:return: the total number of available rotamers
         the average rmsd between each available rotamer
         the total number of viable rotamers
         the average rmsd between each viable rotamer
'''
def rotamers(rot_s, rot_r, s, r, cutoff):
    try:
        rotamer_lib = rotamers.Rotamers(rot_s, list(rot_r.atom)[0])
        a_ls = r.getAtomList()
        r_rmsd_ls = []
        rmsd_ls = []
        counter = 0

        for k, rotamer in enumerate(list(rotamer_lib.rotamers)):
            rotamer.apply()
            rot_a_ls = rot_r.getAtomList()
            no_r_a_ls = [a.index for a in rot_s.atom if a.index not in rot_a_ls and a.chain != 'L']
            clash = steric_clash.clash_volume(rot_s, rot_a_ls, rot_s, no_r_a_ls)

            if 'LEU' in r.pdbres and r.resnum == 167:
                print(k, clash)

            if clash < cutoff:
                counter += 1
                r_rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))

            rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))

        num_rots = len(rotamer_lib.rotamers)
        avg_rot_rmsd = safeAvg(num_rots, rmsd_ls)
        num_r_rots = len(r_rmsd_ls)
        avg_r_rot_rmsd = safeAvg(num_r_rots, r_rmsd_ls)

    except Exception as e:
        if 'ALA' not in r.pdbres and 'GLY' not in r.pdbres and 'PRO' not in r.pdbres:
            print(e)

        num_rots = 0
        avg_rot_rmsd = 0
        num_r_rots = 0
        avg_r_rot_rmsd = 0

    return (num_rots, avg_rot_rmsd, num_r_rots, avg_r_rot_rmsd)


'''
This function first checks if the number of elements is greater than 0, and only  then calculates the average
:param num_elems: number of elements in the list
:param ls: the list
:return: 0 if there are no elements in the list
         otherwise the average of the list
'''
def safeAvg(num_elems, ls):
    if num_elems == 0:
        return 0

    else:
        return statistics.mean(ls)


'''
Get the list of all proteins
:param combind_root: path to the combind root folder
:return: list of protein name strings
'''
def get_proteins(combind_root):
	proteins = sorted(os.listdir(combind_root))
	proteins = [p for p in proteins if p[0] != '.']
	print(proteins)
	return proteins


'''
Get the list of all ligands
:param protein: name of the protein
:param max_ligands: maximum number of ligands to analyze for each protein
:param combind_root: path to the combind root folder
:return: list of ligand name strings
'''
def get_ligands(protein, max_ligands, combind_root):
    ligand_folder = combind_root + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
    return ligands


'''
Get the dictionary of features for every protein structure and pickle it
:param protein: name of the protein
:param ligand: name of the ligand
:param pickle_file: path to the save location of the pickled file
:param combind_root: path to the combind root folder
:return: 
'''
def create_feature_vector(protein, ligand, pickle_file, combind_root):
    ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, ligand, ligand)
    s = list(StructureReader(combind_root + ending_1))[0]
    rot_s = list(StructureReader(combind_root + '/' + ending_1))[0]
    residues = get_all_res(s, rot_s)
    print(residues)
    outfile = open(pickle_file, 'wb')
    pickle.dump(residues, outfile)
    outfile.close()


if __name__ == '__main__':
    task = sys.argv[1]
    combind_root = '/oak/stanford/groups/rondror/projects/combind/bpp_data/'
    result_folder = '/home/users/sidhikab/flexibility_project/flexibility_prediction/Data'
    save_folder = result_folder + '/feature_vectors/'
    partition = 'rondror'
    max_ligands = 25

    if task == 'all':
        proteins = get_proteins(combind_root)
        print(proteins)
        proteins = ['MAPK14', 'MEK1'] #'HSP90AA1',
        #submit jobs for each protein
        cmd = 'sbatch -p {} -t 0:20:00 -o {}_rmsd.out --wrap="$SCHRODINGER/run python3 residue_feature_vector_creator.py  protein {} ligand {}"'
        for prot_name in proteins:
            protein_save_folder = save_folder + prot_name
            os.system('mkdir -p {}'.format(protein_save_folder))
            print(prot_name)
            ligands = get_ligands(prot_name, max_ligands, combind_root)
            for lig_name in ligands:
                print(lig_name)
                if not os.path.exists(protein_save_folder + '/' + lig_name):
                    print(cmd.format(partition, protein_save_folder + '/' + lig_name, prot_name, lig_name))
                    os.system(cmd.format(partition, protein_save_folder + '/' + lig_name, prot_name, lig_name))
                    time.sleep(0.5)
                else:
                    print("Exists")

    if task == 'protein':
        protein = sys.argv[2]
        ligand = sys.argv[4]
        pickle_file = save_folder + protein + '/' + ligand + '.pkl'
        print(protein, ligand, pickle_file)
        create_feature_vector(protein, ligand, pickle_file, combind_root)


# """
# The purpose of this code is to collect the features of each residues of each structure of each protein
# It can be run on sherlock using
# $ ml load chemistry
# $ ml load schrodinger
# $ $SCHRODINGER/run python3 residue_feature_vector_creator.py
# """
#
# import os
# from schrodinger.structure import StructureReader, StructureWriter
# import pickle
# import statistics
# import schrodinger.structutils.analyze as analyze
# from schrodinger.protein import rotamers
# import schrodinger.structutils.rmsd as rmsd
# from schrodinger.structutils.interactions import steric_clash
# from schrodinger.structutils import measure
# import sys
# import time
#
#
# '''
# This function gets the mean and standard deviation of all of the bfactors of a structure
# :param s: the protein structure
# :return: the mean and the standard deviation of the list of bfactors associated with the protein structure
# '''
# def bfactor_stats(s):
#     bfactors = []
#     for m in list(s.molecule):
#         if len(m.residue) != 1:
#             for r in list(m.residue):
#                 bfactors.append(r.temperature_factor)
#     return (statistics.mean(bfactors), statistics.stdev(bfactors))
#
#
# '''
# This function gets all of the residues features
# :param s: the protein structure
# :param rot_s: the protein structure that will be mutated to different rotamers
# :return: the dictionary of each residue's ASL to that residue's features
# '''
# def get_all_res(s, rot_s):
#     cutoff = 50
#     (avg, sdev) = bfactor_stats(s)
#
#     if sdev == 0:
#         return None
#
#     r_dict = {}
#     exception_counter = 0
#
#     for i in range(len(list(s.molecule))):
#         m = list(s.molecule)[i]
#         rot_m = list(rot_s.molecule)[i]
#         residues = list(m.residue)
#
#         for j in range(len(list(m.residue))):
#             r = list(m.residue)[j]
#             rot_r = list(rot_m.residue)[j]
#
#             if r.secondary_structure == -1:
#                 continue
#
#             name = r.pdbres
#             num = r.resnum
#             bfactor = r.temperature_factor
#             normalized_bfactor = normalizedBFactor(residues, i, avg, sdev)
#             prevBfactor = normalizedBFactor(residues, i - 1, avg, sdev)
#             nextBfactor = normalizedBFactor(residues, i + 1, avg, sdev)
#             prev2Bfactor = normalizedBFactor(residues, i - 2, avg, sdev)
#             next2Bfactor = normalizedBFactor(residues, i + 2, avg, sdev)
#             mol_weight = molecularWeight(r)
#             (num_rots, avg_rot_rmsd, num_r_rots, avg_r_rot_rmsd) = rotamers(rot_s, rot_r, s, r, cutoff)
#             packingVal = packing(s, r)
#             sasa = analyze.calculate_sasa(s, r.atom)
#             secondary_structure = r.secondary_structure
#
#             r_dict[r.getAsl()] = (name, num, bfactor, normalized_bfactor, prevBfactor, prev2Bfactor,
#                                   nextBfactor, next2Bfactor, mol_weight, num_rots, avg_rot_rmsd, num_r_rots,
#                                   avg_r_rot_rmsd, packingVal, sasa, secondary_structure)
#
#     print("Num exceptions =", exception_counter)
#     return r_dict
#
#
# #questions: should ligand atoms be excluded?
# def packing(struct, residue):
#     cutoff = 5
#     packingVal = 0
#     for a in residue.getAtomList():
#         for neighbor in measure.get_atoms_close_to_subset(struct, [a], cutoff):
#             if struct.atom[neighbor].chain == 'L' or (struct.atom[a].pdbcode == struct.atom[neighbor].pdbcode and
#                 struct.atom[a].chain == struct.atom[neighbor].chain and
#                 struct.atom[a].resnum == struct.atom[neighbor].resnum and
#                 struct.atom[a].inscode == struct.atom[neighbor].inscode):
#                 continue
#             packingVal += 1 / rmsd.calculate_in_place_rmsd(struct, [a], struct, [neighbor])
#
#     return packingVal
#
#
# '''
# This function finds the normalized bfactor for a particular residue
# :param residues: a list of all residues in the protein structure
# :param index: the index of the particular residue in question
# :param avg: the average bfactor over all residues in the protein structure
# :param sdev: the standard deviation calculated over all residues in the protein structure
# :return: the normalized bfactor value
# '''
# def normalizedBFactor(residues, index, avg, sdev):
#     return (residues[index % len(residues)].temperature_factor - avg) / sdev
#
# '''
# This function finds the molecular weight of a particular residue
# :param residue: a residue object corresponding to the residue in question
# :return: the molecular weight of the residue
# '''
# def molecularWeight(residue):
#     return sum(map(lambda x: x.atomic_weight, list(residue.atom)))
#
#
# '''
# This function finds the rotamer information of a particular residue
# :param rot_s: a structure object corresponding to the residue in question that will be mutated to the different rotamers
# :param rot_r: a structure object corresponding to the residue in question that will not be mutated to the different rotamers
# :param s: a structure object corresponding to the protein structure that will be mutated to the different rotamers
# :param r: a structure object corresponding to the protein structure that will not be mutated to the different rotamers
# :param cutoff: the cutoff for what is considered an acceptable rotamer option
#                if the clash of the rotamer and protein structure is less than cutoff then it is considered a viable rotamer option
# :return: the total number of available rotamers
#          the average rmsd between each available rotamer
#          the total number of viable rotamers
#          the average rmsd between each viable rotamer
# '''
# def rotamers(rot_s, rot_r, s, r, cutoff):
#     try:
#         rotamer_lib = rotamers.Rotamers(rot_s, list(rot_r.atom)[0])
#         a_ls = r.getAtomList()
#         r_rmsd_ls = []
#         rmsd_ls = []
#         counter = 0
#
#         for k, rotamer in enumerate(list(rotamer_lib.rotamers)):
#             rotamer.apply()
#             rot_a_ls = rot_r.getAtomList()
#             no_r_a_ls = [a.index for a in rot_s.atom if a.index not in rot_a_ls and a.chain != 'L']
#             clash = steric_clash.clash_volume(rot_s, rot_a_ls, rot_s, no_r_a_ls)
#
#             if 'LEU' in r.pdbres and r.resnum == 167:
#                 print(k, clash)
#
#             if clash < cutoff:
#                 counter += 1
#                 r_rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))
#
#             rmsd_ls.append(rmsd.calculate_in_place_rmsd(s, a_ls, rot_s, rot_a_ls))
#
#         num_rots = len(rotamer_lib.rotamers)
#         avg_rot_rmsd = safeAvg(num_rots, rmsd_ls)
#         num_r_rots = len(r_rmsd_ls)
#         avg_r_rot_rmsd = safeAvg(num_r_rots, r_rmsd_ls)
#
#     except Exception as e:
#         if 'ALA' not in r.pdbres and 'GLY' not in r.pdbres and 'PRO' not in r.pdbres:
#             print(e)
#
#         num_rots = 0
#         avg_rot_rmsd = 0
#         num_r_rots = 0
#         avg_r_rot_rmsd = 0
#
#     return (num_rots, avg_rot_rmsd, num_r_rots, avg_r_rot_rmsd)
#
#
# '''
# This function first checks if the number of elements is greater than 0, and only  then calculates the average
# :param num_elems: number of elements in the list
# :param ls: the list
# :return: 0 if there are no elements in the list
#          otherwise the average of the list
# '''
# def safeAvg(num_elems, ls):
#     if num_elems == 0:
#         return 0
#
#     else:
#         return statistics.mean(ls)
#
#
# '''
# Get the list of all proteins
# :param combind_root: path to the combind root folder
# :return: list of protein name strings
# '''
# def get_proteins(combind_root):
# 	proteins = sorted(os.listdir(combind_root))
# 	proteins = [p for p in proteins if p[0] != '.']
# 	print(proteins)
# 	return proteins
#
#
# '''
# Get the list of all ligands
# :param protein: name of the protein
# :param max_ligands: maximum number of ligands to analyze for each protein
# :param combind_root: path to the combind root folder
# :return: list of ligand name strings
# '''
# def get_ligands(protein, max_ligands, combind_root):
#     ligand_folder = combind_root + protein + "/docking/grids"
#     ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
#     return ligands
#
#
# '''
# Get the dictionary of features for every protein structure and pickle it
# :param protein: name of the protein
# :param ligand: name of the ligand
# :param pickle_file: path to the save location of the pickled file
# :param combind_root: path to the combind root folder
# :return:
# '''
# def create_feature_vector(protein, ligand, pickle_file, combind_root):
#     ending_1 = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, ligand, ligand)
#     s = list(StructureReader(combind_root + ending_1))[0]
#     rot_s = list(StructureReader(combind_root + '/' + ending_1))[0]
#     residues = get_all_res(s, rot_s)
#     outfile = open(pickle_file, 'wb')
#     pickle.dump(residues, outfile)
#     outfile.close()
#
#
# if __name__ == '__main__':
#     task = sys.argv[1]
#     #save folder: scratch/groups/rondror/combind/flexibility
#     combind_root = '/oak/stanford/groups/rondror/projects/combind/bpp_data/'
#     result_folder = '/home/users/sidhikab/flexibility_project/flexibility_prediction/Data'
#     save_folder = result_folder + '/feature_vectors/'
#     #change partition to rondror
#     partition = 'owners'
#     max_ligands = 25
#
#     if task == 'all':
#         proteins = get_proteins(combind_root)
#         #submit jobs for each protein
#         cmd = 'sbatch -p {} -t 1:00:00 -o {}_rmsd.out --wrap="$SCHRODINGER/run python3 residue_feature_vector_creator.py  protein {} ligand {}"'
#         for prot_name in proteins:
#             protein_save_folder = save_folder + prot_name
#             os.system('mkdir -p {}'.format(protein_save_folder))
#             print(prot_name)
#             ligands = get_ligands(prot_name, max_ligands, combind_root)
#             for lig_name in ligands:
#                 print(lig_name)
#                 if not os.path.exists(protein_save_folder + '/' + lig_name):
#                     os.system(cmd.format(partition, protein_save_folder + '/' + lig_name, prot_name, lig_name))
#                     time.sleep(0.5)
#                 else:
#                     print("Exists")
#
#     if task == 'protein':
#         protein = sys.argv[2]
#         ligand = sys.argv[4]
#         pickle_file = save_folder + protein + '/' + ligand + '.pkl'
#         print(protein, ligand, pickle_file)
#         create_feature_vector(protein, ligand, pickle_file, combind_root)