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
def get_all_res(s):
    (avg, sdev) = bfactor_stats(s)
    if sdev == 0:
        return None
    r_dict = {}
    for m in list(s.molecule):
        if len(m.residue) != 1:
            residues = list(m.residue)
            for i in range(len(residues)):
                if residues[i].secondary_structure == -1:
                    continue
                name = residues[i].pdbres
                num = residues[i].resnum
                bfactor = residues[i].temperature_factor
                normalized_bfactor = (residues[i].temperature_factor - avg) / sdev

                if i == 0:
                    prevBfactor = (residues[len(residues) - 1].temperature_factor - avg) / sdev
                    nextBfactor = (residues[i + 1].temperature_factor - avg) / sdev
                elif i == len(residues) - 1:
                    prevBfactor = (residues[i - 1].temperature_factor - avg) / sdev
                    nextBfactor = (residues[0].temperature_factor - avg) / sdev
                else:
                    prevBfactor = (residues[i - 1].temperature_factor - avg) / sdev
                    nextBfactor = (residues[i + 1].temperature_factor - avg) / sdev

                mol_weight = sum(map(lambda x:x.atomic_weight, list(residues[i].atom)))
                sasa = analyze.calculate_sasa(s, residues[i].atom)
                secondary_structure =  residues[i].secondary_structure

                r_dict[residues[i].getAsl()] = (name, num, bfactor, normalized_bfactor, prevBfactor, nextBfactor, mol_weight, sasa, secondary_structure)
    return r_dict


if __name__ == '__main__':
    folder = "/scratch/PI/rondror/combind/bpp_data/"
    proteins = os.listdir(folder)
    protein_dict = {}
    none_counter = 0

    for i, protein in enumerate(proteins):
        #to monitor progress
        print(i, protein)
        if protein[0] != '.':
            ligand_file = folder + protein + "/structures/aligned_files/"
            ligands = os.listdir(ligand_file)
            ligand_dict = {}
            for ligand in ligands:
                struc = list(StructureReader('{}{}/{}_out.mae'.format(ligand_file, ligand, ligand)))[0]
                residues = get_all_res(struc)
                if residues != None:
                    ligand_dict[ligand] = residues
                else:
                    none_counter += 1

            protein_dict[protein] = ligand_dict

    outfile = open('/home/users/sidhikab/flexibility_project/flexibility_prediction/Data/ASL_to_resinfo_dict', 'wb')
    pickle.dump(protein_dict, outfile)
    outfile.close()
    print(none_counter, 'ligands had a standard deviation of 0')