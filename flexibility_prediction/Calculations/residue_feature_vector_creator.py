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
            for r in list(m.residue):
                r_dict[r.getAsl()] = (r.temperature_factor, r.temperature_factor / avg, (r.temperature_factor - avg) / sdev, list(r.atom)[0].pdbcode, r.secondary_structure)
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