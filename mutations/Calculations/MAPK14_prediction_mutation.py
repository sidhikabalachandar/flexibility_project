"""
The purpose of this code is to select the top 6 residues above a cutoff RMSD from the list of predicted residue RMSDS
These residues are then mutated to alanine
This is repeated for each pair of structures of MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 MAPK14_prediction_mutation.py
"""

import os
import schrodinger.structure as structure
import schrodinger.structutils.build as build
import pickle
import pandas as pd


'''
Get the list of all ligands
:param protein: name of the protein
:param combind_root: path to the combind root folder
:return: list of ligand name strings
'''
def get_ligands(protein, combind_root):
    ligand_folder = combind_root + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))# sorted
    return ligands

if __name__ == '__main__':
    y_pred_data = pd.read_csv('../../flexibility_prediction/Data/MAPK14_y_pred')
    cutoff = 2
    num_mut = 6
    combind_root = '/scratch/PI/rondror/combind/bpp_data/'
    save_location = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/MAPK14_mut_strucs'
    ligands = get_ligands('MAPK14', combind_root)

    protein = 'MAPK14'
    mut_dict = {}
    for start in ligands:
        print('Start', start)
        mut_dict[start] = {}
        for target in ligands:
            if start != target:
                pair_data = y_pred_data[(y_pred_data['start ligand'] == start) & (y_pred_data['target ligand'] == target)]
                pair_data = pair_data.sort_values('y_pred')[(-1 * num_mut):].sort_values('num')
                ending = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, start, start)
                s = list(structure.StructureReader(combind_root + ending))[0]
                counter = 0
                mutate_list = []
                mutate_info = []
                for m in list(s.molecule):
                    if len(m.residue) != 1:
                        for r in list(m.residue):
                            if counter >= len(pair_data):
                                break
                            res_data = pair_data.iloc[counter]
                            if list(r.atom)[0].chain == "A" and r.pdbres == res_data['name'] and r.resnum == res_data['num']:
                                if res_data['y_pred'] > cutoff:
                                    mutate_info.append({'Name':res_data['name'], 'Number':res_data['num'], 'RMSD':res_data['complete rmsd'], 'Y_pred':res_data['y_pred']})
                                    mutate_list.append(list(r.atom)[0])
                                counter += 1

                mut_dict[start][target] = mutate_info

                for atom in mutate_list:
                    build.mutate(s, atom, 'ALA')

                prot_st = s.extract([a.index for a in s.atom if a.chain != 'L'])
                prot_wr = structure.StructureWriter('{}/{}_to_{}.mae'.format(save_location, target, start))
                prot_wr.append(prot_st)
                prot_wr.close()

    with open('/home/users/sidhikab/flexibility_project/mutations/Data/pred', 'wb') as outfile:
        pickle.dump(mut_dict, outfile)


