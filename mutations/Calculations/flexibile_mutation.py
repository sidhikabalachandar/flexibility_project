"""
The purpose of this code is to mutate the flexibile residues of each structure of MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 flexible_mutation.py
"""

import os
import pandas as pd
import schrodinger.structure as structure
import schrodinger.structutils.build as build

def get_ligands(protein, combind_root):
    ligand_folder = combind_root + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))# sorted
    return ligands

if __name__ == '__main__':
    combind_root = '/scratch/PI/rondror/combind/bpp_data/'
    save_location = '/home/users/sidhikab/flexibility_project/mutations/Data/MAPK14'
    ligands = get_ligands('MAPK14', combind_root)

    data = pd.read_csv("../../flexibility_prediction/Data/rmsds/MAPK14_rmsds.csv")
    protein = 'MAPK14'
    # for start in ligands:
    start = ligands[0]
    print('Start', start)
    #     for target in ligands:
    target = ligands[1]
    if start != target:
        with structure.StructureWriter('{}/{}_to_{}.mae'.format(save_location, start, target)) as all:
            ending = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, start, start)
            s = list(structure.StructureReader(combind_root + ending))[0]
            pair_data = data[(data['start ligand'] == start) & (data['target ligand'] == target)]
            counter = 0
            mutate_list = []
            for m in list(s.molecule):
                if len(m.residue) != 1:
                    for r in list(m.residue):
                        if counter >= len(pair_data):
                            break
                        res_data = pair_data.iloc[counter]
                        if list(r.atom)[0].chain == "A" and r.pdbres == res_data['name'] and list(r.atom)[0].resnum == res_data['num'] and res_data['rmsd'] > 2:
                            print(list(r.atom)[0].pdbcode, list(r.atom)[0].resnum)
                            mutate_list.append(list(r.atom)[0])
                            counter += 1

            for atom in mutate_list:
                build.mutate(s, atom, 'ALA')


            all.append(s)



