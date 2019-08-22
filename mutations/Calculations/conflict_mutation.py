"""
The purpose of this code is to mutate the flexibile residues of each structure of MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 conflict_mutation.py
"""

import os
import pandas as pd
import schrodinger.structure as structure
import schrodinger.structutils.build as build
import schrodinger.structutils.interactions.steric_clash as steric_clash
import sys
import time
import pickle

def get_ligands(protein, combind_root):
    ligand_folder = combind_root + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))# sorted
    return ligands

def create_muts(ligand, ligands):
    cutoff = 2
    conflict_dict = {}
    conflict_dict[ligand] = {}
    for target in ligands:
        if ligand != target:
            start_ending = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, ligand, ligand)
            start_s = list(structure.StructureReader(combind_root + start_ending))[0]
            target_ending = '{}/structures/aligned_files/{}/{}_out.mae'.format(protein, target, target)
            target_s = list(structure.StructureReader(combind_root + target_ending))[0]
            mutate_list = []
            mutate_info = []
            for m in list(start_s.molecule):
                if len(m.residue) != 1:
                    for r in list(m.residue):
                        r_atoms = [a.index for a in list(r.atom)]
                        target_lig = [a.index for a in target_s.atom if a.chain == 'L']
                        clash = steric_clash.clash_volume(start_s, r_atoms, target_s, target_lig)
                        if clash > cutoff:
                            mutate_info.append({'Name':r.pdbres, 'Number':r.resnum, 'Clash volume':clash})
                            mutate_list.append(list(r.atom)[0])

            conflict_dict[ligand][target] = mutate_info
    with open('/home/users/sidhikab/flexibility_project/mutations/Data/conflict/' + ligand, 'wb') as outfile:
        pickle.dump(conflict_dict, outfile)

            # for atom in mutate_list:
            #     build.mutate(start_s, atom, 'ALA')
            #
            # prot_st = start_s.extract([a.index for a in start_s.atom if a.chain != 'L'])
            # prot_wr = structure.StructureWriter('{}/{}_to_{}.mae'.format(save_location, target, ligand))
            # prot_wr.append(prot_st)
            # prot_wr.close()


if __name__ == '__main__':
    task = sys.argv[1]
    combind_root = '/scratch/PI/rondror/combind/bpp_data/'
    save_location = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_conflict/MAPK14_mut_strucs'

    partition = 'rondror'
    data = pd.read_csv("/home/users/sidhikab/flexibility_project/flexibility_prediction/Data/rmsds/MAPK14_rmsds.csv")
    protein = 'MAPK14'
    ligands = get_ligands('MAPK14', combind_root)

    if task == 'all':
        #submit jobs for each ligand
        cmd = 'sbatch -p {} -t 1:00:00 -o {}_rmsd.out --wrap="$SCHRODINGER/run python3 conflict_mutation.py ligand {}"'
        for lig_name in ligands:
            print(lig_name)
            os.system(cmd.format(partition, save_location + '/' + lig_name, lig_name))
            time.sleep(0.5)

    if task == 'ligand':
        ligand = sys.argv[2]
        print(ligand)
        create_muts(ligand, ligands)


