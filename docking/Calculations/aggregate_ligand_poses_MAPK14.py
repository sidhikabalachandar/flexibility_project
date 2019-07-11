"""
The purpose of this code is to collect the poses from crystal structures, glide docking, and combind docking
and aggregate into a single .mae file.
This file uses schrodinger structure library to read and save structure files.
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 aggregate_ligand_poses.py
"""

import os
import score.scores as scores
from schrodinger.structure import StructureReader, StructureWriter

def aggregate_structures(data_path, save_location, docking_version, scores_version):
    """
    :param data_path: (string) location of folder data
    :param save_location: (string) where to save the new .mae files, directory path
    :param docking_version: (string) what docking version to use
    :param scores_version: (string) what score version to sue
    :return:
    """

    protein = 'MAPK14'
    print(protein)

    #try to load the score file
    scores_file = "{}/{}/scores/{}/pdb.sc".format(data_path, protein, scores_version)

    if os.path.isfile(scores_file):
        score_data = scores.load_score(scores_file)

        #figure out the structure that was used for docking
        ligand_folder = '{}/{}/structures/aligned_files'.format(data_path, protein)
        struct = sorted(os.listdir(ligand_folder))[0][:4]
        ligands = os.listdir(ligand_folder)
        #pv_template = "{0:}/{1:}/docking/{2:}/{3:}-to-{4:}/{3:}-to-{4:}_pv.maegz".format(data_path, protein, docking_version, '{0:}', struct)

        # write one mae files for all poses
        with StructureWriter('{}/all_{}_pv.mae'.format(save_location, protein)) as all:

            #iterate over each scored ligand
            for i, ligand in enumerate(ligands):
                #ligand = ligand + '_lig'
                #pose_file = pv_template.format(ligand)
                #if not os.path.isfile(pose_file):
                    #pose_file = "{0:}/{1:}/docking/{2:}/{3:}-to-{4:}/{3:}-to-{4:}_raw.maegz".format(data_path, protein, docking_version, ligand, struct)
                pose_file = "{}/{}/structures/aligned_files/{}/{}_out.mae".format(data_path, protein, ligand, ligand)

                #load the structure file containing docking results
                pv = list(StructureReader(pose_file))
                print(pv)

                #add the protein structure but only for the first file
                all.append(pv[0])

if __name__ == '__main__':
    data_path = '/scratch/PI/rondror/combind/bpp_data'
    save_location = '/home/users/sidhikab'
    docking_version = 'confgen_es4'
    scores_version = 'stats7/pdb/standard/1.0-mcss_contact_hbond_sb'
    aggregate_structures(data_path, save_location, docking_version, scores_version)