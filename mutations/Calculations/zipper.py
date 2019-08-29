"""
The purpose of this code is to mutate the flexibile residues of each structure of MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 zipper.py
"""

import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid


if __name__ == '__main__':
    root = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/MAPK14_mut_strucs'
    save = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/grids'
    ligands = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/'

    files = sorted(os.listdir(root))
    grouped_files = []
    n = 8
    os.system('mkdir -p {}'.format(save))
    os.system('mkdir -p {}/run'.format(save))

    for i in range(0, len(files), n):
        grouped_files += [files[i: i+n]]

    for i, group in enumerate(grouped_files):
        print('making grid', i)

        with open('{}/run/grid{}_in.sh'.format(save, i), 'w') as f:

            for s_file in group:
                out_f = s_file[:12]
                os.system('mkdir -p {}/{}'.format(save, out_f))

                with open('{}/{}/{}.in'.format(save, out_f, out_f), 'w') as f_in:

                    if len(s_file) != 16:
                        continue

                    s = next(StructureReader(ligands + s_file[:4] + '_lig.mae'))
                    c = get_centroid(s)
                    x, y, z = c[:3]

                    f_in.write('GRID_CENTER {},{},{}\n'.format(x, y, z))
                    f_in.write('GRIDFILE {}.zip\n'.format(out_f))
                    f_in.write('INNERBOX 15,15,15\n')
                    f_in.write('OUTERBOX 30,30,30\n')
                    f_in.write('RECEP_FILE {}/{}\n'.format(root, s_file))
                    f.write('#!/bin/bash\n')
                    f.write('cd {}/{}\n'.format(save, out_f))
                    f.write('$SCHRODINGER/glide -WAIT {}.in\n'.format(out_f))

        os.chdir('{}/run'.format(save))
        os.system('sbatch -p owners -t 02:00:00 -o grid{}.out grid{}_in.sh'.format(i, i))

