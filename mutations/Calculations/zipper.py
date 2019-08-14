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
    os.system('mkdir -p ../Data/grids')
    root = '/home/users/sidhikab/flexibility_project/mutations/Data/MAPK14/'
    save = '/home/users/sidhikab/flexibility_project/mutations/Data/grids/'
    ligands = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/'

    files = ['4DLI_to_3O8P.mae']
    for s_file in files:
        print(s_file)
        if os.path.exists(save + s_file[:12]):
            print("Exists")
        else:
            s = next(StructureReader(ligands + s_file[:4] + '_lig.mae'))
            c = get_centroid(s)
            x, y, z = c[:3]

            out_f = s_file[:12]
            os.system('mkdir -p {}/{}'.format(save, out_f))
            with open('{}/{}/{}.in'.format(save, out_f, out_f), 'w') as f:
                f.write('GRID_CENTER {},{},{}\n'.format(x, y, z))
                f.write('GRIDFILE {}.zip\n'.format(out_f))
                f.write('INNERBOX 15,15,15\n')
                f.write('OUTERBOX 30,30,30\n')
                f.write('RECEP_FILE {}/{}\n'.format(root, s_file))

            with open('{}/{}/grid_in.sh'.format(save, out_f), 'w') as f:
                f.write('#!/bin/bash\n')
                f.write('$SCHRODINGER/glide -WAIT {}.in'.format(out_f))

            print('making grid', out_f)
            os.chdir('{}/{}'.format(save, out_f))
            os.system('sbatch -p owners -t 00:30:00 -o grid.out grid_in.sh')

