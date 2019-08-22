import os
import sys
import re

from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid


def make_grids():
    root = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/MAPK14_mut_strucs'
    save = '/scratch/PI/rondror/combind/flexibility/MAPK14_mut_pred/grids'
    ligands = '/scratch/PI/rondror/combind/bpp_data/MAPK14/structures/ligands/'

    files = sorted(os.listdir(root))
    grouped_files = []
    n = 15
    os.system('mkdir -p {}'.format(save))

    for s_file in files:
        out_f = s_file[:12]
        s = next(StructureReader(ligands + s_file[:4] + '_lig.mae'))
        c = get_centroid(s)
        x, y, z = c[:3]

        os.system('mkdir -p docking/grids/{}'.format(out_f))
        with open('docking/grids/{}/{}.in'.format(out_f, out_f), 'w') as f:
            f.write('GRID_CENTER {},{},{}\n'.format(x, y, z))
            f.write('GRIDFILE {}.zip\n'.format(out_f))
            f.write('INNERBOX 15,15,15\n')
            f.write('OUTERBOX 30,30,30\n')
            f.write('RECEP_FILE ../../../structures/proteins/{}_prot.mae\n'.format(prot_prefix))

        with open('docking/grids/{}/grid_in.sh'.format(out_f), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('$SCHRODINGER/glide -WAIT {}.in'.format(out_f))

        print('making grid', out_f)
        os.chdir('docking/grids/{}'.format(out_f))
        os.system('sbatch -p owners -t 00:30:00 -o grid.out grid_in.sh')