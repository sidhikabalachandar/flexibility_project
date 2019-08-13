import os
import sys
import re

from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid


def make_grids():
    os.system('mkdir -p ../Data/grids')
    root = '../Data/MAPK14/'
    #for s_file in sorted(os.listdir(root)):
    s_file = sorted(os.listdir(root))[0]
    print(s_file)
    s = next(StructureReader('{}/{}'.format(root, s_file)))
    c = get_centroid(s)
    x, y, z = c[:3]

    out_f = s_file[:13]
    os.system('mkdir -p ../Data/grids/{}'.format(out_f))
    with open('../Data/grids/{}/{}.in'.format(out_f, out_f), 'w') as f:
        f.write('GRID_CENTER {},{},{}\n'.format(x, y, z))
        f.write('GRIDFILE {}.zip\n'.format(out_f))
        f.write('INNERBOX 15,15,15\n')
        f.write('OUTERBOX 30,30,30\n')
        f.write('RECEP_FILE ../Data/MAPK14/{}\n'.format(s_file))

    with open('docking/grids/{}/grid_in.sh'.format(out_f), 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('$SCHRODINGER/glide -WAIT {}.in'.format(out_f))

    print('making grid', out_f)
    os.chdir('docking/grids/{}'.format(out_f))
    os.system('sbatch -p owners -t 00:30:00 -o grid.out grid_in.sh')
    os.chdir('../../..')