"""
The purpose of this code is to mutate the flexibile residues of each structure of MAPK14
It can be run on sherlock using
$ ml load chemistry
$ ml load schrodinger
$ $SCHRODINGER/run python3 zipper.py
"""


import os
from schrodinger.structure import StructureReader, StructureWriter


if __name__ == '__main__':
    root = '/home/users/sidhikab/flexibility_project/mutations/Data/MAPK14/'

    files = sorted(os.listdir(root))
    for s_file in files:
        if len(s_file) == 16:
            print(s_file)
            prot_path = s_file[:12] + '_nolig.mae'
            if os.path.exists(root + prot_path):
                print("Exists")
            else:
                st = next(StructureReader(root + s_file))
                prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
                prot_wr = StructureWriter(root  + prot_path)
                prot_wr.append(prot_st)
                prot_wr.close()
                os.remove(root + s_file)

