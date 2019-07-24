"""
The purpose of this code is to delete all current pose viewer files to increase the amount of available space in the home directory
It can be run on sherlock using
$ ~/miniconda/bin/python3 delete_pose_viewer.py
"""

import os
from docking.docking_class import Docking_Set


folder = "/home/users/sidhikab/all_protein_docking"
proteins = os.listdir(folder)
proteins

for protein in proteins:
    if protein[0] != '.':
        dock_set = Docking_Set()
        dock_set.delete_pose_viewer_files(protein)

print("All pose viewer files removed")