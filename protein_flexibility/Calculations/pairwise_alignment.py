'''
This protocol can be used to find the pairwise alignemnt between the amino acid strings of each pair of proteins
'''
#how to run this file:
# ~/miniconda/bin/python3 pairwise_alignment.py

from Bio import pairwise2
import os
import sys
import pickle
import time

def compute_protein_alignments(protein, seq_file, save_folder):
	infile = open(seq_file,'rb')
	strs = pickle.load(infile)
	infile.close()
	ligands = sorted(strs.keys())
	paired_strs = []
	for i in range(len(ligands)):
		str_s1 = strs[i]
		arr = []
		for j in range(i, len(ligands)):
			str_s2 = strs[j]
			alignments = pairwise2.align.globalxx(str_s1, str_s2)
			arr.append((alignments[0][0], alignments[0][1]))
		paired_strs.append(arr)

	save_file = save_folder + '/{}_alignment'.format(protein)
	outfile = open(save_file, 'wb')
	pickle.dump(paired_strs, outfile)
	outfile.close()

def get_proteins(combind_root):
	'''
	Get the list of all proteins
	:param combind_root: path to the combind root folder
	:return: list of protein name strings
	'''
	proteins = sorted(os.listdir(combind_root))
	proteins = [p for p in proteins if p[0] != '.']
	print(proteins)
	return proteins

if __name__ == '__main__':
	task = sys.argv[1]
	combind_root = '/scratch/PI/rondror/combind/bpp_data'
	result_folder = '/home/users/lxpowers/projects/combind/flexibility/flexibility_project/protein_flexibility/Data'
	save_folder = result_folder+'/alignments'
	partition = 'owners'

	if task == 'all':
		proteins = get_proteins(combind_root)
		#submit jobs for each protein
		cmd = 'sbatch -p {} -t 1:00:00 -o {}.out --wrap="~/miniconda3/bin/python3.4  pairwise_alignment.py  protein {}"'
		for prot_name in proteins:
			os.system(cmd.format(partition, save_folder+'/'+prot_name, prot_name))
			time.sleep(0.5)

	if task == 'protein':
		protein = sys.argv[2]
		seq_file = result_folder + '/sequences/{}_sequences.pkl'.format(protein)
		print(protein, seq_file)
		compute_protein_alignments(protein, seq_file, save_folder)


