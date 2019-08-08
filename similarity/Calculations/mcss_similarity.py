import os
from schrodinger.structure import StructureReader, StructureWriter
import sys
import time

class MCSS:
    """
    Reads and writes MCSS features for a ligand pair.

    There are two key phases of the computation:
        (1) Identification of maximum common substructure(s)
        (2) Computation of RMSDs between the substructures in
            docking results.

    Task (1) is accomplished using Schrodinger's canvasMCSS utility.
    Task (2) is accomplished by identifying all matches of the substructure(s)
    from (1) and finding the pair with the mimimum RMSD. This is a subtlely
    difficult task because of symmetry concerns and details of extracting
    substructures.

    MCSS must be at least half the size of the smaller of the ligands
    or no RMSDs are computed.

    A key design decision is to not specify any file names in this class
    (other than those associated with temp files). The implication of this
    is that MCSSController will be completely in control of this task, while
    this class can be dedicated to actually computing the MCSS feature.
    """

    mcss_cmd = ("$SCHRODINGER/utilities/canvasMCS -imae {} -ocsv {}"
                " -stop {} -atomtype C {}")

    def __init__(self, l1, l2):
        """
        l1, l2: string, ligand names
        """
        if l1 > l2: l1, l2 = l2, l1

        self.l1 = l1
        self.l2 = l2
        self.name = "{}-{}".format(l1, l2)

        self.n_l1_atoms = 0
        self.n_l2_atoms = 0
        self.n_mcss_atoms = 0
        self.smarts_l1 = []
        self.smarts_l2 = []
        self.rmsds = {}

        self.tried_small = False

        # Deprecated.
        self.n_mcss_bonds = 0

    def __str__(self):
        return ','.join(map(str,
                            [self.l1, self.l2,
                             self.n_l1_atoms, self.n_l2_atoms, self.n_mcss_atoms, self.n_mcss_bonds,
                             ';'.join(self.smarts_l1), ';'.join(self.smarts_l2), self.tried_small]
                            ))

    # MCSS Computation Methods.
    def compute_mcss(self, ligands, init_file, mcss_types_file, small=False):
        """
        Compute the MCSS file by calling Schrodinger canvasMCSS.

        Updates instance with MCSSs present in the file
        """
        structure_file = '{}.ligands.mae'.format(init_file)
        mcss_file = '{}.mcss.csv'.format(init_file)
        stwr = StructureWriter(structure_file)
        stwr.append(ligands[self.l1])
        stwr.append(ligands[self.l2])
        stwr.close()
        #set the sizes in atoms of each of the ligands
        self._set_ligand_sizes(structure_file)

        if os.system(self.mcss_cmd.format(structure_file,
                                          mcss_file,
                                          5 if small else 10,
                                          mcss_types_file)):
            assert False, 'MCSS computation failed'
        self._set_mcss(mcss_file)
        self.tried_small = small

        with open(init_file, 'a+') as fp:
            fp.write(str(self) + '\n')

        os.system('rm {} {}'.format(structure_file, mcss_file))

    def _set_ligand_sizes(self, structure_file):
        try:
            refs = [st for st in StructureReader(structure_file)]
        except:
            print('Unable to read MCSS structure file for', self.l1, self.l2)
            return None
        if len(refs) != 2:
            print('Wrong number of structures', self.l1, self.l2)
            return None
        ref1, ref2 = refs
        n_l1_atoms = len([a for a in ref1.atom if a.element != 'H'])
        n_l2_atoms = len([a for a in ref2.atom if a.element != 'H'])

        self.n_l1_atoms = n_l1_atoms
        self.n_l2_atoms = n_l2_atoms

    def _set_mcss(self, mcss_file):
        """
        Updates MCS from the direct output of canvasMCSS.

        Note that there can be multiple maximum common substructures
        of the same size.
        """
        ligs = {}
        n_mcss_atoms = None
        with open(mcss_file) as fp:
            fp.readline()  # Header
            for line in fp:
                smiles, lig, _, _, _, _n_mcss_atoms, _n_mcss_bonds = line.strip().split(',')[:7]
                smarts = line.strip().split(',')[-1]  # There are commas in some of the fields
                _n_mcss_atoms = int(_n_mcss_atoms)

                assert n_mcss_atoms is None or n_mcss_atoms == _n_mcss_atoms, self.name

                if lig not in ligs: ligs[lig] = []
                ligs[lig] += [smarts]
                n_mcss_atoms = _n_mcss_atoms

        if len(ligs) != 2:
            print('Wrong number of ligands in MCSS file', ligs)
            return None
        assert all(smarts for smarts in ligs.values()), ligs

        # MCSS size can change when tautomers change. One particularly prevalent
        # case is when oxyanions are neutralized. Oxyanions are sometimes specified
        # by the smiles string, but nevertheless glide neutralizes them.
        # Can consider initially considering oxyanions and ketones interchangable
        # (see mcss15.typ).
        if self.n_mcss_atoms:
            assert self.n_mcss_atoms <= n_mcss_atoms + 1, 'MCSS size decreased by more than 1.'
            if self.n_mcss_atoms < n_mcss_atoms:
                print(self.name, 'MCSS size increased.')
            if self.n_mcss_atoms > n_mcss_atoms:
                print(self.name, 'MCSS size dencreased by one.')

        self.n_mcss_atoms = n_mcss_atoms


def compute_protein_mcss(protein, ligands, save_folder):
    init_file = '{}/{}_mcss.csv'.format(save_folder, protein)
    for i in range(len(ligands)):
        for j in range(i + 1, len(ligands)):
            l1, l2 = ligands[i], ligands[j]
            l1_path = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(l1, l1)
            l2_path = '/scratch/PI/rondror/combind/bpp_data/MAPK14/ligands/prepared_ligands/{}_lig/{}_lig.mae'.format(l2, l2)
            mcss_types_file = 'mcss_type_file.typ'
            mcss = MCSS(l1, l2)
            with StructureReader(l1_path) as ligand1, StructureReader(l2_path) as ligand2:
                ligands = {l1: next(ligand1), l2: next(ligand2)}
                mcss.compute_mcss(ligands, init_file, mcss_types_file)

def get_ligands(protein, max_ligands, combind_root):
    ligand_folder = combind_root + "/" + protein + "/docking/grids"
    ligands = sorted(os.listdir(ligand_folder))[:max_ligands]  # sorted
    return ligands

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
    result_folder = '/home/users/lxpowers/projects/combind/flexibility/flexibility_project/similarity/Data'
    save_folder = result_folder + '/mcss'
    partition = 'owners'
    max_ligands = 25

    if task == 'all':
        proteins = get_proteins(combind_root)
        # submit jobs for each protein
        cmd = 'sbatch -p {} -t 0:30:00 -o {}.out --wrap="~/miniconda3/bin/python3.4  mcss_similarity.py  protein {}"'
        for prot_name in proteins:
            print(prot_name)
            os.system(cmd.format(partition, save_folder + '/' + prot_name, prot_name))
            time.sleep(0.5)

    if task == 'protein':
        protein = sys.argv[2]
        ligands = get_ligands(protein, max_ligands, combind_root)
        compute_protein_mcss(protein, ligands, save_folder)



