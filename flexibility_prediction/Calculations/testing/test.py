from unittest import TestCase
import flexibility_prediction.Calculations.residue_feature_vector_creator as residue_feature_vector_creator
from schrodinger.structure import StructureReader, StructureWriter


class TestNormalizedBFactor(TestCase):
    def test_normalizedBFactor(self):
        s = \
            list(StructureReader(
                '/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[
                0]
        m = list(s.molecule)[0]
        residues = list(m.residue)
        normalized_bfactor = residue_feature_vector_creator.normalizedBFactor(residues, 0, 1, 2)
        prevBfactor = residue_feature_vector_creator.normalizedBFactor(residues, -1, 1, 2)
        nextBfactor = residue_feature_vector_creator.normalizedBFactor(residues, 1, 1, 2)
        prev2Bfactor = residue_feature_vector_creator.normalizedBFactor(residues, -2, 1, 2)
        next2Bfactor = residue_feature_vector_creator.normalizedBFactor(residues, 2, 1, 2)

        self.assertEqual(normalized_bfactor, 28.928181818181827)
        self.assertEqual(prevBfactor, 31.760714285714286)
        self.assertEqual(nextBfactor, 28.41357142857143)
        self.assertEqual(prev2Bfactor, 32.150625)
        self.assertEqual(next2Bfactor, 25.675714285714285)


class TestMolecularWeight(TestCase):
    def test_molecular_weight(self):
        s = \
        list(StructureReader('/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[0]
        m = list(s.molecule)[0]
        r = list(m.residue)[0]
        mol_weight = residue_feature_vector_creator.molecularWeight(r)

        self.assertEqual(mol_weight, 159.21265)


class TestRotamers(TestCase):
    def test_rotamers(self):
        # find example where rotamer is high or low, and check that thats true
        # rotamers on surface should have more available rotamers than one on the inside
        rot_s = list(StructureReader('/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[0]
        rot_m = list(rot_s.molecule)[0]
        rot_r = list(rot_m.residue)[0]
        s = list(StructureReader('/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[0]
        m = list(s.molecule)[0]
        r = list(m.residue)[0]
        (num_rots, avg_rot_rmsd, num_r_rots, avg_r_rot_rmsd) = residue_feature_vector_creator.rotamers(rot_s, rot_r, s, r, 50)
        rot_s2 = list(StructureReader(
            '/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[
            0]
        rot_m2 = list(rot_s2.molecule)[18]
        rot_r2 = list(rot_m2.residue)[18]
        s2 = list(StructureReader('/oak/stanford/groups/rondror/projects/combind/bpp_data/MAPK14/structures/aligned_files/1KV1/1KV1_out.mae'))[0]
        m2 = list(s2.molecule)[18]
        r2 = list(m2.residue)[18]
        (num_rots2, avg_rot_rmsd2, num_r_rots2, avg_r_rot_rmsd2) = residue_feature_vector_creator.rotamers(rot_s2, rot_r2, s2, r2, 50)

        self.assertEqual(num_rots, 34)
        self.assertEqual(num_r_rots > num_r_rots2)
