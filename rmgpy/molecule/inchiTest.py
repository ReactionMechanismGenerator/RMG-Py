import unittest

from .molecule import Molecule

from .inchi import *

class InChITest(unittest.TestCase):

    def test_constructor(self):
        inchi1 = InChI('InChI=1S/foo')
        self.assertTrue( inchi1 is not None)

        with self.assertRaises(InchiException):
            inchi1 = InChI('foo')

    def test_compare(self):
        inchi1 = InChI('InChI=1S/foo')
        inchi2 = InChI('InChI=1/foo')
        inchi3 = InChI('InChI=1/bar')
        
        self.assertTrue( inchi1 == inchi2)
        self.assertTrue( not inchi1 != inchi2)

        self.assertTrue( (inchi1 < inchi2) is (inchi1 > inchi2))
        self.assertTrue( (inchi1 < inchi3) is not (inchi1 > inchi3))

class AugmentedInChITest(unittest.TestCase):

    def test_constructor(self):
        aug_inchi1 = AugmentedInChI('InChI=1S/foo')

        self.assertTrue( aug_inchi1 == 'foo', aug_inchi1)
        self.assertTrue( aug_inchi1.mult == -1)
        self.assertTrue( aug_inchi1.u_indices is None)

        aug_inchi2 = AugmentedInChI('InChI=1S/foo/mult100')

        self.assertTrue( aug_inchi2 == 'foo/mult100', aug_inchi2)
        self.assertTrue( aug_inchi2.mult == 100)
        self.assertTrue( aug_inchi2.u_indices is None)

        aug_inchi3 = AugmentedInChI('InChI=1S/foo/mult200/u1,3')

        self.assertTrue( aug_inchi3 == 'foo/mult200/u1,3', aug_inchi3)
        self.assertTrue( aug_inchi3.mult == 200)
        self.assertTrue( aug_inchi3.u_indices == [1,3])

    def test_compare(self):
        aug_inchi1 = AugmentedInChI('InChI=1S/foo')
        aug_inchi2 = AugmentedInChI('InChI=1S/foo/mult100')
        aug_inchi3 = AugmentedInChI('InChI=1S/foo/mult200/u1,3')
        aug_inchi4 = AugmentedInChI('InChI=1/foo')
        
        self.assertTrue( aug_inchi1 == aug_inchi4)
        self.assertTrue( not aug_inchi1 != aug_inchi4)

        self.assertTrue( aug_inchi1 != aug_inchi2)
        self.assertTrue( not aug_inchi1 == aug_inchi2        )

        self.assertTrue( aug_inchi3 != aug_inchi2)
        self.assertTrue( not aug_inchi3 == aug_inchi2        )

    def testReduce(self):
        import pickle

        aug_inchi = AugmentedInChI('InChI=1S/foo/mult200/u1,3')
        aug_inchi2 = pickle.loads(pickle.dumps(aug_inchi))

        self.assertTrue( aug_inchi == aug_inchi2)
        self.assertTrue( aug_inchi.mult == aug_inchi2.mult)
        self.assertTrue( aug_inchi.u_indices == aug_inchi2.u_indices)
        

    def testDefaultMultiplicity(self):
        """Test that default multiplicity equals 1."""
        aug_inchi = AugmentedInChI('InChI=1S/CH4/h1H4')

        self.assertTrue( aug_inchi.mult == -1)

class IgnorePrefixTest(unittest.TestCase):

    def test_ignore(self):
        string = 'InChI=1S/foo'
        self.assertTrue( ignore_prefix(string) == 'foo')

        with self.assertRaises(InchiException):
            ignore_prefix('foo')

class ComposeTest(unittest.TestCase):

    def test_compose_aug_inchi(self):
        inchi = 'C2H5/c1-2/h1H2,2H3'
        mult = 2

        aug_inchi = compose_aug_inchi(inchi, U_LAYER_PREFIX + str(mult))
        self.assertTrue( aug_inchi  == INCHI_PREFIX + '/' + inchi + U_LAYER_PREFIX + str(mult), aug_inchi)

class Parse_H_layerTest(unittest.TestCase):

    def test_OCO(self):

        smi = 'O=C-O'
        inchi = Molecule().fromSMILES(smi).toInChI()
        mobile_hs = parse_H_layer(inchi)
        expected = [[2,3]]
        self.assertTrue(mobile_hs == expected)

class Parse_E_LayerTest(unittest.TestCase):
    def test_no_equivalence_layer(self):
        """Test that the absence of an E-layer results in an empty list."""

        auxinfo = "AuxInfo=1/0/N:1/rA:1C/rB:/rC:;"
        e_layer = parse_E_layer(auxinfo)
        self.assertFalse(e_layer)

    def test_C8H22(self):

        auxinfo = "AuxInfo=1/0/N:1,8,4,6,2,7,3,5/E:(1,2)(3,4)(5,6)(7,8)/rA:8C.2C.2CCCCCC/rB:s1;s2;s3;s3;s5;s5;d7;/rC:;;;;;;;;"
        e_layer = parse_E_layer(auxinfo)
        expected = [[1, 2], [3, 4], [5, 6], [7, 8]]
        self.assertTrue(e_layer == expected)

    def test_C7H17(self):

        auxinfo = "AuxInfo=1/0/N:3,5,7,2,4,6,1/E:(1,2,3)(4,5,6)/rA:7CCCCCCC/rB:s1;d2;s1;d4;s1;d6;/rC:;;;;;;;"
        e_layer = parse_E_layer(auxinfo)
        expected = [[1, 2, 3], [4, 5, 6]]
        self.assertTrue(e_layer == expected)


class ParseNLayerTest(unittest.TestCase):
    def test_OCCC(self):
       auxinfo = "AuxInfo=1/0/N:4,3,2,1/rA:4OCCC/rB:s1;s2;s3;/rC:;;;;"
       n_layer = parse_N_layer(auxinfo)
       expected = [4,3,2,1]
       self.assertTrue(n_layer == expected)

if __name__ == '__main__':
    unittest.main()