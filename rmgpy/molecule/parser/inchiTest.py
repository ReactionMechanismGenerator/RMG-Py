import unittest

from rmgpy.molecule.parser.inchi import AugmentedInChI, InChI, ignore_prefix, compose_aug_inchi, INCHI_PREFIX, U_LAYER_PREFIX


class InChITest(unittest.TestCase):

    def test_constructor(self):
        inchi1 = InChI('InChI=1S/foo')
        assert inchi1 is not None

        with self.assertRaises(AssertionError):
            inchi1 = InChI('foo')

    def test_compare(self):
        inchi1 = InChI('InChI=1S/foo')
        inchi2 = InChI('InChI=1/foo')
        inchi3 = InChI('InChI=1/bar')
        
        assert inchi1 == inchi2
        assert not inchi1 != inchi2

        assert (inchi1 < inchi2) is (inchi1 > inchi2)
        assert (inchi1 < inchi3) is not (inchi1 > inchi3)

class AugmentedInChITest(unittest.TestCase):

    def test_constructor(self):
        aug_inchi1 = AugmentedInChI('InChI=1S/foo')

        assert aug_inchi1 == 'foo', aug_inchi1
        assert aug_inchi1.mult == 1
        assert aug_inchi1.u_indices is None

        aug_inchi2 = AugmentedInChI('InChI=1S/foo/mult100')

        assert aug_inchi2 == 'foo/mult100', aug_inchi2
        assert aug_inchi2.mult == 100
        assert aug_inchi2.u_indices is None

        aug_inchi3 = AugmentedInChI('InChI=1S/foo/mult200/u1,3')

        assert aug_inchi3 == 'foo/mult200/u1,3', aug_inchi3
        assert aug_inchi3.mult == 200
        assert aug_inchi3.u_indices == [1,3]

    def test_compare(self):
        aug_inchi1 = AugmentedInChI('InChI=1S/foo')
        aug_inchi2 = AugmentedInChI('InChI=1S/foo/mult100')
        aug_inchi3 = AugmentedInChI('InChI=1S/foo/mult200/u1,3')
        aug_inchi4 = AugmentedInChI('InChI=1/foo')
        
        assert aug_inchi1 == aug_inchi4
        assert not aug_inchi1 != aug_inchi4

        assert aug_inchi1 != aug_inchi2
        assert not aug_inchi1 == aug_inchi2        

        assert aug_inchi3 != aug_inchi2
        assert not aug_inchi3 == aug_inchi2        

    def testReduce(self):
        import pickle

        aug_inchi = AugmentedInChI('InChI=1S/foo/mult200/u1,3')
        aug_inchi2 = pickle.loads(pickle.dumps(aug_inchi))

        assert aug_inchi == aug_inchi2
        assert aug_inchi.mult == aug_inchi2.mult
        assert aug_inchi.u_indices == aug_inchi2.u_indices
        
class IgnorePrefixTest(unittest.TestCase):

    def test_ignore(self):
        string = 'InChI=1S/foo'
        assert ignore_prefix(string) == 'foo'

        with self.assertRaises(AssertionError):
            ignore_prefix('foo')

class ComposeTest(unittest.TestCase):

    def test_compose_aug_inchi(self):
        inchi = 'C2H5/c1-2/h1H2,2H3'
        mult = 2

        aug_inchi = compose_aug_inchi(inchi, U_LAYER_PREFIX + str(mult))
        assert aug_inchi  == INCHI_PREFIX + '/' + inchi + U_LAYER_PREFIX + str(mult), aug_inchi


if __name__ == '__main__':
    unittest.main()