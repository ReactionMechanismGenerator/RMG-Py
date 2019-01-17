"""
This scripts runs tests on the database
"""
from rmgpy import settings

import unittest
import rmgpy.molecule.group

import nose
import nose.tools


class TestMoleculeGeneration(unittest.TestCase):
    """
    Contains unit tests for the database for rigorous error checking.
    """
    @classmethod
    def setUpClass(cls):
        """
        Load the database before running the tests.
        """
        pass

    def test_aworking(self):
        g = rmgpy.molecule.group.Group()

        # Cb-Cs from thermo groups
        g.fromAdjacencyList("""
        1 * Cb u0 {2,S}
        2   Cs u0 {1,S}
        """)

        g.makeSampleMolecule()

    def test_broken(self):
        g = rmgpy.molecule.group.Group()

        # Cb-Cl from thermo groups
        g.fromAdjacencyList("""
        1 * Cb u0 {2,S}
        2   Cl1s u0 {1,S}
        """)

        g.makeSampleMolecule()


if __name__ == '__main__':
    nose.run(argv=[__file__, '-v', '--nologcapture'], defaultTest=__name__)
