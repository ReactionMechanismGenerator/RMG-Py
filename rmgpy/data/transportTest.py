#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.quantity import DipoleMoment, Length, Volume, Energy
from rmgpy.data.transport import CriticalPointGroupContribution, TransportDatabase
from rmgpy.transport import TransportData
from rmgpy.species import Species

################################################################################

class TestCriticalPointGroupContribution(unittest.TestCase):
    """
    Contains unit test of the :class: 'criticalPointGroupContribution' class
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tc = 0.0141
        self.Pc = -.0012
        self.Vc = 65
        self.Tb = 23.58
        self.structureIndex = 1

        self.criticalPointContribution = CriticalPointGroupContribution(
            Tc=self.Tc,
            Pc=self.Pc,
            Vc=self.Vc,
            Tb=self.Tb,
            structureIndex=self.structureIndex,
        )

    def test_Tc(self):
        """
        Test that the CriticalPointGroupContribution Tc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Tc, self.Tc, 6)

    def test_Pc(self):
        """
        Test that the CriticalPointGroupContribution Pc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Pc, self.Pc, 6)

    def test_Vc(self):
        """
        Test that the CriticalPointGroupContribution Vc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Vc, self.Vc, 6)

    def test_Tb(self):
        """
        Test that the CriticalPointGroupContribution Tb property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Tb, self.Tb, 6)

    def test_structureIndex(self):
        """
        Test that the CriticalPointGroupContribution structureIndex property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, self.structureIndex, 6)

    def test_pickle(self):
        """
        Test that a CriticalPointGroupContribution object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        criticalPointContribution = cPickle.loads(cPickle.dumps(self.criticalPointContribution,-1))
        self.assertAlmostEqual(self.criticalPointContribution.Tc, criticalPointContribution.Tc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Pc, criticalPointContribution.Pc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Vc, criticalPointContribution.Vc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Tb, criticalPointContribution.Tb, 4)
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, criticalPointContribution.structureIndex, 4)

    def test_repr(self):
        """
        Test that a CriticalPointGroupContribution object can be reconstructed from its repr() output with no loss of information
        """
        criticalPointContribution = None
        exec('criticalPointContribution = {0!r}'.format(self.criticalPointContribution))
        self.assertAlmostEqual(self.criticalPointContribution.Tc, criticalPointContribution.Tc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Pc, criticalPointContribution.Pc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Vc, criticalPointContribution.Vc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Tb, criticalPointContribution.Tb, 4)
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, criticalPointContribution.structureIndex, 4)


class TestTransportDatabase(unittest.TestCase):
    """
    Contains unit tests of the :class:`TransportDatabase` class.
    """
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        self.database = TransportDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'transport'),
                           ['GRI-Mech', 'PrimaryTransportLibrary'])

        self.speciesList = [
            Species().fromSMILES('C'),
            Species().fromSMILES('CCCC'),
            Species().fromSMILES('O'),
            Species().fromSMILES('[CH3]'),
            Species().fromSMILES('[OH]'),
            Species().fromSMILES('c1ccccc1'),
        ]

    def testJoback(self):
        """Test transport property estimation via Joback groups."""
        self.testCases = [
            ['acetone', 'CC(=O)C', Length(5.36421, 'angstroms'), Energy(3.20446, 'kJ/mol'), "Epsilon & sigma estimated with Tc=500.53 K, Pc=47.11 bar (from Joback method)"],
            ['cyclopenta-1,2-diene', 'C1=C=CCC1', None, None, None],  # not sure what to expect, we just want to make sure it doesn't crash
            ['benzene', 'c1ccccc1', None, None, None],
            ]

        #values calculate from joback's estimations
        for name, smiles, sigma, epsilon, comment in self.testCases:
            species = Species().fromSMILES(smiles)
            transportData, blank, blank2 = self.database.getTransportPropertiesViaGroupEstimates(species)
            # check Joback worked.
            # If we don't know what to expect, don't check (just make sure we didn't crash)
            if comment:
                self.assertTrue(transportData.comment == comment)
            if sigma:
                self.assertAlmostEqual(transportData.sigma.value_si * 1e10, sigma.value_si * 1e10, 4)
            if epsilon:
                self.assertAlmostEqual(transportData.epsilon.value_si, epsilon.value_si, 1)

    @work_in_progress
    def testJobackOnBenzeneBonds(self):
        """Test Joback doesn't crash on Cb desription of benzene"""
        species = Species().fromAdjacencyList("""
                                              1  C u0 p0 {2,B} {6,B} {7,S}
                                              2  C u0 p0 {1,B} {3,B} {8,S}
                                              3  C u0 p0 {2,B} {4,B} {9,S}
                                              4  C u0 p0 {3,B} {5,B} {10,S}
                                              5  C u0 p0 {4,B} {6,B} {11,S}
                                              6  C u0 p0 {1,B} {5,B} {12,S}
                                              7  H u0 p0 {1,S}
                                              8  H u0 p0 {2,S}
                                              9  H u0 p0 {3,S}
                                              10 H u0 p0 {4,S}
                                              11 H u0 p0 {5,S}
                                              12 H u0 p0 {6,S}
                                              """)
        transportData, blank, blank2 = self.database.getTransportPropertiesViaGroupEstimates(species)
        self.assertIsNotNone(transportData)

    def testGetTransportProperties(self):
        """Test that we can retrieve best transport properties for a species."""

        for species in self.speciesList:
            transport = self.database.getTransportProperties(species)
            self.assertIsNotNone(transport)
            self.assertTrue(isinstance(transport, tuple))
            self.assertTrue(isinstance(transport[0], TransportData))

    def testGetAllTransportProperties(self):
        """Test that we can retrieve transport properties from all sources for a species.

        Used for transport search on website."""

        for species in self.speciesList:
            transport = self.database.getAllTransportProperties(species)
            self.assertIsNotNone(transport)
            for result in transport:
                self.assertTrue(isinstance(result, tuple))
                self.assertTrue(isinstance(result[0], TransportData))


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

