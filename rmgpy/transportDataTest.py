#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This script contains unit test of the :mod: 'rmgpy.transport' module and :mod: 'rmgpy.data.transport' module
"""

import unittest
import os
import rmgpy.constants as constants
from rmgpy import settings
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.quantity import DipoleMoment, Length, Volume, Energy
from rmgpy.transport import TransportData
from rmgpy.data.transport import CriticalPointGroupContribution, TransportDatabase
#################################################################################

class TestTransportData(unittest.TestCase):
    """
    Contains unit test of the :class: 'transportData' class
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.shapeIndex = 1
        self.epsilon = Energy(2.104, 'kJ/mol')
        self.sigma = Length(3.402, 'angstroms')
        self.dipoleMoment = DipoleMoment(1.000, 'C*m')
        self.polarizability = Volume(0.134, 'angstroms^3')
        self.rotrelaxcollnum = 0.000
        self.comment = 'test'

        self.transport = TransportData(
            shapeIndex=self.shapeIndex,
            epsilon=self.epsilon,
            sigma=self.sigma,
            dipoleMoment=self.dipoleMoment,
            polarizability=self.polarizability,
            rotrelaxcollnum=self.rotrelaxcollnum,
            comment=self.comment,
        )

    def test_shapeIndex(self):
        """
        Test that the TransportData shapeIndex property was properly set.
        """
        self.assertAlmostEqual(self.transport.shapeIndex, self.shapeIndex, 6)

    def test_epsilon(self):
        """
        Test that the TransportData epsilon property was properly set.
        """
        self.assertAlmostEqual(self.transport.epsilon.value_si, self.epsilon.value_si, 6)

    def test_sigma(self):
        """
        Test that the TransportData sigma property was properly set.
        """
        self.assertAlmostEqual(self.transport.sigma.value_si * 1e10, self.sigma.value_si * 1e10, 6)

    def test_dipoleMoment(self):
        """
        Test that the TransportData dipoleMoment property was properly set.
        """
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, self.dipoleMoment.value_si, 6)

    def test_polarizability(self):
        """
        Test that the TransportData polarizability property was properly set.
        """
        self.assertAlmostEqual(self.transport.polarizability.value_si, self.polarizability.value_si, 6)

    def test_rotrelaxcollnum(self):
        """
        Test that the TransportData rotrelaxcollnum property was properly set.
        """
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, self.rotrelaxcollnum, 6)

    def test_comment(self):
        """
        Test that the TransportData comment property was properly set.
        """
        self.assertEqual(self.transport.comment, self.comment)

    def test_getCollisionFrequency(self):
        """
        Test the LennardJones.getCollisionFrequency() method.
        """
        T = 1000; P = 1.0e5
        M = P / constants.R / T
        mu = 1.0
        omega = self.transport.getCollisionFrequency(T, M, mu)
        self.assertAlmostEqual(omega / 1.17737e10, 1.0, 4)
 
    def test_pickle(self):
        """
        Test that a TransportData object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        transport = cPickle.loads(cPickle.dumps(self.transport,-1))
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)

    def test_repr(self):
        """
        Test that a TransportData object can be reconstructed from its repr() output with no loss of information
        """
        transport = None
        exec('transport = {0!r}'.format(self.transport))
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)

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
    Contains unit test of the :class: 'TransportDatabase' class
    """

    def setUp(self):
        """
        a function run before each unit test in this class
        """
        self.libraries = ['GRI-Mech', 'PrimaryTransportLibrary']
        self.groups = ['ring', 'nonring']
        self.libraryOrder = []

        path = os.path.join(settings['database.directory'], 'transport')
        self.transportdb = TransportDatabase()
        self.transportdb.load(path, self.libraries)

    def testJoback(self):
        #values calculate from joback's estimations
        self.testCases = [
            ['acetone', 'CC(=O)C', Length(5.36421, 'angstroms'), Energy(3.20446, 'kJ/mol'), "Epsilon & sigma estimated with Tc=500.53 K, Pc=47.11 bar (from Joback method)"],
            ['cyclopenta-1,2-diene', 'C1=C=CCC1', None, None, None],  # not sure what to expect, we just want to make sure it doesn't crash
            ['benzene', 'c1ccccc1', None, None, None],
            ]
        for name, smiles, sigma, epsilon, comment in self.testCases:
            molecule=Molecule(SMILES=smiles)
            species = Species(molecule=[molecule])
            transportData, blank, blank2 = self.transportdb.getTransportPropertiesViaGroupEstimates(species)
            # check Joback worked.
            # If we don't know what to expect, don't check (just make sure we didn't crash)
            if comment:
                self.assertTrue(transportData.comment == comment)
            if sigma:
                self.assertAlmostEqual(transportData.sigma.value_si * 1e10, sigma.value_si * 1e10, 4)
            if epsilon:
                self.assertAlmostEqual(transportData.epsilon.value_si, epsilon.value_si, 1)

    def testJobackOnBenzeneBonds(self):
        "Test Joback doesn't crash on Cb desription of beneze"
        adjlist = """
                    1  C u0 p0 {2,D} {6,S} {7,S}
                    2  C u0 p0 {1,D} {3,S} {8,S}
                    3  C u0 p0 {2,S} {4,D} {9,S}
                    4  C u0 p0 {3,D} {5,S} {10,S}
                    5  C u0 p0 {4,S} {6,D} {11,S}
                    6  C u0 p0 {1,S} {5,D} {12,S}
                    7  H u0 p0 {1,S}
                    8  H u0 p0 {2,S}
                    9  H u0 p0 {3,S}
                    10 H u0 p0 {4,S}
                    11 H u0 p0 {5,S}
                    12 H u0 p0 {6,S}
                    """
        m = Molecule().fromAdjacencyList(adjlist)
        species = Species(molecule=[m])
        transportData, blank, blank2 = self.transportdb.getTransportPropertiesViaGroupEstimates(species)
        self.assertIsNotNone(transportData)

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
