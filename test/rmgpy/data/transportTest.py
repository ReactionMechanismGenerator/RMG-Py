#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import os
import unittest

from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.transport import CriticalPointGroupContribution, TransportDatabase
from rmgpy.molecule import Molecule
from rmgpy.quantity import Energy, Length
from rmgpy.species import Species
from rmgpy.transport import TransportData


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

    def test__tc(self):
        """
        Test that the CriticalPointGroupContribution Tc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Tc, self.Tc, 6)

    def test__pc(self):
        """
        Test that the CriticalPointGroupContribution Pc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Pc, self.Pc, 6)

    def test__vc(self):
        """
        Test that the CriticalPointGroupContribution Vc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Vc, self.Vc, 6)

    def test__tb(self):
        """
        Test that the CriticalPointGroupContribution Tb property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Tb, self.Tb, 6)

    def test_structure_index(self):
        """
        Test that the CriticalPointGroupContribution structureIndex property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, self.structureIndex, 6)

    def test_pickle(self):
        """
        Test that a CriticalPointGroupContribution object can be pickled and unpickled with no loss of information.
        """
        import pickle
        criticalPointContribution = pickle.loads(pickle.dumps(self.criticalPointContribution, -1))
        self.assertAlmostEqual(self.criticalPointContribution.Tc, criticalPointContribution.Tc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Pc, criticalPointContribution.Pc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Vc, criticalPointContribution.Vc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Tb, criticalPointContribution.Tb, 4)
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, criticalPointContribution.structureIndex, 4)

    def test_repr(self):
        """
        Test that a CriticalPointGroupContribution object can be reconstructed from its repr() output with no loss of information
        """
        namespace = {}
        exec('criticalPointContribution = {0!r}'.format(self.criticalPointContribution), globals(), namespace)
        self.assertIn('criticalPointContribution', namespace)
        criticalPointContribution = namespace['criticalPointContribution']
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
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.database = TransportDatabase()
        cls.database.load(os.path.join(settings['database.directory'], 'transport'),
                          ['GRI-Mech', 'PrimaryTransportLibrary'])

        cls.speciesList = [
            Species().from_smiles('C'),
            Species().from_smiles('CCCC'),
            Species().from_smiles('O'),
            Species().from_smiles('[CH3]'),
            Species().from_smiles('[OH]'),
            Species().from_smiles('c1ccccc1'),
        ]

    def test_joback(self):
        """Test transport property estimation via Joback groups."""
        self.testCases = [
            ['acetone', 'CC(=O)C', Length(5.36421, 'angstroms'), Energy(3.20446, 'kJ/mol'), "Epsilon & sigma estimated with Tc=500.53 K, Pc=47.11 bar (from Joback method)"],
            ['cyclopenta-1,2-diene', 'C1=C=CCC1', None, None, None],  # not sure what to expect, we just want to make sure it doesn't crash
            ['benzene', 'c1ccccc1', None, None, None],
            ['N-methylmethanamine', 'CNC', None, None, None],
            ['imidazole', 'c1ncc[nH]1', None, None, None],
        ]

        # values calculate from joback's estimations
        for name, smiles, sigma, epsilon, comment in self.testCases:
            species = Species().from_smiles(smiles)
            transport_data, blank, blank2 = self.database.get_transport_properties_via_group_estimates(species)
            # check Joback worked.
            # If we don't know what to expect, don't check (just make sure we didn't crash)
            if comment:
                self.assertTrue(transport_data.comment == comment)
            if sigma:
                self.assertAlmostEqual(transport_data.sigma.value_si * 1e10, sigma.value_si * 1e10, 4)
            if epsilon:
                self.assertAlmostEqual(transport_data.epsilon.value_si, epsilon.value_si, 1)

    @work_in_progress
    def test_joback_on_benzene_bonds(self):
        """
        Test Joback doesn't crash on Cb desription of benzene

        Notes:
            If this test fails as WIP, then the method ``get_transport_properties_via_group_estimates`` needs to be
            updated so that aromatic molecules are not kekulized. See RMG-Py PR#1936
        """
        molecule = Molecule().from_adjacency_list("""
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
        
        critical_point = self.database.estimate_critical_properties_via_group_additivity(molecule)
        self.assertIsNotNone(critical_point)

    def test_Tb_correction_for_halogens(self):
        """
        Test that the halogen `Tb` correction is applied to the critical point estimated from 
        group additivity
        """
        partial_F_mol1 = Molecule(smiles='CCF') # partially fluorinated without other halogens
        partial_F_mol2 = Molecule(smiles='ClCCF') # partially fluorinated with other halogens
        per_F_mol = Molecule(smiles='FC(F)(F)C(F)(F)F') # perfluorinated
        partial_hal_mol = Molecule(smiles='BrCCCl') # partially halogenated without fluorine
        per_hal_mol1 = Molecule(smiles='BrC(F)(Cl)C(Br)(F)Cl') # perhalogenated with fluorine
        per_hal_mol2 = Molecule(smiles='BrC(Cl)(Cl)C(Cl)(Br)Cl') # perhalogenated without fluorine

        for mol,comment in [
            (partial_F_mol1,'with partial fluorination Tb correction (-25 K)'),
            (partial_F_mol2,'with partial fluorination Tb correction (-25 K)'),
            (per_F_mol,'with perfluorinated Tb correction (-45.57 K)'),
            (partial_hal_mol,'with partial halogenation Tb correction (+11.43 K)'),
            (per_hal_mol1,'with perhalogenated Tb correction (-53.55 K)'),
            (per_hal_mol2,'with perhalogenated Tb correction (-53.55 K)')
        ]:
            critical_point = self.database.estimate_critical_properties_via_group_additivity(mol)
            self.assertEqual(critical_point.comment,comment)

    def test_get_transport_properties(self):
        """Test that we can retrieve best transport properties for a species."""

        for species in self.speciesList:
            transport = self.database.get_transport_properties(species)
            self.assertIsNotNone(transport)
            self.assertTrue(isinstance(transport, tuple))
            self.assertTrue(isinstance(transport[0], TransportData))

    def test_get_all_transport_properties(self):
        """Test that we can retrieve transport properties from all sources for a species.

        Used for transport search on website."""

        for species in self.speciesList:
            transport = self.database.get_all_transport_properties(species)
            self.assertIsNotNone(transport)
            for result in transport:
                self.assertTrue(isinstance(result, tuple))
                self.assertTrue(isinstance(result[0], TransportData))


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
