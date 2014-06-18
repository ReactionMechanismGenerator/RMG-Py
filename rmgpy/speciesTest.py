#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.species module.
"""

import unittest

from rmgpy.species import Species
from rmgpy.transport import TransportData
from rmgpy.molecule import Molecule
from rmgpy.thermo import ThermoData
from rmgpy.statmech import Conformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator

################################################################################

class TestSpecies(unittest.TestCase):
    """
    Contains unit tests for the Species class.
    """
    
    def setUp(self):
        """
        A method that is run before each unit test in this class.
        """
        self.species = Species(
            index=1,
            label='C2H4',
            thermo=ThermoData(
                Tdata=([300.0,400.0,500.0,600.0,800.0,1000.0,1500.0],'K'),
                Cpdata=([3.0,4.0,5.0,6.0,8.0,10.0,15.0],'cal/(mol*K)'),
                H298=(-20.0,'kcal/mol'),
                S298=(50.0,'cal/(mol*K)'),
                Tmin=(300.0,'K'),
                Tmax=(2000.0,'K'),
            ),
            conformer=Conformer(
                E0=(0.0,'kJ/mol'),
                modes=[
                    IdealGasTranslation(mass=(28.03,'amu')),
                    NonlinearRotor(inertia=([5.6952e-47, 2.7758e-46, 3.3454e-46],'kg*m^2'), symmetry=1),
                    HarmonicOscillator(frequencies=([834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0],'cm^-1')),
                ],
                spinMultiplicity=1,
                opticalIsomers=1,
            ),
            molecule=[Molecule().fromSMILES('C=C')],
            transportData=TransportData(sigma=(1, 'angstrom'), epsilon=(100, 'K')),
            molecularWeight=(28.03,'amu'),
            reactive=True,
        )
        
    def testPickle(self):
        """
        Test that a Species object can be pickled and unpickled.
        
        ...with no loss of information.
        """
        import cPickle
        species = cPickle.loads(cPickle.dumps(self.species,-1))
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.molecule[0].multiplicity, species.molecule[0].multiplicity)
        self.assertEqual(self.species.thermo.H298.value_si, species.thermo.H298.value_si)
        self.assertEqual(self.species.thermo.H298.units, species.thermo.H298.units)
        self.assertEqual(len(self.species.conformer.modes), len(species.conformer.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].isIsomorphic(species.molecule[0]))
        self.assertEqual(self.species.conformer.E0.value_si, species.conformer.E0.value_si)
        self.assertEqual(self.species.conformer.E0.units, species.conformer.E0.units)
        self.assertEqual(self.species.transportData.sigma.value_si, species.transportData.sigma.value_si)
        self.assertEqual(self.species.transportData.sigma.units, species.transportData.sigma.units)
        self.assertAlmostEqual(self.species.transportData.epsilon.value_si / 1.381e-23, species.transportData.epsilon.value_si / 1.381e-23, 4)
        self.assertEqual(self.species.transportData.epsilon.units, species.transportData.epsilon.units)
        self.assertEqual(self.species.molecularWeight.value_si, species.molecularWeight.value_si)
        self.assertEqual(self.species.molecularWeight.units, species.molecularWeight.units)
        self.assertEqual(self.species.reactive, species.reactive)

    def testOutput(self):
        """
        Test that a Species object can be reconstructed from its repr().
        
        ...with no loss of information.
        """
        species = None
        exec('species = {0!r}'.format(self.species))
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.molecule[0].multiplicity, species.molecule[0].multiplicity)
        self.assertEqual(self.species.thermo.H298.value_si, species.thermo.H298.value_si)
        self.assertEqual(self.species.thermo.H298.units, species.thermo.H298.units)
        self.assertEqual(len(self.species.conformer.modes), len(species.conformer.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].isIsomorphic(species.molecule[0]))
        self.assertEqual(self.species.conformer.E0.value_si, species.conformer.E0.value_si)
        self.assertEqual(self.species.conformer.E0.units, species.conformer.E0.units)
        self.assertEqual(self.species.transportData.sigma.value_si, species.transportData.sigma.value_si)
        self.assertEqual(self.species.transportData.sigma.units, species.transportData.sigma.units)
        self.assertAlmostEqual(self.species.transportData.epsilon.value_si, species.transportData.epsilon.value_si, 3)
        self.assertEqual(self.species.transportData.epsilon.units, species.transportData.epsilon.units)
        self.assertEqual(self.species.molecularWeight.value_si, species.molecularWeight.value_si)
        self.assertEqual(self.species.molecularWeight.units, species.molecularWeight.units)
        self.assertEqual(self.species.reactive, species.reactive)
        
    def testToAdjacencyList(self):
        """
        Test that toAdjacencyList() works as expected.
        """
        string = self.species.toAdjacencyList()
        self.assertTrue(string.startswith(self.species.molecule[0].toAdjacencyList(label=self.species.label,removeH=False)),string)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
