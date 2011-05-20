#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.species module.
"""

import unittest

from rmgpy.molecule import Molecule
from rmgpy.species import *
from rmgpy.thermo import ThermoData
from rmgpy.statmech import *
from rmgpy.quantity import constants

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
                Tdata=[300.0,400.0,500.0,600.0,800.0,1000.0,1500.0],
                Cpdata=[3.0,4.0,5.0,6.0,8.0,10.0,15.0],
                H298=-20.0*4184,
                S298=50.0*4.184,
                Tmin=300.0,
                Tmax=2000.0,
            ),
            states=StatesModel(modes=[
                Translation(mass=0.02803),
                RigidRotor(linear=False, inertia=[5.6952e-47, 2.7758e-46, 3.3454e-46], symmetry=1),
                HarmonicOscillator(frequencies=[834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0]),
            ]),
            molecule=[Molecule().fromSMILES('C=C')],
            E0=0.0,
            lennardJones=LennardJones(sigma=1e-10, epsilon=constants.kB),
            molecularWeight=0.02803,
            reactive=True
        )
        
    def testPickle(self):
        """
        Test that a Species object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        species = cPickle.loads(cPickle.dumps(self.species))
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.thermo.H298.value, species.thermo.H298.value)
        self.assertEqual(self.species.thermo.H298.units, species.thermo.H298.units)
        self.assertEqual(len(self.species.states.modes), len(species.states.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].isIsomorphic(species.molecule[0]))
        self.assertEqual(self.species.E0.value, species.E0.value)
        self.assertEqual(self.species.E0.units, species.E0.units)
        self.assertEqual(self.species.lennardJones.sigma.value, species.lennardJones.sigma.value)
        self.assertEqual(self.species.lennardJones.sigma.units, species.lennardJones.sigma.units)
        self.assertAlmostEqual(self.species.lennardJones.epsilon.value / 1.381e-23, species.lennardJones.epsilon.value / 1.381e-23, 4)
        self.assertEqual(self.species.lennardJones.epsilon.units, species.lennardJones.epsilon.units)
        self.assertEqual(self.species.molecularWeight.value, species.molecularWeight.value)
        self.assertEqual(self.species.molecularWeight.units, species.molecularWeight.units)
        self.assertEqual(self.species.reactive, species.reactive)

    def testOutput(self):
        """
        Test that a Species object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('species = %r' % self.species)
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.thermo.H298.value, species.thermo.H298.value)
        self.assertEqual(self.species.thermo.H298.units, species.thermo.H298.units)
        self.assertEqual(len(self.species.states.modes), len(species.states.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].isIsomorphic(species.molecule[0]))
        self.assertEqual(self.species.E0.value, species.E0.value)
        self.assertEqual(self.species.E0.units, species.E0.units)
        self.assertEqual(self.species.lennardJones.sigma.value, species.lennardJones.sigma.value)
        self.assertEqual(self.species.lennardJones.sigma.units, species.lennardJones.sigma.units)
        self.assertAlmostEqual(self.species.lennardJones.epsilon.value / 1.381e-23, species.lennardJones.epsilon.value / 1.381e-23, 4)
        self.assertEqual(self.species.lennardJones.epsilon.units, species.lennardJones.epsilon.units)
        self.assertEqual(self.species.molecularWeight.value, species.molecularWeight.value)
        self.assertEqual(self.species.molecularWeight.units, species.molecularWeight.units)
        self.assertEqual(self.species.reactive, species.reactive)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
