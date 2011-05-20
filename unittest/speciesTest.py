#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from rmgpy.chem.molecule import Molecule
from rmgpy.chem.species import *
from rmgpy.chem.thermo import ThermoData
from rmgpy.chem.states import *
import rmgpy.chem.constants as constants

################################################################################

class SpeciesTest(unittest.TestCase):
    """
    Contains unit tests for the rmgpy.chem.species module, used for working
    with chemical species objects.
    """
    
    def testPickle(self):
        """
        Test that a Species object can be successfully pickled and
        unpickled with no loss of information.
        """
        spec0 = Species(
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
        import cPickle
        spec = cPickle.loads(cPickle.dumps(spec0))

        self.assertEqual(spec0.index, spec.index)
        self.assertEqual(spec0.label, spec.label)
        self.assertEqual(spec0.thermo.H298, spec.thermo.H298)
        self.assertEqual(len(spec0.states.modes), len(spec.states.modes))
        self.assertEqual(len(spec0.molecule), len(spec.molecule))
        self.assertTrue(spec0.molecule[0].isIsomorphic(spec.molecule[0]))
        self.assertEqual(spec0.E0, spec.E0)
        self.assertEqual(spec0.lennardJones.sigma, spec.lennardJones.sigma)
        self.assertEqual(spec0.lennardJones.epsilon, spec.lennardJones.epsilon)
        self.assertEqual(spec0.molecularWeight, spec.molecularWeight)
        self.assertEqual(spec0.reactive, spec.reactive)

    def testOutput(self):
        """
        Test that a Species object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        spec0 = Species(
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
        exec('spec = %r' % spec0)
        
        self.assertEqual(spec0.index, spec.index)
        self.assertEqual(spec0.label, spec.label)
        self.assertEqual(spec0.thermo.H298, spec.thermo.H298)
        self.assertEqual(len(spec0.states.modes), len(spec.states.modes))
        self.assertEqual(len(spec0.molecule), len(spec.molecule))
        self.assertTrue(spec0.molecule[0].isIsomorphic(spec.molecule[0]))
        self.assertEqual(spec0.E0, spec.E0)
        self.assertAlmostEqual(spec0.lennardJones.sigma, spec.lennardJones.sigma, 6)
        self.assertAlmostEqual(spec0.lennardJones.epsilon, spec.lennardJones.epsilon, 6)
        self.assertAlmostEqual(spec0.molecularWeight, spec.molecularWeight, 6)
        self.assertEqual(spec0.reactive, spec.reactive)


################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
