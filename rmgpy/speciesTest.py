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
    
    def testSpeciesProps(self):
        """
        Test a key-value pair is added to the props attribute of Species.
        """
        self.species.props['foo'] = 'bar'
        self.assertIsInstance(self.species.props, dict)
        self.assertEquals(self.species.props['foo'], 'bar')
        
    def testSpeciesProps_object_attribute(self):
        """
        Test that Species's props dictionaries are independent of each other.
        
        Create a test in which is checked whether props is an object attribute rather
        than a class attribute
        """
        spc2 = Species()
        self.species.props['foo'] = 'bar'
        spc3 = Species()
        spc3.props['foo'] = 'bla'
        self.assertEquals(self.species.props['foo'], 'bar')
        self.assertDictEqual(spc2.props, {})
        self.assertDictEqual(spc3.props, {'foo': 'bla'})

    def testResonanceIsomersGenerated(self):
        "Test that 1-penten-3-yl makes 2-penten-1-yl resonance isomer"
        spec = Species().fromSMILES('C=C[CH]CC')
        spec.generateResonanceIsomers()
        self.assertEquals(len(spec.molecule), 2)
        self.assertEquals(spec.molecule[1].toSMILES(), "[CH2]C=CCC")

    def testResonaceIsomersRepresented(self):
        "Test that both resonance forms of 1-penten-3-yl are printed by __repr__"
        spec = Species().fromSMILES('C=C[CH]CC')
        spec.generateResonanceIsomers()
        exec('spec2 = {0!r}'.format(spec))
        self.assertEqual(len(spec.molecule), len(spec2.molecule))
        for i, j in zip(spec.molecule, spec2.molecule):
            self.assertTrue(i.isIsomorphic(j))

    def testCopy(self):
        """Test that we can make a copy of a Species object."""

        spc_cp = self.species.copy()

        self.assertTrue(id(self.species) != id(spc_cp))
        self.assertTrue(self.species.isIsomorphic(spc_cp))
        self.assertEquals(self.species.label, spc_cp.label)
        self.assertEquals(self.species.index, spc_cp.index)

        self.assertTrue(self.species.molecularWeight.equals(spc_cp.molecularWeight))
        self.assertEquals(self.species.reactive, spc_cp.reactive)
        
    def testCantera(self):
        """
        Test that a Cantera Species object is created correctly.
        """
        from rmgpy.thermo import NASA, NASAPolynomial
        import cantera as ct
        rmgSpecies = Species(label="Ar", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(SMILES="[Ar]")], transportData=TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstrom'), dipoleMoment=(2,'De'), polarizability=(1,'angstrom^3'), rotrelaxcollnum=15.0, comment="""GRI-Mech"""))
        
        rmg_ctSpecies = rmgSpecies.toCantera()
        
        ctSpecies = ct.Species.fromCti("""species(name=u'Ar',
        atoms='Ar:1',
        thermo=(NASA([200.00, 1000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       4.37967000E+00]),
                NASA([1000.00, 6000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       4.37967000E+00])),
        transport=gas_transport(geom='atom',
                                diam=3.33,
                                well_depth=136.501,
                                dipole=2.0,
                                polar=1.0,
                                rot_relax=15.0))""")
        self.assertEqual(type(rmg_ctSpecies),type(ctSpecies))
        self.assertEqual(rmg_ctSpecies.name, ctSpecies.name)
        self.assertEqual(rmg_ctSpecies.composition, ctSpecies.composition)
        self.assertEqual(rmg_ctSpecies.size, ctSpecies.size)
        self.assertEqual(type(rmg_ctSpecies.thermo), type(ctSpecies.thermo))
        self.assertEqual(type(rmg_ctSpecies.transport), type(ctSpecies.transport))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
