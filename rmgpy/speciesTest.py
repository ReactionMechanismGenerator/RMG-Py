#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
        self.assertEqual(self.species.getThermoData().H298.value_si, species.getThermoData().H298.value_si)
        self.assertEqual(self.species.getThermoData().H298.units, species.getThermoData().H298.units)
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
        self.assertEqual(self.species.getThermoData().H298.value_si, species.getThermoData().H298.value_si)
        self.assertEqual(self.species.getThermoData().H298.units, species.getThermoData().H298.units)
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
        spec.generate_resonance_structures()
        self.assertEquals(len(spec.molecule), 2)
        self.assertEquals(spec.molecule[1].toSMILES(), "[CH2]C=CCC")

    def testResonaceIsomersRepresented(self):
        "Test that both resonance forms of 1-penten-3-yl are printed by __repr__"
        spec = Species().fromSMILES('C=C[CH]CC')
        spec.generate_resonance_structures()
        exec('spec2 = {0!r}'.format(spec))
        self.assertEqual(len(spec.molecule), len(spec2.molecule))
        for i, j in zip(spec.molecule, spec2.molecule):
            self.assertTrue(j.isIsomorphic(i), msg='i is not isomorphic with j, where i is {} and j is {}'.format(i.toSMILES(), j.toSMILES()))

    def test_is_isomorphic_to_filtered_resonance_structure(self):
        """
        Test that a Species containing a non-representative resonance structure is isomorphic
        with the "correct" Species containing only representative structures (which were not filtered out)

        When generating resonance isomers for N/O/S atoms, a large number of resonance structures per species could
        potentially be generated, yet most are filtered out and only the "correct" / "representative" structures
        are kept. This test makes sure that if a non-representative structure (i.e., a structure that was filtered out)
        is generated, RMG finds the Species it belongs to, if the last exists.
        """

        spc1_correct = Species().fromSMILES('[O]N=O')  # check charge separation with higher octet deviation
        spc1_nonrepresentative = Species().fromAdjacencyList("""multiplicity 2
                                                                1 N u1 p1 c0 {2,S} {3,S}
                                                                2 O u0 p3 c-1 {1,S}
                                                                3 O u0 p2 c+1 {1,S}""")
        spc2_correct = Species().fromSMILES('[N]=NON=O')  # check atoms with val 6
        spc2_nonrepresentative = Species().fromAdjacencyList("""multiplicity 2
                                                                1 O u0 p2 c0 {2,S} {3,S}
                                                                2 N u1 p1 c0 {1,S} {4,S}
                                                                3 N u0 p2 c-1 {1,S} {5,S}
                                                                4 N u0 p2 c0 {2,S}
                                                                5 O u0 p2 c+1 {3,S}""")
        spc3_correct = Species().fromSMILES('[O]S(O)=O')  # check O4tc penalty
        spc3_nonrepresentative = Species().fromAdjacencyList("""multiplicity 2
                                                                1 S u0 p1 c-1 {2,S} {3,S} {4,T}
                                                                2 O u0 p2 c0 {1,S} {5,S}
                                                                3 O u1 p2 c0 {1,S}
                                                                4 O u0 p1 c+1 {1,T}
                                                                5 H u0 p0 c0 {2,S}""")
        spc4_correct = Species().fromSMILES('N#[N+][S-](O)O')  # check O4dc penalty
        spc4_nonrepresentative = Species().fromAdjacencyList("""1 S u0 p0 c+1 {2,D} {3,D} {4,S}
                                                                2 N u0 p1 c0 {1,D} {5,S}
                                                                3 O u0 p1 c+1 {1,D} {6,S}
                                                                4 O u0 p2 c0 {1,S} {7,S}
                                                                5 N u0 p3 c-2 {2,S}
                                                                6 H u0 p0 c0 {3,S}
                                                                7 H u0 p0 c0 {4,S}""")
        spc5_correct = Species().fromSMILES('[N]=NO')  # check val N > 8 penalty
        spc5_nonrepresentative = Species().fromAdjacencyList("""multiplicity 2
                                                                1 N u1 p0 c0 {2,S} {3,T}
                                                                2 O u0 p2 c0 {1,S} {4,S}
                                                                3 N u0 p1 c0 {1,T}
                                                                4 H u0 p0 c0 {2,S}""")
        spc6_correct = Species().fromSMILES('[O][S]')  # checks birad penalty
        spc6_nonrepresentative = Species().fromAdjacencyList("""multiplicity 3
                                                                1 O u0 p2 c0 {2,D}
                                                                2 S u2 p1 c0 {1,D}""")
        spc7_correct = Species().fromSMILES('N#[N+]SS[O-]')  # checks the S#S case
        spc7_nonrepresentative = Species().fromAdjacencyList("""1 S u0 p1 c0 {2,S} {3,T}
                                                                2 N u0 p0 c+1 {1,S} {4,T}
                                                                3 S u0 p1 c0 {1,T} {5,S}
                                                                4 N u0 p1 c0 {2,T}
                                                                5 O u0 p3 c-1 {3,S}""")

        # check that the structures are not isomorphic if resonance structures are not generated:
        self.assertFalse(spc1_correct.isIsomorphic(spc1_nonrepresentative, generate_res=False))

        # check that the nonrepresentative structure is isomorphic by generating resonance structures:
        self.assertTrue(spc1_correct.isIsomorphic(spc1_nonrepresentative, generate_res=True))
        self.assertTrue(spc2_correct.isIsomorphic(spc2_nonrepresentative, generate_res=True))
        self.assertTrue(spc3_correct.isIsomorphic(spc3_nonrepresentative, generate_res=True))
        self.assertTrue(spc4_correct.isIsomorphic(spc4_nonrepresentative, generate_res=True))
        self.assertTrue(spc5_correct.isIsomorphic(spc5_nonrepresentative, generate_res=True))
        self.assertTrue(spc6_correct.isIsomorphic(spc6_nonrepresentative, generate_res=True))
        self.assertTrue(spc7_correct.isIsomorphic(spc7_nonrepresentative, generate_res=True))

    def testGetResonanceHybrid(self):
        """
        Tests that getResonanceHybrid returns an isomorphic structure
        which has intermediate bond orders.
        
        This check is for C=C[CH]CC which has another resonance structure,
        [CH2]C=CC. When these structures are merged, the bond structure should be,
        C~C~CC, where '~' is a hybrid bond of order 1.5. 
        """
        spec = Species().fromSMILES('C=C[CH]CC')
        hybridMol = spec.getResonanceHybrid()
        
        self.assertTrue(hybridMol.toSingleBonds().isIsomorphic(spec.molecule[0].toSingleBonds()))
        
        # a rough check for intermediate bond orders
        expected_orders = [1,1.5]
        bonds = []
        # ensure all bond orders are expected
        for atom in hybridMol.atoms:
            for atom2 in atom.bonds:
                bond = hybridMol.getBond(atom,atom2)
                self.assertTrue(any([bond.isOrder(otherOrder) for otherOrder in expected_orders]), 'Unexpected bond order {}'.format(bond.getOrderNum()))
                bonds.append(bond)
                
        # ensure all expected orders are present
        for expected_order in expected_orders:
            self.assertTrue(any([bond.isOrder(expected_order) for bond in bonds]),'No bond of order {} found'.format(expected_order))
            
            
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
        
        rmg_ctSpecies = rmgSpecies.toCantera(useChemkinIdentifier = True)
        
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

    def testGetTransportData(self):
        """
        Test that transport data can be retrieved correctly via the getTransportData method.
        """

        spc = Species(label="Ar", molecule=[Molecule(SMILES="[Ar]")], transportData=TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstrom'), dipoleMoment=(2,'De'), polarizability=(1,'angstrom^3'), rotrelaxcollnum=15.0, comment="""GRI-Mech"""))

        self.assertTrue(spc.getTransportData() is spc.transportData)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
