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
                Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
                Cpdata=([3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0], 'cal/(mol*K)'),
                H298=(-20.0, 'kcal/mol'),
                S298=(50.0, 'cal/(mol*K)'),
                Tmin=(300.0, 'K'),
                Tmax=(2000.0, 'K'),
            ),
            conformer=Conformer(
                E0=(0.0, 'kJ/mol'),
                modes=[
                    IdealGasTranslation(mass=(28.03, 'amu')),
                    NonlinearRotor(inertia=([5.6952e-47, 2.7758e-46, 3.3454e-46], 'kg*m^2'), symmetry=1),
                    HarmonicOscillator(frequencies=(
                    [834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3, 3121.6, 3136.7, 3192.5, 3221.0],
                    'cm^-1')),
                ],
                spin_multiplicity=1,
                optical_isomers=1,
            ),
            molecule=[Molecule().from_smiles('C=C')],
            transport_data=TransportData(sigma=(1, 'angstrom'), epsilon=(100, 'K')),
            molecular_weight=(28.03, 'amu'),
            reactive=True,
        )

        self.species2 = Species().from_adjacency_list(
            """
            1  C u0 p0 c0 {2,D} {6,S} {7,S}
            2  C u0 p0 c0 {1,D} {3,S} {8,S}
            3  C u0 p0 c0 {2,S} {4,D} {9,S}
            4  C u0 p0 c0 {3,D} {5,S} {10,S}
            5  C u0 p0 c0 {4,S} {6,D} {11,S}
            6  C u0 p0 c0 {1,S} {5,D} {12,S}
            7  H u0 p0 c0 {1,S}
            8  H u0 p0 c0 {2,S}
            9  H u0 p0 c0 {3,S}
            10 H u0 p0 c0 {4,S}
            11 H u0 p0 c0 {5,S}
            12 H u0 p0 c0 {6,S}
            """)

        self.species3 = Species().from_adjacency_list(
            """
            multiplicity 2
            1 O u1 p2 c0 {3,S}
            2 O u0 p2 c0 {3,D}
            3 N u0 p1 c0 {1,S} {2,D}
            """)

        self.species4 = Species().from_adjacency_list(
            """
            Propane     
            multiplicity 1
            1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
            2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
            3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
            4  H u0 p0 c0 {1,S}
            5  H u0 p0 c0 {1,S}
            6  H u0 p0 c0 {1,S}
            7  H u0 p0 c0 {2,S}
            8  H u0 p0 c0 {2,S}
            9  H u0 p0 c0 {3,S}
            10 H u0 p0 c0 {3,S}
            11 H u0 p0 c0 {3,S}
            """)

    def test_pickle(self):
        """
        Test that a Species object can be pickled and unpickled.
        
        ...with no loss of information.
        """
        import pickle
        species = pickle.loads(pickle.dumps(self.species, -1))
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.molecule[0].multiplicity, species.molecule[0].multiplicity)
        self.assertEqual(self.species.get_thermo_data().H298.value_si, species.get_thermo_data().H298.value_si)
        self.assertEqual(self.species.get_thermo_data().H298.units, species.get_thermo_data().H298.units)
        self.assertEqual(len(self.species.conformer.modes), len(species.conformer.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].is_isomorphic(species.molecule[0]))
        self.assertEqual(self.species.conformer.E0.value_si, species.conformer.E0.value_si)
        self.assertEqual(self.species.conformer.E0.units, species.conformer.E0.units)
        self.assertEqual(self.species.transport_data.sigma.value_si, species.transport_data.sigma.value_si)
        self.assertEqual(self.species.transport_data.sigma.units, species.transport_data.sigma.units)
        self.assertAlmostEqual(self.species.transport_data.epsilon.value_si / 1.381e-23,
                               species.transport_data.epsilon.value_si / 1.381e-23, 4)
        self.assertEqual(self.species.transport_data.epsilon.units, species.transport_data.epsilon.units)
        self.assertEqual(self.species.molecular_weight.value_si, species.molecular_weight.value_si)
        self.assertEqual(self.species.molecular_weight.units, species.molecular_weight.units)
        self.assertEqual(self.species.reactive, species.reactive)

    def test_output(self):
        """
        Test that a Species object can be reconstructed from its repr().
        
        ...with no loss of information.
        """
        namespace = {}
        exec('species = {0!r}'.format(self.species), globals(), namespace)
        self.assertIn('species', namespace)
        species = namespace['species']
        self.assertEqual(self.species.index, species.index)
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.molecule[0].multiplicity, species.molecule[0].multiplicity)
        self.assertEqual(self.species.get_thermo_data().H298.value_si, species.get_thermo_data().H298.value_si)
        self.assertEqual(self.species.get_thermo_data().H298.units, species.get_thermo_data().H298.units)
        self.assertEqual(len(self.species.conformer.modes), len(species.conformer.modes))
        self.assertEqual(len(self.species.molecule), len(species.molecule))
        self.assertTrue(self.species.molecule[0].is_isomorphic(species.molecule[0]))
        self.assertEqual(self.species.conformer.E0.value_si, species.conformer.E0.value_si)
        self.assertEqual(self.species.conformer.E0.units, species.conformer.E0.units)
        self.assertEqual(self.species.transport_data.sigma.value_si, species.transport_data.sigma.value_si)
        self.assertEqual(self.species.transport_data.sigma.units, species.transport_data.sigma.units)
        self.assertAlmostEqual(self.species.transport_data.epsilon.value_si, species.transport_data.epsilon.value_si, 3)
        self.assertEqual(self.species.transport_data.epsilon.units, species.transport_data.epsilon.units)
        self.assertEqual(self.species.molecular_weight.value_si, species.molecular_weight.value_si)
        self.assertEqual(self.species.molecular_weight.units, species.molecular_weight.units)
        self.assertEqual(self.species.reactive, species.reactive)

    def test_equality(self):
        """Test that we can perform equality comparison with Species objects"""
        spc1 = Species(smiles='C')
        spc2 = Species(smiles='C')

        self.assertNotEqual(spc1, spc2)
        self.assertEqual(spc1, spc1)
        self.assertEqual(spc2, spc2)

    def test_less_than(self):
        """Test that we can perform less than comparison with Species objects"""
        spc1 = Species(index=1, label='a', smiles='C')
        spc2 = Species(index=2, label='a', smiles='C')
        spc3 = Species(index=2, label='b', smiles='C')
        spc4 = Species(index=1, label='a', smiles='CC')

        self.assertLess(spc1, spc2)
        self.assertLess(spc1, spc3)
        self.assertLess(spc2, spc3)
        self.assertLess(spc1, spc4)
        self.assertLess(spc2, spc4)
        self.assertLess(spc3, spc4)

    def test_greater_than(self):
        """Test that we can perform greater than comparison with Species objects"""
        spc1 = Species(index=1, label='a', smiles='C')
        spc2 = Species(index=2, label='a', smiles='C')
        spc3 = Species(index=2, label='b', smiles='C')
        spc4 = Species(index=1, label='a', smiles='CC')

        self.assertGreater(spc2, spc1)
        self.assertGreater(spc3, spc1)
        self.assertGreater(spc3, spc2)
        self.assertGreater(spc4, spc1)
        self.assertGreater(spc4, spc2)
        self.assertGreater(spc4, spc3)

    def test_hash(self):
        """Test behavior of Species hashing using dictionaries and sets"""
        spc1 = Species(index=1, label='a', smiles='C')
        spc2 = Species(index=2, label='a', smiles='C')
        spc3 = Species(index=2, label='b', smiles='C')
        spc4 = Species(index=1, label='a', smiles='CC')

        # Test dictionary behavior
        self.assertEqual(len(dict.fromkeys([spc1, spc2, spc3, spc4])), 4)

        # Test set behavior
        self.assertEqual(len({spc1, spc2, spc3, spc4}), 4)

    def test_to_adjacency_list(self):
        """
        Test that to_adjacency_list() works as expected.
        """
        string = self.species.to_adjacency_list()
        self.assertTrue(
            string.startswith(self.species.molecule[0].to_adjacency_list(label=self.species.label, remove_h=False)),
            string)

    def test_species_props(self):
        """
        Test a key-value pair is added to the props attribute of Species.
        """
        self.species.props['foo'] = 'bar'
        self.assertIsInstance(self.species.props, dict)
        self.assertEquals(self.species.props['foo'], 'bar')

    def test_species_props_object_attribute(self):
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

    def test_resonance_isomers_generated(self):
        """Test that 1-penten-3-yl makes 2-penten-1-yl resonance isomer"""
        spec = Species().from_smiles('C=C[CH]CC')
        spec.generate_resonance_structures()
        self.assertEquals(len(spec.molecule), 2)
        self.assertEquals(spec.molecule[1].to_smiles(), "[CH2]C=CCC")

    def test_resonace_isomers_represented(self):
        """Test that both resonance forms of 1-penten-3-yl are printed by __repr__"""
        spec = Species().from_smiles('C=C[CH]CC')
        spec.generate_resonance_structures()
        namespace = {}
        exec('spec2 = {0!r}'.format(spec), globals(), namespace)
        self.assertIn('spec2', namespace)
        spec2 = namespace['spec2']
        self.assertEqual(len(spec.molecule), len(spec2.molecule))
        for i, j in zip(spec.molecule, spec2.molecule):
            self.assertTrue(j.is_isomorphic(i),
                            msg='i is not isomorphic with j, where i is {} and j is {}'.format(i.to_smiles(),
                                                                                               j.to_smiles()))

    def test_is_isomorphic_to_filtered_resonance_structure(self):
        """
        Test that a Species containing a non-representative resonance structure is isomorphic
        with the "correct" Species containing only representative structures (which were not filtered out)

        When generating resonance isomers for N/O/S atoms, a large number of resonance structures per species could
        potentially be generated, yet most are filtered out and only the "correct" / "representative" structures
        are kept. This test makes sure that if a non-representative structure (i.e., a structure that was filtered out)
        is generated, RMG finds the Species it belongs to, if the last exists.
        """

        spc1_correct = Species().from_smiles('[O]N=O')  # check charge separation with higher octet deviation
        spc1_nonrepresentative = Species().from_adjacency_list("""multiplicity 2
                                                                1 N u1 p1 c0 {2,S} {3,S}
                                                                2 O u0 p3 c-1 {1,S}
                                                                3 O u0 p2 c+1 {1,S}""")
        spc2_correct = Species().from_smiles('[N]=NON=O')  # check atoms with val 6
        spc2_nonrepresentative = Species().from_adjacency_list("""multiplicity 2
                                                                1 O u0 p2 c0 {2,S} {3,S}
                                                                2 N u1 p1 c0 {1,S} {4,S}
                                                                3 N u0 p2 c-1 {1,S} {5,S}
                                                                4 N u0 p2 c0 {2,S}
                                                                5 O u0 p2 c+1 {3,S}""")
        spc3_correct = Species().from_smiles('[O]S(O)=O')  # check O4tc penalty
        spc3_nonrepresentative = Species().from_adjacency_list("""multiplicity 2
                                                                1 S u0 p1 c-1 {2,S} {3,S} {4,T}
                                                                2 O u0 p2 c0 {1,S} {5,S}
                                                                3 O u1 p2 c0 {1,S}
                                                                4 O u0 p1 c+1 {1,T}
                                                                5 H u0 p0 c0 {2,S}""")
        spc4_correct = Species().from_smiles('OS(=[N+]=[N-])O')  # check O4dc penalty
        spc4_nonrepresentative = Species().from_adjacency_list("""1 S u0 p0 c+1 {2,D} {3,D} {4,S}
                                                                2 N u0 p1 c0 {1,D} {5,S}
                                                                3 O u0 p1 c+1 {1,D} {6,S}
                                                                4 O u0 p2 c0 {1,S} {7,S}
                                                                5 N u0 p3 c-2 {2,S}
                                                                6 H u0 p0 c0 {3,S}
                                                                7 H u0 p0 c0 {4,S}""")
        spc5_correct = Species().from_smiles('[O][S]')  # checks birad penalty
        spc5_nonrepresentative = Species().from_adjacency_list("""multiplicity 3
                                                                1 O u0 p2 c0 {2,D}
                                                                2 S u2 p1 c0 {1,D}""")
        spc6_correct = Species().from_smiles('[N-]=[N+]=S=S=O')  # checks the S#S case
        spc6_nonrepresentative = Species().from_adjacency_list("""1 S u0 p1 c0 {2,S} {3,T}
                                                                2 N u0 p0 c+1 {1,S} {4,T}
                                                                3 S u0 p1 c0 {1,T} {5,S}
                                                                4 N u0 p1 c0 {2,T}
                                                                5 O u0 p3 c-1 {3,S}""")

        # check that the structures are not isomorphic if resonance structures are not generated:
        self.assertFalse(spc1_correct.is_isomorphic(spc1_nonrepresentative, strict=True))

        # check that the nonrepresentative structure is isomorphic by generating resonance structures:
        self.assertTrue(spc1_correct.is_isomorphic(spc1_nonrepresentative, strict=False))
        self.assertTrue(spc2_correct.is_isomorphic(spc2_nonrepresentative, strict=False))
        self.assertTrue(spc3_correct.is_isomorphic(spc3_nonrepresentative, strict=False))
        self.assertTrue(spc4_correct.is_isomorphic(spc4_nonrepresentative, strict=False))
        self.assertTrue(spc5_correct.is_isomorphic(spc5_nonrepresentative, strict=False))
        self.assertTrue(spc6_correct.is_isomorphic(spc6_nonrepresentative, strict=False))

    def test_get_resonance_hybrid(self):
        """
        Tests that get_resonance_hybrid returns an isomorphic structure
        which has intermediate bond orders.
        
        This check is for C=C[CH]CC which has another resonance structure,
        [CH2]C=CC. When these structures are merged, the bond structure should be,
        C~C~CC, where '~' is a hybrid bond of order 1.5. 
        """
        spec = Species().from_smiles('C=C[CH]CC')
        hybrid_mol = spec.get_resonance_hybrid()

        self.assertTrue(hybrid_mol.to_single_bonds().is_isomorphic(spec.molecule[0].to_single_bonds()))

        # a rough check for intermediate bond orders
        expected_orders = [1, 1.5]
        bonds = []
        # ensure all bond orders are expected
        for atom in hybrid_mol.atoms:
            for atom2 in atom.bonds:
                bond = hybrid_mol.get_bond(atom, atom2)
                self.assertTrue(any([bond.is_order(otherOrder) for otherOrder in expected_orders]),
                                'Unexpected bond order {}'.format(bond.get_order_num()))
                bonds.append(bond)

        # ensure all expected orders are present
        for expected_order in expected_orders:
            self.assertTrue(any([bond.is_order(expected_order) for bond in bonds]),
                            'No bond of order {} found'.format(expected_order))

    def test_copy(self):
        """Test that we can make a copy of a Species object."""

        spc_cp = self.species.copy()

        self.assertTrue(id(self.species) != id(spc_cp))
        self.assertTrue(self.species.is_isomorphic(spc_cp))
        self.assertEquals(self.species.label, spc_cp.label)
        self.assertEquals(self.species.index, spc_cp.index)

        self.assertTrue(self.species.molecular_weight.equals(spc_cp.molecular_weight))
        self.assertEquals(self.species.reactive, spc_cp.reactive)

    def test_cantera(self):
        """
        Test that a Cantera Species object is created correctly.
        """
        from rmgpy.thermo import NASA, NASAPolynomial
        import cantera as ct
        rmg_species = Species(label="Ar", thermo=NASA(
            polynomials=[NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 4.37967], Tmin=(200, 'K'), Tmax=(1000, 'K')),
                         NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, -745.375, 4.37967], Tmin=(1000, 'K'),
                                        Tmax=(6000, 'K'))], Tmin=(200, 'K'), Tmax=(6000, 'K'), comment="""
Thermo library: primaryThermoLibrary
"""), molecule=[Molecule(smiles="[Ar]")], transport_data=TransportData(shapeIndex=0, epsilon=(1134.93, 'J/mol'),
                                                                       sigma=(3.33, 'angstrom'), dipoleMoment=(2, 'De'),
                                                                       polarizability=(1, 'angstrom^3'),
                                                                       rotrelaxcollnum=15.0, comment="""GRI-Mech"""))

        rmg_ct_species = rmg_species.to_cantera(use_chemkin_identifier=True)

        ct_species = ct.Species.fromCti("""species(name=u'Ar',
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
        self.assertEqual(type(rmg_ct_species), type(ct_species))
        self.assertEqual(rmg_ct_species.name, ct_species.name)
        self.assertEqual(rmg_ct_species.composition, ct_species.composition)
        self.assertEqual(rmg_ct_species.size, ct_species.size)
        self.assertEqual(type(rmg_ct_species.thermo), type(ct_species.thermo))
        self.assertEqual(type(rmg_ct_species.transport), type(ct_species.transport))

    def test_get_transport_data(self):
        """
        Test that transport data can be retrieved correctly via the get_transport_data method.
        """

        spc = Species(label="Ar", molecule=[Molecule(smiles="[Ar]")],
                      transport_data=TransportData(shapeIndex=0, epsilon=(1134.93, 'J/mol'), sigma=(3.33, 'angstrom'),
                                                   dipoleMoment=(2, 'De'), polarizability=(1, 'angstrom^3'),
                                                   rotrelaxcollnum=15.0, comment="""GRI-Mech"""))

        self.assertTrue(spc.get_transport_data() is spc.transport_data)

    def test_fingerprint_property(self):
        """Test that the fingerprint property works"""
        self.assertEqual(self.species2.fingerprint, 'C06H06N00O00S00')

    def test_inchi_property(self):
        """Test that the InChI property works"""
        self.assertEqual(self.species2.inchi, 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H')

    def test_multiplicity_property(self):
        """Test that the fingerprint property works"""
        self.assertEqual(self.species2.multiplicity, 1)

    def test_smiles_property(self):
        """Test that the InChI property works"""
        self.assertEqual(self.species2.smiles, 'C1=CC=CC=C1')

    def test_inchi_instantiation(self):
        """Test that we can create a species using the InChI argument"""
        test = Species(inchi='InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H')

        self.assertTrue(test.is_isomorphic(self.species2))

    def test_smiles_instantiation(self):
        """Test that we can create a species using the SMILES argument"""
        test = Species(smiles='C1=CC=CC=C1')

        self.assertTrue(test.is_isomorphic(self.species2))
        
    def test_is_isomorphic_strict(self):
        """Test that the strict argument to Species.is_isomorphic works"""
        spc1 = Species(smiles='[CH2]C1=CC=CC2=C1C=CC1=C2C=CC=C1')
        spc2 = Species(smiles='C=C1C=CC=C2C1=C[CH]C1=C2C=CC=C1')
        spc3 = Species(smiles='[CH2]C1=CC2=C(C=C1)C1=C(C=CC=C1)C=C2')

        self.assertFalse(spc1.is_isomorphic(spc2, strict=True))
        self.assertTrue(spc1.is_isomorphic(spc2, strict=False))
        self.assertFalse(spc1.is_isomorphic(spc3, strict=True))
        self.assertFalse(spc1.is_isomorphic(spc3, strict=False))

        spc1.generate_resonance_structures()
        spc2.generate_resonance_structures()
        spc3.generate_resonance_structures()

        self.assertTrue(spc1.is_isomorphic(spc2, strict=True))
        self.assertTrue(spc1.is_isomorphic(spc2, strict=False))
        self.assertFalse(spc1.is_isomorphic(spc3, strict=True))
        self.assertFalse(spc1.is_isomorphic(spc3, strict=False))

    def test_species_label(self):
        """Test that the species label is not being assigned with the multiplicity string"""
        self.assertEqual(self.species3.label, '')
        self.assertEqual(self.species4.label, 'Propane')


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
