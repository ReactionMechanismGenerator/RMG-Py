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
This script contains unit tests of the :mod:`arkane.isodesmic` module.
"""

import unittest

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species

from arkane.encorr.isodesmic import AtomConstraint, BondConstraint, Connection, ErrorCancelingScheme, \
    ErrorCancelingSpecies, ErrorCancelingReaction, IsodesmicScheme, SpeciesConstraints, bond_centric_constraints
from arkane.modelchem import LevelOfTheory

################################################################################


class TestErrorCancelingReactionAndSpecies(unittest.TestCase):
    """
    Tests that ErrorCancelingReaction objects and ErrorCancelingSpecies object are properly implemented
    """

    @classmethod
    def setUpClass(cls):
        """
        A method called before each unit test in this class.
        """
        cls.molecule1 = Molecule(smiles='CC')
        cls.molecule2 = Molecule(smiles='[CH3]')
        cls.species = Species(smiles='CC')

    def test_error_canceling_species(self):
        """
        Test that ErrorCancelingSpecies can be created properly
        """
        lot = LevelOfTheory('test')
        error_canceling_species = ErrorCancelingSpecies(self.molecule1, (123.4, 'kcal/mol'), lot, (100.0, 'kJ/mol'))
        self.assertIsInstance(error_canceling_species, ErrorCancelingSpecies)
        self.assertAlmostEqual(error_canceling_species.low_level_hf298.value_si, 123.4*4184)

        # For target species the high level data is not given
        target_species = ErrorCancelingSpecies(self.molecule2, (10.1, 'J/mol'), lot)
        self.assertIs(target_species.high_level_hf298, None)

    def test_molecule_input_in_error_canceling_species(self):
        """
        Test that an exception is raised if an rmgpy Molecule object is not passed to an ErrorCancelingSpecies
        """
        with self.assertRaises(ValueError):
            ErrorCancelingSpecies(self.species, (100.0, 'J/mol'), LevelOfTheory('test'))

    def test_error_canceling_reactions(self):
        """
        Test that ErrorCancelingReaction object can be created and that hf298 can be calculated for the target
        """
        # Take ethane as the target
        lot = LevelOfTheory('test')
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, 'kJ/mol'), lot)
        methyl = ErrorCancelingSpecies(self.molecule2, (20.0, 'kcal/mol'), lot, (21000.0, 'cal/mol'))

        # This reaction is not an isodesmic reaction, but that does not matter for the unit test
        rxn = ErrorCancelingReaction(ethane, {methyl: 2})
        self.assertAlmostEqual(rxn.calculate_target_thermo().value_si, 2*21000.0*4.184-(2*20.0*4184-100.0*1000))

    def test_level_of_theory_consistency(self):
        """
        Test that ErrorCancelingReaction objects properly check that all species use the same level of theory
        """
        # Take ethane as the target
        ethane = ErrorCancelingSpecies(self.molecule1, (100.0, 'kJ/mol'), LevelOfTheory('test_A'))
        methyl = ErrorCancelingSpecies(self.molecule2, (20.0, 'kcal/mol'), LevelOfTheory('test_B'),
                                       (21000.0, 'cal/mol'))

        # This should throw an exception because the model chemistry is different
        with self.assertRaises(ValueError):
            ErrorCancelingReaction(ethane, {methyl: 2})


class TestConstraintObjects(unittest.TestCase):
    """
    A class for testing AtomConstraint, BondConstraint, and Connection objects
    """
    @classmethod
    def setUpClass(cls):
        # Simple Atoms
        cls.c4_1 = AtomConstraint(label='C4')
        cls.c4_2 = AtomConstraint(label='C4')
        cls.c3 = AtomConstraint(label='C3')
        cls.h1 = AtomConstraint(label='H1')
        cls.o2 = AtomConstraint(label='O2')
        cls.c = AtomConstraint(label='C')
        cls.h = AtomConstraint(label='H')
        cls.o = AtomConstraint(label='O')

        # Simple Bonds
        cls.ch_bond1 = BondConstraint(atom1=cls.c4_1, atom2=cls.h1, bond_order=1)
        cls.ch_bond2 = BondConstraint(atom1=cls.h1, atom2=cls.c4_1, bond_order=1)
        cls.ch_bond3 = BondConstraint(atom1=cls.c4_2, atom2=cls.h1, bond_order=1)
        cls.hh_bond = BondConstraint(atom1=cls.h1, atom2=cls.h1, bond_order=1)

        # Simple Connections
        cls.cs_connect_1 = Connection(atom=cls.c4_1, bond_order=1)
        cls.cs_connect_2 = Connection(atom=cls.c4_2, bond_order=1)
        cls.cd_connect = Connection(atom=cls.c4_1, bond_order=2)
        cls.h_connect = Connection(atom=cls.h1, bond_order=1)
        cls.os_connect = Connection(atom=cls.o2, bond_order=1)

        # Complex Atoms
        cls.c_chho_1 = AtomConstraint(label='C4', connections=[cls.cs_connect_1, cls.h_connect,
                                                               cls.h_connect, cls.os_connect])
        cls.c_chho_2 = AtomConstraint(label='C4', connections=[cls.cs_connect_2, cls.h_connect,
                                                               cls.h_connect, cls.os_connect])
        cls.c_choo = AtomConstraint(label='C4', connections=[cls.cs_connect_1, cls.h_connect,
                                                             cls.os_connect, cls.os_connect])

        # Complex Bonds
        cls.chho_chho = BondConstraint(atom1=cls.c_chho_1, atom2=cls.c_chho_2, bond_order=1)
        cls.chho_choo = BondConstraint(atom1=cls.c_chho_1, atom2=cls.c_choo, bond_order=1)
        cls.choo_chho = BondConstraint(atom1=cls.c_choo, atom2=cls.c_chho_2, bond_order=1)
        cls.chho_choo_d = BondConstraint(atom1=cls.c_chho_1, atom2=cls.c_choo, bond_order=2)

        # Trial molecules
        hf = (100.0, 'J/mol')
        lot = LevelOfTheory('test')
        cls.ethane = ErrorCancelingSpecies(Molecule(smiles='CC'), hf, lot)
        cls.ethanol = ErrorCancelingSpecies(Molecule(smiles='CCO'), hf, lot)
        cls.benzene = ErrorCancelingSpecies(Molecule(smiles='c1ccccc1'), hf, lot)

        # RC2 constraints
        cls.rc2_cc = BondConstraint(atom1=cls.c, atom2=cls.c, bond_order=1)
        cls.rc2_ch = BondConstraint(atom1=cls.c, atom2=cls.h, bond_order=1)
        cls.rc2_co = BondConstraint(atom1=cls.c, atom2=cls.o, bond_order=1)
        cls.rc2_oh = BondConstraint(atom1=cls.o, atom2=cls.h, bond_order=1)
        cls.rc2_ccd = BondConstraint(atom1=cls.c, atom2=cls.c, bond_order=2)

        # RC3 constraints
        cls.rc3_cc = BondConstraint(atom1=cls.c4_1, atom2=cls.c4_1, bond_order=1)
        cls.rc3_c3c3 = BondConstraint(atom1=cls.c3, atom2=cls.c3, bond_order=1)
        cls.rc3_ch = BondConstraint(atom1=cls.c4_1, atom2=cls.h1, bond_order=1)
        cls.rc3_c3h = BondConstraint(atom1=cls.c3, atom2=cls.h1, bond_order=1)
        cls.rc3_co = BondConstraint(atom1=cls.c4_1, atom2=cls.o2, bond_order=1)
        cls.rc3_oh = BondConstraint(atom1=cls.o2, atom2=cls.h1, bond_order=1)
        cls.rc3_ccd = BondConstraint(atom1=cls.c3, atom2=cls.c3, bond_order=2)

        # RC4 constraints
        cls.c_chhh = AtomConstraint(label='C4', connections=[cls.cs_connect_1, cls.h_connect,
                                                             cls.h_connect, cls.h_connect])

        cls.chho_chhh = BondConstraint(atom1=cls.c_chho_1, atom2=cls.c_chhh, bond_order=1)
        cls.chhh_chhh = BondConstraint(atom1=cls.c_chhh, atom2=cls.c_chhh, bond_order=1)

    def test_simple_atom_constraints(self):
        self.assertEqual(self.c4_1, self.c4_2)
        self.assertNotEqual(self.c4_1, self.h1)

    def test_simple_bond_constraints(self):
        self.assertEqual(self.ch_bond1, self.ch_bond2)
        self.assertEqual(self.ch_bond1, self.ch_bond3)
        self.assertNotEqual(self.ch_bond1, self.hh_bond)

    def test_connection(self):
        self.assertEqual(self.cs_connect_1, self.cs_connect_2)
        self.assertNotEqual(self.cs_connect_1, self.cd_connect)
        self.assertNotEqual(self.cs_connect_1, self.h_connect)

    def test_multiple_connection_constraints(self):
        self.assertEqual(self.c_chho_1, self.c_chho_2)
        self.assertNotEqual(self.c_chho_1, self.c_choo)

        self.assertEqual(self.chho_choo, self.choo_chho)
        self.assertNotEqual(self.chho_choo, self.chho_chho)
        self.assertNotEqual(self.chho_choo, self.chho_choo_d)

    def test_constraint_repr(self):
        self.assertEqual(repr(self.chho_choo), 'C4(-C4)(-H1)(-H1)(-O2)-C4(-C4)(-H1)(-O2)(-O2)')

    def test_bond_centric_constraints(self):
        constraints = bond_centric_constraints(species=self.ethane, constraint_class='rc2')
        self.assertEqual(len(constraints), 7)
        self.assertIn(self.rc2_cc, constraints)
        self.assertIn(self.rc2_ch, constraints)

        constraints = bond_centric_constraints(species=self.ethanol, constraint_class='rc2')
        self.assertEqual(len(constraints), 8)
        self.assertIn(self.rc2_cc, constraints)
        self.assertIn(self.rc2_ch, constraints)
        self.assertIn(self.rc2_co, constraints)
        self.assertIn(self.rc2_oh, constraints)

        constraints = bond_centric_constraints(species=self.benzene, constraint_class='rc2')
        self.assertEqual(len(constraints), 12)
        self.assertIn(self.rc2_cc, constraints)
        self.assertIn(self.rc2_ch, constraints)
        self.assertIn(self.rc2_ccd, constraints)

    def test_buerger_rc3(self):
        constraints = bond_centric_constraints(species=self.ethane, constraint_class='rc3')
        self.assertEqual(len(constraints), 7)
        self.assertIn(self.rc3_cc, constraints)
        self.assertIn(self.rc3_ch, constraints)

        constraints = bond_centric_constraints(species=self.ethanol, constraint_class='rc3')
        self.assertEqual(len(constraints), 8)
        self.assertIn(self.rc3_cc, constraints)
        self.assertIn(self.rc3_ch, constraints)
        self.assertIn(self.rc3_co, constraints)
        self.assertIn(self.rc3_oh, constraints)

        constraints = bond_centric_constraints(species=self.benzene, constraint_class='rc3')
        self.assertEqual(len(constraints), 12)
        self.assertIn(self.rc3_c3c3, constraints)
        self.assertIn(self.rc3_c3h, constraints)
        self.assertIn(self.rc3_ccd, constraints)

    def test_buerger_rc4(self):
        constraints = bond_centric_constraints(species=self.ethane, constraint_class='rc4')
        self.assertEqual(len(constraints), 7)
        self.assertIn(self.chhh_chhh, constraints)

        constraints = bond_centric_constraints(species=self.ethanol, constraint_class='rc4')
        self.assertEqual(len(constraints), 8)
        self.assertIn(self.chho_chhh, constraints)


class TestErrorCancelingScheme(unittest.TestCase):
    """
    A class for testing that the ErrorCancelingScheme class functions properly
    """

    @classmethod
    def setUpClass(cls):
        try:
            #import pyomo as pyo
            pass
        except ImportError:
            pyo = None
        cls.pyo = None

        lot = LevelOfTheory('test')
        cls.propene = ErrorCancelingSpecies(Molecule(smiles='CC=C'), (100, 'kJ/mol'), lot, (105, 'kJ/mol'))
        cls.propane = ErrorCancelingSpecies(Molecule(smiles='CCC'), (75, 'kJ/mol'), lot, (80, 'kJ/mol'))
        cls.butane = ErrorCancelingSpecies(Molecule(smiles='CCCC'), (150, 'kJ/mol'), lot, (145, 'kJ/mol'))
        cls.butene = ErrorCancelingSpecies(Molecule(smiles='C=CCC'), (175, 'kJ/mol'), lot, (180, 'kJ/mol'))
        cls.pentane = ErrorCancelingSpecies(Molecule(smiles='CCCCC'), (200, 'kJ/mol'), lot, (190, 'kJ/mol'))
        cls.pentene = ErrorCancelingSpecies(Molecule(smiles='C=CCCC'), (225, 'kJ/mol'), lot, (220, 'kJ/mol'))
        cls.hexane = ErrorCancelingSpecies(Molecule(smiles='CCCCCC'), (250, 'kJ/mol'), lot, (260, 'kJ/mol'))
        cls.hexene = ErrorCancelingSpecies(Molecule(smiles='C=CCCCC'), (275, 'kJ/mol'), lot, (275, 'kJ/mol'))
        cls.benzene = ErrorCancelingSpecies(Molecule(smiles='c1ccccc1'), (-50, 'kJ/mol'), lot, (-80, 'kJ/mol'))
        cls.caffeine = ErrorCancelingSpecies(Molecule(smiles='CN1C=NC2=C1C(=O)N(C(=O)N2C)C'), (300, 'kJ/mol'), lot)
        cls.ethyne = ErrorCancelingSpecies(Molecule(smiles='C#C'), (200, 'kJ/mol'), lot)

    def test_creating_error_canceling_schemes(self):
        scheme = ErrorCancelingScheme(self.propene, [self.butane, self.benzene, self.caffeine, self.ethyne],
                                      isodesmic_class='rc2', conserve_ring_size=True, limit_charges=True,
                                      limit_scope=True)

        self.assertEqual(scheme.reference_species, [self.butane])

        isodesmic_scheme = IsodesmicScheme(self.propene, [self.butane, self.benzene, self.caffeine, self.ethyne])

        self.assertEqual(isodesmic_scheme.reference_species, [self.butane, self.benzene])

    def test_find_error_canceling_reaction(self):
        """
        Test that the MILP problem can be solved to find a single isodesmic reaction
        """
        scheme = IsodesmicScheme(self.propene, [self.propane, self.butane, self.butene, self.caffeine, self.ethyne])

        # Note that caffeine and ethyne will not be allowed, so for the full set the indices are [0, 1, 2]
        rxn, _ = scheme._find_error_canceling_reaction([0, 1, 2], milp_software=['lpsolve'])
        self.assertEqual(rxn.species[self.butane], -1)
        self.assertEqual(rxn.species[self.propane], 1)
        self.assertEqual(rxn.species[self.butene], 1)

        if self.pyo is not None:
            rxn, _ = scheme._find_error_canceling_reaction([0, 1, 2], milp_software=['pyomo'])
            self.assertEqual(rxn.species[self.butane], -1)
            self.assertEqual(rxn.species[self.propane], 1)
            self.assertEqual(rxn.species[self.butene], 1)

    def test_multiple_error_canceling_reactions(self):
        """
        Test that multiple error canceling reactions can be found
        """
        scheme = IsodesmicScheme(self.propene, [self.propane, self.butane, self.butene, self.pentane, self.pentene,
                                                self.hexane, self.hexene, self.benzene])

        reaction_list = scheme.multiple_error_canceling_reaction_search(n_reactions_max=20)
        self.assertEqual(len(reaction_list), 20)
        reaction_string = reaction_list.__repr__()
        # Consider both permutations of the products in the reaction string
        rxn_str1 = '<ErrorCancelingReaction 1*C=CC + 1*CCCC <=> 1*CCC + 1*C=CCC >'
        rxn_str2 = '<ErrorCancelingReaction 1*C=CC + 1*CCCC <=> 1*C=CCC + 1*CCC >'
        self.assertTrue(any(rxn_string in reaction_string for rxn_string in [rxn_str1, rxn_str2]))

        if self.pyo is not None:
            # pyomo is slower, so don't test as many
            reaction_list = scheme.multiple_error_canceling_reaction_search(n_reactions_max=5, milp_software=['pyomo'])
            self.assertEqual(len(reaction_list), 5)
            reaction_string = reaction_list.__repr__()
            self.assertTrue(any(rxn_string in reaction_string for rxn_string in [rxn_str1, rxn_str2]))

    def test_calculate_target_enthalpy(self):
        """
        Test that ErrorCancelingScheme is able to calculate thermochemistry for the target species
        """
        scheme = IsodesmicScheme(self.propene, [self.propane, self.butane, self.butene, self.pentane, self.pentene,
                                                self.hexane, self.hexene, self.benzene])

        target_thermo, rxn_list = scheme.calculate_target_enthalpy(n_reactions_max=3, milp_software=['lpsolve'])
        self.assertEqual(target_thermo.value_si, 115000.0)
        self.assertIsInstance(rxn_list[0], ErrorCancelingReaction)

        if self.pyo is not None:
            target_thermo, _ = scheme.calculate_target_enthalpy(n_reactions_max=3, milp_software=['pyomo'])
            self.assertEqual(target_thermo.value_si, 115000.0)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
