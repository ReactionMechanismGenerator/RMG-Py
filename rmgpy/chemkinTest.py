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
from unittest import mock

import rmgpy
from rmgpy.chemkin import get_species_identifier, load_chemkin_file, load_transport_file, mark_duplicate_reactions, \
    read_kinetics_entry, read_reaction_comments, read_thermo_entry, save_chemkin_file, save_chemkin_surface_file, \
    save_species_dictionary, save_transport_file, write_thermo_entry
from rmgpy.chemkin import _remove_line_breaks, _process_duplicate_reactions
from rmgpy.data.kinetics import LibraryReaction
from rmgpy.exceptions import ChemkinError
from rmgpy.kinetics.arrhenius import Arrhenius, MultiArrhenius
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.quantity import Quantity
from rmgpy.kinetics.surface import SurfaceArrhenius, StickingCoefficient


###################################################

class ChemkinTest(unittest.TestCase):
    @mock.patch('rmgpy.chemkin.logging')
    def test_read_thermo_entry_bad_element_count(self, mock_logging):
        """
        Test that invalid element count logs the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   XC   X          L   100.000  5000.000  827.28      1
                2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
                -1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
                -8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
                """
        with self.assertRaises(ValueError):
            read_thermo_entry(entry)

        mock_logging.info.assert_called_with(
            "Trouble reading line 'C2H6                    H   XC   X          L   100.000  5000.000  827.28      1' element segment 'H   X'")

    @mock.patch('rmgpy.chemkin.logging')
    def test_read_thermo_entry_not_gas_phase(self, mock_logging):
        """
        Test that non gas phase data logs the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   6C   2          L   100.000  5000.000  827.28      1
                2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
                -1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
                -8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
                """
        species, thermo, formula = read_thermo_entry(entry)

        mock_logging.warning.assert_called_with("Was expecting gas phase thermo data for C2H6. Skipping thermo data.")
        self.assertEqual(species, 'C2H6')
        self.assertIsNone(formula)
        self.assertIsNone(thermo)

    @mock.patch('rmgpy.chemkin.logging')
    def test_read_thermo_entry_not_float(self, mock_logging):
        """
        Test that non-float parameters log the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   6C   2          G   100.000  5000.000  827.28      1
 X.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
-1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
-8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
"""
        species, thermo, formula = read_thermo_entry(entry)

        mock_logging.warning.assert_called_with("could not convert string to float: 'X.44813916E+00'")
        self.assertEqual(species, 'C2H6')
        self.assertIsNone(formula)
        self.assertIsNone(thermo)

    def test_read_thermo_entry_no_temperature_range(self):
        """Test that missing temperature range can be handled for thermo entry."""
        entry = """C2H6                    H   6C   2          G                                  1
 2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
-1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
-8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
"""
        species, thermo, formula = read_thermo_entry(entry, Tmin=100.0, Tint=827.28, Tmax=5000.0)

        self.assertEqual(species, 'C2H6')
        self.assertEqual(formula, {'H': 6, 'C': 2})
        self.assertTrue(isinstance(thermo, NASA))

    def test_read_and_write_and_read_template_reaction_family_for_minimal_example(self):
        """
        This example tests if family and templates info can be correctly
        parsed from comments like '!Template reaction: R_Recombination'.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkin_path = os.path.join(folder, 'minimal', 'chem.inp')
        dictionary_path = os.path.join(folder, 'minimal', 'species_dictionary.txt')

        # read original chemkin file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        # ensure correct reading
        reaction1 = reactions[0]
        self.assertEqual(reaction1.family, "R_Recombination")
        self.assertEqual(frozenset('C_methyl;C_methyl'.split(';')), frozenset(reaction1.template))
        reaction2 = reactions[1]
        self.assertEqual(reaction2.family, "H_Abstraction")
        self.assertEqual(frozenset('C/H3/Cs\H3;C_methyl'.split(';')), frozenset(reaction2.template))
        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, 'minimal', 'chem_new.inp')
        dictionary_save_path = os.path.join(folder, 'minimal', 'species_dictionary_new.txt')

        save_chemkin_file(chemkin_save_path, species, reactions, verbose=True, check_for_duplicates=True)
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        self.assertTrue(os.path.isfile(chemkin_save_path))
        self.assertTrue(os.path.isfile(dictionary_save_path))

        # read newly written chemkin file to make sure the entire cycle works
        _, reactions2 = load_chemkin_file(chemkin_save_path, dictionary_save_path)

        reaction1_new = reactions2[0]
        self.assertEqual(reaction1_new.family, reaction1_new.family)
        self.assertEqual(reaction1_new.template, reaction1_new.template)
        self.assertEqual(reaction1_new.degeneracy, reaction1_new.degeneracy)

        reaction2_new = reactions2[1]
        self.assertEqual(reaction2_new.family, reaction2_new.family)
        self.assertEqual(reaction2_new.template, reaction2_new.template)
        self.assertEqual(reaction2_new.degeneracy, reaction2_new.degeneracy)

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)

    def test_read_and_write_template_reaction_family_for_pdd_example(self):
        """
        This example is mainly to ensure comments like
        '! Kinetics were estimated in this direction instead
        of the reverse because:' or '! This direction matched
        an entry in H_Abstraction, the other was just an estimate.'
        won't interfere reaction family info retrival.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkin_path = os.path.join(folder, 'pdd', 'chem.inp')
        dictionary_path = os.path.join(folder, 'pdd', 'species_dictionary.txt')

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        reaction1 = reactions[0]
        self.assertEqual(reaction1.family, "H_Abstraction")

        reaction2 = reactions[1]
        self.assertEqual(reaction2.family, "H_Abstraction")

        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, 'minimal', 'chem_new.inp')
        dictionary_save_path = os.path.join(folder, 'minimal', 'species_dictionary_new.txt')

        save_chemkin_file(chemkin_save_path, species, reactions, verbose=False, check_for_duplicates=False)
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        self.assertTrue(os.path.isfile(chemkin_save_path))
        self.assertTrue(os.path.isfile(dictionary_save_path))

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)

    def test_transport_data_read_and_write(self):
        """
        Test that we can write to chemkin and recreate the same transport object
        """
        Ar = Species(label="Ar",
                     transport_data=TransportData(shapeIndex=0, epsilon=(1134.93, 'J/mol'), sigma=(3.33, 'angstrom'),
                                                  dipoleMoment=(0, 'De'), polarizability=(0, 'angstrom^3'),
                                                  rotrelaxcollnum=0.0, comment="""GRI-Mech"""))

        Ar_write = Species(label="Ar")
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data')

        temp_transport_path = os.path.join(folder, 'tran_temp.dat')

        save_transport_file(temp_transport_path, [Ar])
        species_dict = {'Ar': Ar_write}
        load_transport_file(temp_transport_path, species_dict)
        self.assertEqual(repr(Ar), repr(Ar_write))

        os.remove(temp_transport_path)

    def test_use_chemkin_names(self):
        """
        Test that the official chemkin names are used as labels for the created Species objects.
        """

        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkin_path = os.path.join(folder, 'minimal', 'chem.inp')
        dictionary_path = os.path.join(folder, 'minimal', 'species_dictionary.txt')

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, use_chemkin_names=True)

        expected = [
            'Ar',
            'He',
            'Ne',
            'N2',
            'ethane',
            'CH3',
            'C2H5',
            'C'
        ]

        for spc, label in zip(species, expected):
            self.assertEqual(spc.label, label)

    def test_reactant_n2_is_reactive_and_gets_right_species_identifier(self):
        """
        Test that after loading chemkin files, species such as N2, which is in the default
        inert list of RMG, should be treated as reactive species and given right species
        Identifier when it's reacting in reactions.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkin_path = os.path.join(folder, 'NC', 'chem.inp')
        dictionary_path = os.path.join(folder, 'NC', 'species_dictionary.txt')

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, use_chemkin_names=True)

        for n2 in species:
            if n2.label == 'N2':
                break
        self.assertTrue(n2.reactive)

        self.assertEqual(get_species_identifier(n2), 'N2(35)')

    def test_read_specific_collider(self):
        """
        Test that a Chemkin reaction with a specific species as a third body collider can be properly read
        even if the species name contains parenthesis
        """
        entry = """O2(4)+H(5)(+N2(5))<=>HO2(10)(+N2(5))                          4.651e+12 0.440     0.000"""
        species_dict = {}
        s1 = Species().from_adjacency_list("""O2(4)
                                         multiplicity 3
                                         1 O u1 p2 c0 {2,S}
                                         2 O u1 p2 c0 {1,S}""")
        s2 = Species().from_adjacency_list("""H(5)
                                         multiplicity 2
                                         1 H u1 p0 c0""")
        s3 = Species().from_adjacency_list("""N2(5)
                                         1 N u0 p1 c0 {2,T}
                                         2 N u0 p1 c0 {1,T}""")
        s4 = Species().from_adjacency_list("""HO2(10)
                                         multiplicity 2
                                         1 O u0 p2 c0 {2,S} {3,S}
                                         2 O u1 p2 c0 {1,S}
                                         3 H u0 p0 c0 {1,S}""")
        species_dict['O2(4)'] = s1
        species_dict['H(5)'] = s2
        species_dict['N2(5)'] = s3
        species_dict['HO2(10)'] = s4
        A_units = ['', 's^-1', 'cm^3/(mol*s)', 'cm^6/(mol^2*s)', 'cm^9/(mol^3*s)']
        A_units_surf = []
        E_units = 'kcal/mol'
        reaction = read_kinetics_entry(entry, species_dict, A_units, A_units_surf, E_units)

        self.assertEqual(reaction.specific_collider.label, 'N2(5)')

    def test_process_duplicate_reactions(self):
        """
        Test that duplicate reactions are handled correctly when
        loading a Chemkin file.
        """
        s1 = Species().from_smiles('CC')
        s2 = Species().from_smiles('[CH3]')
        s3 = Species().from_smiles('[OH]')
        s4 = Species().from_smiles('C[CH2]')
        s5 = Species().from_smiles('O')
        r1 = Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius())
        r2 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Arrhenius())
        r3 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Arrhenius())
        r4 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=False, kinetics=Arrhenius())
        r5 = LibraryReaction(reactants=[s1, s3], products=[s4, s5], duplicate=True,
                             kinetics=Arrhenius(), library='lib1')
        r6 = LibraryReaction(reactants=[s1, s3], products=[s4, s5], duplicate=True,
                             kinetics=Arrhenius(), library='lib2')
        r7 = LibraryReaction(reactants=[s1, s3], products=[s4, s5], duplicate=True,
                             kinetics=Chebyshev(), library='lib1')
        r8 = LibraryReaction(reactants=[s1, s3], products=[s4, s5], duplicate=True,
                             kinetics=Arrhenius(), library='lib1')
        r9 = LibraryReaction(reactants=[s1, s3], products=[s4, s5], duplicate=False,
                             kinetics=MultiArrhenius(arrhenius=[Arrhenius(), Arrhenius()]), library='lib1')
        reaction_list_with_duplicate = [r1, r2, r3]
        reaction_list_with_duplicate2 = [r1, r2, r3]
        reaction_list_unmarked_duplicate = [r1, r2, r4]
        reaction_list_unequal_libraries = [r1, r5, r6]
        reaction_list_mixed_kinetics = [r1, r5, r7]
        reaction_list_mergeable = [r1, r5, r8]
        reaction_list_merged = [r1, r9]

        # Test that duplicates are not removed for non-library reactions
        _process_duplicate_reactions(reaction_list_with_duplicate)
        self.assertEqual(reaction_list_with_duplicate, reaction_list_with_duplicate2)

        # Test that unmarked duplicate reactions are detected if both
        # reactions are p-dep or p-indep
        self.assertRaisesRegexp(ChemkinError,
                                'Encountered unmarked duplicate reaction',
                                _process_duplicate_reactions,
                                reaction_list_unmarked_duplicate)

        # Test that unequal libraries are recognized
        self.assertRaisesRegexp(ChemkinError,
                                'from different libraries',
                                _process_duplicate_reactions,
                                reaction_list_unequal_libraries)

        # Test that an error is raised for reactions with kinetics
        # that cannot be merged
        self.assertRaisesRegexp(ChemkinError,
                                'Mixed kinetics for duplicate reaction',
                                _process_duplicate_reactions,
                                reaction_list_mixed_kinetics)

        # Test that duplicate library reactions are merged successfully
        _process_duplicate_reactions(reaction_list_mergeable)
        self.assertEqual(len(reaction_list_mergeable), len(reaction_list_merged))
        self.assertEqual(reaction_list_mergeable[0], reaction_list_merged[0])
        rtest = reaction_list_mergeable[1]
        rtrue = reaction_list_merged[1]
        self.assertEqual(rtest.reactants, rtrue.reactants)
        self.assertEqual(rtest.products, rtrue.products)
        self.assertEqual(rtest.duplicate, rtrue.duplicate)
        self.assertEqual(rtest.library, rtrue.library)
        self.assertTrue(isinstance(rtest.kinetics, MultiArrhenius))
        self.assertTrue(all(isinstance(k, Arrhenius) for k in rtest.kinetics.arrhenius))

    def test_mark_duplicate_reactions(self):
        """Test that we can properly mark duplicate reactions for Chemkin."""
        s1 = Species().from_smiles('CC')
        s2 = Species().from_smiles('[CH3]')
        s3 = Species().from_smiles('[OH]')
        s4 = Species().from_smiles('C[CH2]')
        s5 = Species().from_smiles('O')
        s6 = Species().from_smiles('[H]')

        # Try initializing with duplicate=False
        reaction_list = [
            Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius()),
            Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius()),
            Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=False, kinetics=Arrhenius()),
            Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=False, kinetics=Chebyshev()),
            Reaction(reactants=[s1], products=[s4, s6], duplicate=False, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s1], products=[s4, s6], duplicate=False, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s5], products=[s3, s6], duplicate=False, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s3, s6], products=[s5], duplicate=False, kinetics=Arrhenius(), reversible=False),
        ]

        expected_flags = [True, True, False, False, True, True, False, False]

        mark_duplicate_reactions(reaction_list)
        duplicate_flags = [rxn.duplicate for rxn in reaction_list]

        self.assertEqual(duplicate_flags, expected_flags)

        # Try initializing with duplicate=True
        reaction_list = [
            Reaction(reactants=[s1], products=[s2, s2], duplicate=True, kinetics=Arrhenius()),
            Reaction(reactants=[s1], products=[s2, s2], duplicate=True, kinetics=Arrhenius()),
            Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Arrhenius()),
            Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Chebyshev()),
            Reaction(reactants=[s1], products=[s4, s6], duplicate=True, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s1], products=[s4, s6], duplicate=True, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s5], products=[s3, s6], duplicate=True, kinetics=Arrhenius(), reversible=False),
            Reaction(reactants=[s3, s6], products=[s5], duplicate=True, kinetics=Arrhenius(), reversible=False),
        ]

        mark_duplicate_reactions(reaction_list)
        duplicate_flags = [rxn.duplicate for rxn in reaction_list]

        self.assertEqual(duplicate_flags, expected_flags)

    def test_read_write_coverage_dependence(self):
        """Test that we can properly read and write coverage dependent parameters"""

        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        s_x_entry = """X                       X   1               G   100.000  5000.000 1554.80      1
 1.60299900E-01-2.52235409E-04 1.14181275E-07-1.21471653E-11 3.85790025E-16    2
-7.08100885E+01-9.09527530E-01 7.10139498E-03-4.25619522E-05 8.98533016E-08    3
-7.80193649E-11 2.32465471E-14-8.76101712E-01-3.11211229E-02                   4"""
        s_hx_entry = """H*                      H   1X   1          G   100.000  5000.000  952.91      1
 2.80339655E+00-5.41047017E-04 4.99507978E-07-7.54963647E-11 3.06772366E-15    2
-2.34636021E+03-1.59436787E+01-3.80965452E-01 5.47228709E-03 2.60912778E-06    3
-9.64961980E-09 4.63946753E-12-1.40561079E+03 1.01725550E+00                   4"""
        s_h2_entry = """H2                      H   2               G   100.000  5000.000 1959.08      1
 2.78816619E+00 5.87640475E-04 1.59010635E-07-5.52739465E-11 4.34311304E-15    2
-5.96144481E+02 1.12730527E-01 3.43536411E+00 2.12710383E-04-2.78625110E-07    3
 3.40267219E-10-7.76032129E-14-1.03135984E+03-3.90841731E+00                   4"""
        s_ohx_entry = """OH*                     H   1O   1X   1     G   100.000  5000.000  914.54      1
 2.43541504E+00 4.64599933E-03-2.39987608E-06 4.26351871E-10-2.60607840E-14    2
-2.17972457E+04-1.03879495E+01-1.29522792E+00 3.36487597E-02-7.07760603E-05    3
 6.54375105E-08-2.19437885E-11-2.16453879E+04 4.37657050E+00                   4"""
        s_ox_entry = """O*                      O   1X   1          G   100.000  5000.000  888.26      1
 1.89893631E+00 2.03295404E-03-1.19976562E-06 2.32680628E-10-1.53508256E-14    2
-2.21111565E+04-9.64104147E+00-7.59013115E-01 1.89868504E-02-3.82473769E-05    3
 3.43558429E-08-1.13974388E-11-2.18356104E+04 1.76017413E+00                   4"""
 
        s_h2 = Species().from_smiles("[H][H]")
        s_x = Species().from_adjacency_list("1 X u0 p0")
        s_x.label = 'X'
        s_hx = Species().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")
        s_hx.label = 'HX'
        s_ohx = Species().from_smiles("O[*]")
        s_ohx.label = 'OHX'
        s_ox = Species().from_smiles("O=[*]")
        s_ox.label = 'OX'

        species, thermo, formula = read_thermo_entry(s_x_entry)
        s_x.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_hx_entry)
        s_hx.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_h2_entry)
        s_h2.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_ohx_entry)
        s_ohx.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_ox_entry)
        s_ox.thermo = thermo


        species = [s_h2, s_x, s_hx, s_ohx, s_ox]

        # test a coverage dependent arrhenius and a coverage dependent sticking coefficient
        rxn1_entry = """H2+X+X=H*+H*                                        4.000e-02 0.000     7.098
    STICK
    COV / H*                                        0         0         1.099     /"""

        rxn2_entry = """X+OH*=O*+H*                                         7.390000e+19 0.000     18.475   
    COV / O*                                        0         0         -17.500   /"""

        species_dict  = {}
        species_dict['H2'] = s_h2
        species_dict['X'] = s_x
        species_dict['H*'] = s_hx
        species_dict['OH*'] = s_ohx
        species_dict['O*'] = s_ox

        A_units = ['', 's^-1', 'cm^3/(mol*s)',
                   'cm^6/(mol^2*s)', 'cm^9/(mol^3*s)']
        A_units_surf = ['', 's^-1', 'cm^2/(mol*s)',
                        'cm^4/(mol^2*s)', 'cm^6/(mol^3*s)']

        E_units = 'kcal/mol'
        self.rxn1 = read_kinetics_entry(
            rxn1_entry, species_dict, A_units, A_units_surf, E_units)

        self.rxn2 = read_kinetics_entry(
            rxn2_entry, species_dict, A_units, A_units_surf, E_units)

        reactions = [self.rxn1, self.rxn2]

        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, 'surface', 'chem-surface_new.inp')
        dictionary_save_path = os.path.join(folder, 'surface', 'species_dictionary_new.txt')
        save_chemkin_file(chemkin_save_path, species, reactions, verbose=True, check_for_duplicates=True)
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        # load chemkin file
        species_load, reactions_load = load_chemkin_file(chemkin_save_path, dictionary_save_path)

        # compare only labels because the objects are different
        species_label = []
        species_load_label = []
        for s in range(len(species)):
            species_label.append(species[s].label)
            species_load_label.append(species_load[s].label)

        # check the written chemkin file matches what is expected
        self.assertTrue(all(i in species_load_label for i in species_label))

        for r in range(len(reactions)):
            for s in range(len(reactions[r].products)):
                self.assertEqual(reactions[r].products[s].label,reactions_load[r].products[s].label)
            for s in range(len(reactions[r].reactants)):
                self.assertEqual(reactions[r].reactants[s].label, reactions_load[r].reactants[s].label)
            self.assertIsNotNone(reactions[r].kinetics.coverage_dependence)
            self.assertIsNotNone(reactions_load[r].kinetics.coverage_dependence)

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)

class TestThermoReadWrite(unittest.TestCase):

    def setUp(self):
        """This method is run once before each test."""
        coeffs_low = [4.03055, -0.00214171, 4.90611e-05, -5.99027e-08, 2.38945e-11, -11257.6, 3.5613]
        coeffs_high = [-0.307954, 0.0245269, -1.2413e-05, 3.07724e-09, -3.01467e-13, -10693, 22.628]
        Tmin = 300.
        Tmax = 3000.
        Tint = 650.73
        E0 = -782292.  # J/mol.
        comment = "C2H6"
        self.nasa = NASA(
            polynomials=[
                NASAPolynomial(coeffs=coeffs_low, Tmin=(Tmin, "K"), Tmax=(Tint, "K")),
                NASAPolynomial(coeffs=coeffs_high, Tmin=(Tint, "K"), Tmax=(Tmax, "K")),
            ],
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            E0=(E0, "J/mol"),
            comment=comment,
        )

        # Chemkin entries for testing - note that the values are all the same
        self.entry1 = """C2H6                    C   2H   6          G   300.000  3000.000  650.73      1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

        self.entry2 = """CH3NO2X                                     G   300.000  3000.000  650.73      1&
C 1 H 3 N 1 O 2 X 1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

        self.entry3 = """CH3NO2SX                                    G   300.000  3000.000  650.73      1&
C 1 H 3 N 1 O 2 S 1 X 1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

    def test_write_thermo_block(self):
        """Test that we can write a normal thermo block"""
        species = Species(smiles='CC')
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        self.assertEqual(result, self.entry1)

    def test_read_thermo_block(self):
        """Test that we can read a normal thermo block"""
        species, thermo, formula = read_thermo_entry(self.entry1)

        self.assertEqual(species, 'C2H6')
        self.assertEqual(formula, {'H': 6, 'C': 2})
        self.assertTrue(self.nasa.is_identical_to(thermo))

    def test_write_thermo_block_5_elem(self):
        """Test that we can write a thermo block for a species with 5 elements"""
        species = Species().from_adjacency_list("""
1 O u0 p3 c-1 {3,S}
2 O u0 p2 c0 {3,D}
3 N u0 p0 c+1 {1,S} {2,D} {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 X u0 p0 c0
""")
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        self.assertEqual(result, self.entry2)

    def test_read_thermo_block_5_elem(self):
        """Test that we can read a thermo block with 5 elements"""
        species, thermo, formula = read_thermo_entry(self.entry2)

        self.assertEqual(species, 'CH3NO2X')
        self.assertEqual(formula, {'X': 1, 'C': 1, 'O': 2, 'H': 3, 'N': 1})
        self.assertTrue(self.nasa.is_identical_to(thermo))

    def test_write_thermo_block_6_elem(self):
        """Test that we can write a thermo block for a species with 6 elements"""
        species = Species().from_adjacency_list("""
1 O u0 p3 c-1 {2,S}
2 N u0 p0 c+1 {1,S} {3,D} {4,S}
3 O u0 p2 c0 {2,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 S u0 p2 c0 {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 X u0 p0 c0
""")
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        self.assertEqual(result, self.entry3)

    def test_read_thermo_block_6_elem(self):
        """Test that we can read a thermo block with 6 elements"""
        species, thermo, formula = read_thermo_entry(self.entry3)

        self.assertEqual(species, 'CH3NO2SX')
        self.assertEqual(formula, {'X': 1, 'C': 1, 'O': 2, 'H': 3, 'N': 1, 'S': 1})
        self.assertTrue(self.nasa.is_identical_to(thermo))
    
    def test_write_bidentate_species(self):
        """Test that species with 2 or more surface sites get proper formatting"""

        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')
        chemkin_path = os.path.join(folder, 'surface', 'chem-surface.inp')
        dictionary_path = os.path.join(folder, 'surface', 'species_dictionary.txt')
        chemkin_save_path = os.path.join(folder, 'surface', 'chem-surface-test.inp')
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        surface_atom_count = species[3].molecule[0].get_num_atoms('X')
        self.assertEqual(surface_atom_count, 3)
        save_chemkin_surface_file(chemkin_save_path, species, reactions, verbose=False, check_for_duplicates=False)
        
        bidentate_test = "    CH2OX2(52)/2/             \n"
        tridentate_test = "    CHOX3(61)/3/             \n" 
        with open(chemkin_save_path, "r") as f:
            for i, line in enumerate(f):
                if i == 3:
                    bidentate_read = line
                if i == 4:
                    tridentate_read = line

        self.assertEqual(bidentate_test.strip(), bidentate_read.strip())
        self.assertEqual(tridentate_test.strip(), tridentate_read.strip())
        
        os.remove(chemkin_save_path)

class TestReadReactionComments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r = Species().from_smiles('[CH3]')
        r.label = '[CH3]'
        p = Species().from_smiles('CC')
        p.label = 'CC'

        cls.reaction = Reaction(reactants=[r, r],
                                products=[p],
                                kinetics=Arrhenius(A=(8.26e+17, 'cm^3/(mol*s)'),
                                                   n=-1.4,
                                                   Ea=(1, 'kcal/mol'), T0=(1, 'K'))
                                )
        cls.comments_list = ["""
                              Reaction index: Chemkin #1; RMG #1
                              Template reaction: R_Recombination
                              Exact match found for rate rule (C_methyl;C_methyl)
                              Multiplied by reaction path degeneracy 0.5
                              """,
                             """
                             Reaction index: Chemkin #2; RMG #4
                             Template reaction: H_Abstraction
                             Estimated using template (C/H3/Cs;C_methyl) for rate rule (C/H3/Cs\H3;C_methyl)
                             Multiplied by reaction path degeneracy 6
                             """,
                             """
                              Reaction index: Chemkin #13; RMG #8
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using an average for rate rule [C/H3/Cs\H3;C_rad/H2/Cs]
                              Multiplied by reaction path degeneracy 6.0
                              """,
                             """
                              Reaction index: Chemkin #17; RMG #31
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using average of templates [C/H3/Cs;H_rad] + [C/H3/Cs\H3;Y_rad] for rate rule [C/H3/Cs\H3;H_rad]
                              Multiplied by reaction path degeneracy 6.0
                              """,
                             """
                              Reaction index: Chemkin #69; RMG #171
                              Template reaction: intra_H_migration
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] for rate rule
                              [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
                              Multiplied by reaction path degeneracy 3.0
                              """,
                             """
                              Reaction index: Chemkin #3; RMG #243
                              Template reaction: Disproportionation
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Average of [Average of [O2b;O_Csrad] + Average of [O_atom_triplet;O_Csrad + CH2_triplet;O_Csrad] + Average of [Average of [Ct_rad/Ct;O_Csrad from
                              training reaction 0] + Average of [O_pri_rad;O_Csrad + Average of [O_rad/NonDeC;O_Csrad + O_rad/NonDeO;O_Csrad]] + Average of [Cd_pri_rad;O_Csrad] +
                              Average of [CO_pri_rad;O_Csrad] + Average of [C_methyl;O_Csrad + Average of [C_rad/H2/Cs;O_Csrad + C_rad/H2/Cd;O_Csrad + C_rad/H2/O;O_Csrad] + Average
                              of [C_rad/H/NonDeC;O_Csrad] + Average of [Average of [C_rad/Cs3;O_Csrad]]] + H_rad;O_Csrad]]
                              Estimated using template [Y_rad_birad_trirad_quadrad;O_Csrad] for rate rule [CH_quartet;O_Csrad]
                              """,
                             """
                              Reaction index: Chemkin #4; RMG #303
                              Template reaction: Disproportionation
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Matched reaction 0 C2H + CH3O <=> C2H2 + CH2O in Disproportionation/training
                              """,
                             """
                              Reaction index: Chemkin #51; RMG #136
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using an average for rate rule [C/H3/Cd\H_Cd\H2;C_rad/H2/Cs]
                              Euclidian distance = 0
                              Multiplied by reaction path degeneracy 3.0
                              """,
                             """
                             Reaction index: Chemkin #32; RMG #27
                             Template reaction: R_Recombination
                             Flux pairs: [CH3], CC; [CH3], CC; 
                             Matched reaction 20 CH3 + CH3 <=> C2H6 in R_Recombination/training
                             This reaction matched rate rule [C_methyl;C_methyl]
                             """,
                             """
                             Reaction index: Chemkin #2; RMG #4
                             Template reaction: R_Recombination
                             Flux pairs: [CH3], CC; [CH3], CC; 
                             From training reaction 21 used for C_rad/H2/Cs;C_methyl
                             Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
                             Euclidian distance = 0
                             """]
        cls.template_list = [['C_methyl', 'C_methyl'],
                             ['C/H3/Cs\H3', 'C_methyl'],
                             ['C/H3/Cs\H3', 'C_rad/H2/Cs'],
                             ['C/H3/Cs\H3', 'H_rad'],
                             ['R3H_SS_Cs', 'O_rad_out', 'Cs_H_out_2H'],
                             ['CH_quartet', 'O_Csrad'],
                             None,
                             ['C/H3/Cd\H_Cd\H2', 'C_rad/H2/Cs'],
                             ['C_methyl', 'C_methyl'],
                             ['C_rad/H2/Cs', 'C_methyl']]
        cls.family_list = ['R_Recombination',
                           'H_Abstraction',
                           'H_Abstraction',
                           'H_Abstraction',
                           'intra_H_migration',
                           'Disproportionation',
                           'Disproportionation',
                           'H_Abstraction',
                           'R_Recombination',
                           'R_Recombination', ]
        cls.degeneracy_list = [0.5,
                               6,
                               6,
                               6,
                               3,
                               1,
                               1,
                               3,
                               1,
                               1]
        cls.expected_lines = [4, 4, 5, 5, 5, 5, 4, 6, 5, 6]

    def test_read_reaction_comments_template(self):
        """
        Test that the template is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = read_reaction_comments(self.reaction, comment)

            # only check template if meant to find one
            if self.template_list[index]:
                self.assertTrue(new_rxn.template,
                                'The template was not saved from the reaction comment {}'.format(comment))
                self.assertEqual(frozenset(new_rxn.template), frozenset(self.template_list[index]),
                                 'The reaction template does not match')
            else:
                self.assertFalse(new_rxn.template)

    def test_read_reaction_comments_family(self):
        """
        Test that the family is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = read_reaction_comments(self.reaction, comment)

            self.assertEqual(new_rxn.family, self.family_list[index], 'wrong reaction family stored')

    def test_read_reaction_comments_degeneracy(self):
        """
        Test that the degeneracy is picked up from reading reaction comments.

        Also checks that reaction rate was not modified in the process.
        """
        for index, comment in enumerate(self.comments_list):
            # Clear any leftover kinetics comments
            self.reaction.kinetics.comment = ''
            previous_rate = self.reaction.kinetics.A.value_si
            new_rxn = read_reaction_comments(self.reaction, comment)
            new_rate = new_rxn.kinetics.A.value_si

            self.assertEqual(new_rxn.degeneracy, self.degeneracy_list[index], 'wrong degeneracy was stored')
            self.assertEqual(previous_rate, new_rate)

            # Check that the comment only appears once in the kinetics comment
            if new_rxn.degeneracy != 1:
                self.assertEqual(new_rxn.kinetics.comment.count(
                    'Multiplied by reaction path degeneracy {}'.format(new_rxn.degeneracy)), 1,
                    'Reaction degeneracy comment duplicated while reading Chemkin comments')
            else:
                self.assertTrue('Multiplied by reaction path degeneracy' not in new_rxn.kinetics.comment)

    def test_remove_line_breaks(self):
        """
        tests that _remove_line_breaks functions properly
        """
        for index, comment in enumerate(self.comments_list):
            new_comment = _remove_line_breaks(comment)
            new_comment_lines = len(new_comment.strip().splitlines())
            self.assertEqual(new_comment_lines, self.expected_lines[index],
                             'Found {} more lines than expected for comment \n\n""{}""\n\n which converted to \n\n""{}""'.format(
                                 new_comment_lines - self.expected_lines[index], comment.strip(), new_comment.strip()))
