################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest
import mock
import os
from chemkin import *
from chemkin import _removeLineBreaks
import rmgpy
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics.arrhenius import Arrhenius


###################################################

class ChemkinTest(unittest.TestCase):
    @mock.patch('rmgpy.chemkin.logging')
    def test_readThermoEntry_BadElementCount(self, mock_logging):
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
            readThermoEntry(entry)

        mock_logging.info.assert_called_with("Trouble reading line 'C2H6                    H   XC   X          L   100.000  5000.000  827.28      1' element segment 'H   X'")

    @mock.patch('rmgpy.chemkin.logging')
    def test_readThermoEntry_NotGasPhase(self, mock_logging):
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
        species, thermo, formula = readThermoEntry(entry)

        mock_logging.warning.assert_called_with("Was expecting gas phase thermo data for C2H6. Skipping thermo data.")
        self.assertEqual(species, 'C2H6')
        self.assertIsNone(formula)
        self.assertIsNone(thermo)

    @mock.patch('rmgpy.chemkin.logging')
    def test_readThermoEntry_NotFloat(self, mock_logging):
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
        species, thermo, formula = readThermoEntry(entry)

        mock_logging.warning.assert_called_with("could not convert string to float: X.44813916E+00")
        self.assertEqual(species, 'C2H6')
        self.assertIsNone(formula)
        self.assertIsNone(thermo)

    def test_readThermoEntry_NoTRange(self):
        """Test that missing temperature range can be handled for thermo entry."""
        entry = """C2H6                    H   6C   2          G                                  1
 2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
-1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
-8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
"""
        species, thermo, formula = readThermoEntry(entry, Tmin=100.0, Tint=827.28, Tmax=5000.0)

        self.assertEqual(species, 'C2H6')
        self.assertEqual(formula, {'H': 6, 'C': 2})
        self.assertTrue(isinstance(thermo, NASA))

    def testReadAndWriteAndReadTemplateReactionFamilyForMinimalExample(self):
        """
        This example tests if family and templates info can be correctly
        parsed from comments like '!Template reaction: R_Recombination'.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkinPath = os.path.join(folder, 'minimal', 'chem.inp')
        dictionaryPath = os.path.join(folder, 'minimal', 'species_dictionary.txt')

        # read original chemkin file
        species, reactions = loadChemkinFile(chemkinPath, dictionaryPath)

        # ensure correct reading
        reaction1 = reactions[0]
        self.assertEqual(reaction1.family, "R_Recombination")
        self.assertEqual(frozenset('C_methyl;C_methyl'.split(';')), frozenset(reaction1.template))
        reaction2 = reactions[1]
        self.assertEqual(reaction2.family, "H_Abstraction")
        self.assertEqual(frozenset('C/H3/Cs\H3;C_methyl'.split(';')), frozenset(reaction2.template))
        # saveChemkinFile
        chemkinSavePath = os.path.join(folder, 'minimal', 'chem_new.inp')
        dictionarySavePath = os.path.join(folder, 'minimal', 'species_dictionary_new.txt')

        saveChemkinFile(chemkinSavePath, species, reactions, verbose=True, checkForDuplicates=True)
        saveSpeciesDictionary(dictionarySavePath, species, oldStyle=False)

        self.assertTrue(os.path.isfile(chemkinSavePath))
        self.assertTrue(os.path.isfile(dictionarySavePath))

        # read newly written chemkin file to make sure the entire cycle works
        _, reactions2 =loadChemkinFile(chemkinSavePath, dictionarySavePath)

        reaction1_new = reactions2[0]
        self.assertEqual(reaction1_new.family, reaction1_new.family)
        self.assertEqual(reaction1_new.template, reaction1_new.template)
        self.assertEqual(reaction1_new.degeneracy, reaction1_new.degeneracy)

        reaction2_new = reactions2[1]
        self.assertEqual(reaction2_new.family, reaction2_new.family)
        self.assertEqual(reaction2_new.template, reaction2_new.template)
        self.assertEqual(reaction2_new.degeneracy, reaction2_new.degeneracy)

        # clean up
        os.remove(chemkinSavePath)
        os.remove(dictionarySavePath)

    def testReadAndWriteTemplateReactionFamilyForPDDExample(self):
        """
        This example is mainly to ensure comments like
        '! Kinetics were estimated in this direction instead
        of the reverse because:' or '! This direction matched
        an entry in H_Abstraction, the other was just an estimate.'
        won't interfere reaction family info retrival.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkinPath = os.path.join(folder, 'pdd', 'chem.inp')
        dictionaryPath = os.path.join(folder, 'pdd', 'species_dictionary.txt')

        # loadChemkinFile
        species, reactions = loadChemkinFile(chemkinPath, dictionaryPath)

        reaction1 = reactions[0]
        self.assertEqual(reaction1.family, "H_Abstraction")

        reaction2 = reactions[1]
        self.assertEqual(reaction2.family, "H_Abstraction")

        # saveChemkinFile
        chemkinSavePath = os.path.join(folder, 'minimal', 'chem_new.inp')
        dictionarySavePath = os.path.join(folder, 'minimal', 'species_dictionary_new.txt')

        saveChemkinFile(chemkinSavePath, species, reactions, verbose=False, checkForDuplicates=False)
        saveSpeciesDictionary(dictionarySavePath, species, oldStyle=False)

        self.assertTrue(os.path.isfile(chemkinSavePath))
        self.assertTrue(os.path.isfile(dictionarySavePath))

        # clean up
        os.remove(chemkinSavePath)
        os.remove(dictionarySavePath)

    def testTransportDataReadAndWrite(self):
        """
        Test that we can write to chemkin and recreate the same transport object
        """
        from rmgpy.species import Species
        from rmgpy.molecule import Molecule
        from rmgpy.transport import TransportData

        Ar = Species(label="Ar",
                     transportData=TransportData(shapeIndex=0, epsilon=(1134.93, 'J/mol'), sigma=(3.33, 'angstrom'),
                                                 dipoleMoment=(0, 'De'), polarizability=(0, 'angstrom^3'),
                                                 rotrelaxcollnum=0.0, comment="""GRI-Mech"""))

        Ar_write = Species(label="Ar")
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data')

        tempTransportPath = os.path.join(folder, 'tran_temp.dat')

        saveTransportFile(tempTransportPath, [Ar])
        speciesDict = {'Ar': Ar_write}
        loadTransportFile(tempTransportPath, speciesDict)
        self.assertEqual(repr(Ar), repr(Ar_write))

        os.remove(tempTransportPath)

    def testUseChemkinNames(self):
        """
        Test that the official chemkin names are used as labels for the created Species objects.
        """

        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkinPath = os.path.join(folder, 'minimal', 'chem.inp')
        dictionaryPath = os.path.join(folder, 'minimal', 'species_dictionary.txt')

        # loadChemkinFile
        species, reactions = loadChemkinFile(chemkinPath, dictionaryPath, useChemkinNames=True)

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

    def testReactantN2IsReactiveAndGetsRightSpeciesIdentifier(self):
        """
        Test that after loading chemkin files, species such as N2, which is in the default
        inert list of RMG, should be treated as reactive species and given right species
        Identifier when it's reacting in reactions.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'test_data/chemkin/chemkin_py')

        chemkinPath = os.path.join(folder, 'NC', 'chem.inp')
        dictionaryPath = os.path.join(folder, 'NC', 'species_dictionary.txt')

        # loadChemkinFile
        species, reactions = loadChemkinFile(chemkinPath, dictionaryPath, useChemkinNames=True)

        for n2 in species:
            if n2.label == 'N2':
                break
        self.assertTrue(n2.reactive)

        self.assertEqual(getSpeciesIdentifier(n2), 'N2(35)')

    def testReadSpecificCollider(self):
        """
        Test that a Chemkin reaction with a specific species as a third body collider can be properly read
        even if the species name contains parenthesis
        """
        entry = """O2(4)+H(5)(+N2(5))<=>HO2(10)(+N2(5))                          4.651e+12 0.440     0.000"""
        speciesDict = {}
        s1 = Species().fromAdjacencyList("""O2(4)
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}""")
        s2 = Species().fromAdjacencyList("""H(5)
multiplicity 2
1 H u1 p0 c0""")
        s3 = Species().fromAdjacencyList("""N2(5)
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}""")
        s4 = Species().fromAdjacencyList("""HO2(10)
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}""")
        speciesDict['O2(4)'] = s1
        speciesDict['H(5)'] = s2
        speciesDict['N2(5)'] = s3
        speciesDict['HO2(10)'] = s4
        Aunits = ['','s^-1','cm^3/(mol*s)','cm^6/(mol^2*s)','cm^9/(mol^3*s)']
        Eunits = 'kcal/mol'
        reaction = readKineticsEntry(entry, speciesDict, Aunits, Eunits)

        self.assertEqual(reaction.specificCollider.label, 'N2(5)')

class TestReadReactionComments(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        r = Species().fromSMILES('[CH3]')
        r.label = '[CH3]'
        p = Species().fromSMILES('CC')
        p.label = 'CC'

        self.reaction = Reaction(reactants=[r,r],
                                 products=[p],
                                 kinetics = Arrhenius(A=(8.26e+17,'cm^3/(mol*s)'),
                                                      n=-1.4,
                                                      Ea=(1,'kcal/mol'), T0=(1,'K'))
                                )
        self.comments_list = ["""
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
From training reaction 21 for rate rule [C_rad/H2/Cs;C_methyl]
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
"""]
        self.template_list = [['C_methyl','C_methyl'],
                              ['C/H3/Cs\H3','C_methyl'],
                              ['C/H3/Cs\H3','C_rad/H2/Cs'],
                              ['C/H3/Cs\H3','H_rad'],
                              ['R3H_SS_Cs','O_rad_out','Cs_H_out_2H'],
                              ['CH_quartet','O_Csrad'],
                              None,
                              ['C/H3/Cd\H_Cd\H2','C_rad/H2/Cs'],
                              ['C_methyl','C_methyl'],
                              ['C_rad/H2/Cs','C_methyl']]
        self.family_list = ['R_Recombination',
                            'H_Abstraction',
                            'H_Abstraction',
                            'H_Abstraction',
                            'intra_H_migration',
                            'Disproportionation',
                            'Disproportionation',
                            'H_Abstraction',
                            'R_Recombination',
                            'R_Recombination',]
        self.degeneracy_list = [0.5,
                                6,
                                6,
                                6,
                                3,
                                1,
                                1,
                                3,
                                1,
                                1]
        self.expected_lines = [4,4,5,5,5,5,4,6,5,6]

    def testReadReactionCommentsTemplate(self):
        """
        Test that the template is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = readReactionComments(self.reaction, comment)

            # only check template if meant to find one
            if self.template_list[index]:
                self.assertTrue(new_rxn.template,'The template was not saved from the reaction comment {}'.format(comment))
                self.assertEqual(frozenset(new_rxn.template),frozenset(self.template_list[index]),'The reaction template does not match')
            else:
                self.assertFalse(new_rxn.template)

    def testReadReactionCommentsFamily(self):
        """
        Test that the family is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = readReactionComments(self.reaction, comment)

            self.assertEqual(new_rxn.family, self.family_list[index], 'wrong reaction family stored')

    def testReadReactionCommentsDegeneracy(self):
        """
        Test that the degeneracy is picked up from reading reaction comments.

        Also checks that reaction rate was not modified in the process.
        """
        for index, comment in enumerate(self.comments_list):
            previous_rate = self.reaction.kinetics.A.value_si
            new_rxn = readReactionComments(self.reaction, comment)
            new_rate = new_rxn.kinetics.A.value_si

            self.assertEqual(new_rxn.degeneracy, self.degeneracy_list[index], 'wrong degeneracy was stored')
            self.assertEqual(previous_rate, new_rate)

    def testRemoveLineBreaks(self):
        """
        tests that _removeLineBreaks functions properly
        """
        for index, comment in enumerate(self.comments_list):
            new_comment = _removeLineBreaks(comment)
            new_comment_lines = len(new_comment.strip().splitlines())
            self.assertEqual(new_comment_lines, self.expected_lines[index],
                             'Found {} more lines than expected for comment \n\n""{}""\n\n which converted to \n\n""{}""'.format(new_comment_lines - self.expected_lines[index],comment.strip(), new_comment.strip()))
