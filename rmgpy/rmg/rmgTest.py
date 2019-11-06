#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.rmg.react import react_species
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species
from rmgpy.util import parse_command_line_arguments


###################################################

class TestRMGWorkFlow(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        # set-up RMG object
        cls.rmg = RMG()
        cls.rmg.reaction_model = CoreEdgeReactionModel()

        # load kinetic database and forbidden structures
        cls.rmg.database = RMGDatabase()
        path = os.path.join(settings['test_data.directory'], 'testing_database')

        # kinetics family Disproportionation loading
        cls.rmg.database.load_kinetics(os.path.join(path, 'kinetics'),
                                       kinetics_families=['H_Abstraction', 'R_Addition_MultipleBond'],
                                       reaction_libraries=[])

        # load empty forbidden structures
        for family in cls.rmg.database.kinetics.families.values():
            family.forbidden = ForbiddenStructures()
        cls.rmg.database.forbidden_structures = ForbiddenStructures()

    @classmethod
    def tearDownClass(cls):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    @work_in_progress
    def test_deterministic_reaction_template_matching(self):
        """
        Test RMG work flow can match reaction template for kinetics estimation 
        deterministically. 

        In this test, a change of molecules order in a reacting species should 
        not change the reaction template matched.

        However, this is inherently impossible with the existing reaction
        generation algorithm. Currently, the first reaction will be the one
        that is kept if the reactions are identical. If different templates
        are a result of different transition states, all are kept.
        
        {O=C-[C]=C, [O]-C=C=C} -> H + C=C=C=O
        """

        # react
        spc = Species().from_smiles("O=C[C]=C")
        spc.generate_resonance_structures()
        new_reactions = react_species((spc,))

        # try to pick out the target reaction 
        mol_H = Molecule().from_smiles("[H]")
        mol_C3H2O = Molecule().from_smiles("C=C=C=O")

        target_rxns = find_target_rxns_containing(mol_H, mol_C3H2O, new_reactions)
        self.assertEqual(len(target_rxns), 2)

        # reverse the order of molecules in spc
        spc.molecule = list(reversed(spc.molecule))

        # react again
        new_reactions_reverse = []
        new_reactions_reverse.extend(react_species((spc,)))

        # try to pick out the target reaction 
        target_rxns_reverse = find_target_rxns_containing(mol_H, mol_C3H2O, new_reactions_reverse)
        self.assertEqual(len(target_rxns_reverse), 2)

        # whatever order of molecules in spc, the reaction template matched should be same
        self.assertEqual(target_rxns[0].template, target_rxns_reverse[0].template)

    def test_check_for_existing_species_for_bi_aromatics(self):
        """
        Test RMG check_for_existing_species can correctly check isomorphism for biaromatics.
        In this test, DPP is a species already stored in rmg species_dict, mol_test is a newly
        created molecule which has one kekulized benzene ring and one double_bond-single_bond
        benzene ring.
        """

        rmg_test = RMG()
        rmg_test.reaction_model = CoreEdgeReactionModel()
        DPP = Species().from_smiles('C1=CC=C(C=C1)CCCC1C=CC=CC=1')
        DPP.generate_resonance_structures()
        formula = DPP.molecule[0].get_formula()
        if formula in rmg_test.reaction_model.species_dict:
            rmg_test.reaction_model.species_dict[formula].append(DPP)
        else:
            rmg_test.reaction_model.species_dict[formula] = [DPP]

        mol_test = Molecule().from_adjacency_list(
"""
1     C u0 p0 c0 {2,S} {3,S} {16,S} {17,S}
2     C u0 p0 c0 {1,S} {4,S} {18,S} {19,S}
3     C u0 p0 c0 {1,S} {5,S} {20,S} {21,S}
4     C u0 p0 c0 {2,S} {6,B} {7,B}
5     C u0 p0 c0 {3,S} {8,D} {9,S}
6     C u0 p0 c0 {4,B} {10,B} {22,S}
7     C u0 p0 c0 {4,B} {12,B} {24,S}
8     C u0 p0 c0 {5,D} {14,S} {27,S}
9     C u0 p0 c0 {5,S} {15,D} {28,S}
10    C u0 p0 c0 {6,B} {11,B} {23,S}
11    C u0 p0 c0 {10,B} {12,B} {25,S}
12    C u0 p0 c0 {7,B} {11,B} {26,S}
13    C u0 p0 c0 {14,D} {15,S} {29,S}
14    C u0 p0 c0 {8,S} {13,D} {30,S}
15    C u0 p0 c0 {9,D} {13,S} {31,S}
16    H u0 p0 c0 {1,S}
17    H u0 p0 c0 {1,S}
18    H u0 p0 c0 {2,S}
19    H u0 p0 c0 {2,S}
20    H u0 p0 c0 {3,S}
21    H u0 p0 c0 {3,S}
22    H u0 p0 c0 {6,S}
23    H u0 p0 c0 {10,S}
24    H u0 p0 c0 {7,S}
25    H u0 p0 c0 {11,S}
26    H u0 p0 c0 {12,S}
27    H u0 p0 c0 {8,S}
28    H u0 p0 c0 {9,S}
29    H u0 p0 c0 {13,S}
30    H u0 p0 c0 {14,S}
31    H u0 p0 c0 {15,S}
""")
        spec = rmg_test.reaction_model.check_for_existing_species(mol_test)
        self.assertIsNotNone(spec)


def find_target_rxns_containing(mol1, mol2, reactions):
    target_rxns = []
    for rxn in reactions:
        reactants = rxn.reactants
        products = rxn.products
        rxn_specs = reactants + products
        for rxn_spec in rxn_specs:
            if rxn_spec.is_isomorphic(mol1):
                for rxn_spec1 in rxn_specs:
                    if rxn_spec1.is_isomorphic(mol2):
                        target_rxns.append(rxn)
    return target_rxns


class TestRMGScript(unittest.TestCase):
    """
    Contains unit tests for rmg.py
    """

    def test_parse_command_line_arguments_defaults(self):
        """
        Test the default values for the parseCommandLineArguments module
        """

        # Acquire default arguments
        args = parse_command_line_arguments(['input.py'])

        # Test default values
        self.assertEqual(args.walltime, '00:00:00:00')
        self.assertEqual(args.output_directory, os.path.abspath(os.path.dirname('./')))
        self.assertEqual(args.debug, False)
        self.assertEqual(args.file, 'input.py')
        self.assertEqual(args.maxiter, None)
        self.assertEqual(args.kineticsdatastore, False)
        self.assertEqual(args.postprocess, False)
        self.assertEqual(args.profile, False)
        self.assertEqual(args.quiet, False)
        self.assertEqual(args.restart, '')
        self.assertEqual(args.verbose, False)

    def test_parse_command_line_non_defaults(self):
        """
        Test user command line inputs into rmg.py
        """

        # Acquire arguments
        args = parse_command_line_arguments(['other_name.py', '-d', '-o', '/test/output/dir/', '-r', 'test/seed/', '-P',
                                             '-t', '01:20:33:45', '-k', '-i', '100'])

        # Test expected values
        self.assertEqual(args.walltime, '01:20:33:45')
        self.assertEqual(args.output_directory, '/test/output/dir/')
        self.assertEqual(args.debug, True)
        self.assertEqual(args.file, 'other_name.py')
        self.assertEqual(args.maxiter, 100)
        self.assertEqual(args.kineticsdatastore, True)
        self.assertEqual(args.postprocess, True)
        self.assertEqual(args.profile, True)
        self.assertEqual(args.restart, 'test/seed/')
