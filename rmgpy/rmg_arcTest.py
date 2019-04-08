#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
This script contains unit tests of the `rmg_arc` module.
"""

import os
import shutil
import unittest

import rmg_arc
from rmgpy import settings
from rmgpy.species import Species
from rmgpy.thermo.model import HeatCapacityModel

################################################################################


class TestRMGARC(unittest.TestCase):
    """
    Contains unit tests that ensure that the functions in rmg_arc work as expected
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run ONCE before all unit tests in this class.
        """
        input_file = """
# Data sources
 database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.9,
    },
    terminationTime=(1e6,'s'),
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)


arc(
SA species: 10
SA reactions: 10
SA pdep: True
all core species: False
SA observables: True
collision violators: True
kwargs: {ess_settings: {'gaussian': 'c3ddb', 'molpro': 'pharos', 'qchem': 'pharos'},
          level_of_theory: 'b3lyp/6-311+g(3df,2p)//b3lyp/6-31+g(d,p)',
           scan_rotors: False, fine: False}
)

"""
        base_path = os.path.join(settings['test_data.directory'], 'arc_rmg')
        if not os.path.exists(base_path):
            os.mkdir(base_path)
        cls.input_file_path = os.path.join(base_path, 'rmg_input_1.py')
        with open(cls.input_file_path, 'w') as f:
            f.write(input_file)

    def test_parse_input_file(self):
        """Test parsing the arc argument from an input file"""
        rmg_input_file, arc_arguments = rmg_arc.parse_input_file(self.input_file_path)
        self.assertFalse('arc' in rmg_input_file)
        self.assertEqual(arc_arguments['SA species'], 10)
        self.assertEqual(arc_arguments['SA reactions'], 10)
        self.assertTrue(arc_arguments['SA pdep'])
        self.assertFalse(arc_arguments['all core species'])
        self.assertTrue(arc_arguments['SA observables'])
        self.assertTrue(arc_arguments['collision violators'])
        self.assertEqual(arc_arguments['kwargs']['ess_settings']['gaussian'], 'c3ddb')
        self.assertEqual(arc_arguments['kwargs']['ess_settings']['molpro'], 'pharos')
        self.assertEqual(arc_arguments['kwargs']['ess_settings']['qchem'], 'pharos')
        self.assertEqual(arc_arguments['kwargs']['level_of_theory'], 'b3lyp/6-311+g(3df,2p)//b3lyp/6-31+g(d,p)')
        self.assertFalse(arc_arguments['kwargs']['scan_rotors'])
        self.assertFalse(arc_arguments['kwargs']['fine'])

    def test_determine_species_to_calculate(self):
        """Test that we correctly determine the species to be calculated from an RMG job"""
        # TODO: add an actual case with SA and PDep and coll violators
        pass

    def test_print_summary(self):
        pass

    def test_add_rmg_libraries(self):
        pass

    def testget_unconverged_species(self):
        pass

    def test_calc_based_on_thermo_comment(self):
        """Test which species are selected for calculation based on their thermo comment"""
        spc1 = Species().fromSMILES('CC')
        spc1.thermo = HeatCapacityModel()
        spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH)' \
                              ' + group(Cs-CsHHH) + radical(RCCJ)'
        spc2 = Species().fromSMILES('CC')
        spc2.thermo = HeatCapacityModel()
        spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'
        spc3 = Species().fromSMILES('CC')
        spc3.thermo = HeatCapacityModel()
        spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'

        self.assertTrue(rmg_arc.calc_based_on_thermo_comment(spc1))
        self.assertTrue(rmg_arc.calc_based_on_thermo_comment(spc2))
        self.assertFalse(rmg_arc.calc_based_on_thermo_comment(spc3))

    def test_add_new_library_to_line(self):
        """Test that a new library label is correctly added to an existing RMG input file"""
        lines = ["thermoLibraries=['BurkeH2O2','primaryNS'],\n",
                 "thermoLibraries=['BurkeH2O2','primaryNS'], reactionLibraries=['BurkeH2O2inN2'],",
                 "seedMechanisms=[], thermoLibraries=['BurkeH2O2','primaryNS'], reactionLibraries=['BurkeH2O2inN2']"]
        library = 'arc_thermo_library'
        new_lines = []
        for line in lines:
            new_lines.append(rmg_arc.add_new_library_to_line(line, library))
        self.assertEqual(new_lines[0], "thermoLibraries=['BurkeH2O2','primaryNS','arc_thermo_library'],\n")
        self.assertEqual(new_lines[1], "thermoLibraries=['BurkeH2O2','primaryNS','arc_thermo_library'],"
                                       " reactionLibraries=['BurkeH2O2inN2'],")
        self.assertEqual(new_lines[2], "seedMechanisms=[], thermoLibraries=['BurkeH2O2','primaryNS',"
                                       "'arc_thermo_library'], reactionLibraries=['BurkeH2O2inN2']")

    @classmethod
    def tearDownClass(cls):
        """
        A function that is run ONCE after all unit tests in this class.
        """
        base_path = os.path.join(settings['test_data.directory'], 'arc_rmg')
        shutil.rmtree(base_path)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
