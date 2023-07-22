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
This module contains unit tests of the :mod:`arkane.ess.gaussian` module.
"""

import os
import shutil
import unittest

from nose.plugins.attrib import attr

import rmgpy

from arkane.ess.gaussian import GaussianLog
from arkane.main import Arkane
from arkane.output import prettify, get_str_xyz
from rmgpy.species import Species

################################################################################


@attr('functional')
class OutputTest(unittest.TestCase):
    """
    Contains functional tests for Arkane's output module.
    """
    def test_prettify(self):
        """Test that the prettify function works for an Arkane job"""
        benzyl_path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)),
                                   'examples', 'arkane', 'species', 'Benzyl')

        arkane = Arkane(input_file=os.path.join(benzyl_path, 'input.py'), output_directory=benzyl_path)
        arkane.plot = False
        arkane.execute()
        with open(os.path.join(benzyl_path, 'output.py'), 'r') as f:
            lines = f.readlines()
        self.assertIn('conformer(\n', lines)
        self.assertIn("    E0 = (193.749, 'kJ/mol'),\n", lines)
        self.assertIn('thermo(\n', lines)
        self.assertIn("        Cp0 = (33.2579, 'J/(mol*K)'),\n", lines)

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        benzyl_path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)),
                                   'examples', 'arkane', 'species', 'Benzyl')
        extensions_to_delete = ['pdf', 'csv', 'txt', 'inp']
        files_to_delete = ['arkane.log', 'output.py']
        for name in os.listdir(benzyl_path):
            item_path = os.path.join(benzyl_path, name)
            if os.path.isfile(item_path):
                extension = name.split('.')[-1]
                if name in files_to_delete or extension in extensions_to_delete:
                    os.remove(item_path)
            else:
                if os.path.split(item_path)[-1] in ['r0']:
                    continue
                # This is a sub-directory. remove.
                shutil.rmtree(item_path)


class OutputUnitTest(unittest.TestCase):
    """
    Contains unit tests for the Arkane's output module.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.data_path = os.path.join(os.path.dirname(__file__), 'data')

    def test_prettify(self):
        """Test that ``prettify`` returns the expected result"""
        input_str = ("conformer(label='C7H7', E0=(193.749,'kJ/mol'), modes=[IdealGasTranslation(mass=(91.0548,'amu')), "
                     "NonlinearRotor(inertia=([91.0567,186.675,277.733],'amu*angstrom^2'), symmetry=2), "
                     "HarmonicOscillator(frequencies=([199.381,360.536,413.795,480.347,536.285,630.723,687.118,709.613,"
                     "776.662,830.404,834.386,901.841,973.498,975.148,993.349,998.606,1040.14,1120.69,1179.22,1189.07,"
                     "1292.86,1332.91,1357.18,1479.46,1495.36,1507.91,1583.14,1604.63,3156.85,3170.22,3172.78,3185.05,"
                     "3189.8,3203.23,3253.99],'cm^-1')), HinderedRotor(inertia=(1.70013,'amu*angstrom^2'), symmetry=2, "
                     "fourier=([[-0.315923,-27.7665,0.177678,3.2437,0.0509515],[-0.00164953,-0.0021925,-0.00386396,"
                     "-0.000912068,0.00274206]],'kJ/mol'), quantum=True, semiclassical=False)], spin_multiplicity=2, "
                     "optical_isomers=1)")
        expected_output = """conformer(
    label = 'C7H7',
    E0 = (193.749, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(91.0548, 'amu')),
        NonlinearRotor(
            inertia = ([91.0567, 186.675, 277.733], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(
            frequencies = ([199.381, 360.536, 413.795, 480.347, 536.285, 630.723, 687.118, 709.613, 776.662, 830.404, 834.386, 901.841, 973.498, 975.148, 993.349, 998.606, 1040.14, 1120.69, 1179.22, 1189.07, 1292.86, 1332.91, 1357.18, 1479.46, 1495.36, 1507.91, 1583.14, 1604.63, 3156.85, 3170.22, 3172.78, 3185.05, 3189.8, 3203.23, 3253.99], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (1.70013, 'amu*angstrom^2'),
            symmetry = 2,
            fourier = (
                [
                    [-0.315923, -27.7665, 0.177678, 3.2437, 0.0509515],
                    [-0.00164953, -0.0021925, -0.00386396, -0.000912068, 0.00274206],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)"""
        self.assertEqual(prettify(input_str), expected_output)

    def test_get_str_xyz(self):
        """Test generating an xyz string from the species.conformer object"""
        log = GaussianLog(os.path.join(self.data_path, 'gaussian', 'ethylene_G3.log'))
        conformer = log.load_conformer()[0]
        coords, number, mass = log.load_geometry()
        conformer.coordinates, conformer.number, conformer.mass = (coords, "angstroms"), number, (mass, "amu")
        spc1 = Species(smiles='C=C')
        spc1.conformer = conformer
        xyz_str = get_str_xyz(spc1)
        expected_xyz_str = """C       0.00545100    0.00000000    0.00339700
H       0.00118700    0.00000000    1.08823200
H       0.97742900    0.00000000   -0.47841600
C      -1.12745800    0.00000000   -0.70256500
H      -1.12319800    0.00000000   -1.78740100
H      -2.09943900    0.00000000   -0.22075700"""
        self.assertEqual(xyz_str, expected_xyz_str)

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
