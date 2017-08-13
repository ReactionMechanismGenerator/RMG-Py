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
This script contains unit tests of the :mod:`rmgpy.pdep.network` module.
"""

import unittest

from rmgpy.pdep.configuration import Configuration
from rmgpy.transport import TransportData
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.rotation import NonlinearRotor
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.statmech.torsion import HinderedRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.species import Species
from rmgpy.pdep.collision import SingleExponentialDown

################################################################################

class TestConfiguration(unittest.TestCase):
    """
    Contains unit tests of the :class:`Network` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.nC4H10O = Species(
            label = 'n-C4H10O',
            conformer = Conformer(
                E0 = (-317.807,'kJ/mol'),
                modes = [
                    IdealGasTranslation(mass=(74.07,"g/mol")),
                    NonlinearRotor(inertia=([41.5091,215.751,233.258],"amu*angstrom^2"), symmetry=1),
                    HarmonicOscillator(frequencies=([240.915,341.933,500.066,728.41,809.987,833.93,926.308,948.571,1009.3,1031.46,1076,1118.4,1184.66,1251.36,1314.36,1321.42,1381.17,1396.5,1400.54,1448.08,1480.18,1485.34,1492.24,1494.99,1586.16,2949.01,2963.03,2986.19,2988.1,2995.27,3026.03,3049.05,3053.47,3054.83,3778.88],"cm^-1")),
                    HinderedRotor(inertia=(0.854054,"amu*angstrom^2"), symmetry=1, fourier=([[0.25183,-1.37378,-2.8379,0.0305112,0.0028088], [0.458307,0.542121,-0.599366,-0.00283925,0.0398529]],"kJ/mol")),
                    HinderedRotor(inertia=(8.79408,"amu*angstrom^2"), symmetry=1, fourier=([[0.26871,-0.59533,-8.15002,-0.294325,-0.145357], [1.1884,0.99479,-0.940416,-0.186538,0.0309834]],"kJ/mol")),
                    HinderedRotor(inertia=(7.88153,"amu*angstrom^2"), symmetry=1, fourier=([[-4.67373,2.03735,-6.25993,-0.27325,-0.048748], [-0.982845,1.76637,-1.57619,0.474364,-0.000681718]],"kJ/mol")),
                    HinderedRotor(inertia=(2.81525,"amu*angstrom^2"), symmetry=3, barrier=(2.96807,"kcal/mol")),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
            molecularWeight = (74.07,"g/mol"),
            transportData=TransportData(sigma=(5.94, 'angstrom'), epsilon=(559, 'K')),
            energyTransferModel = SingleExponentialDown(alpha0=(447.5*0.011962,"kJ/mol"), T0=(300,"K"), n=0.85),
        )
        
        self.nC4H8 = Species(
            label = 'n-C4H8',
            conformer = Conformer(
                E0 = (-17.8832,'kJ/mol'),
                modes = [
                    IdealGasTranslation(mass=(56.06,"g/mol")),
                    NonlinearRotor(inertia=([22.2748,122.4,125.198],"amu*angstrom^2"), symmetry=1),
                    HarmonicOscillator(frequencies=([308.537,418.67,636.246,788.665,848.906,936.762,979.97,1009.48,1024.22,1082.96,1186.38,1277.55,1307.65,1332.87,1396.67,1439.09,1469.71,1484.45,1493.19,1691.49,2972.12,2994.31,3018.48,3056.87,3062.76,3079.38,3093.54,3174.52],"cm^-1")),
                    HinderedRotor(inertia=(5.28338,"amu*angstrom^2"), symmetry=1, fourier=([[-0.579364,-0.28241,-4.46469,0.143368,0.126756], [1.01804,-0.494628,-0.00318651,-0.245289,0.193728]],"kJ/mol")),
                    HinderedRotor(inertia=(2.60818,"amu*angstrom^2"), symmetry=3, fourier=([[0.0400372,0.0301986,-6.4787,-0.0248675,-0.0324753], [0.0312541,0.0538,-0.493785,0.0965968,0.125292]],"kJ/mol")),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
        )
        
        self.H2O = Species(
            label = 'H2O',
            conformer = Conformer(
                E0 = (-269.598,'kJ/mol'),
                modes = [
                    IdealGasTranslation(mass=(18.01,"g/mol")),
                    NonlinearRotor(inertia=([0.630578,1.15529,1.78586],"amu*angstrom^2"), symmetry=2),
                    HarmonicOscillator(frequencies=([1622.09,3771.85,3867.85],"cm^-1")),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
        )

        self.configuration = Configuration(self.nC4H8, self.H2O)
        self.configuration_unimolecular = Configuration(self.nC4H10O)

    def test_repr(self):
        """
        Test that the `repr` representation contains desired properties.
        """
        output= repr(self.configuration)
        # ensure species strings
        labels = ['H2O','n-C4H8']
        for label in labels:
            self.assertIn(label, output)

        # ensure classes are used as well
        attributes = ['Configuration','Species','Conformer','NonlinearRotor',
                   'HarmonicOscillator','frequencies','IdealGasTranslation',
                   'HinderedRotor','E0','mass','symmetry','fourier']
        for label in attributes:
            self.assertIn(label, output)

    def test_str(self):
        """
        Test that the string representation contains desired properties.
        """
        output = str(self.configuration)
        # ensure species strings
        labels = ['H2O','n-C4H8']
        for label in labels:
            self.assertIn(label, output)

        # ensure this extra fluff is not in Network string
        attributes = ['Species','Conformer','Molecule','NonlinearRotor',
                   'HarmonicOscillator','frequencies','spinMultiplicity','TransportData',
                   'molecularWeight','SingleExponentialDown']
        for label in attributes:
            self.assertNotIn(label, output)


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
