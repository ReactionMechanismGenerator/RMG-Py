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

import unittest
import logging
from copy import deepcopy

from rmgpy.pdep.network import Network
from rmgpy.pdep.configuration import Configuration
from rmgpy.reaction import Reaction
from rmgpy.statmech.vibration import  HarmonicOscillator
from rmgpy.statmech.torsion import  HinderedRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.species import Species, TransitionState
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.transport import TransportData
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.rotation import  NonlinearRotor
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.kinetics.arrhenius import Arrhenius

###################################################

class TestPdep(unittest.TestCase):

    def setUp(self):
        """
        A method that is run before each unit test in this class.
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
        
        self.nC4H10O.fromSMILES('CCCCO')
        
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
            
        self.nC4H8.fromSMILES('CCC=C')
        
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
        
        self.H2O.fromSMILES('O')

        self.N2 = Species(
            label = 'N2',
            molecularWeight = (28.04,"g/mol"),
            transportData=TransportData(sigma=(3.41, "angstrom"), epsilon=(124, "K")),
            energyTransferModel = None,
        )
        
        self.N2.fromSMILES('N#N')
        
        logging.error('to TS')
                           
        self.TS = TransitionState(
            label = 'TS',
            conformer = Conformer(
                E0 = (-42.4373,"kJ/mol"),
                modes = [
                    IdealGasTranslation(mass=(74.07,"g/mol")),
                    NonlinearRotor(inertia=([40.518,232.666,246.092],"u*angstrom**2"), symmetry=1, quantum=False),
                    HarmonicOscillator(frequencies=([134.289,302.326,351.792,407.986,443.419,583.988,699.001,766.1,777.969,829.671,949.753,994.731,1013.59,1073.98,1103.79,1171.89,1225.91,1280.67,1335.08,1373.9,1392.32,1417.43,1469.51,1481.61,1490.16,1503.73,1573.16,2972.85,2984.3,3003.67,3045.78,3051.77,3082.37,3090.44,3190.73,3708.52],"kayser")),
                    HinderedRotor(inertia=(2.68206,"amu*angstrom^2"), symmetry=3, barrier=(3.35244,"kcal/mol")),
                    HinderedRotor(inertia=(9.77669,"amu*angstrom^2"), symmetry=1, fourier=([[0.208938,-1.55291,-4.05398,-0.105798,-0.104752], [2.00518,-0.020767,-0.333595,0.137791,-0.274578]],"kJ/mol")),
                ],
                spinMultiplicity = 1,
                opticalIsomers = 1,
            ),
            frequency=(-2038.34,'cm^-1'),
        )
        
        self.reaction = Reaction(
            label = 'dehydration',
            reactants = [self.nC4H10O],
            products = [self.nC4H8, self.H2O],
            transitionState = self.TS,
            kinetics = Arrhenius(A=(0.0387, 'm^3/(mol*s)'), n=2.7, Ea=(2.6192e4, 'J/mol'), T0=(1, 'K'))
        )
        
        self.network = Network(
            label = 'n-butanol',
            isomers = [Configuration(self.nC4H10O)],
            reactants = [],
            products = [Configuration(self.nC4H8, self.H2O)],
            pathReactions = [self.reaction],
            bathGas = {self.N2: 1.0},
        )
        
        self.pdepnetwork = deepcopy(self.network)
        self.pdepnetwork.__class__ = PDepNetwork
        self.pdepnetwork.source = [self.pdepnetwork.isomers[0].species[0]]
        self.pdepnetwork.index = 1
        self.pdepnetwork.explored = []

    def test_energy_filter(self):
        rxns = self.pdepnetwork.get_energy_filtered_reactions(1000.0,0.0)
        self.assertEquals(len(rxns),1)
        self.assertEquals(rxns[0],self.pdepnetwork.pathReactions[0])

    def test_flux_filter(self):
        prods = self.pdepnetwork.get_rate_filtered_products(1000.0,100000.0,1.0)
        self.assertEquals(len(prods),0)

if __name__ == '__main__':
    unittest.main()
