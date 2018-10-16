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
This script contains unit tests of the :mod:`rmgpy.kinetics.diffusionLimited` module.
"""

import unittest
import os

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.kinetics import Arrhenius
from rmgpy.data.solvation import SolvationDatabase
from rmgpy.data.thermo import ThermoData
from rmgpy.kinetics.diffusionLimited import diffusionLimiter
################################################################################

class TestDiffusionLimited(unittest.TestCase):
    """
    Contains unit tests of the DiffusionLimited class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        octyl_pri = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs=[-0.772759,0.093255,-5.84447e-05,1.8557e-08,-2.37127e-12,-3926.9,37.6131], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[25.051,0.036948,-1.25765e-05,1.94628e-09,-1.12669e-13,-13330.1,-102.557], Tmin=(1390,'K'), Tmax=(5000,'K'))
            ],
            Tmin=(298,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'),
            comment="""Thermo library: JetSurF0.2"""), molecule=[Molecule(SMILES="[CH2]CCCCCCC")])
        octyl_sec = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs=[-0.304233,0.0880077,-4.90743e-05,1.21858e-08,-8.87773e-13,-5237.93,36.6583], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[24.9044,0.0366394,-1.2385e-05,1.90835e-09,-1.10161e-13,-14713.5,-101.345], Tmin=(1383,'K'), Tmax=(5000,'K'))
            ],
            Tmin=(298,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'),
            comment="""Thermo library: JetSurF0.2"""), molecule=[Molecule(SMILES="CC[CH]CCCCC")])
        ethane = Species(label="", thermo=ThermoData(
            Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([10.294,12.643,14.933,16.932,20.033,22.438,26.281],'cal/(mol*K)'), H298=(12.549,'kcal/mol'), S298=(52.379,'cal/(mol*K)'),
            Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CH"""),
            molecule=[Molecule(SMILES="C=C")])
        decyl = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs=[-1.31358,0.117973,-7.51843e-05,2.43331e-08,-3.17523e-12,-9689.68,43.501], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[31.5697,0.0455818,-1.54995e-05,2.39711e-09,-1.3871e-13,-21573.8,-134.709], Tmin=(1390,'K'), Tmax=(5000,'K'))
            ],
            Tmin=(298,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(719.202,'J/(mol*K)'),
            comment="""Thermo library: JetSurF0.2"""), molecule=[Molecule(SMILES="[CH2]CCCCCCCCC")])
        acetone = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs = [3.75568, 0.0264934, -6.55661e-05, 1.94971e-07, -1.82059e-10, -27905.3, 9.0162], Tmin = (10, 'K'), Tmax = (422.477, 'K')),
            NASAPolynomial(coeffs = [0.701289, 0.0344988, -1.9736e-05, 5.48052e-09, -5.92612e-13, -27460.6, 23.329],Tmin = (422.477, 'K'),Tmax = (3000, 'K'))
            ],
            Tmin = (10, 'K'), Tmax = (3000, 'K'), E0 = (-232.025, 'kJ/mol'), Cp0 = (33.2579, 'J/(mol*K)'), CpInf = (232.805, 'J/(mol*K)')),
            molecule=[Molecule(SMILES="CC(=O)C")])
        peracetic_acid = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs = [3.81786, 0.016419, 3.32204e-05, -8.98403e-08, 6.63474e-11, -42057.8, 9.65245], Tmin = (10, 'K'), Tmax = (354.579, 'K')),
            NASAPolynomial(coeffs = [2.75993, 0.0283534, -1.72659e-05, 5.08158e-09, -5.77773e-13, -41982.8, 13.6595], Tmin = (354.579, 'K'), Tmax = (3000, 'K'))
            ],
            Tmin = (10, 'K'), Tmax = (3000, 'K'), E0 = (-349.698, 'kJ/mol'),Cp0 = (33.2579, 'J/(mol*K)'), CpInf = (199.547, 'J/(mol*K)')),
            molecule=[Molecule(SMILES="CC(=O)OO")])
        acetic_acid = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs = [3.97665, 0.00159915, 8.5542e-05, -1.76486e-07, 1.20201e-10, -53911.5, 8.99309], Tmin = (10, 'K'), Tmax = (375.616, 'K')),
            NASAPolynomial(coeffs = [1.57088, 0.0272146, -1.67357e-05, 5.01453e-09, -5.82273e-13, -53730.7, 18.2442], Tmin = (375.616, 'K'), Tmax = (3000, 'K'))
            ],
            Tmin = (10, 'K'), Tmax = (3000, 'K'), E0 = (-448.245, 'kJ/mol'), Cp0 = (33.2579, 'J/(mol*K)'), CpInf = (182.918, 'J/(mol*K)')),
            molecule=[Molecule(SMILES="CC(=O)O")])
        criegee = Species(label="", thermo=NASA(polynomials=[
            NASAPolynomial(coeffs = [3.23876, 0.0679583, -3.35611e-05, 7.91519e-10, 3.13038e-12, -77986, 13.6438], Tmin = (10, 'K'), Tmax = (1053.46, 'K')),
            NASAPolynomial(coeffs = [9.84525, 0.0536795, -2.86165e-05, 7.39945e-09, -7.48482e-13, -79977.6, -21.4187], Tmin = (1053.46, 'K'), Tmax = (3000, 'K'))
            ],
            Tmin = (10, 'K'), Tmax = (3000, 'K'), E0 = (-648.47, 'kJ/mol'),Cp0 = (33.2579, 'J/(mol*K)'), CpInf = (457.296, 'J/(mol*K)')),
            molecule=[Molecule(SMILES="CC(=O)OOC(C)(O)C")])
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
        self.solvent = 'octane'
        diffusionLimiter.enable(self.database.getSolventData(self.solvent), self.database)
        self.T = 298
        self.uni_reaction = Reaction(reactants=[octyl_pri], products=[octyl_sec])
        self.uni_reaction.kinetics = Arrhenius(A=(2.0, '1/s'), n=0, Ea=(0,'kJ/mol'))
        self.bi_uni_reaction = Reaction(reactants=[octyl_pri, ethane], products=[decyl])
        self.bi_uni_reaction.kinetics = Arrhenius(A=(1.0E-22, 'cm^3/molecule/s'), n=0, Ea=(0,'kJ/mol'))
        self.tri_bi_reaction = Reaction(reactants=[acetone, peracetic_acid, acetic_acid],
                                        products=[criegee, acetic_acid])
        self.tri_bi_reaction.kinetics = Arrhenius(A=(1.07543e-11, 'cm^6/(mol^2*s)'), n=5.47295, Ea=(-38.5379, 'kJ/mol'))
        self.intrinsic_rates = {
            self.uni_reaction: self.uni_reaction.kinetics.getRateCoefficient(self.T, P=100e5),
            self.bi_uni_reaction: self.bi_uni_reaction.kinetics.getRateCoefficient(self.T, P=100e5),
            self.tri_bi_reaction: self.tri_bi_reaction.kinetics.getRateCoefficient(self.T, P=100e5),
            }

    def tearDown(self):
        diffusionLimiter.disable()

    def testGetEffectiveRateUnimolecular(self):
        """
        Tests that the effective rate is the same as the intrinsic rate for
        unimiolecular reactions.
        """
        effective_rate = diffusionLimiter.getEffectiveRate(self.uni_reaction, self.T)
        self.assertEqual(effective_rate, self.intrinsic_rates[self.uni_reaction])

    def testGetEffectiveRate2to1(self):
        """
        Tests that the effective rate is limited in the forward direction for
        a 2 -> 1 reaction
        """
        effective_rate = diffusionLimiter.getEffectiveRate(self.bi_uni_reaction, self.T)
        self.assertTrue(effective_rate < self.intrinsic_rates[self.bi_uni_reaction])
        self.assertTrue(effective_rate >= 0.2 * self.intrinsic_rates[self.bi_uni_reaction])

    def testGetEffectiveRate3to2(self):
        """
        Tests that the effective rate is limited for a 3 -> 2 reaction
        """
        effective_rate = diffusionLimiter.getEffectiveRate(self.tri_bi_reaction, self.T)
        self.assertTrue(effective_rate < self.intrinsic_rates[self.tri_bi_reaction])
        self.assertTrue(effective_rate >= 0.2 * self.intrinsic_rates[self.tri_bi_reaction])

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
