#!/usr/bin/env python
# encoding: utf-8

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
from rmgpy.kinetics.diffusionLimited import DiffusionLimited, diffusionLimiter
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
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
        self.solvent = 'octane'
        diffusionLimiter.enable(self.database.getSolventData(self.solvent), self.database)
        self.T = 298
        self.uni_reaction = Reaction(reactants=[octyl_pri], products=[octyl_sec])
        self.uni_reaction.kinetics = Arrhenius(A=(2.0, '1/s'), n=0, Ea=(0,'kJ/mol'))
        self.bi_uni_reaction = Reaction(reactants=[octyl_pri, ethane], products=[decyl])
        self.bi_uni_reaction.kinetics = Arrhenius(A=(1.0E-22, 'cm^3/molecule/s'), n=0, Ea=(0,'kJ/mol'))
        self.intrinsic_rates = {
            self.uni_reaction: self.uni_reaction.kinetics.getRateCoefficient(self.T, P=100e5),
            self.bi_uni_reaction: self.bi_uni_reaction.kinetics.getRateCoefficient(self.T, P=100e5),
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

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
