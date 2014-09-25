#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This script contains unit tests of the :mod:`rmgpy.quantity` module.
"""
import unittest
import numpy
import os
from rmgpy.cantherm import CanTherm
import rmgpy.constants as constants
################################################################################

class CommonTest(unittest.TestCase):
    """
    Contains unit tests of the Cantherm common functions.
    """

    def test_checkConformerEnergy(self):
        """
        test the checkConformerEnergy function with an list of energies.
        """
        Vlist = [-272.2779012225, -272.2774933703, -272.2768397635, -272.2778432059, -272.278645477, -272.2789602654, -272.2788749196, -272.278496709, -272.2779350675, -272.2777008843, -272.2777167286, -272.2780937643, -272.2784838846, -272.2788050464, -272.2787865352, -272.2785091607, -272.2779977452, -272.2777957743, -272.2779134906, -272.2781827547, -272.278443339, -272.2788244214, -272.2787748749]
        Vlist = numpy.array(Vlist, numpy.float64)
        Vdiff = (Vlist[0] - numpy.min(Vlist)) * constants.E_h * constants.Na / 1000
        self.assertAlmostEqual(Vdiff / 2.7805169838282797, 1, 5)


class testCanthermJob(unittest.TestCase):
    """
    Contains unit tests of the Cantherm module and its interactions with other RMG modules.
    """
    def setUp(self):

        cantherm = CanTherm()
        
        jobList = cantherm.loadInputFile(os.path.join(os.path.dirname(os.path.abspath(__file__)),r'files/methoxy.py'))
        pdepjob = jobList[-1]
        self.kineticsjob = jobList[0]
        pdepjob.activeJRotor = True
        network = pdepjob.network
        self.Nisom = len(network.isomers)
        self.Nreac = len(network.reactants)
        self.Nprod = len(network.products)
        self.Npath = len(network.pathReactions)
        self.PathReaction2 = network.pathReactions[2]
        self.TminValue = pdepjob.Tmin.value
        self.Tmaxvalue = pdepjob.Tmax.value
        self.TmaxUnits = pdepjob.Tmax.units
        self.TlistValue = pdepjob.Tlist.value
        self.PminValue = pdepjob.Pmin.value
        self.Pcount = pdepjob.Pcount
        self.Tcount = pdepjob.Tcount
        self.GenTlist = pdepjob.generateTemperatureList()
        self.PlistValue = pdepjob.Plist.value
        self.maximumGrainSizeValue = pdepjob.maximumGrainSize.value
        self.method = pdepjob.method
        self.rmgmode = pdepjob.rmgmode

# test Cantherm's interactions with the network module
    def testNisom(self):
        """
        Test the number of isomers identified.
        """
        self.assertEqual(self.Nisom, 2, msg=None)

    def testNreac(self):
        """
        Test the number of reactants identified.
        """
        self.assertEqual(self.Nreac, 1, msg=None)

    def testNprod(self):
        """
        Test the number of products identified.
        """
        self.assertEqual(self.Nprod, 1, msg=None)

    def testNpathReactions(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.Npath, 3, msg=None)

    def testPathReactions(self):
        """
        Test a path reaction label
        """
        self.assertEqual(str(self.PathReaction2), 'CH2OH <=> methoxy', msg=None)

# test Cantherm's interactions with the pdep module
    def testTemperaturesUnits(self):
        """
        Test the Temperature Units.
        """
        self.assertEqual(str(self.TmaxUnits), 'K', msg=None)

    def testTemperaturesValue(self):
        """
        Test the temperature value.
        """
        self.assertEqual(self.TminValue, 450.0, msg=None)

    def testTemperaturesList(self):
        """
        Test the temperature list.
        """
        self.assertEqual(numpy.array_equal(self.TlistValue, numpy.array([450, 500, 678, 700])), True, msg=None)

    def testPminValue(self):
        """
        Test the minimum pressure value.
        """
        self.assertEqual("%0.7f" % self.PminValue, str(0.0101325), msg=None)

    def testPcount(self):
        """
        Test the number pressures specified.
        """
        self.assertEqual(self.Pcount, 7, msg=None)

    def testTcount(self):
        """
        Test the number temperatures specified.
        """
        self.assertEqual(self.Tcount, 4, msg=None)

    def testPressureList(self):
        """
        Test the pressure list.
        """
        self.assertEqual(numpy.array_equal(self.PlistValue, numpy.array([0.01, 0.1, 1, 3, 10, 100, 1000])), True, msg=None)

    def testGenerateTemperatureList(self):
        """
        Test the generated temperature list.
        """
        self.assertEqual(list(self.GenTlist), [450.0, 500.0, 678.0, 700.0], msg=None)

    def testmaximumGrainSizeValue(self):
        """
        Test the max grain size value.
        """
        self.assertEqual(self.maximumGrainSizeValue, 0.5, msg=None)

    def testMethod(self):
        """
        Test the master equation solution method chosen.
        """
        self.assertEqual(self.method, 'modified strong collision', msg=None)

    def testRmgmode(self):
        """
        Test the whether or not RMG mode is turned on.
        """
        self.assertEqual(self.rmgmode, False, msg=None)

# Test cantherms interactions with the kinetics module
    def testCalculateTSTRateCoefficient(self):
        """
        Test the calculation of the high-pressure limit rate coef for one of the kinetics jobs at Tmin and Tmax.
        """
        self.assertEqual("%0.7f" % self.kineticsjob.reaction.calculateTSTRateCoefficient(self.TminValue), str(46608.5904933), msg=None)
        self.assertEqual("%0.5f" % self.kineticsjob.reaction.calculateTSTRateCoefficient(self.Tmaxvalue), str(498796.64535), msg=None)

    def testTunneling(self):
        """
        Test the whether or not tunneling has been included in a specific kinetics job.
        """
        self.assertEqual(self.kineticsjob.reaction.transitionState.tunneling, None, msg=None)

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
