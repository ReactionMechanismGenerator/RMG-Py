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
This script contains unit tests of the :mod:`arkane.multidimensionalTorsions` module.
"""
from __future__ import division

import unittest
import os
import zipfile
import shutil

from rmgpy.statmech.ndTorsions import HinderedRotor2D, HinderedRotorClassicalND
from arkane.util import determine_qm_software

class TestHinderedRotor2D(unittest.TestCase):
    """
    Contains unit tests of the StatmechJob class.
    """
    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        cls.path = os.path.join(os.path.dirname(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]),
                                                         'arkane', 'data', 'CH2CHOOH', 'CH2CHOOHscans')
        if not os.path.exists(cls.path):
            zippath = os.path.join(os.path.dirname(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]),
                                                             'arkane', 'data', 'CH2CHOOH', 'CH2CHOOHscans.zip')
            with zipfile.ZipFile(zippath,'r') as zip_ref:
                zip_ref.extractall(os.path.split(cls.path)[0])
        
        cls.hd2d= HinderedRotor2D(calcPath=cls.path,name="r0",torsigma1=1,
                      torsigma2=1,symmetry='b',pivots1=[6,7],pivots2=[1,6],top1=[7,8],top2=[6,7,8])
    
    @unittest.skipIf(not os.path.isfile(os.path.join(os.path.split(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0])[0], 'external', 'Q2DTor', 'src', 'Q2DTor.py')),
                     "Q2DTor not installed")
    def test_q2dtor_setup(self):
        self.hd2d.readScan()
        self.assertAlmostEquals(self.hd2d.Es[0]/10**9,-594373977.268/10**9,3)
        self.hd2d.getTorsions()
        self.assertEqual(self.hd2d.torsion1, [2, 1, 6, 7])
        self.hd2d.writeInp()
        self.hd2d.writePes()
        self.hd2d.getIcsFile()

    def test_partition_function_calc(self):
        self.hd2d.readEigvals()
        self.assertAlmostEqual(self.hd2d.getPartitionFunction(300.0),3.29752, 4)
        
    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        if os.path.exists(cls.path):
            shutil.rmtree(cls.path) #delete unzipped and created files
            os.remove(os.path.join(os.path.dirname(cls.path),'r0','IOfiles','r0.pes'))
            os.remove(os.path.join(os.path.dirname(cls.path),'r0','r0.out'))
        
class TestHinderedRotorClassicalND(unittest.TestCase):
    """
    Contains unit tests of the StatmechJob class.
    """
    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        freqpath = os.path.join(os.path.dirname(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]),
                                                         'arkane', 'data', 'TolueneFreq.log')
        rotpath = os.path.join(os.path.dirname(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]),
                                                         'arkane', 'data', 'TolueneRot1.log')
        lg = determine_qm_software(freqpath)
        
        conf,unscaled_freqs = lg.loadConformer(symmetry=1,spinMultiplicity=1,
                                                                    opticalIsomers=1,
                                                                    label='Toulene')
        coordinates, number, mass = lg.loadGeometry()
        conf.coordinates = (coordinates, "angstroms")
        conf.number = number
        conf.mass = (mass, "amu")
        
        F = lg.loadForceConstantMatrix()
        
        cls.hdnd= HinderedRotorClassicalND(pivots=[[3,12]],tops=[[12,13,14,15]],sigmas=[6.0],
                                           calcPath=rotpath,conformer=conf,F=F,semiclassical=True)
        
    def test_hindered_rotor_ND(self):
        self.hdnd.readScan()
        self.assertAlmostEqual(self.hdnd.Es[0],20.048316823666962,4)
        self.hdnd.fit()
        self.assertAlmostEqual(self.hdnd.calcPartitionFunction(300.0),2.85254214434672,5)
        
if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
