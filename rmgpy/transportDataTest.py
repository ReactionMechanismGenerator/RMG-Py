#!/usr/bin/env python
# encoding: utf-8

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
This script contains unit test of the :mod: 'rmgpy.transport' module and :mod: 'rmgpy.data.transport' module
"""

import unittest
import numpy

from rmgpy.data.transport import CriticalPointGroupContribution
#################################################################################

class TestTransportData(unittest.TestCase):
    """
    Contains unit test of the :class: 'transportData' class
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.shapeIndex = 1
        self.epsilon = 2.104
        self.sigma = 3.402
        self.dipoleMoment = 1.000
        self.polarizability = 0.134
        self.rotrelaxcollnum = 0.000
        self.comment = 'test'
        
        self.transport = transportData(
            shapeIndex = self.shapeIndex,
            epsilon = self.epsilon,
            sigma = self.sigma,
            dipoleMoment = self.dipoleMoment,
            polarizability = self.polarizability,
            rotrelaxcollnum = self.rotrelaxcollnum,
            comment = self.comment,
        )
        
    def test_shapeIndex(self):
        """
        Test that the TransportData shapeIndex property was properly set.
        """
        self.assertAlmostEqual(self.transport.shapeIndex.value_si, self.shapeIndex, 6)
    
    def test_epsilon(self):
        """
        Test that the TransportData epsilon property was properly set.
        """
        self.assertAlmostEqual(self.transport.epsilon.value_si, self.epsilon, 6)
        
    def test_sigma(self):
        """
        Test that the TransportData sigma property was properly set.
        """
        self.assertAlmostEqual(self.transport.sigma.value_si, self.sigma, 6)
        
    def test_dipoleMoment(self):
        """
        Test that the TransportData dipoleMoment property was properly set.
        """
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, self.dipoleMoment, 6)    
    
    def test_polarizability(self):
        """
        Test that the TransportData polarizability property was properly set.
        """
        self.assertAlmostEqual(self.transport.polarizability.value_si, self.polarizability, 6)
    
    def test_rotrelaxcollnum(self):
        """
        Test that the TransportData rotrelaxcollnum property was properly set.
        """
        self.assertAlmostEqual(self.transport.rotrelaxcollnum.value_si, self.rotrelaxcollnum, 6)
    
    def test_comment(self):
        """
        Test that the TransportData comment property was properly set.
        """
        self.assertEqual(self.transport.comment, self.comment)
    
    def test_pickle(self):
        """
        Test that a TransportData object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        transport = cPickle.loads(cPickle.dumps(self.transport))
        self.assertAlmostEqual(self.transport.shapeIndex.value, transport.shapeIndex.value, 4)
        self.assertAlmostEqual(self.transport.epsilon.value, transport.epsilon.value, 4)
        self.assertAlmostEqual(self.transport.sigma.value, transport.sigma.value, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value, transport.dipoleMoment.value, 4)
        self.assertAlmostEqual(self.transport.polarizability.value, transport.polarizability.value, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum.value, transport.rotrelaxcollnum.value, 4)
        self.assertEqal(self.transport.comment, transport.comment)
    
    def test_repr(self):
        """
        Test that a TransportData object can be reconstructed from its repr() output with no loss of information
        """
        transport = None
        exec('transport = {0!r}'.format(self.transport))
        self.assertAlmostEqual(self.transport.shapeIndex.value, transport.shapeIndex.value, 4)
        self.assertAlmostEqual(self.transport.epsilon.value, transport.epsilon.value, 4)
        self.assertAlmostEqual(self.transport.sigma.value, transport.sigma.value, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value, transport.dipoleMoment.value, 4)
        self.assertAlmostEqual(self.transport.polarizability.value, transport.polarizability.value, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum.value, transport.rotrelaxcollnum.value, 4)
        self.assertEqal(self.transport.comment, transport.comment)
        
class TestCriticalPointGroupContribution(unittest.TestCase):
    """
    Contains unit test of the :class: 'criticalPointGroupContribution' class
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tc = 0.0141
        self.Pc = -.0012
        self.Vc = 65
        self.Tb = 23.58
        self.structureIndex = 1
        
        self.criticalPointContribution = CriticalPointGroupContribution(
            Tc = self.Tc,
            Pc = self.Pc,
            Vc = self.Vc,
            Tb = self.Tb,
            structureIndex = self.structureIndex,
        )
        
    def test_Tc(self):
        """
        Test that the CriticalPointGroupContribution Tc property was properly set.
        """
        self.assertAlmostEqual(self.transport.Tc.value_si, self.Tc, 6)
    
    def test_Pc(self):
        """
        Test that the CriticalPointGroupContribution Pc property was properly set.
        """
        self.assertAlmostEqual(self.transport.Pc.value_si, self.Pc, 6)
        
    def test_Vc(self):
        """
        Test that the CriticalPointGroupContribution Vc property was properly set.
        """
        self.assertAlmostEqual(self.transport.Vc.value_si, self.Vc, 6)
        
    def test_Tb(self):
        """
        Test that the CriticalPointGroupContribution Tb property was properly set.
        """
        self.assertAlmostEqual(self.transport.Tb.value_si, self.Tb, 6)    
    
    def test_structureIndex(self):
        """
        Test that the CriticalPointGroupContribution structureIndex property was properly set.
        """
        self.assertAlmostEqual(self.transport.structureIndex.value_si, self.structureIndex, 6)
    
    def test_pickle(self):
        """
        Test that a CriticalPointGroupContribution object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        criticalPointContribution = cPickle.loads(cPickle.dumps(self.criticalPointContribution))
        self.assertAlmostEqual(self.transport.Tc.value, transport.Tc.value, 4)
        self.assertAlmostEqual(self.transport.Pc.value, transport.Pc.value, 4)
        self.assertAlmostEqual(self.transport.Vc.value, transport.Vc.value, 4)
        self.assertAlmostEqual(self.transport.Tb.value, transport.Tb.value, 4)
        self.assertAlmostEqual(self.transport.structureIndex.value, transport.structureIndex.value, 4)
    
    def test_repr(self):
        """
        Test that a CriticalPointGroupContribution object can be reconstructed from its repr() output with no loss of information
        """
        criticalPointContribution = None
        exec('criticalPointGroupContribution = {0!r}'.format(self.criticalPointContribution))
        self.assertAlmostEqual(self.transport.Tc.value, transport.Tc.value, 4)
        self.assertAlmostEqual(self.transport.Pc.value, transport.Pc.value, 4)
        self.assertAlmostEqual(self.transport.Vc.value, transport.Vc.value, 4)
        self.assertAlmostEqual(self.transport.Tb.value, transport.Tb.value, 4)
        self.assertAlmostEqual(self.transport.structureIndex.value, transport.structureIndex.value, 4)
    
class TestTransportDatabase(unittest.TestCase):
    """
    Contains unit test of the :class: 'TransportDatabase' class
    """
    
    def setUp(self):
        """
        a function run before each unit test in this class
        """ 
        self.libraries = ['GRI-Mech', 'PrimaryTransportLibrary']
        self.groups = ['ring', 'nonring']
        self.libraryOrder = []
        
        transportdb = TransportDatabase(
            libraries = self.libraries,
            groups = self.groups,
            libraryOrder = self.libraryOrder)
        