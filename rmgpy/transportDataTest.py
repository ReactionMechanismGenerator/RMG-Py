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
import os
from rmgpy.species import Species
from rmgpy.molecule.molecule import Molecule
from rmgpy.quantity import DipoleMoment, Temperature, Length, Polarizability
from rmgpy.transport import TransportData
from rmgpy.data.transport import CriticalPointGroupContribution, TransportDatabase
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
        self.dipoleMoment = 1.000,
        self.polarizability = 0.134,
        self.rotrelaxcollnum = 0.000
        self.comment = 'test'
        
        self.transport = TransportData(
            shapeIndex = self.shapeIndex,
            epsilon = (self.epsilon, 'K'),
            sigma = (self.sigma, 'angstroms'),
            dipoleMoment = (self.dipoleMoment, 'C*m'),
            polarizability = (self.polarizability, 'C*m^2*V^-1'),
            rotrelaxcollnum = self.rotrelaxcollnum,
            comment = self.comment,
        )
        
    def test_shapeIndex(self):
        """
        Test that the TransportData shapeIndex property was properly set.
        """
        self.assertAlmostEqual(self.transport.shapeIndex, self.shapeIndex, 6)
    
    def test_epsilon(self):
        """
        Test that the TransportData epsilon property was properly set.
        """
        self.assertAlmostEqual(self.transport.epsilon.value_si, self.epsilon, 6)
        
    def test_sigma(self):
        """
        Test that the TransportData sigma property was properly set.
        """
        self.assertAlmostEqual(self.transport.sigma.value_si*1e10, self.sigma, 6)
        
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
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, self.rotrelaxcollnum, 6)
    
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
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)
    
    def test_repr(self):
        """
        Test that a TransportData object can be reconstructed from its repr() output with no loss of information
        """
        transport = None
        exec('transport = {0!r}'.format(self.transport))
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)
        
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
        self.assertAlmostEqual(self.criticalPointContribution.Tc, self.Tc, 6)
    
    def test_Pc(self):
        """
        Test that the CriticalPointGroupContribution Pc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Pc, self.Pc, 6)
        
    def test_Vc(self):
        """
        Test that the CriticalPointGroupContribution Vc property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Vc, self.Vc, 6)
        
    def test_Tb(self):
        """
        Test that the CriticalPointGroupContribution Tb property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.Tb, self.Tb, 6)    
    
    def test_structureIndex(self):
        """
        Test that the CriticalPointGroupContribution structureIndex property was properly set.
        """
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, self.structureIndex, 6)
    
    def test_pickle(self):
        """
        Test that a CriticalPointGroupContribution object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        criticalPointContribution = cPickle.loads(cPickle.dumps(self.criticalPointContribution))
        self.assertAlmostEqual(self.criticalPointContribution.Tc, criticalPointContribution.Tc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Pc, criticalPointContribution.Pc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Vc, criticalPointContribution.Vc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Tb, criticalPointContribution.Tb, 4)
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, criticalPointContribution.structureIndex, 4)
    
    def test_repr(self):
        """
        Test that a CriticalPointGroupContribution object can be reconstructed from its repr() output with no loss of information
        """
        criticalPointContribution = None
        exec('criticalPointContribution = {0!r}'.format(self.criticalPointContribution))
        self.assertAlmostEqual(self.criticalPointContribution.Tc, criticalPointContribution.Tc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Pc, criticalPointContribution.Pc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Vc, criticalPointContribution.Vc, 4)
        self.assertAlmostEqual(self.criticalPointContribution.Tb, criticalPointContribution.Tb, 4)
        self.assertAlmostEqual(self.criticalPointContribution.structureIndex, criticalPointContribution.structureIndex, 4)
    
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
        
        path = os.path.normpath(os.path.join( os.path.dirname(os.path.abspath(__file__)), '../../../RMG Database/RMG-database/input/transport'))
        self.transportdb = TransportDatabase()
        self.transportdb.load(path,self.libraries)
    def testJoback(self):    
        
        #values calculate from joback's estimations
        self.testCases = [
            ['acetone', 'CC(=O)C', 2.536255455e-10, 5.43244513e-21]
            ]
        
        for name, smiles, sigma, epsilon in self.testCases:
            species = Species(molecule=[Molecule(SMILES=smiles)])
            transportData = self.transportdb.getTransportPropertiesViaGroupEstimates(species)
            print name, transportData
            print self.assertAlmostEqual(transportData.sigma.value, sigma*1e10,4)
            print self.assertAlmostEqual(transportData.epsilon.value, epsilon)
            