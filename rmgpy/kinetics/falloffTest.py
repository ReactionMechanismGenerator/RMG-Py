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
This script contains unit tests of the :mod:`rmgpy.kinetics.falloff` module.
"""

import unittest
import numpy

from rmgpy.kinetics.falloff import ThirdBody, Lindemann, Troe
from rmgpy.kinetics.arrhenius import Arrhenius

################################################################################

class TestThirdBody(unittest.TestCase):
    """
    Contains unit tests of the ThirdBody class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        self.efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        self.Tmin = 300.
        self.Tmax = 2000.
        self.Pmin = 0.01
        self.Pmax = 100.
        self.comment = """H + CH3 -> CH4"""
        self.thirdBody = ThirdBody(
            arrheniusLow = self.arrheniusLow,
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            Pmin = (self.Pmin,"bar"),
            Pmax = (self.Pmax,"bar"),
            efficiencies = self.efficiencies,
            comment = self.comment,
        )
        
    def test_arrheniusLow(self):
        """
        Test that the ThirdBody arrhenius property was properly set.
        """
        self.assertTrue(self.thirdBody.arrheniusLow is self.arrheniusLow)
        
    def test_Tmin(self):
        """
        Test that the ThirdBody Tmin property was properly set.
        """
        self.assertAlmostEqual(self.thirdBody.Tmin.value_si, self.Tmin, 6)
        
    def test_Tmax(self):
        """
        Test that the ThirdBody Tmax property was properly set.
        """
        self.assertAlmostEqual(self.thirdBody.Tmax.value_si, self.Tmax, 6)

    def test_Pmin(self):
        """
        Test that the ThirdBody Pmin property was properly set.
        """
        self.assertAlmostEqual(self.thirdBody.Pmin.value_si*1e-5, self.Pmin, 6)
        
    def test_Pmax(self):
        """
        Test that the ThirdBody Pmax property was properly set.
        """
        self.assertAlmostEqual(self.thirdBody.Pmax.value_si*1e-5, self.Pmax, 6)
        
    def test_comment(self):
        """
        Test that the ThirdBody comment property was properly set.
        """
        self.assertEqual(self.thirdBody.comment, self.comment)

    def test_isPressureDependent(self):
        """
        Test the ThirdBody.isPressureDependent() method.
        """
        self.assertTrue(self.thirdBody.isPressureDependent())
    
    def test_getEffectivePressure(self):
        """
        Test the ThirdBody.getEffectivePressure() method.
        """
        P = 1.0
        # Test that each pure bath gas gives the correct effective pressure
        species = self.efficiencies.keys()
        for spec, eff in self.efficiencies.items():
            fractions = numpy.zeros(len(species))
            i = species.index(spec)
            fractions[i] = 1.0
            Peff = self.thirdBody.getEffectivePressure(P, species, fractions)
            self.assertEqual(P * eff, Peff)
        # Also test a mixture of bath gases
        fractions = numpy.zeros(len(species))
        fractions[0] = 0.5
        fractions[1] = 0.5
        eff = 0.5 * self.efficiencies[species[0]] + 0.5 * self.efficiencies[species[1]]
        Peff = self.thirdBody.getEffectivePressure(P, species, fractions)
        self.assertEqual(P * eff, Peff)
            
    def test_getRateCoefficient(self):
        """
        Test the ThirdBody.getRateCoefficient() method.
        """
        Tlist = numpy.array([300,500,1000,1500])
        Plist = numpy.array([1e4,1e5,1e6])
        Kexp = numpy.array([
            [2.83508e+08, 2.83508e+09, 2.83508e+10],
            [7.68759e+07, 7.68759e+08, 7.68759e+09],
            [4.84353e+06, 4.84353e+07, 4.84353e+08],
            [7.05740e+05, 7.05740e+06, 7.05740e+07],
        ])
        for t in range(Tlist.shape[0]):
            for p in range(Plist.shape[0]):
                Kact = self.thirdBody.getRateCoefficient(Tlist[t], Plist[p])
                self.assertAlmostEqual(Kact, Kexp[t,p], delta=1e-4*Kexp[t,p])

    def test_pickle(self):
        """
        Test that a ThirdBody object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thirdBody = cPickle.loads(cPickle.dumps(self.thirdBody,-1))
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.A.value, thirdBody.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.thirdBody.arrheniusLow.A.units, thirdBody.arrheniusLow.A.units)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.n.value, thirdBody.arrheniusLow.n.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.n.units, thirdBody.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.Ea.value, thirdBody.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.Ea.units, thirdBody.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.T0.value, thirdBody.arrheniusLow.T0.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.T0.units, thirdBody.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.thirdBody.Tmin.value, thirdBody.Tmin.value, 4)
        self.assertEqual(self.thirdBody.Tmin.units, thirdBody.Tmin.units)
        self.assertAlmostEqual(self.thirdBody.Tmax.value, thirdBody.Tmax.value, 4)
        self.assertEqual(self.thirdBody.Tmax.units, thirdBody.Tmax.units)
        self.assertAlmostEqual(self.thirdBody.Pmin.value, thirdBody.Pmin.value, 4)
        self.assertEqual(self.thirdBody.Pmin.units, thirdBody.Pmin.units)
        self.assertAlmostEqual(self.thirdBody.Pmax.value, thirdBody.Pmax.value, 4)
        self.assertEqual(self.thirdBody.Pmax.units, thirdBody.Pmax.units)
        self.assertEqual(self.thirdBody.efficiencies, thirdBody.efficiencies)
        self.assertEqual(self.thirdBody.comment, thirdBody.comment)

    def test_repr(self):
        """
        Test that a ThirdBody object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        thirdBody = None
        exec('thirdBody = {0!r}'.format(self.thirdBody))
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.A.value, thirdBody.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.thirdBody.arrheniusLow.A.units, thirdBody.arrheniusLow.A.units)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.n.value, thirdBody.arrheniusLow.n.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.n.units, thirdBody.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.Ea.value, thirdBody.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.Ea.units, thirdBody.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.thirdBody.arrheniusLow.T0.value, thirdBody.arrheniusLow.T0.value, 4)
        self.assertEqual(self.thirdBody.arrheniusLow.T0.units, thirdBody.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.thirdBody.Tmin.value, thirdBody.Tmin.value, 4)
        self.assertEqual(self.thirdBody.Tmin.units, thirdBody.Tmin.units)
        self.assertAlmostEqual(self.thirdBody.Tmax.value, thirdBody.Tmax.value, 4)
        self.assertEqual(self.thirdBody.Tmax.units, thirdBody.Tmax.units)
        self.assertAlmostEqual(self.thirdBody.Pmin.value, thirdBody.Pmin.value, 4)
        self.assertEqual(self.thirdBody.Pmin.units, thirdBody.Pmin.units)
        self.assertAlmostEqual(self.thirdBody.Pmax.value, thirdBody.Pmax.value, 4)
        self.assertEqual(self.thirdBody.Pmax.units, thirdBody.Pmax.units)
        self.assertEqual(self.thirdBody.efficiencies, thirdBody.efficiencies)
        self.assertEqual(self.thirdBody.comment, thirdBody.comment)
        
    def test_changeRate(self):
        """
        Test the ThirdBody.changeRate() method.
        """
        Tlist = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500])
        k0list = numpy.array([self.thirdBody.getRateCoefficient(T,1e5) for T in Tlist])
        self.thirdBody.changeRate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.thirdBody.getRateCoefficient(T,1e5)
            self.assertAlmostEqual(2*kexp, kact, delta=1e-6*kexp)

################################################################################

class TestLindemann(unittest.TestCase):
    """
    Contains unit tests of the Lindemann class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.arrheniusHigh = Arrhenius(
            A = (1.39e+16,"cm^3/(mol*s)"), 
            n = -0.534, 
            Ea = (2.243,"kJ/mol"), 
            T0 = (1,"K"),
        )
        self.arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        self.efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        self.Tmin = 300.
        self.Tmax = 2000.
        self.Pmin = 0.01
        self.Pmax = 100.
        self.comment = """H + CH3 -> CH4"""
        self.lindemann = Lindemann(
            arrheniusHigh = self.arrheniusHigh,
            arrheniusLow = self.arrheniusLow,
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            Pmin = (self.Pmin,"bar"),
            Pmax = (self.Pmax,"bar"),
            efficiencies = self.efficiencies,
            comment = self.comment,
        )
        
    def test_arrheniusHigh(self):
        """
        Test that the Lindemann arrheniusHigh property was properly set.
        """
        self.assertTrue(self.lindemann.arrheniusHigh is self.arrheniusHigh)
        
    def test_arrheniusLow(self):
        """
        Test that the Lindemann arrheniusLow property was properly set.
        """
        self.assertTrue(self.lindemann.arrheniusLow is self.arrheniusLow)
        
    def test_Tmin(self):
        """
        Test that the Lindemann Tmin property was properly set.
        """
        self.assertAlmostEqual(self.lindemann.Tmin.value_si, self.Tmin, 6)
        
    def test_Tmax(self):
        """
        Test that the Lindemann Tmax property was properly set.
        """
        self.assertAlmostEqual(self.lindemann.Tmax.value_si, self.Tmax, 6)

    def test_Pmin(self):
        """
        Test that the Lindemann Pmin property was properly set.
        """
        self.assertAlmostEqual(self.lindemann.Pmin.value_si*1e-5, self.Pmin, 6)
        
    def test_Pmax(self):
        """
        Test that the Lindemann Pmax property was properly set.
        """
        self.assertAlmostEqual(self.lindemann.Pmax.value_si*1e-5, self.Pmax, 6)
        
    def test_comment(self):
        """
        Test that the Lindemann comment property was properly set.
        """
        self.assertEqual(self.lindemann.comment, self.comment)

    def test_isPressureDependent(self):
        """
        Test the Lindemann.isPressureDependent() method.
        """
        self.assertTrue(self.lindemann.isPressureDependent())
    
    def test_getRateCoefficient(self):
        """
        Test the Lindemann.getRateCoefficient() method.
        """
        Tlist = numpy.array([300,500,1000,1500])
        Plist = numpy.array([1e4,1e5,1e6])
        Kexp = numpy.array([
            [1.38023e+08, 2.45661e+08, 2.66439e+08], 
            [6.09146e+07, 2.12349e+08, 2.82604e+08],
            [4.75671e+06, 4.09594e+07, 1.71441e+08],
            [7.03616e+05, 6.85062e+06, 5.42111e+07],
        ])
        for t in range(Tlist.shape[0]):
            for p in range(Plist.shape[0]):
                Kact = self.lindemann.getRateCoefficient(Tlist[t], Plist[p])
                self.assertAlmostEqual(Kact, Kexp[t,p], delta=1e-4*Kexp[t,p])

    def test_pickle(self):
        """
        Test that a Lindemann object can be pickled and unpickled with no loss
        of information.
        """
        import cPickle
        lindemann = cPickle.loads(cPickle.dumps(self.lindemann,-1))
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.A.value, lindemann.arrheniusHigh.A.value, delta=1e0)
        self.assertEqual(self.lindemann.arrheniusHigh.A.units, lindemann.arrheniusHigh.A.units)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.n.value, lindemann.arrheniusHigh.n.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.n.units, lindemann.arrheniusHigh.n.units, 4)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.Ea.value, lindemann.arrheniusHigh.Ea.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.Ea.units, lindemann.arrheniusHigh.Ea.units)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.T0.value, lindemann.arrheniusHigh.T0.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.T0.units, lindemann.arrheniusHigh.T0.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.A.value, lindemann.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.lindemann.arrheniusLow.A.units, lindemann.arrheniusLow.A.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.n.value, lindemann.arrheniusLow.n.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.n.units, lindemann.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.Ea.value, lindemann.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.Ea.units, lindemann.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.T0.value, lindemann.arrheniusLow.T0.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.T0.units, lindemann.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.lindemann.Tmin.value, lindemann.Tmin.value, 4)
        self.assertEqual(self.lindemann.Tmin.units, lindemann.Tmin.units)
        self.assertAlmostEqual(self.lindemann.Tmax.value, lindemann.Tmax.value, 4)
        self.assertEqual(self.lindemann.Tmax.units, lindemann.Tmax.units)
        self.assertAlmostEqual(self.lindemann.Pmin.value, lindemann.Pmin.value, 4)
        self.assertEqual(self.lindemann.Pmin.units, lindemann.Pmin.units)
        self.assertAlmostEqual(self.lindemann.Pmax.value, lindemann.Pmax.value, 4)
        self.assertEqual(self.lindemann.Pmax.units, lindemann.Pmax.units)
        self.assertEqual(self.lindemann.efficiencies, lindemann.efficiencies)
        self.assertEqual(self.lindemann.comment, lindemann.comment)

    def test_repr(self):
        """
        Test that a Lindemann object can be reconstructed from its repr()
        output with no loss of information.
        """
        lindemann = None
        exec('lindemann = {0!r}'.format(self.lindemann))
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.A.value, lindemann.arrheniusHigh.A.value, delta=1e0)
        self.assertEqual(self.lindemann.arrheniusHigh.A.units, lindemann.arrheniusHigh.A.units)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.n.value, lindemann.arrheniusHigh.n.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.n.units, lindemann.arrheniusHigh.n.units, 4)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.Ea.value, lindemann.arrheniusHigh.Ea.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.Ea.units, lindemann.arrheniusHigh.Ea.units)
        self.assertAlmostEqual(self.lindemann.arrheniusHigh.T0.value, lindemann.arrheniusHigh.T0.value, 4)
        self.assertEqual(self.lindemann.arrheniusHigh.T0.units, lindemann.arrheniusHigh.T0.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.A.value, lindemann.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.lindemann.arrheniusLow.A.units, lindemann.arrheniusLow.A.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.n.value, lindemann.arrheniusLow.n.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.n.units, lindemann.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.Ea.value, lindemann.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.Ea.units, lindemann.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.lindemann.arrheniusLow.T0.value, lindemann.arrheniusLow.T0.value, 4)
        self.assertEqual(self.lindemann.arrheniusLow.T0.units, lindemann.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.lindemann.Tmin.value, lindemann.Tmin.value, 4)
        self.assertEqual(self.lindemann.Tmin.units, lindemann.Tmin.units)
        self.assertAlmostEqual(self.lindemann.Tmax.value, lindemann.Tmax.value, 4)
        self.assertEqual(self.lindemann.Tmax.units, lindemann.Tmax.units)
        self.assertAlmostEqual(self.lindemann.Pmin.value, lindemann.Pmin.value, 4)
        self.assertEqual(self.lindemann.Pmin.units, lindemann.Pmin.units)
        self.assertAlmostEqual(self.lindemann.Pmax.value, lindemann.Pmax.value, 4)
        self.assertEqual(self.lindemann.Pmax.units, lindemann.Pmax.units)
        self.assertEqual(self.lindemann.efficiencies, lindemann.efficiencies)
        self.assertEqual(self.lindemann.comment, lindemann.comment)
        
    def test_changeRate(self):
        """
        Test the Lindemann.changeRate() method.
        """
        Tlist = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500])
        k0list = numpy.array([self.lindemann.getRateCoefficient(T,1e5) for T in Tlist])
        self.lindemann.changeRate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.lindemann.getRateCoefficient(T,1e5)
            self.assertAlmostEqual(2*kexp, kact, delta=1e-6*kexp)

################################################################################

class TestTroe(unittest.TestCase):
    """
    Contains unit tests of the Troe class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.arrheniusHigh = Arrhenius(
            A = (1.39e+16,"cm^3/(mol*s)"), 
            n = -0.534, 
            Ea = (2.243,"kJ/mol"), 
            T0 = (1,"K"),
        )
        self.arrheniusLow = Arrhenius(
            A = (2.62e+33,"cm^6/(mol^2*s)"), 
            n = -4.76, 
            Ea = (10.21,"kJ/mol"), 
            T0 = (1,"K"),
        )
        self.alpha = 0.783
        self.T3 = 74
        self.T1 = 2941
        self.T2 = 6964
        self.efficiencies = {"C": 3, "C(=O)=O": 2, "CC": 3, "O": 6, "[Ar]": 0.7, "[C]=O": 1.5, "[H][H]": 2}
        self.Tmin = 300.
        self.Tmax = 2000.
        self.Pmin = 0.01
        self.Pmax = 100.
        self.comment = """H + CH3 -> CH4"""
        self.troe = Troe(
            arrheniusHigh = self.arrheniusHigh,
            arrheniusLow = self.arrheniusLow,
            alpha = self.alpha,
            T3 = (self.T3,"K"),
            T1 = (self.T1,"K"),
            T2 = (self.T2,"K"),
            Tmin = (self.Tmin,"K"),
            Tmax = (self.Tmax,"K"),
            Pmin = (self.Pmin,"bar"),
            Pmax = (self.Pmax,"bar"),
            efficiencies = self.efficiencies,
            comment = self.comment,
        )
        
    def test_arrheniusHigh(self):
        """
        Test that the Troe arrheniusHigh property was properly set.
        """
        self.assertTrue(self.troe.arrheniusHigh is self.arrheniusHigh)
        
    def test_arrheniusLow(self):
        """
        Test that the Troe arrheniusLow property was properly set.
        """
        self.assertTrue(self.troe.arrheniusLow is self.arrheniusLow)
        
    def test_alpha(self):
        """
        Test that the Troe alpha property was properly set.
        """
        self.assertEqual(self.troe.alpha, self.alpha)
        
    def test_T3(self):
        """
        Test that the Troe T3 property was properly set.
        """
        self.assertAlmostEqual(self.troe.T3.value_si, self.T3, 6)
        
    def test_T1(self):
        """
        Test that the Troe T1 property was properly set.
        """
        self.assertAlmostEqual(self.troe.T1.value_si, self.T1, 6)
        
    def test_T2(self):
        """
        Test that the Troe T2 property was properly set.
        """
        self.assertAlmostEqual(self.troe.T2.value_si, self.T2, 6)
        
    def test_Tmin(self):
        """
        Test that the Troe Tmin property was properly set.
        """
        self.assertAlmostEqual(self.troe.Tmin.value_si, self.Tmin, 6)
        
    def test_Tmax(self):
        """
        Test that the Troe Tmax property was properly set.
        """
        self.assertAlmostEqual(self.troe.Tmax.value_si, self.Tmax, 6)

    def test_Pmin(self):
        """
        Test that the Troe Pmin property was properly set.
        """
        self.assertAlmostEqual(self.troe.Pmin.value_si*1e-5, self.Pmin, 6)
        
    def test_Pmax(self):
        """
        Test that the Troe Pmax property was properly set.
        """
        self.assertAlmostEqual(self.troe.Pmax.value_si*1e-5, self.Pmax, 6)
        
    def test_comment(self):
        """
        Test that the Troe comment property was properly set.
        """
        self.assertEqual(self.troe.comment, self.comment)

    def test_isPressureDependent(self):
        """
        Test the Troe.isPressureDependent() method.
        """
        self.assertTrue(self.troe.isPressureDependent())
    
    def test_getRateCoefficient(self):
        """
        Test the Troe.getRateCoefficient() method.
        """
        Tlist = numpy.array([300,500,1000,1500])
        Plist = numpy.array([1e4,1e5,1e6])
        Kexp = numpy.array([
            [1.00866e+08, 2.03759e+08, 2.55190e+08],
            [4.74623e+07, 1.41629e+08, 2.47597e+08],
            [3.97397e+06, 2.89521e+07, 9.57569e+07],
            [5.91277e+05, 5.14013e+06, 3.12239e+07],
        ])
        for t in range(Tlist.shape[0]):
            for p in range(Plist.shape[0]):
                Kact = self.troe.getRateCoefficient(Tlist[t], Plist[p])
                self.assertAlmostEqual(Kact, Kexp[t,p], delta=1e-4*Kexp[t,p])

    def test_pickle(self):
        """
        Test that a Troe object can be pickled and unpickled with no loss of
        information.
        """
        import cPickle
        troe = cPickle.loads(cPickle.dumps(self.troe,-1))
        self.assertAlmostEqual(self.troe.arrheniusHigh.A.value, troe.arrheniusHigh.A.value, delta=1e0)
        self.assertEqual(self.troe.arrheniusHigh.A.units, troe.arrheniusHigh.A.units)
        self.assertAlmostEqual(self.troe.arrheniusHigh.n.value, troe.arrheniusHigh.n.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.n.units, troe.arrheniusHigh.n.units, 4)
        self.assertAlmostEqual(self.troe.arrheniusHigh.Ea.value, troe.arrheniusHigh.Ea.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.Ea.units, troe.arrheniusHigh.Ea.units)
        self.assertAlmostEqual(self.troe.arrheniusHigh.T0.value, troe.arrheniusHigh.T0.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.T0.units, troe.arrheniusHigh.T0.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.A.value, troe.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.troe.arrheniusLow.A.units, troe.arrheniusLow.A.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.n.value, troe.arrheniusLow.n.value, 4)
        self.assertEqual(self.troe.arrheniusLow.n.units, troe.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.troe.arrheniusLow.Ea.value, troe.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.troe.arrheniusLow.Ea.units, troe.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.T0.value, troe.arrheniusLow.T0.value, 4)
        self.assertEqual(self.troe.arrheniusLow.T0.units, troe.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.troe.alpha, troe.alpha, 6)
        self.assertAlmostEqual(self.troe.T3.value, troe.T3.value, 6)
        self.assertEqual(self.troe.T3.units, troe.T3.units)
        self.assertAlmostEqual(self.troe.T1.value, troe.T1.value, 6)
        self.assertEqual(self.troe.T1.units, troe.T1.units)
        self.assertAlmostEqual(self.troe.T2.value, troe.T2.value, 6)
        self.assertEqual(self.troe.T2.units, troe.T2.units)
        self.assertAlmostEqual(self.troe.Tmin.value, troe.Tmin.value, 4)
        self.assertEqual(self.troe.Tmin.units, troe.Tmin.units)
        self.assertAlmostEqual(self.troe.Tmax.value, troe.Tmax.value, 4)
        self.assertEqual(self.troe.Tmax.units, troe.Tmax.units)
        self.assertAlmostEqual(self.troe.Pmin.value, troe.Pmin.value, 4)
        self.assertEqual(self.troe.Pmin.units, troe.Pmin.units)
        self.assertAlmostEqual(self.troe.Pmax.value, troe.Pmax.value, 4)
        self.assertEqual(self.troe.Pmax.units, troe.Pmax.units)
        self.assertEqual(self.troe.efficiencies, troe.efficiencies)
        self.assertEqual(self.troe.comment, troe.comment)

    def test_repr(self):
        """
        Test that a Troe object can be reconstructed from its repr() output
        with no loss of information.
        """
        troe = None
        exec('troe = {0!r}'.format(self.troe))
        self.assertAlmostEqual(self.troe.arrheniusHigh.A.value, troe.arrheniusHigh.A.value, delta=1e0)
        self.assertEqual(self.troe.arrheniusHigh.A.units, troe.arrheniusHigh.A.units)
        self.assertAlmostEqual(self.troe.arrheniusHigh.n.value, troe.arrheniusHigh.n.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.n.units, troe.arrheniusHigh.n.units, 4)
        self.assertAlmostEqual(self.troe.arrheniusHigh.Ea.value, troe.arrheniusHigh.Ea.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.Ea.units, troe.arrheniusHigh.Ea.units)
        self.assertAlmostEqual(self.troe.arrheniusHigh.T0.value, troe.arrheniusHigh.T0.value, 4)
        self.assertEqual(self.troe.arrheniusHigh.T0.units, troe.arrheniusHigh.T0.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.A.value, troe.arrheniusLow.A.value, delta=1e0)
        self.assertEqual(self.troe.arrheniusLow.A.units, troe.arrheniusLow.A.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.n.value, troe.arrheniusLow.n.value, 4)
        self.assertEqual(self.troe.arrheniusLow.n.units, troe.arrheniusLow.n.units, 4)
        self.assertAlmostEqual(self.troe.arrheniusLow.Ea.value, troe.arrheniusLow.Ea.value, 4)
        self.assertEqual(self.troe.arrheniusLow.Ea.units, troe.arrheniusLow.Ea.units)
        self.assertAlmostEqual(self.troe.arrheniusLow.T0.value, troe.arrheniusLow.T0.value, 4)
        self.assertEqual(self.troe.arrheniusLow.T0.units, troe.arrheniusLow.T0.units)
        self.assertAlmostEqual(self.troe.alpha, troe.alpha, 6)
        self.assertAlmostEqual(self.troe.T3.value, troe.T3.value, 6)
        self.assertEqual(self.troe.T3.units, troe.T3.units)
        self.assertAlmostEqual(self.troe.T1.value, troe.T1.value, 6)
        self.assertEqual(self.troe.T1.units, troe.T1.units)
        self.assertAlmostEqual(self.troe.T2.value, troe.T2.value, 6)
        self.assertEqual(self.troe.T2.units, troe.T2.units)
        self.assertAlmostEqual(self.troe.Tmin.value, troe.Tmin.value, 4)
        self.assertEqual(self.troe.Tmin.units, troe.Tmin.units)
        self.assertAlmostEqual(self.troe.Tmax.value, troe.Tmax.value, 4)
        self.assertEqual(self.troe.Tmax.units, troe.Tmax.units)
        self.assertAlmostEqual(self.troe.Pmin.value, troe.Pmin.value, 4)
        self.assertEqual(self.troe.Pmin.units, troe.Pmin.units)
        self.assertAlmostEqual(self.troe.Pmax.value, troe.Pmax.value, 4)
        self.assertEqual(self.troe.Pmax.units, troe.Pmax.units)
        self.assertEqual(self.troe.efficiencies, troe.efficiencies)
        self.assertEqual(self.troe.comment, troe.comment)
        
    def test_changeRate(self):
        """
        Test the Troe.changeRate() method.
        """
        Tlist = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500])
        k0list = numpy.array([self.troe.getRateCoefficient(T,1e5) for T in Tlist])
        self.troe.changeRate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.troe.getRateCoefficient(T,1e5)
            self.assertAlmostEqual(2*kexp, kact, delta=1e-6*kexp)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
