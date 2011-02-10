#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest

import rmgpy.chem.constants as constants
from rmgpy.chem.thermo import *

################################################################################

class ThermoTest(unittest.TestCase):
    """
    Contains unit tests for the rmgpy.chem.thermo module, used for working with
    thermodynamics models.
    """
    
    def testWilhoit(self):
        """
        Tests the Wilhoit thermodynamics model functions.
        """
        
        # CC(=O)O[O]
        wilhoit = WilhoitModel(cp0=4.0*constants.R, cpInf=21.0*constants.R, a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=500.0, H0=-6.151e+04, S0=-790.2)
        
        Tlist = numpy.arange(200.0, 2001.0, 200.0, numpy.float64)
        Cplist0 = [ 64.398,  94.765, 116.464, 131.392, 141.658, 148.830, 153.948, 157.683, 160.469, 162.589]
        Hlist0 = [-166312., -150244., -128990., -104110., -76742.9, -47652.6, -17347.1, 13834.8, 45663.0, 77978.1]
        Slist0 = [287.421, 341.892, 384.685, 420.369, 450.861, 477.360, 500.708, 521.521, 540.262, 557.284]
        Glist0 = [-223797., -287002., -359801., -440406., -527604., -620485., -718338., -820599., -926809., -1036590.]

        Cplist = wilhoit.getHeatCapacities(Tlist)
        Hlist = wilhoit.getEnthalpies(Tlist)
        Slist = wilhoit.getEntropies(Tlist)
        Glist = wilhoit.getFreeEnergies(Tlist)

        for i in range(len(Tlist)):
            self.assertAlmostEqual(Cplist[i] / Cplist0[i], 1.0, 4)
            self.assertAlmostEqual( Hlist[i] /  Hlist0[i], 1.0, 4)
            self.assertAlmostEqual( Slist[i] /  Slist0[i], 1.0, 4)
            self.assertAlmostEqual( Glist[i] /  Glist0[i], 1.0, 4)

    def testPickleThermoGA(self):
        """
        Test that a ThermoGAModel object can be successfully pickled and
        unpickled with no loss of information.
        """
        Tdata = [300.0,400.0,500.0,600.0,800.0,1000.0,1500.0]
        Cpdata = [3.0,4.0,5.0,6.0,8.0,10.0,15.0]
        thermo0 = ThermoGAModel(Tdata, Cpdata, H298=-2000.0, S298=50.0, dCp=None, dH=100.0, dS=2.0, Tmin=300.0, Tmax=2000.0, comment='This data is completely made up')
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(thermo0))

        self.assertEqual(len(thermo0.Tdata), len(thermo.Tdata))
        self.assertEqual(len(thermo0.Cpdata), len(thermo.Cpdata))
        self.assertEqual(len(thermo0.dCp), len(thermo.dCp))
        for i in range(len(thermo0.Tdata)):
            self.assertEqual(thermo0.Tdata[i], thermo.Tdata[i])
            self.assertEqual(thermo0.Cpdata[i], thermo.Cpdata[i])
            self.assertEqual(thermo0.dCp[i], thermo.dCp[i])
        self.assertEqual(thermo0.H298, thermo.H298)
        self.assertEqual(thermo0.S298, thermo.S298)
        self.assertEqual(thermo0.dH, thermo.dH)
        self.assertEqual(thermo0.dS, thermo.dS)

        self.assertEqual(thermo0.Tmin, thermo.Tmin)
        self.assertEqual(thermo0.Tmax, thermo.Tmax)
        self.assertEqual(thermo0.comment, thermo.comment)

    def testPickleWilhoit(self):
        """
        Test that a WilhoitModel object can be successfully pickled and
        unpickled with no loss of information.
        """
        thermo0 = WilhoitModel(cp0=4.0*constants.R, cpInf=21.0*constants.R, a0=-3.95, a1=9.26, a2=-15.6, a3=8.55, B=500.0, H0=-6.151e+04, S0=-790.2, Tmin=300.0, Tmax=2000.0, comment='CC(=O)O[O]')
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(thermo0))

        self.assertEqual(thermo0.cp0, thermo.cp0)
        self.assertEqual(thermo0.cpInf, thermo.cpInf)
        self.assertEqual(thermo0.a0, thermo.a0)
        self.assertEqual(thermo0.a1, thermo.a1)
        self.assertEqual(thermo0.a2, thermo.a2)
        self.assertEqual(thermo0.a3, thermo.a3)
        self.assertEqual(thermo0.H0, thermo.H0)
        self.assertEqual(thermo0.S0, thermo.S0)
        self.assertEqual(thermo0.B, thermo.B)

        self.assertEqual(thermo0.Tmin, thermo.Tmin)
        self.assertEqual(thermo0.Tmax, thermo.Tmax)
        self.assertEqual(thermo0.comment, thermo.comment)

    def testPickleNASA(self):
        """
        Test that a NASAModel object can be successfully pickled and
        unpickled with no loss of information.
        """

        nasa0 = NASAPolynomial(coeffs=[11.0,12.0,13.0,14.0,15.0,16.0,17.0], Tmin=300.0, Tmax=1000.0, comment='This data is completely made up and unphysical')
        nasa1 = NASAPolynomial(coeffs=[21.0,22.0,23.0,24.0,25.0,26.0,27.0], Tmin=1000.0, Tmax=6000.0, comment='This data is also completely made up and unphysical')

        thermo0 = NASAModel(polynomials=[nasa0, nasa1], Tmin=300.0, Tmax=6000.0, comment='This data is completely made up and unphysical')
        import cPickle
        thermo = cPickle.loads(cPickle.dumps(thermo0))

        self.assertEqual(len(thermo0.polynomials), len(thermo.polynomials))
        for poly0, poly in zip(thermo0.polynomials, thermo.polynomials):
            self.assertEqual(poly0.cm2, poly.cm2)
            self.assertEqual(poly0.cm1, poly.cm1)
            self.assertEqual(poly0.c0, poly.c0)
            self.assertEqual(poly0.c1, poly.c1)
            self.assertEqual(poly0.c2, poly.c2)
            self.assertEqual(poly0.c3, poly.c3)
            self.assertEqual(poly0.c4, poly.c4)
            self.assertEqual(poly0.c5, poly.c5)
            self.assertEqual(poly0.c6, poly.c6)
            self.assertEqual(poly0.Tmin, poly.Tmin)
            self.assertEqual(poly0.Tmax, poly.Tmax)
            self.assertEqual(poly0.comment, poly.comment)

        self.assertEqual(thermo0.Tmin, thermo.Tmin)
        self.assertEqual(thermo0.Tmax, thermo.Tmax)
        self.assertEqual(thermo0.comment, thermo.comment)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
