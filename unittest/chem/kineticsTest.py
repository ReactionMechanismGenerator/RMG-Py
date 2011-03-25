#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import unittest

from rmgpy.chem.kinetics import *

################################################################################

class KineticsTest(unittest.TestCase):
    """
    Contains unit tests for the rmgpy.chem.kinetics module, used for working with
    kinetics models.
    """

    def testPickleArrhenius(self):
        """
        Test that an Arrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """
        
        kinetics0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        self.assertEqual(kinetics0.A.value, kinetics.A.value)
        self.assertEqual(kinetics0.n.value, kinetics.n.value)
        self.assertEqual(kinetics0.T0.value, kinetics.T0.value)
        self.assertEqual(kinetics0.Ea.value, kinetics.Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleArrheniusEP(self):
        """
        Test that an ArrheniusEP object can be successfully pickled and
        unpickled with no loss of information.
        """

        kinetics0 = ArrheniusEP(A=1.0e6, n=1.0, alpha=0.5, E0=10000.0, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        self.assertAlmostEqual(kinetics0.A.value, kinetics.A.value, 4)
        self.assertAlmostEqual(kinetics0.n.value, kinetics.n.value, 4)
        self.assertAlmostEqual(kinetics0.alpha.value, kinetics.alpha.value, 4)
        self.assertAlmostEqual(kinetics0.E0.value, kinetics.E0.value, 4)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleMultiArrhenius(self):
        """
        Test that a MultiArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """

        arrh0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrh1 = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')

        kinetics0 = MultiArrhenius(arrheniusList=[arrh0, arrh1], Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        Narrh = 2
        self.assertEqual(len(kinetics0.arrheniusList), Narrh)
        self.assertEqual(len(kinetics.arrheniusList), Narrh)
        for i in range(Narrh):
            self.assertEqual(kinetics0.arrheniusList[i].A.value, kinetics.arrheniusList[i].A.value)
            self.assertEqual(kinetics0.arrheniusList[i].n.value, kinetics.arrheniusList[i].n.value)
            self.assertEqual(kinetics0.arrheniusList[i].T0.value, kinetics.arrheniusList[i].T0.value)
            self.assertEqual(kinetics0.arrheniusList[i].Ea.value, kinetics.arrheniusList[i].Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPicklePDepArrhenius(self):
        """
        Test that a PDepArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """

        P0 = 1e3; P1 = 1e5
        arrh0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrh1 = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')

        kinetics0 = PDepArrhenius(pressures=[P0,P1], arrhenius=[arrh0, arrh1], Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        Narrh = 2
        self.assertEqual(len(kinetics0.pressures.values), Narrh)
        self.assertEqual(len(kinetics.pressures.values), Narrh)
        self.assertEqual(len(kinetics0.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(kinetics0.pressures.values[i], kinetics.pressures.values[i])
            self.assertEqual(kinetics0.arrhenius[i].A.value, kinetics.arrhenius[i].A.value)
            self.assertEqual(kinetics0.arrhenius[i].n.value, kinetics.arrhenius[i].n.value)
            self.assertEqual(kinetics0.arrhenius[i].T0.value, kinetics.arrhenius[i].T0.value)
            self.assertEqual(kinetics0.arrhenius[i].Ea.value, kinetics.arrhenius[i].Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleChebyshev(self):
        """
        Test that a Chebyshev object can be successfully pickled and
        unpickled with no loss of information.
        """

        coeffs = numpy.array([[1.0,2.0,3.0,4.0], [5.0,6.0,7.0,8.0], [9.0,10.0,11.0,12.0]], numpy.float64)

        kinetics0 = Chebyshev(coeffs=coeffs, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        degreeT = 3; degreeP = 4
        self.assertEqual(kinetics0.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics0.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics0.degreeT, degreeT)
        self.assertEqual(kinetics0.degreeP, degreeP)
        self.assertEqual(kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics.degreeT, degreeT)
        self.assertEqual(kinetics.degreeP, degreeP)
        for i in range(degreeT):
            for j in range(degreeP):
                self.assertEqual(kinetics0.coeffs[i,j], kinetics.coeffs[i,j])
        
        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleThirdBody(self):
        """
        Test that a ThirdBody object can be successfully pickled and
        unpickled with no loss of information.
        """

        arrhHigh = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = ThirdBody(arrheniusHigh=arrhHigh, efficiencies=efficiencies, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleLindemann(self):
        """
        Test that a Lindemann object can be successfully pickled and
        unpickled with no loss of information.
        """

        arrhLow = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrhHigh = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = Lindemann(arrheniusLow=arrhLow, arrheniusHigh=arrhHigh, efficiencies=efficiencies, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        self.assertEqual(kinetics0.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(kinetics0.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(kinetics0.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(kinetics0.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testPickleTroe(self):
        """
        Test that a Troe object can be successfully pickled and
        unpickled with no loss of information.
        """

        arrhLow = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrhHigh = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = Troe(arrheniusLow=arrhLow, arrheniusHigh=arrhHigh, efficiencies=efficiencies, alpha=0.6, T3=1000.0, T1=500.0, T2=300.0, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(kinetics0))

        self.assertEqual(kinetics0.alpha.value, kinetics.alpha.value)
        self.assertEqual(kinetics0.T3.value, kinetics.T3.value)
        self.assertEqual(kinetics0.T1.value, kinetics.T1.value)
        self.assertEqual(kinetics0.T2.value, kinetics.T2.value)
        self.assertEqual(kinetics0.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(kinetics0.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(kinetics0.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(kinetics0.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputArrhenius(self):
        """
        Test that an Arrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        kinetics0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        exec('kinetics = %r' % kinetics0)

        self.assertEqual(kinetics0.A.value, kinetics.A.value)
        self.assertEqual(kinetics0.n.value, kinetics.n.value)
        self.assertEqual(kinetics0.T0.value, kinetics.T0.value)
        self.assertEqual(kinetics0.Ea.value, kinetics.Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputArrheniusEP(self):
        """
        Test that an ArrheniusEP object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        kinetics0 = ArrheniusEP(A=1.0e6, n=1.0, alpha=0.5, E0=10000.0, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        exec('kinetics = %r' % kinetics0)

        self.assertAlmostEqual(kinetics0.A.value, kinetics.A.value, 4)
        self.assertAlmostEqual(kinetics0.n.value, kinetics.n.value, 4)
        self.assertAlmostEqual(kinetics0.alpha.value, kinetics.alpha.value, 4)
        self.assertAlmostEqual(kinetics0.E0.value, kinetics.E0.value, 4)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputMultiArrhenius(self):
        """
        Test that a MultiArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        arrh0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrh1 = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')

        kinetics0 = MultiArrhenius(arrheniusList=[arrh0, arrh1], Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        exec('kinetics = %r' % kinetics0)

        Narrh = 2
        self.assertEqual(len(kinetics0.arrheniusList), Narrh)
        self.assertEqual(len(kinetics.arrheniusList), Narrh)
        for i in range(Narrh):
            self.assertEqual(kinetics0.arrheniusList[i].A.value, kinetics.arrheniusList[i].A.value)
            self.assertEqual(kinetics0.arrheniusList[i].n.value, kinetics.arrheniusList[i].n.value)
            self.assertEqual(kinetics0.arrheniusList[i].T0.value, kinetics.arrheniusList[i].T0.value)
            self.assertEqual(kinetics0.arrheniusList[i].Ea.value, kinetics.arrheniusList[i].Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputPDepArrhenius(self):
        """
        Test that a PDepArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        P0 = 1e3; P1 = 1e5
        arrh0 = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrh1 = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')

        kinetics0 = PDepArrhenius(pressures=[P0,P1], arrhenius=[arrh0, arrh1], Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up')
        exec('kinetics = %r' % kinetics0)

        Narrh = 2
        self.assertEqual(len(kinetics0.pressures.values), Narrh)
        self.assertEqual(len(kinetics.pressures.values), Narrh)
        self.assertEqual(len(kinetics0.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(kinetics0.pressures.values[i], kinetics.pressures.values[i])
            self.assertEqual(kinetics0.arrhenius[i].A.value, kinetics.arrhenius[i].A.value)
            self.assertEqual(kinetics0.arrhenius[i].n.value, kinetics.arrhenius[i].n.value)
            self.assertEqual(kinetics0.arrhenius[i].T0.value, kinetics.arrhenius[i].T0.value)
            self.assertEqual(kinetics0.arrhenius[i].Ea.value, kinetics.arrhenius[i].Ea.value)

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputChebyshev(self):
        """
        Test that a Chebyshev object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        coeffs = numpy.array([[1.0,2.0,3.0,4.0], [5.0,6.0,7.0,8.0], [9.0,10.0,11.0,12.0]], numpy.float64)

        kinetics0 = Chebyshev(coeffs=coeffs, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        exec('kinetics = %r' % kinetics0)

        degreeT = 3; degreeP = 4
        self.assertEqual(kinetics0.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics0.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics0.degreeT, degreeT)
        self.assertEqual(kinetics0.degreeP, degreeP)
        self.assertEqual(kinetics.coeffs.shape[0], degreeT)
        self.assertEqual(kinetics.coeffs.shape[1], degreeP)
        self.assertEqual(kinetics.degreeT, degreeT)
        self.assertEqual(kinetics.degreeP, degreeP)
        for i in range(degreeT):
            for j in range(degreeP):
                self.assertEqual(kinetics0.coeffs[i,j], kinetics.coeffs[i,j])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputThirdBody(self):
        """
        Test that a ThirdBody object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        arrhHigh = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = ThirdBody(arrheniusHigh=arrhHigh, efficiencies=efficiencies, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        exec('kinetics = %r' % kinetics0)

        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputLindemann(self):
        """
        Test that a Lindemann object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        arrhLow = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrhHigh = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = Lindemann(arrheniusLow=arrhLow, arrheniusHigh=arrhHigh, efficiencies=efficiencies, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        exec('kinetics = %r' % kinetics0)

        self.assertEqual(kinetics0.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(kinetics0.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(kinetics0.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(kinetics0.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

    def testOutputTroe(self):
        """
        Test that a Troe object can be successfully reconstructed
        from its repr() output with no loss of information.
        """

        arrhLow = Arrhenius(A=1.0e6, n=1.0, Ea=10000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are completely made up')
        arrhHigh = Arrhenius(A=1.0e12, n=0.0, Ea=20000.0, T0=298.15, Tmin=300, Tmax=2000, numReactants=2, comment='These parameters are also completely made up')
        efficiencies = {'N2': 0.5, 'Ar': 1.5}

        kinetics0 = Troe(arrheniusLow=arrhLow, arrheniusHigh=arrhHigh, efficiencies=efficiencies, alpha=0.6, T3=1000.0, T1=500.0, T2=300.0, Tmin=300, Tmax=2000, Pmin=1e3, Pmax=1e5, numReactants=2, comment='These parameters are completely made up and unrealistic')
        exec('kinetics = %r' % kinetics0)

        self.assertEqual(kinetics0.alpha.value, kinetics.alpha.value)
        self.assertEqual(kinetics0.T3.value, kinetics.T3.value)
        self.assertEqual(kinetics0.T1.value, kinetics.T1.value)
        self.assertEqual(kinetics0.T2.value, kinetics.T2.value)
        self.assertEqual(kinetics0.arrheniusLow.A.value, kinetics.arrheniusLow.A.value)
        self.assertEqual(kinetics0.arrheniusLow.n.value, kinetics.arrheniusLow.n.value)
        self.assertEqual(kinetics0.arrheniusLow.T0.value, kinetics.arrheniusLow.T0.value)
        self.assertEqual(kinetics0.arrheniusLow.Ea.value, kinetics.arrheniusLow.Ea.value)
        self.assertEqual(kinetics0.arrheniusHigh.A.value, kinetics.arrheniusHigh.A.value)
        self.assertEqual(kinetics0.arrheniusHigh.n.value, kinetics.arrheniusHigh.n.value)
        self.assertEqual(kinetics0.arrheniusHigh.T0.value, kinetics.arrheniusHigh.T0.value)
        self.assertEqual(kinetics0.arrheniusHigh.Ea.value, kinetics.arrheniusHigh.Ea.value)
        for collider, efficiency in kinetics0.efficiencies.iteritems():
            self.assertTrue(collider in kinetics.efficiencies)
            self.assertEqual(efficiency, kinetics.efficiencies[collider])

        self.assertEqual(kinetics0.Tmin.value, kinetics.Tmin.value)
        self.assertEqual(kinetics0.Tmax.value, kinetics.Tmax.value)
        self.assertEqual(kinetics0.Pmin.value, kinetics.Pmin.value)
        self.assertEqual(kinetics0.Pmax.value, kinetics.Pmax.value)
        self.assertEqual(kinetics0.numReactants, kinetics.numReactants)
        self.assertEqual(kinetics0.comment, kinetics.comment)

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
