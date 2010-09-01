#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module contains the kinetics models that are available in ChemPy.
All such models derive from the :class:`KineticsModel` base class.
"""

################################################################################

import math
import numpy
import numpy.linalg
import cython

import constants
from exception import InvalidKineticsModelError

################################################################################

class KineticsModel:
    """
    Represent a set of kinetic data. The details of the form of the kinetic
    data are left to a derived class. The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `Tmin`          :class:`float`  The minimum absolute temperature in K at which the model is valid
    `Tmax`          :class:`float`  The maximum absolute temperature in K at which the model is valid
    `Pmin`          :class:`float`  The minimum absolute pressure in Pa at which the model is valid
    `Pmax`          :class:`float`  The maximum absolute pressure in Pa at which the model is valid
    `numReactants`  :class:`int`    The number of reactants (used to determine the units of the kinetics)
    `comment`       :class:`str`    A string containing information about the model (e.g. its source)
    =============== =============== ============================================
    
    """

    def __init__(self, Tmin=0.0, Tmax=1.0e10, Pmin=0.0, Pmax=1.0e100, numReactants=-1, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.numReactants = numReactants
        self.comment = comment

    def isTemperatureValid(self, T):
        """
        Return :data:`True` if temperature `T` in K is within the valid 
        temperature range and :data:`False` if not. 
        """
        return (self.Tmin <= T and T <= self.Tmax)

    def isPressureValid(self, P):
        """
        Return :data:`True` if pressure `P` in Pa is within the valid pressure
        range, and :data:`False` if not.
        """
        return (self.Pmin <= P and P <= self.Pmax)

    def getRateCoefficients(self, Tlist):
        """
        Return the rate coefficient k(T) in SI units at temperatures
        `Tlist` in K.
        """
        return numpy.array([self.getRateCoefficient(T) for T in Tlist], numpy.float64)

################################################################################

class ArrheniusModel(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics. The kinetic expression has
    the form

    .. math:: k(T) = A \\left( \\frac{T}{T_0} \\right)^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

    where :math:`A`, :math:`n`, :math:`E_\\mathrm{a}`, and :math:`T_0` are the
    parameters to be set, :math:`T` is absolute temperature, and :math:`R` is
    the gas law constant. The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `A`             :class:`float`  The preexponential factor in s^-1, m^3/mol*s, etc.
    `T0`            :class:`float`  The reference temperature in K
    `n`             :class:`float`  The temperature exponent
    `Ea`            :class:`float`  The activation energy in J/mol
    =============== =============== ============================================
    
    """
    
    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=298.15):
        KineticsModel.__init__(self)
        self.A = A
        self.T0 = T0
        self.n = n
        self.Ea = Ea
    
    def __str__(self):
        return 'k(T) = %g * (T / %g) ** %g * exp(-%g / RT)    %g < T < %g' % (self.A, self.T0, self.n, self.Ea, self.Tmin, self.Tmax)
    
    def __repr__(self):
        return '<ArrheniusModel A=%g Ea=%g kJ/mol n=%g T0=%g K>' % (self.A,self.Ea/1000.0, self.n, self.T0)
    
    def getRateCoefficient(self, T, P=1e5):
        """
        Return the rate coefficient k(T) in SI units at temperature 
        `T` in K.
        """
        return self.A * (T / self.T0)** self.n * math.exp(-self.Ea / constants.R / T)

    def changeT0(self, T0):
        """
        Changes the reference temperature used in the exponent to `T0`, and
        adjusts the preexponential accordingly.
        """
        self.A = (self.T0 / T0)**self.n
        self.T0 = T0

    def fitToData(self, Tlist, klist, T0=298.15):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        corresponding to a set of temperatures `Tlist` in K. A linear least-
        squares fit is used, which guarantees that the resulting parameters
        provide the best possible approximation to the data.
        """
        import numpy.linalg
        A = numpy.zeros((len(Tlist),3), numpy.float64)
        A[:,0] = numpy.ones_like(Tlist)
        A[:,1] = numpy.log(Tlist / T0)
        A[:,2] = -1.0 / constants.R / Tlist
        b = numpy.log(klist)
        x = numpy.linalg.lstsq(A,b)[0]
        
        self.A = math.exp(x[0])
        self.n = x[1]
        self.Ea = x[2]
        self.T0 = T0
        return self
    
################################################################################

class ArrheniusEPModel(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
    kinetic expression has the form

    .. math:: k(T) = A T^n \\exp \\left( - \\frac{E_0 + \\alpha \\Delta H_\\mathrm{rxn}}{RT} \\right)

    The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `A`             :class:`float`  The preexponential factor in s^-1, m^3/mol*s, etc.
    `n`             :class:`float`  The temperature exponent
    `E0`            :class:`float`  The activation energy at zero enthalpy of reaction in J/mol
    `alpha`         :class:`float`  The linear dependence of activation energy on enthalpy of reaction
    =============== =============== ============================================
    
    """

    def __init__(self, A=0.0, E0=0.0, n=0.0, alpha=0.0):
        KineticsModel.__init__(self)
        self.A = A
        self.E0 = E0
        self.n = n
        self.alpha = alpha

    def __str__(self):
        return 'k(T) = %g * T ** %g * exp(-(%g + %g * dHrxn) / RT)    %g < T < %g' % (self.A, self.n, self.E0, self.alpha, self.Tmin, self.Tmax)
        
    def __repr__(self):
        return '<ArrheniusEPModel A=%g E0=%g kJ/mol n=%g alpha=%.1g>' % (self.A, self.E0/1000.0, self.n, self.alpha)
    
    def getActivationEnergy(self, dHrxn):
        """
        Return the activation energy in J/mol using the enthalpy of reaction 
        `dHrxn` in J/mol.
        """
        return self.E0 + self.alpha * dHrxn
    
    def getRateCoefficient(self, T, dHrxn):
        """
        Return the rate coefficient k(T, P) in SI units at a 
        temperature `T` in K for a reaction having an enthalpy of reaction 
        `dHrxn` in J/mol.
        """
        Ea = cython.declare(cython.double)
        Ea = self.getActivationEnergy(dHrxn)
        return self.A * (T ** self.n) * math.exp(-self.Ea / constants.R / T)

    def toArrhenius(self, dHrxn):
        """
        Return an :class:`ArrheniusModel` object corresponding to this object
        by using the provided enthalpy of reaction `dHrxn` in J/mol to calculate
        the activation energy.
        """
        return ArrheniusModel(A=self.A, n=self.n, Ea=self.getActivationEnergy(dHrxn), T0=1.0)

################################################################################

class PDepArrheniusModel(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = A(P) T^{n(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

    where the modified Arrhenius parameters are stored at a variety of pressures
    and interpolated between on a logarithmic scale. The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `pressures`     :class:`list`   The list of pressures in Pa
    `arrhenius`     :class:`list`   The list of :class:`ArrheniusModel` objects at each pressure
    =============== =============== ============================================
    
    """

    def __init__(self, pressures=None, arrhenius=None):
        KineticsModel.__init__(self)
        self.pressures = pressures or []
        self.arrhenius = arrhenius or []

    def __getAdjacentExpressions(self, P):
        """
        Returns the pressures and ArrheniusModel expressions for the pressures that
        most closely bound the specified pressure `P` in Pa.
        """
        cython.declare(Plow=cython.double, Phigh=cython.double)
        cython.declare(arrh=ArrheniusModel)
        cython.declare(i=cython.int, ilow=cython.int, ihigh=cython.int)
        
        if P in self.pressures:
            arrh = self.arrhenius[self.pressures.index(P)]
            return P, P, arrh, arrh
        else:
            ilow = 0; ihigh = -1; Plow = self.pressures[0]; Phigh = 0.0
            for i in range(1, len(self.pressures)):
                if self.pressures[i] <= P:
                    ilow = i; Plow = P
                if self.pressures[i] > P and ihigh is None:
                    ihigh = i; Phigh = P
            
            return Plow, Phigh, self.arrhenius[ilow], self.arrhenius[ihigh]
    
    def getRateCoefficient(self, T, P):
        """
        Return the rate constant k(T, P) in SI units at a temperature 
        `Tlist` in K and pressure `P` in Pa by evaluating the pressure-
        dependent Arrhenius expression.
        """
        cython.declare(Plow=cython.double, Phigh=cython.double)
        cython.declare(alow=ArrheniusModel, ahigh=ArrheniusModel)
        cython.declare(j=cython.int, klist=cython.double, klow=cython.double, khigh=cython.double)
        
        k = 0.0
        Plow, Phigh, alow, ahigh = self.__getAdjacentExpressions(P)
        if Plow == Phigh:
            k = alow.getRateCoefficient(T)
        else:
            klow = alow.getRateCoefficient(T)
            khigh = ahigh.getRateCoefficient(T)
            k = 10**(math.log10(P/Plow)/math.log10(Phigh/Plow)*math.log10(khigh/klow))
        return k

    def fitToData(self, Tlist, Plist, K, T0=298.0):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` corresponding to a set of temperatures `Tlist` in
        K and pressures `Plist` in Pa. An Arrhenius model is fit at each
        pressure.
        """
        cython.declare(i=cython.int)
        self.pressures = list(Plist)
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = ArrheniusModel()
            arrhenius.fitToData(Tlist, K[:,i], T0)
            self.arrhenius.append(arrhenius)

################################################################################

class ChebyshevModel(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression
    
    .. math:: \\log k(T,P) = \\sum_{t=1}^{N_T} \\sum_{p=1}^{N_P} \\alpha_{tp} \\phi_t(\\tilde{T}) \\phi_p(\\tilde{P})
    
    where :math:`\\alpha_{tp}` is a constant, :math:`\\phi_n(x)` is the
    Chebyshev polynomial of degree :math:`n` evaluated at :math:`x`, and
    
    .. math:: \\tilde{T} \\equiv \\frac{2T^{-1} - T_\\mathrm{min}^{-1} - T_\\mathrm{max}^{-1}}{T_\\mathrm{max}^{-1} - T_\\mathrm{min}^{-1}}
    
    .. math:: \\tilde{P} \\equiv \\frac{2 \\log P - \\log P_\\mathrm{min} - \\log P_\\mathrm{max}}{\\log P_\\mathrm{max} - \\log P_\\mathrm{min}}
    
    are reduced temperature and reduced pressures designed to map the ranges
    :math:`(T_\\mathrm{min}, T_\\mathrm{max})` and
    :math:`(P_\\mathrm{min}, P_\\mathrm{max})` to :math:`(-1, 1)`.
    The attributes are:
    
    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `coeffs`        :class:`list`   Matrix of Chebyshev coefficients
    `degreeT`       :class:`int`    The number of terms in the inverse temperature direction
    `degreeP`       :class:`int`    The number of terms in the log pressure direction
    =============== =============== ============================================
    
    """

    def __init__(self, Tmin=0.0, Tmax=0.0, Pmin=0.0, Pmax=0.0, coeffs=None):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax)
        self.coeffs = coeffs
        self.degreeT = 0
        self.degreeP = 0

    def __chebyshev(self, n, x):
        if n == 0:
            return 1
        elif n == 1:
            return x
        elif n == 2:
            return -1 + 2*x*x
        elif n == 3:
            return x * (-3 + 4*x*x)
        elif n == 4:
            return 1 + x*x*(-8 + 8*x*x)
        elif n == 5:
            return x * (5 + x*x*(-20 + 16*x*x))
        elif n == 6:
            return -1 + x*x*(18 + x*x*(-48 + 32*x*x))
        elif n == 7:
            return x * (-7 + x*x*(56 + x*x*(-112 + 64*x*x)))
        elif n == 8:
            return 1 + x*x*(-32 + x*x*(160 + x*x*(-256 + 128*x*x)))
        elif n == 9:
            return x * (9 + x*x*(-120 + x*x*(432 + x*x*(-576 + 256*x*x))))
        elif cython.compiled:
            return cos(n * acos(x))
        else:
            return math.cos(n * math.acos(x))

    def __getReducedTemperature(self, T):
        return (2.0/T - 1.0/self.Tmin - 1.0/self.Tmax) / (1.0/self.Tmax - 1.0/self.Tmin)
    
    def __getReducedPressure(self, P):
        if cython.compiled:
            return (2.0*log10(P) - log10(self.Pmin) - log10(self.Pmax)) / (log10(self.Pmax) - log10(self.Pmin))
        else:
            return (2.0*math.log(P) - math.log(self.Pmin) - math.log(self.Pmax)) / (math.log(self.Pmax) - math.log(self.Pmin))
    
    def getRateCoefficient(self, T, P):
        """
        Return the rate constant k(T, P) in SI units at a temperature 
        `Tlist` in K and pressure `P` in Pa by evaluating the Chebyshev 
        expression.
        """
        
        cython.declare(Tred=cython.double, Pred=cython.double, k=cython.double)
        cython.declare(i=cython.int, j=cython.int, t=cython.int, p=cython.int)
        
        k = 0.0
        Tred = self.__getReducedTemperature(T)
        Pred = self.__getReducedPressure(P)
        for t in range(self.degreeT):
            for p in range(self.degreeP):
                k += self.coeffs[t,p] * self.__chebyshev(t, Tred) * self.__chebyshev(p, Pred)
        return 10.0**k

    def fitToData(self, Tlist, Plist, K, degreeT, degreeP, Tmin, Tmax, Pmin, Pmax):
        """
        Fit a Chebyshev kinetic model to a set of rate coefficients `K`, which
        is a matrix corresponding to the temperatures `Tlist` in K and pressures
        `Plist` in Pa. `degreeT` and `degreeP` are the degree of the polynomials
        in temperature and pressure, while `Tmin`, `Tmax`, `Pmin`, and `Pmax`
        set the edges of the valid temperature and pressure ranges in K and Pa,
        respectively.
        """

        cython.declare(nT=cython.int, nP=cython.int, Tred=list, Pred=list)
        cython.declare(A=numpy.ndarray, b=numpy.ndarray)
        cython.declare(t1=cython.int, p1=cython.int, t2=cython.int, p2=cython.int)
        cython.declare(T=cython.double, P=cython.double)

        nT = len(Tlist); nP = len(Plist)

        self.degreeT = degreeT; self.degreeP = degreeP

        # Set temperature and pressure ranges
        self.Tmin = Tmin; self.Tmax = Tmax
        self.Pmin = Pmin; self.Pmax = Pmax

        # Calculate reduced temperatures and pressures
        Tred = [self.__getReducedTemperature(T) for T in Tlist]
        Pred = [self.__getReducedPressure(P) for P in Plist]

        # Create matrix and vector for coefficient fit (linear least-squares)
        A = numpy.zeros((nT*nP, degreeT*degreeP), numpy.float64)
        b = numpy.zeros((nT*nP), numpy.float64)
        for t1, T in enumerate(Tred):
            for p1, P in enumerate(Pred):
                for t2 in range(degreeT):
                    for p2 in range(degreeP):
                        A[p1*nT+t1, p2*degreeT+t2] = self.__chebyshev(t2, T) * self.__chebyshev(p2, P)
                b[p1*nT+t1] = math.log10(K[t1,p1])

        # Do linear least-squares fit to get coefficients
        x, residues, rank, s = numpy.linalg.lstsq(A, b)

        # Extract coefficients
        self.coeffs = numpy.zeros((degreeT,degreeP), numpy.float64)
        for t2 in range(degreeT):
            for p2 in range(degreeP):
                self.coeffs[t2,p2] = x[p2*degreeT+t2]
    
