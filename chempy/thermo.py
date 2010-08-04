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
This module contains the thermodynamics models that are available in ChemPy.
All such models derive from the :class:`ThermoModel` base class.
"""

################################################################################

import math
import numpy
import cython

import constants
from exception import InvalidThermoModelError

################################################################################

class ThermoError:
    """
    An exception class for errors that occur while working with thermodynamics
    models. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class ThermoModel:
    """
    A base class for thermodynamics models, containing several attributes
    common to all models:
    
    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `Tmin`          :class:`float`  The minimum temperature in K at which the model is valid
    `Tmax`          :class:`float`  The maximum temperature in K at which the model is valid
    `comment`       :class:`str`    A string containing information about the model (e.g. its source)
    =============== =============== ============================================
    
    """
    
    def __init__(self, Tmin=0.0, Tmax=1.0e10, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comment = comment
    
    def isTemperatureValid(self, T):
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the thermodynamic data, or ``False`` if not.
        """
        return self.Tmin <= T and T <= self.Tmax

    def getHeatCapacity(self, Tlist):
        raise InvalidThermoModelError('Unexpected call to ThermoModel.getHeatCapacity(); you should be using a class derived from ThermoModel.')

    def getEnthalpy(self, Tlist):
        raise InvalidThermoModelError('Unexpected call to ThermoModel.getEnthalpy(); you should be using a class derived from ThermoModel.')

    def getEntropy(self, Tlist):
        raise InvalidThermoModelError('Unexpected call to ThermoModel.getEntropy(); you should be using a class derived from ThermoModel.')

    def getFreeEnergy(self, Tlist):
        raise InvalidThermoModelError('Unexpected call to ThermoModel.getFreeEnergy(); you should be using a class derived from ThermoModel.')

################################################################################

class ThermoGAModel(ThermoModel):
    """
    A thermodynamic model defined by a set of heat capacities. The attributes
    are:

    =========== =================== ============================================
    Attribute   Type                Description
    =========== =================== ============================================
    `Tdata`     ``numpy.ndarray``   The temperatures at which the heat capacity data is provided in K
    `Cpdata`    ``numpy.ndarray``   The standard heat capacity in J/mol*K at each temperature in `Tdata`
    `H298`      ``double``          The standard enthalpy of formation at 298 K in J/mol
    `S298`      ``double``          The standard entropy of formation at 298 K in J/mol*K
    =========== =================== ============================================
    """

    def __init__(self, Tdata=None, Cpdata=None, H298=0.0, S298=0.0, Tmin=0.0, Tmax=99999.9, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Tdata = Tdata
        self.Cpdata = Cpdata
        self.H298 = H298
        self.S298 = S298
    
    def __repr__(self):
        string = 'ThermoGAModel(Tdata=%s, Cpdata=%s, H298=%s, S298=%s)' % (self.Tdata, self.Cpdata, self.H298, self.S298)
        return string

    def __str__(self):
        """
        Return a string summarizing the thermodynamic data.
        """
        string = ''
        string += 'Enthalpy of formation: %g kJ/mol\n' % (self.H298 / 1000.0)
        string += 'Entropy of formation: %g J/mol*K\n' % (self.S298)
        string += 'Heat capacity (J/mol*K): '
        for T, Cp in zip(self.Tdata, self.Cpdata):
            string += '%.1f(%g K) ' % (Cp,T)
        string += '\n'
        string += 'Comment: %s' % (self.comment)
        return string

    def __add__(self, other):
        """
        Add two sets of thermodynamic data together. All parameters are
        considered additive. Returns a new :class:`ThermoGAModel` object that is
        the sum of the two sets of thermodynamic data.
        """
        cython.declare(i=int, new=ThermoGAModel)
        if len(self.Tdata) != len(other.Tdata) or any([T1 != T2 for T1, T2 in zip(self.Tdata, other.Tdata)]):
            raise Exception('Cannot add these ThermoGAModel objects due to their having different temperature points.')
        new = ThermoGAModel()
        new.H298 = self.H298 + other.H298
        new.S298 = self.S298 + other.S298
        new.Tdata = self.Tdata
        new.Cpdata = self.Cpdata + other.Cpdata
        if self.comment == '': new.comment = other.comment
        elif other.comment == '': new.comment = self.comment
        else: new.comment = self.comment + ' + ' + other.comment
        return new

    def getHeatCapacity(self, Tlist):
        """
        Return the constant-pressure heat capacity (Cp) in J/mol*K at temperature `T` in K.
        """
        cython.declare(Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
        cython.declare(Cplist=numpy.ndarray)
        Cplist = numpy.zeros_like(Tlist)
        for i, T in enumerate(Tlist):
            if not self.isTemperatureValid(T):
                raise ThermoError('Invalid temperature "%g K" for heat capacity estimation.' % T)
            if T < numpy.min(Tlist):
                Cplist[i] = self.Cpdata[0]
            elif T >= numpy.max(Tlist):
                Cplist[i] = self.Cpdata[0]
            else:
                for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata[:-1], self.Tdata[1:], self.Cpdata[:-1], self.Cpdata[1:]):
                    if Tmin <= T and T < Tmax:
                        Cplist[i] = (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin
        return Cplist
    
    def getEnthalpy(self, Tlist):
        """
        Return the enthalpy in J/mol at temperature `T` in K.
        """
        cython.declare(H=cython.double, slope=cython.double, intercept=cython.double,
             Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
        cython.declare(Hlist=numpy.ndarray)
        Hlist = self.H298 * numpy.ones_like(Tlist)
        for i, T in enumerate(Tlist):
            if not self.isTemperatureValid(T):
                raise ThermoError('Invalid temperature "%g K" for enthalpy estimation.' % T)
            for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata[:-1], self.Tdata[1:], self.Cpdata[:-1], self.Cpdata[1:]):
                if T > Tmin:
                    slope = (Cpmax - Cpmin) / (Tmax - Tmin)
                    intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
                    if T < Tmax:	Hlist[i] += 0.5 * slope * (T*T - Tmin*Tmin) + intercept * (T - Tmin)
                    else:			Hlist[i] += 0.5 * slope * (Tmax*Tmax - Tmin*Tmin) + intercept * (Tmax - Tmin)
            if T > self.Tdata[-1]:
                Hlist[i] += self.Cp[-1] * (T - self.Tdata[-1])
        return Hlist

    def getEntropy(self, Tlist):
        """
        Return the entropy in J/mol*K at temperature `T` in K.
        """

        Slist = self.H298 * numpy.ones_like(Tlist)
        for i, T in enumerate(Tlist):
            if not self.isTemperatureValid(T):
                raise ThermoError('Invalid temperature "%g K" for entropy estimation.' % T)
            for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata[:-1], self.Tdata[1:], self.Cpdata[:-1], self.Cpdata[1:]):
                if T > Tmin:
                    slope = (Cpmax - Cpmin) / (Tmax - Tmin)
                    intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
                    if T < Tmax:	Slist[i] += slope * (T - Tmin) + intercept * math.log(T/Tmin)
                    else:			Slist[i] += slope * (Tmax - Tmin) + intercept * math.log(Tmax/Tmin)
            if T > self.Tdata[-1]:
                Slist[i] += self.Cp[-1] * math.log(T / self.Tdata[-1])
        return Slist

    def getFreeEnergy(self, Tlist):
        """
        Return the Gibbs free energy in J/mol at temperature `T` in K.
        """
        for T in Tlist:
            if not self.isTemperatureValid(T):
                raise ThermoError('Invalid temperature "%g K" for Gibbs free energy estimation.' % T)
        return self.getEnthalpy(Tlist) - T * self.getEntropy(Tlist)

################################################################################

class WilhoitModel(ThermoModel):
    """
    A thermodynamics model based on the Wilhoit equation for heat capacity,

    .. math::
        C_\\mathrm{p}(T) = C_\\mathrm{p}(0) + \\left[ C_\\mathrm{p}(\\infty) -
        C_\\mathrm{p}(0) \\right] y^2 \\left[ 1 + (y - 1) \\sum_{i=0}^3 a_i y^i \\right]

    where :math:`y \\equiv \\frac{T}{T + B}` is a scaled temperature that ranges
    from zero to one. (The characteristic temperature :math:`B` is chosen by
    default to be 500 K.) This formulation has the advantage of correctly
    reproducting the heat capacity behavior as :math:`T \\rightarrow 0` and
    :math:`T \\rightarrow \\infty`. The low-temperature limit 
    :math:`C_\\mathrm{p}(0)` is taken to be :math:`3.5R` for linear molecules
    and :math:`4R` for nonlinear molecules. The high-temperature limit 
    :math:`C_\\mathrm{p}(\\infty)` is taken to be 
    :math:`\\left[ 3 N_\\mathrm{atoms} - 1.5 \\right] R` for linear molecules and
    :math:`\\left[ 3 N_\\mathrm{atoms} - (2 + 0.5 N_\\mathrm{rotors}) \\right] R`
    for nonlinear molecules, for a molecule composed of :math:`N_\\mathrm{atoms}`
    atoms and :math:`N_\\mathrm{rotors}` internal rotors.
    
    The Wilhoit parameters are stored in the attributes `cp0`, `cpInf`, `a0`,
    `a1`, `a2`, `a3`, and `B`. There are also integration constants `H0` and
    `S0` that are needed to evaluate the enthalpy and entropy, respectively.
    """

    def __init__(self, cp0=0.0, cpInf=0.0, a0=0.0, a1=0.0, a2=0.0, a3=0.0, H0=0.0, S0=0.0, comment='', B=500.0):
        ThermoModel.__init__(self, comment=comment)
        self.cp0 = cp0
        self.cpInf = cpInf
        self.B = B
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.H0 = H0
        self.S0 = S0
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        return 'WilhoitModel(cp0=%g, cpInf=%g, a0=%g, a1=%g, a2=%g, a3=%g, H0=%g, S0=%g, B=%g)' % (self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B)
    
    def getHeatCapacity(self, Tlist):
        """
        Return the constant-pressure heat capacity (Cp) in J/mol*K at the
        specified temperatures `Tlist` in K.
        """
        cython.declare(y=numpy.ndarray)
        y = Tlist/(Tlist+self.B)
        return self.cp0+(self.cpInf-self.cp0)*y*y*( 1 +
            (y-1)*(self.a0 + y*(self.a1 + y*(self.a2 + y*self.a3))) )
    
    def getEnthalpy(self, Tlist):
        """
        Return the enthalpy in J/mol at the specified temperatures `Tlist` in
        K. The formula is

        .. math::
            H(T) & = H_0 +
            C_\\mathrm{p}(0) T + \\left[ C_\\mathrm{p}(\\infty) - C_\\mathrm{p}(0) \\right] T \\\\
            & \\left\\{ \\left[ 2 + \\sum_{i=0}^3 a_i \\right]
            \\left[ \\frac{1}{2}y - 1 + \\left( \\frac{1}{y} - 1 \\right) \\ln \\frac{T}{y} \\right]
            + y^2 \\sum_{i=0}^3 \\frac{y^i}{(i+2)(i+3)} \\sum_{j=0}^3 f_{ij} a_j
            \\right\\}

        where :math:`f_{ij} = 3 + j` if :math:`i = j`, :math:`f_{ij} = 1` if
        :math:`i > j`, and :math:`f_{ij} = 0` if :math:`i < j`.
        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
        cython.declare(y=numpy.ndarray, y2=numpy.ndarray, logBplust=numpy.ndarray)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
        y = Tlist/(Tlist+B)
        y2 = y*y
        logBplust = numpy.log(B + Tlist)
        return self.H0 + cp0*Tlist - (cpInf-cp0)*Tlist*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplust))
    
    def getEntropy(self, Tlist):
        """
        Return the entropy in J/mol*K at the specified temperatures `Tlist` in
        K. The formula is

        .. math::
            S(T) = S_0 +
            C_\\mathrm{p}(\\infty) \\ln T - \\left[ C_\\mathrm{p}(\\infty) - C_\\mathrm{p}(0) \\right]
            \\left[ \\ln y + \\left( 1 + y \\sum_{i=0}^3 \\frac{a_i y^i}{2+i} \\right) y
            \\right]

        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
        cython.declare(y=numpy.ndarray, logt=numpy.ndarray, logy=numpy.ndarray)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
        y = Tlist/(Tlist+B)
        logt = numpy.log(Tlist)
        logy = numpy.log(y)
        return self.S0 + cpInf*logt-(cpInf-cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
    
    def getFreeEnergy(self, Tlist):
        """
        Return the Gibbs free energy in J/mol at the specified temperatures
        `Tlist` in K.
        """
        return self.getEnthalpy(Tlist) - Tlist * self.getEntropy(Tlist)
    
    def __residual(self, B, Tlist, Cplist, linear, nFreq, nRotors, H298, S298):
        # The residual corresponding to the fitToData() method
        # Parameters are the same as for that method
        cython.declare(Cp_fit=numpy.ndarray)
        self.fitToDataForConstantB(Tlist, Cplist, linear, nFreq, nRotors, B, H298, S298)
        Cp_fit = self.getHeatCapacity(Tlist)
        # Objective function is linear least-squares
        return numpy.sum( (Cp_fit - Cplist) * (Cp_fit - Cplist) )
    
    def fitToData(self, Tlist, Cplist, linear, nFreq, nRotors, H298, S298, B0=500.0):
        """
        Fit a Wilhoit model to the data points provided, allowing the 
        characteristic temperature `B` to vary so as to improve the fit. This
        procedure requires an optimization, using the ``fminbound`` function
        in the ``scipy.optimize`` module. The data consists of a set
        of dimensionless heat capacity points `Cplist` at a given set of
        temperatures `Tlist` in K. The linearity of the molecule, number of
        vibrational frequencies, and number of internal rotors (`linear`,
        `nFreq`, and `nRotors`, respectively) is used to set the limits at
        zero and infinite temperature.
        """
        self.B = B0
        import scipy.optimize
        scipy.optimize.fminbound(self.__residual, 300.0, 3000.0, args=(Tlist, Cplist, linear, nFreq, nRotors, H298, S298))
        return self
    
    def fitToDataForConstantB(self, Tlist, Cplist, linear, nFreq, nRotors, B, H298, S298):
        """
        Fit a Wilhoit model to the data points provided using a specified value
        of the characteristic temperature `B`. The data consists of a set
        of dimensionless heat capacity points `Cplist` at a given set of
        temperatures `Tlist` in K. The linearity of the molecule, number of
        vibrational frequencies, and number of internal rotors (`linear`,
        `nFreq`, and `nRotors`, respectively) is used to set the limits at
        zero and infinite temperature.
        """
        
        cython.declare(y=numpy.ndarray, A=numpy.ndarray, b=numpy.ndarray, x=numpy.ndarray)
        
        # Set the Cp(T) limits as T -> and T -> infinity
        self.cp0 = 3.5 * constants.R if linear else 4.0 * constants.R
        self.cpInf = self.cp0 + (nFreq + 0.5 * nRotors) * constants.R
        
        # What remains is to fit the polynomial coefficients (a0, a1, a2, a3)
        # This can be done directly - no iteration required
        y = Tlist / (Tlist + B)
        A = numpy.zeros((len(Cplist),4), numpy.float64)
        for j in range(4):
            A[:,j] = (y*y*y - y*y) * y**j
        b = ((Cplist - self.cp0) / (self.cpInf - self.cp0) - y*y)
        x, residues, rank, s = numpy.linalg.lstsq(A, b)
        
        self.B = float(B)
        self.a0 = float(x[0])
        self.a1 = float(x[1])
        self.a2 = float(x[2])
        self.a3 = float(x[3])

        Tlist0 = numpy.array([298.15], numpy.float64)
        self.H0 = 0.0; self.S0 = 0.0
        self.H0 = H298 - self.getEnthalpy(Tlist0)[0]
        self.S0 = S298 - self.getEntropy(Tlist0)[0]

        return self

################################################################################

class NASAPolynomial(ThermoModel):
    """
    A single NASA polynomial for thermodynamic data. The `coeffs` attribute
    stores the seven polynomial coefficients
    :math:`\\mathbf{a} = \\left[a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6\\ a_7 \\right]`
    from which the relevant thermodynamic parameters are evaluated via the
    expressions
    
    .. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4
    
    .. math:: \\frac{H(T)}{RT} = a_1 + \\frac{1}{2} a_2 T + \\frac{1}{3} a_3 T^2 + \\frac{1}{4} a_4 T^3 + \\frac{1}{5} a_5 T^4 + \\frac{a_6}{T}
    
    .. math:: \\frac{S(T)}{R} = a_1 \\ln T + a_2 T + \\frac{1}{2} a_3 T^2 + \\frac{1}{3} a_4 T^3 + \\frac{1}{4} a_5 T^4 + a_7
    
    The above was adapted from `this page <http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html>`_.
    """
    
    def __init__(self, Tmin=0.0, Tmax=0.0, coeffs=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        return 'NASAPolynomial(Tmin=%g, Tmax=%g, coeffs=[%g, %g, %g, %g, %g, %g, %g])' % (self.Tmin, self.Tmax, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
    
    def getHeatCapacity(self, Tlist):
        """
        Return the constant-pressure heat capacity (Cp) in J/mol*K at the
        specified temperatures `Tlist` in K.
        """
        # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        return (self.c0 + Tlist*(self.c1 + Tlist*(self.c2 + Tlist*(self.c3 + self.c4*Tlist)))) * constants.R
    
    def getEnthalpy(self, Tlist):
        """
        Return the enthalpy in J/mol at the specified temperatures `Tlist` in
        K.
        """
        cython.declare(T2=numpy.ndarray, T4=numpy.ndarray)
        T2 = Tlist*Tlist
        T4 = T2*T2
        # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
        return (self.c0 + self.c1*Tlist/2 + self.c2*T2/3 + self.c3*T2*Tlist/4 + self.c4*T4/5 + self.c5/Tlist) * constants.R * Tlist
    
    def getEntropy(self, Tlist):
        """
        Return the entropy in J/mol*K at the specified temperatures `Tlist` in
        K.
        """
        cython.declare(T2=numpy.ndarray, T4=numpy.ndarray)
        T2 = Tlist*Tlist
        T4 = T2*T2
        # S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
        return ( self.c0*numpy.log(Tlist) + self.c1*Tlist + self.c2*T2/2 +
            self.c3*T2*Tlist/3 + self.c4*T4/4 + self.c6 ) * constants.R
    
    def getFreeEnergy(self, Tlist):
        """
        Return the Gibbs free energy in J/mol at the specified temperatures
        `Tlist` in K.
        """
        return self.getEnthalpy(Tlist) - Tlist * self.getEntropy(Tlist)

################################################################################

class NASAModel(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASAPolynomial` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """
    
    def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.polynomials = polynomials or []
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the 
        object.
        """
        return 'NASAPolynomial(Tmin=%g, Tmax=%g, polynomials=%s)' % (self.Tmin, self.Tmax, self.polynomials)
    
    def getHeatCapacity(self, Tlist):
        """
        Return the constant-pressure heat capacity (Cp) in J/mol*K at the
        specified temperatures `Tlist` in K.
        """
        cython.declare(poly=NASAPolynomial, Cp=numpy.ndarray)
        Cp = numpy.zeros_like(Tlist)
        for i in range(len(Tlist)):
            poly = self.__selectPolynomialForTemperature(Tlist[i])
            Cp[i] = poly.getHeatCapacity(Tlist[i])
        return Cp
    
    def getEnthalpy(self, Tlist):
        """
        Return the enthalpy in J/mol at the specified temperatures `Tlist` in
        K.
        """
        cython.declare(poly=NASAPolynomial, H=numpy.ndarray)
        H = numpy.zeros_like(Tlist)
        for i in range(len(Tlist)):
            poly = self.__selectPolynomialForTemperature(Tlist[i])
            H[i] = poly.getEnthalpy(Tlist[i])
        return H
    
    def getEntropy(self, Tlist):
        """
        Return the entropy in J/mol*K at the specified temperatures `Tlist` in
        K.
        """
        cython.declare(poly=NASAPolynomial, S=numpy.ndarray)
        S = numpy.zeros_like(Tlist)
        for i in range(len(Tlist)):
            poly = self.__selectPolynomialForTemperature(Tlist[i])
            S[i] = poly.getEntropy(Tlist[i])
        return S
    
    def getFreeEnergy(self, Tlist):
        """
        Return the Gibbs free energy in J/mol at the specified temperatures
        `Tlist` in K.
        """
        cython.declare(poly=NASAPolynomial, G=numpy.ndarray)
        G = numpy.zeros_like(Tlist)
        for i in range(len(Tlist)):
            poly = self.__selectPolynomialForTemperature(Tlist[i])
            G[i] = poly.getFreeEnergy(Tlist[i])
        return G
    
    def __selectPolynomialForTemperature(self, T):
        poly = cython.declare(NASAPolynomial)
        for poly in self.polynomials:
            if poly.isTemperatureValid(T): return poly
        else:
            raise InvalidThermoModelError("No polynomial found for T=%s" % T)

################################################################################