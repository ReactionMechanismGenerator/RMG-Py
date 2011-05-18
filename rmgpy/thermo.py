#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
This module contains classes and methods for working with thermodynamics 
models. All such models derive from the :class:`ThermoModel` base class, and
generally vary by how the heat capacity data is represented:

* :class:`ThermoData` - A thermodynamics model using a discrete set of heat
  capacity data points.

* :class:`Wilhoit` - A thermodynamics model using the Wilhoit polynomial
  equation for heat capacity.

* :class:`MultiNASA` - A thermodynamics model using a set of
  :class:`NASA` objects, each representing a seven-coefficient or
  nine-coefficient polynomial equation for heat capacity, enthalpy, and entropy.

"""

import math
import numpy
import logging
import cython

from quantity import Quantity, constants

################################################################################

class ThermoError(Exception):
    """
    An exception to be raised when an error occurs while working with 
    thermodynamics data. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

################################################################################

class ThermoModel:
    """
    A base class for thermodynamics models, containing several attributes
    common to all models:
    
    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          :class:`Quantity`   The minimum temperature at which the model is valid, or ``None`` if unknown or undefined
    `Tmax`          :class:`Quantity`   The maximum temperature at which the model is valid, or ``None`` if unknown or undefined
    `comment`       ``str``             Information about the model (e.g. its source)
    =============== =================== ========================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Quantity(Tmin)
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Quantity(Tmax)
        else:
            self.Tmax = None
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoModel object.
        """
        return 'ThermoModel(Tmin={0!r}, Tmax={1!r}, comment="""{2}""")'.format(self.Tmin, self.Tmax, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoModel object.
        """
        return (ThermoModel, (self.Tmin, self.Tmax, self.comment))

    def isTemperatureValid(self, T):
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the thermodynamic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self.Tmin is None or self.Tmin.value <= T) and (self.Tmax is None or T <= self.Tmax.value)

    def getHeatCapacity(self, T):
        """
        Return the constant-pressure heat capacity in J/mol*K at temperature 
        `T` in K. This method must be overloaded in the derived class.
        """
        raise ThermoError('Unexpected call to ThermoModel.getHeatCapacity(); you should be using a class derived from ThermoModel.')

    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise ThermoError('Unexpected call to ThermoModel.getEnthalpy(); you should be using a class derived from ThermoModel.')

    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise ThermoError('Unexpected call to ThermoModel.getEntropy(); you should be using a class derived from ThermoModel.')

    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol at temperature `T` in K. This 
        method must be overloaded in the derived class.
        """
        raise ThermoError('Unexpected call to ThermoModel.getFreeEnergy(); you should be using a class derived from ThermoModel.')

    def getHeatCapacities(self, Tlist):
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperatures `Tlist` in K, as a numpy array.
        """
        return numpy.array([self.getHeatCapacity(T) for T in Tlist], numpy.float)

    def getEnthalpies(self, Tlist):
        """
        Return the enthalpy in J/mol at the specified temperatures `Tlist` in
        K, as a numpy array.
        """
        return numpy.array([self.getEnthalpy(T) for T in Tlist], numpy.float)

    def getEntropies(self, Tlist):
        """
        Return the entropy in J/mol*K at the specified temperatures `Tlist` in
        K, as a numpy array.
        """
        return numpy.array([self.getEntropy(T) for T in Tlist], numpy.float)

    def getFreeEnergies(self, Tlist):
        """
        Return the Gibbs free energy in J/mol at the specified temperatures
        `Tlist` in K, as a numpy array.
        """
        return numpy.array([self.getFreeEnergy(T) for T in Tlist], numpy.float)

################################################################################

class ThermoData(ThermoModel):
    """
    A thermodynamic model defined by a set of heat capacities. The attributes
    are:

    =========== =================== ============================================
    Attribute   Type                Description
    =========== =================== ============================================
    `Tdata`     :class:`Quantity`   The temperatures at which the heat capacity data is provided
    `Cpdata`    :class:`Quantity`   The standard heat capacity at each temperature in `Tdata`
    `H298`      :class:`Quantity`   The standard enthalpy of formation at 298 K
    `S298`      :class:`Quantity`   The standard entropy of formation at 298 K
    =========== =================== ============================================
    
    """

    def __init__(self, Tdata, Cpdata, H298, S298, Tmin=None, Tmax=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Tdata = Quantity(Tdata)
        self.Cpdata = Quantity(Cpdata)
        self.H298 = Quantity(H298)
        self.S298 = Quantity(S298)
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoData object.
        """
        string = 'ThermoData(Tdata={0!r}, Cpdata={1!r}, H298={2!r}, S298={3!r}'.format(self.Tdata, self.Cpdata, self.H298, self.S298)
        if self.Tmin: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoData object.
        """
        return (ThermoData, (self.Tdata, self.Cpdata, self.H298, self.S298, self.Tmin, self.Tmax, self.comment))

    def __add__(self, other):
        """
        Add two sets of thermodynamic data together. All parameters are
        considered additive. Returns a new :class:`ThermoData` object that is
        the sum of the two sets of thermodynamic data.
        """
        cython.declare(i=int, new=ThermoData)
        if len(self.Tdata.values) != len(other.Tdata.values) or any([T1 != T2 for T1, T2 in zip(self.Tdata.values, other.Tdata.values)]):
            raise ThermoError('Cannot add these ThermoData objects due to their having different temperature points.')
        new = ThermoData(
            Tdata = self.Tdata,
            Cpdata = self.Cpdata + other.Cpdata,
            H298 = self.H298 + other.H298,
            S298 = self.S298 + other.S298,
        )
        if self.comment == '': new.comment = other.comment
        elif other.comment == '': new.comment = self.comment
        else: new.comment = self.comment + ' + ' + other.comment
        return new

    def getHeatCapacity(self, T):
        """
        Return the constant-pressure heat capacity in J/mol*K at temperature 
        `T` in K.
        """
        cython.declare(Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
        cython.declare(Cp=cython.double)
        Cp = 0.0
        if not self.isTemperatureValid(T):
            raise ThermoError('Invalid temperature {0:g} K for heat capacity estimation.'.format(T))
        if T < numpy.min(self.Tdata.values):
            Cp = self.Cpdata.values[0]
        elif T >= numpy.max(self.Tdata.values):
            Cp = self.Cpdata.values[-1]
        else:
            for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata.values[:-1], self.Tdata.values[1:], self.Cpdata.values[:-1], self.Cpdata.values[1:]):
                if Tmin <= T and T < Tmax:
                    Cp = (Cpmax - Cpmin) * ((T - Tmin) / (Tmax - Tmin)) + Cpmin
        return Cp
    
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at temperature `T` in K.
        """
        cython.declare(H=cython.double, slope=cython.double, intercept=cython.double,
             Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
        H = self.H298.value
        if self.Tdata.values[0] > 298:
            H += self.Cpdata.values[0] * (self.Tdata.values[0] - 298)
        if not self.isTemperatureValid(T):
            raise ThermoError('Invalid temperature {0:g} K for enthalpy estimation.'.format(T))
        for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata.values[:-1], self.Tdata.values[1:], self.Cpdata.values[:-1], self.Cpdata.values[1:]):
            if T > Tmin:
                slope = (Cpmax - Cpmin) / (Tmax - Tmin)
                intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
                if T < Tmax:	H += 0.5 * slope * (T*T - Tmin*Tmin) + intercept * (T - Tmin)
                else:			H += 0.5 * slope * (Tmax*Tmax - Tmin*Tmin) + intercept * (Tmax - Tmin)
        if T > self.Tdata.values[-1]:
            H += self.Cpdata.values[-1] * (T - self.Tdata.values[-1])
        return H

    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at temperature `T` in K.
        """
        cython.declare(S=cython.double, slope=cython.double, intercept=cython.double,
             Tmin=cython.double, Tmax=cython.double, Cpmin=cython.double, Cpmax=cython.double)
        S = self.S298.value
        if self.Tdata.values[0] > 298:
            S += self.Cpdata.values[0] * math.log(self.Tdata.values[0] / 298)
        if not self.isTemperatureValid(T):
            raise ThermoError('Invalid temperature {0:g} K for entropy estimation.'.format(T))
        for Tmin, Tmax, Cpmin, Cpmax in zip(self.Tdata.values[:-1], self.Tdata.values[1:], self.Cpdata.values[:-1], self.Cpdata.values[1:]):
            if T > Tmin:
                slope = (Cpmax - Cpmin) / (Tmax - Tmin)
                intercept = (Cpmin * Tmax - Cpmax * Tmin) / (Tmax - Tmin)
                if T < Tmax:	S += slope * (T - Tmin) + intercept * math.log(T/Tmin)
                else:			S += slope * (Tmax - Tmin) + intercept * math.log(Tmax/Tmin)
        if T > self.Tdata.values[-1]:
            S += self.Cpdata.values[-1] * math.log(T / self.Tdata.values[-1])
        return S

    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol at temperature `T` in K.
        """
        if not self.isTemperatureValid(T):
            raise ThermoError('Invalid temperature {0:g} K for Gibbs free energy estimation.'.format(T))
        return self.getEnthalpy(T) - T * self.getEntropy(T)

################################################################################

class Wilhoit(ThermoModel):
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

    def __init__(self, cp0=0.0, cpInf=0.0, a0=0.0, a1=0.0, a2=0.0, a3=0.0, H0=0.0, S0=0.0, B=500.0, Tmin=None, Tmax=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.cp0 = Quantity(cp0)
        self.cpInf = Quantity(cpInf)
        self.B = Quantity(B)
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.H0 = Quantity(H0)
        self.S0 = Quantity(S0)
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Wilhoit object.
        """
        string = 'Wilhoit(cp0={0!r}, cpInf={1!r}, a0={2!r}, a1={3!r}, a2={4!r}, a3={5!r}, H0={6!r}, S0={7!r}, B={8!r}'.format(self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B)
        if self.Tmin: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Wilhoit object.
        """
        return (Wilhoit, (self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B, self.Tmin, self.Tmax, self.comment))

    def getHeatCapacity(self, T):
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
        cython.declare(y=cython.double)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0.value, self.cpInf.value, self.B.value, self.a0, self.a1, self.a2, self.a3
        y = T/(T+B)
        return cp0+(cpInf-cp0)*y*y*( 1 + (y-1)*(a0 + y*(a1 + y*(a2 + y*a3))) )
            
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at the specified temperature `T` in
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
        cython.declare(y=cython.double, y2=cython.double, logBplust=cython.double)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0.value, self.cpInf.value, self.B.value, self.a0, self.a1, self.a2, self.a3
        y = T/(T+B)
        y2 = y*y
        logBplust = math.log(B + T)
        return self.H0.value + cp0*T - (cpInf-cp0)*T*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplust))
    
    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at the specified temperature `T` in
        K. The formula is

        .. math::
            S(T) = S_0 +
            C_\\mathrm{p}(\\infty) \\ln T - \\left[ C_\\mathrm{p}(\\infty) - C_\\mathrm{p}(0) \\right]
            \\left[ \\ln y + \\left( 1 + y \\sum_{i=0}^3 \\frac{a_i y^i}{2+i} \\right) y
            \\right]

        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
        cython.declare(y=cython.double, logt=cython.double, logy=cython.double)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0.value, self.cpInf.value, self.B.value, self.a0, self.a1, self.a2, self.a3
        y = T/(T+B)
        logt = math.log(T)
        logy = math.log(y)
        return self.S0.value + cpInf*logt-(cpInf-cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
    
    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol at the specified temperature
        `T` in K.
        """
        return self.getEnthalpy(T) - T * self.getEntropy(T)
    
    def __residual(self, B, Tlist, Cplist, linear, nFreq, nRotors, H298, S298):
        # The residual corresponding to the fitToData() method
        # Parameters are the same as for that method
        cython.declare(Cp_fit=numpy.ndarray)
        self.fitToDataForConstantB(Tlist, Cplist, linear, nFreq, nRotors, B, H298, S298)
        Cp_fit = self.getHeatCapacities(Tlist)
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
        self.B = Quantity(B0,"K")
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
        self.cp0 = Quantity(3.5 * constants.R if linear else 4.0 * constants.R, "J/(mol*K)")
        self.cpInf = Quantity(self.cp0.value + (nFreq + 0.5 * nRotors) * constants.R, "J/(mol*K)")
        
        # What remains is to fit the polynomial coefficients (a0, a1, a2, a3)
        # This can be done directly - no iteration required
        y = Tlist / (Tlist + B)
        A = numpy.zeros((len(Cplist),4), numpy.float64)
        for j in range(4):
            A[:,j] = (y*y*y - y*y) * y**j
        b = ((Cplist - self.cp0.value) / (self.cpInf.value - self.cp0.value) - y*y)
        x, residues, rank, s = numpy.linalg.lstsq(A, b)
        
        self.B = Quantity(float(B), "K")
        self.a0 = float(x[0])
        self.a1 = float(x[1])
        self.a2 = float(x[2])
        self.a3 = float(x[3])

        self.H0 = Quantity(0.0,"J/mol"); self.S0 = Quantity(0.0,"J/(mol*K)")
        self.H0.value = H298 - self.getEnthalpy(298.15)
        self.S0.value = S298 - self.getEntropy(298.15)

        return self

################################################################################

class NASA(ThermoModel):
    """
    A single NASA polynomial for thermodynamic data. The `coeffs` attribute
    stores the seven or nine polynomial coefficients
    :math:`\\mathbf{a} = \\left[a_{-2}\\ a_{-1}\\ a_0\\ a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6 \\right]`
    from which the relevant thermodynamic parameters are evaluated via the
    expressions
    
    .. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_{-2} T^{-2} + a_{-1} T^{-1} + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
    
    .. math:: \\frac{H(T)}{RT} = - a_{-2} T^{-2} + a_{-1} T^{-1} \\ln T + a_0 + \\frac{1}{2} a_1 T + \\frac{1}{3} a_2 T^2 + \\frac{1}{4} a_3 T^3 + \\frac{1}{5} a_4 T^4 + \\frac{a_5}{T}
    
    .. math:: \\frac{S(T)}{R} = -\\frac{1}{2} a_{-2} T^{-2} - a_{-1} T^{-1} + a_0 \\ln T + a_1 T + \\frac{1}{2} a_2 T^2 + \\frac{1}{3} a_3 T^3 + \\frac{1}{4} a_4 T^4 + a_6
    
    The coefficients are stored internally in the nine-coefficient format, even
    when only seven coefficients are provided.
    """
    
    def __init__(self, coeffs, Tmin=None, Tmax=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if len(coeffs) == 7:
            self.cm2 = 0.0; self.cm1 = 0.0
            self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        elif len(coeffs) == 9:
            self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        else:
            raise ThermoError('Invalid number of NASA polynomial coefficients; should be 7 or 9.')
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'NASA(Tmin={0!r}, Tmax={1!r}'.format(self.Tmin, self.Tmax)
        if self.cm2 == 0 and self.cm1 == 0:
            string += ', coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g}]'.format(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        else:
            string += ', coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g}]'.format(self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (NASA, ([self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6], self.Tmin, self.Tmax, self.comment))

    def getHeatCapacity(self, T):
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        # Cp/R = am2 T^-2 + am1 T^-1 + a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4
        return ((self.cm2 / T + self.cm1) / T + self.c0 + T*(self.c1 + T*(self.c2 + T*(self.c3 + self.c4*T)))) * constants.R
    
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at the specified temperature `T` in
        K.
        """
        cython.declare(T2=cython.double, T4=cython.double)
        T2 = T*T
        T4 = T2*T2
        # H/RT = -am2 T^-2 + am1 ln(T)/T + a0 + a1 T /2 + a2 T^2 /3 + a3 T^3 /4 + a4 T^4 /5 + a5/T
        return ((-self.cm2 / T + self.cm1 * math.log(T)) / T + self.c0 + self.c1*T/2. + self.c2*T2/3. + self.c3*T2*T/4. + self.c4*T4/5. + self.c5/T) * constants.R * T
    
    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at the specified temperature `T` in
        K.
        """
        cython.declare(T2=cython.double, T4=cython.double)
        T2 = T*T
        T4 = T2*T2
        # S/R  = -am2/2 T^-2 - am1 T^-1 + a0 ln T + a1 T + a2 T^2 /2 + a3 T^3 /3 + a4 T^4 /4 + a6
        return ((-self.cm2 / T / 2. - self.cm1) / T + self.c0*math.log(T) + self.c1*T + self.c2*T2/2. +
            self.c3*T2*T/3. + self.c4*T4/4. + self.c6 ) * constants.R
    
    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol at the specified temperature
        `T` in K.
        """
        return self.getEnthalpy(T) - T * self.getEntropy(T)

################################################################################

class MultiNASA(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """
    
    def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.polynomials = polynomials or []
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiNASA object.
        """
        string = 'MultiNASA(Tmin={0!r}, Tmax={1!r}'.format(self.Tmin, self.Tmax)
        string += ', polynomials=[{0}]'.format(','.join(['%r' % poly for poly in self.polynomials]))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a MultiNASA object.
        """
        return (MultiNASA, (self.polynomials, self.Tmin, self.Tmax, self.comment))

    def getHeatCapacity(self, T):
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        return self.__selectPolynomialForTemperature(T).getHeatCapacity(T)
    
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at the specified temperature `T` in K.
        """
        return self.__selectPolynomialForTemperature(T).getEnthalpy(T)
    
    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        return self.__selectPolynomialForTemperature(T).getEntropy(T)
    
    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol at the specified temperature
        `T` in K.
        """
        return self.__selectPolynomialForTemperature(T).getFreeEnergy(T)
    
    def __selectPolynomialForTemperature(self, T):
        poly = cython.declare(NASA)
        for poly in self.polynomials:
            if poly.isTemperatureValid(T): return poly
        else:
            raise ThermoError("No valid NASA polynomial found for T={0:g} K".format(T))

################################################################################