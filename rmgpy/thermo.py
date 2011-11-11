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
        # it might be faster to compute this using numpy.linalg.lstsq "residues",
        # computed during fitToDataForConstantB, but if speed here is not critical, this is fine
        return numpy.sum( (Cp_fit - Cplist) * (Cp_fit - Cplist) )
    
    def fitToData(self, Tlist, Cplist, linear, nFreq, nRotors, H298, S298, B0=500.0, Bmin=300.0, Bmax=3000.0):
        """
        Fit a Wilhoit model to the data points provided, allowing the 
        characteristic temperature `B` to vary so as to improve the fit. This
        procedure requires an optimization, using the ``fminbound`` function
        in the ``scipy.optimize`` module. The data consists of a set
        of dimensionless heat capacity points `Cplist` at a given set of
        temperatures `Tlist` in K. The linearity of the molecule, number of
        vibrational frequencies (not including internal rotors), and number
        of internal rotors (`linear`,`nFreq`, and `nRotors`, respectively) is
        used to set the limits at zero and infinite temperature.
        """
        self.B = Quantity(B0,"K")
        import scipy.optimize
        scipy.optimize.fminbound(self.__residual, Bmin, Bmax, args=(Tlist, Cplist, linear, nFreq, nRotors, H298, S298))
        #compute the rmsErr of Cp/R
        #note that IF the "residues" from the least squares fit was used instead, the error is
        #minimized with respect to Cpdata/(CpInf-Cp0) + const, so we would need to
        #scale back by (CpInf-Cp0) to get the error in Cpdata
        rmsErr = math.sqrt(self.__residual(self.B.value, Tlist, Cplist, linear, nFreq, nRotors, H298, S298)/len(Tlist))/constants.R
        self.comment = self.comment + 'Wilhoit polynomial fit to ThermoData with RMS error = %.3f*R;'%(rmsErr)
        #print a warning if the rms fit is worse that 0.25*R
        if(rmsErr > 0.25):
            logging.warning("Poor ThermoData-to-Wilhoit fit quality: RMS error = %.3f*R" % (rmsErr))
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
        
        if nFreq == 0:
            # Monatomic species
            assert nRotors == 0
            self.cp0 = Quantity(2.5 * constants.R, "J/(mol*K)")
            self.cpInf = Quantity(2.5 * constants.R, "J/(mol*K)")
            self.B = Quantity(float(B), "K")
            self.a0 = 0.0
            self.a1 = 0.0
            self.a2 = 0.0
            self.a3 = 0.0
    
        else:
            # Polyatomic species
    
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

def convertThermoModel(model, thermoClass, **kwargs):
    """
    Convert a given thermodynamics `model` to an object of the class 
    `thermoClass`. Depending on the desired conversion, various additional
    keyword arguments may be required:
    
    * In general, if you are converting to a :class:`Wilhoit` model -- even if 
      this conversion is only part of a conversion to another format -- you 
      need to specify whether or not the molecule is linear, the number of 
      vibrational modes (excluding internal hindered rotors),and the number of
      internal hindered rotor modes as the `linear`, `nFreq`, and `nRotors`
      keyword arguments, so that the correct limits at zero and infinite
      temperature can be determined.
    
    * If you are converting to a :class:`ThermoData` model, you must provide 
      the set of temperatures in K to use for the heat capacity points via the
      `Tdata` keyword argument.
    
    A :class:`ThermoError` will be raised if you do not provide the appropriate
    keyword arguments for the conversion you are attempting.
    """
    # Make sure we've requested a valid output format
    if thermoClass not in [ThermoData, Wilhoit, MultiNASA]:
        raise ValueError('Unexpected value "{0}" for thermoClass parameter; should be one of ["ThermoData", "Wilhoit", "MultiNASA"].'.format(thermoClass))
    
    # If the output format is the same as the input format, just return the same model
    if isinstance(model, thermoClass):
        return model
    
    # Do the conversion
    output = None
    if isinstance(model, ThermoData) and thermoClass == Wilhoit:
        try:
            linear = kwargs['linear']
            nFreq = kwargs['nFreq']
            nRotors = kwargs['nRotors']
        except KeyError:
            raise ValueError('To convert ThermoData to Wilhoit or MultiNASA, you must provide the keyword arguments linear, nFreq, and nRotors.')
        output = Wilhoit().fitToData(model.Tdata.values, model.Cpdata.values, linear, nFreq, nRotors, model.H298.value, model.S298.value)
    
    elif isinstance(model, ThermoData) and thermoClass == MultiNASA:
        # First convert it to a Wilhoit, then to a MultiNASA
        output = convertThermoModel(convertThermoModel(model, Wilhoit, **kwargs), MultiNASA, **kwargs)
	# compute error for the overall conversion
        Cp_fit = output.getHeatCapacities(model.Tdata.values)
        rmsErr = math.sqrt(numpy.sum( (Cp_fit - model.Cpdata.values) * (Cp_fit - model.Cpdata.values) )/len(model.Tdata.values))/constants.R
        output.comment = output.comment + 'Overall conversion of ThermoData to MultiNASA with RMS error = %.3f*R;'%(rmsErr)
        #there is already a warning in model.py's generateThermoData
        #if(rmsErr > 0.50):#print a warning if the rms fit is worse that 0.50*R
        #    logging.warning("Poor ThermoData-to-MultiNASA fit quality: RMS error = %.3f*R" % (rmsErr))

    elif isinstance(model, Wilhoit) and thermoClass == ThermoData:
        try:
            Tdata = kwargs['Tdata']
        except KeyError:
            raise ValueError('To convert Wilhoit to ThermoData, you must provide the keyword argument Tdata.')
        output = ThermoData(
            Tdata = (Tdata,"K"),
            Cpdata = ([model.getHeatCapacity(T) for T in Tdata],"J/(mol*K)"),
            H298 = (model.getEnthalpy(298) / 1000.0,"kJ/mol"),
            S298 = (model.getEntropy(298),"J/(mol*K)"),
            Tmin = model.Tmin,
            Tmax = model.Tmax,
            comment = model.comment,
        )
        
    elif isinstance(model, Wilhoit) and thermoClass == MultiNASA:
        try:
            Tmin = kwargs['Tmin']
            Tmax = kwargs['Tmax']
            Tint = kwargs['Tint']
        except KeyError:
            raise ValueError('To convert Wilhoit to MultiNASA, you must provide the keyword arguments Tmin, Tmax, and Tint.')
        output = convertWilhoitToNASA(model, Tmin, Tmax, Tint)
        
    elif isinstance(model, MultiNASA) and thermoClass == ThermoData:
        try:
            Tdata = kwargs['Tdata']
        except KeyError:
            raise ValueError('To convert MultiNASA to ThermoData, you must provide the keyword argument Tdata.')
        output = ThermoData(
            Tdata = (Tdata,"K"),
            Cpdata = ([model.getHeatCapacity(T) for T in Tdata],"J/(mol*K)"),
            H298 = (model.getEnthalpy(298) / 1000.0,"kJ/mol"),
            S298 = (model.getEntropy(298),"J/(mol*K)"),
            Tmin = model.Tmin,
            Tmax = model.Tmax,
            comment = model.comment,
        )
    
    elif isinstance(model, MultiNASA) and thermoClass == Wilhoit:
        # First convert it to a ThermoData, then to a Wilhoit
        Tdata = numpy.arange(model.Tmin.value, model.Tmax.value, 50.0)
        output = convertThermoModel(convertThermoModel(model, ThermoData, Tdata=Tdata), Wilhoit, **kwargs)
    
    return output

def convertWilhoitToNASA(wilhoit, Tmin, Tmax, Tint, fixedTint=False, weighting=True, continuity=3):
    """
    Convert a :class:`Wilhoit` object `Wilhoit` to a :class:`MultiNASA` 
    object. You must specify the minimum and maximum temperatures of the fit
    `Tmin` and `Tmax`, as well as the intermediate temperature `Tint` to use
    as the bridge between the two fitted polynomials. The remaining parameters
    can be used to modify the fitting algorithm used:
    
    * `fixedTint` - ``False`` to allow `Tint` to vary in order to improve the fit, or ``True`` to keep it fixed

    * `weighting` - ``True`` to weight the fit by :math:`T^{-1}` to emphasize good fit at lower temperatures, or ``False`` to not use weighting

    * `continuity` - The number of continuity constraints to enforce at `Tint`:

        - 0: no constraints on continuity of :math:`C_\\mathrm{p}(T)` at `Tint`

        - 1: constrain :math:`C_\\mathrm{p}(T)` to be continous at `Tint`

        - 2: constrain :math:`C_\\mathrm{p}(T)` and :math:`\\frac{d C_\\mathrm{p}}{dT}` to be continuous at `Tint`

        - 3: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, and :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}` to be continuous at `Tint`

        - 4: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, and :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}` to be continuous at `Tint`

        - 5: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}`, and :math:`\\frac{d^4 C_\\mathrm{p}}{dT^4}` to be continuous at `Tint`
        
    Note that values of `continuity` of 5 or higher effectively constrain all
    the coefficients to be equal and should be equivalent to fitting only one
    polynomial (rather than two).

    Returns the fitted :class:`MultiNASA` object containing the two fitted
    :class:`NASA` objects.
    """

    # Scale the temperatures to kK
    Tmin /= 1000.
    Tint /= 1000.
    Tmax /= 1000.
    
    # Make copy of Wilhoit data so we don't modify the original
    wilhoit_scaled = Wilhoit(wilhoit.cp0, wilhoit.cpInf, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3, wilhoit.H0, wilhoit.S0, wilhoit.B, Tmin=wilhoit.Tmin, Tmax=wilhoit.Tmax, comment=wilhoit.comment)
    # Rescale Wilhoit parameters
    wilhoit_scaled.cp0.value /= constants.R
    wilhoit_scaled.cpInf.value /= constants.R
    wilhoit_scaled.B.value /= 1000.

    #if we are using fixed Tint, do not allow Tint to float
    if fixedTint:
        nasa_low, nasa_high = Wilhoit2NASA(wilhoit_scaled, Tmin, Tmax, Tint, weighting, continuity)
    else:
        nasa_low, nasa_high, Tint = Wilhoit2NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting, continuity)
    iseUnw = TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, 0, continuity) #the scaled, unweighted ISE (integral of squared error)
    rmsUnw = math.sqrt(iseUnw/(Tmax-Tmin))
    rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
    if(weighting == 1):
        iseWei= TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, weighting, continuity) #the scaled, weighted ISE
        rmsWei = math.sqrt(iseWei/math.log(Tmax/Tmin))
        rmsStr = 'Weighted RMS error = %.3f*R;'%(rmsWei)+rmsStr

    #print a warning if the rms fit is worse that 0.25*R
    if(rmsUnw > 0.25 or rmsWei > 0.25):
        logging.warning("Poor Wilhoit-to-NASA fit quality: RMS error = %.3f*R" % (rmsWei if weighting == 1 else rmsUnw))

    #restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
    Tint *= 1000.
    Tmin *= 1000.
    Tmax *= 1000.

    nasa_low.c1 /= 1000.
    nasa_low.c2 /= 1000000.
    nasa_low.c3 /= 1000000000.
    nasa_low.c4 /= 1000000000000.

    nasa_high.c1 /= 1000.
    nasa_high.c2 /= 1000000.
    nasa_high.c3 /= 1000000000.
    nasa_high.c4 /= 1000000000000.

    # output comment
    comment = 'NASA function fitted to Wilhoit function with B = ' + wilhoit.B + '. ' + rmsStr + wilhoit.comment
    nasa_low.Tmin = Quantity(Tmin,"K"); nasa_low.Tmax = Quantity(Tint,"K")
    nasa_low.comment = 'Low temperature range polynomial'
    nasa_high.Tmin = Quantity(Tint,"K"); nasa_high.Tmax = Quantity(Tmax,"K")
    nasa_high.comment = 'High temperature range polynomial'

    #for the low polynomial, we want the results to match the Wilhoit value at 298.15K
    #low polynomial enthalpy:
    Hlow = (wilhoit.getEnthalpy(298) - nasa_low.getEnthalpy(298))/constants.R
    #low polynomial entropy:
    Slow = (wilhoit.getEntropy(298) - nasa_low.getEntropy(298))/constants.R

    # update last two coefficients
    nasa_low.c5 = Hlow
    nasa_low.c6 = Slow

    #for the high polynomial, we want the results to match the low polynomial value at tint
    #high polynomial enthalpy:
    Hhigh = (nasa_low.getEnthalpy(Tint) - nasa_high.getEnthalpy(Tint))/constants.R
    #high polynomial entropy:
    Shigh = (nasa_low.getEntropy(Tint) - nasa_high.getEntropy(Tint))/constants.R

    # update last two coefficients
    #polynomial_high.coeffs = (b6,b7,b8,b9,b10,Hhigh,Shigh)
    nasa_high.c5 = Hhigh
    nasa_high.c6 = Shigh

    return MultiNASA(Tmin=Tmin, Tmax=Tmax, polynomials=[nasa_low,nasa_high], comment=comment)

def Wilhoit2NASA(wilhoit, tmin, tmax, tint, weighting, contCons):
    """
    input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3,
           Tmin (minimum temperature (in kiloKelvin),
           Tmax (maximum temperature (in kiloKelvin),
           Tint (intermediate temperature, in kiloKelvin)
           weighting (boolean: should the fit be weighted by 1/T?)
           contCons: a measure of the continutity constraints on the fitted NASA polynomials; possible values are:
            5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
            4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at tint
            3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at tint
            2: constrain Cp and dCp/dT to be continuous at tint
            1: constrain Cp to be continous at tint
            0: no constraints on continuity of Cp(T) at tint
            note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at tint by the form of the Cp(T) function
    output: NASA polynomials (nasa_low, nasa_high) with scaled parameters
    """
    #construct (typically 13*13) symmetric A matrix (in A*x = b); other elements will be zero
    A = numpy.zeros([10+contCons,10+contCons])
    b = numpy.zeros([10+contCons])

    if weighting:
        A[0,0] = 2*math.log(tint/tmin)
        A[0,1] = 2*(tint - tmin)
        A[0,2] = tint*tint - tmin*tmin
        A[0,3] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
        A[0,4] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
        A[1,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
        A[2,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
        A[3,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
        A[4,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
    else:
        A[0,0] = 2*(tint - tmin)
        A[0,1] = tint*tint - tmin*tmin
        A[0,2] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
        A[0,3] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
        A[0,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
        A[1,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
        A[2,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
        A[3,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
        A[4,4] = 2.*(tint*tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/9
    A[1,1] = A[0,2]
    A[1,2] = A[0,3]
    A[1,3] = A[0,4]
    A[2,2] = A[0,4]
    A[2,3] = A[1,4]
    A[3,3] = A[2,4]

    if weighting:
        A[5,5] = 2*math.log(tmax/tint)
        A[5,6] = 2*(tmax - tint)
        A[5,7] = tmax*tmax - tint*tint
        A[5,8] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
        A[5,9] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
        A[6,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
        A[7,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
        A[8,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
        A[9,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
    else:
        A[5,5] = 2*(tmax - tint)
        A[5,6] = tmax*tmax - tint*tint
        A[5,7] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
        A[5,8] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
        A[5,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
        A[6,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
        A[7,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
        A[8,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
        A[9,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint*tint)/9
    A[6,6] = A[5,7]
    A[6,7] = A[5,8]
    A[6,8] = A[5,9]
    A[7,7] = A[5,9]
    A[7,8] = A[6,9]
    A[8,8] = A[7,9]

    if(contCons > 0):#set non-zero elements in the 11th column for Cp(T) continuity contraint
        A[0,10] = 1.
        A[1,10] = tint
        A[2,10] = tint*tint
        A[3,10] = A[2,10]*tint
        A[4,10] = A[3,10]*tint
        A[5,10] = -A[0,10]
        A[6,10] = -A[1,10]
        A[7,10] = -A[2,10]
        A[8,10] = -A[3,10]
        A[9,10] = -A[4,10]
        if(contCons > 1): #set non-zero elements in the 12th column for dCp/dT continuity constraint
            A[1,11] = 1.
            A[2,11] = 2*tint
            A[3,11] = 3*A[2,10]
            A[4,11] = 4*A[3,10]
            A[6,11] = -A[1,11]
            A[7,11] = -A[2,11]
            A[8,11] = -A[3,11]
            A[9,11] = -A[4,11]
            if(contCons > 2): #set non-zero elements in the 13th column for d2Cp/dT2 continuity constraint
                A[2,12] = 2.
                A[3,12] = 6*tint
                A[4,12] = 12*A[2,10]
                A[7,12] = -A[2,12]
                A[8,12] = -A[3,12]
                A[9,12] = -A[4,12]
                if(contCons > 3): #set non-zero elements in the 14th column for d3Cp/dT3 continuity constraint
                    A[3,13] = 6
                    A[4,13] = 24*tint
                    A[8,13] = -A[3,13]
                    A[9,13] = -A[4,13]
                    if(contCons > 4): #set non-zero elements in the 15th column for d4Cp/dT4 continuity constraint
                        A[4,14] = 24
                        A[9,14] = -A[4,14]

    # make the matrix symmetric
    for i in range(1,10+contCons):
        for j in range(0, i):
            A[i,j] = A[j,i]

    #construct b vector
    w0int = Wilhoit_integral_T0(wilhoit, tint)
    w1int = Wilhoit_integral_T1(wilhoit, tint)
    w2int = Wilhoit_integral_T2(wilhoit, tint)
    w3int = Wilhoit_integral_T3(wilhoit, tint)
    w0min = Wilhoit_integral_T0(wilhoit, tmin)
    w1min = Wilhoit_integral_T1(wilhoit, tmin)
    w2min = Wilhoit_integral_T2(wilhoit, tmin)
    w3min = Wilhoit_integral_T3(wilhoit, tmin)
    w0max = Wilhoit_integral_T0(wilhoit, tmax)
    w1max = Wilhoit_integral_T1(wilhoit, tmax)
    w2max = Wilhoit_integral_T2(wilhoit, tmax)
    w3max = Wilhoit_integral_T3(wilhoit, tmax)
    if weighting:
        wM1int = Wilhoit_integral_TM1(wilhoit, tint)
        wM1min = Wilhoit_integral_TM1(wilhoit, tmin)
        wM1max = Wilhoit_integral_TM1(wilhoit, tmax)
    else:
        w4int = Wilhoit_integral_T4(wilhoit, tint)
        w4min = Wilhoit_integral_T4(wilhoit, tmin)
        w4max = Wilhoit_integral_T4(wilhoit, tmax)

    if weighting:
        b[0] = 2*(wM1int - wM1min)
        b[1] = 2*(w0int - w0min)
        b[2] = 2*(w1int - w1min)
        b[3] = 2*(w2int - w2min)
        b[4] = 2*(w3int - w3min)
        b[5] = 2*(wM1max - wM1int)
        b[6] = 2*(w0max - w0int)
        b[7] = 2*(w1max - w1int)
        b[8] = 2*(w2max - w2int)
        b[9] = 2*(w3max - w3int)
    else:
        b[0] = 2*(w0int - w0min)
        b[1] = 2*(w1int - w1min)
        b[2] = 2*(w2int - w2min)
        b[3] = 2*(w3int - w3min)
        b[4] = 2*(w4int - w4min)
        b[5] = 2*(w0max - w0int)
        b[6] = 2*(w1max - w1int)
        b[7] = 2*(w2max - w2int)
        b[8] = 2*(w3max - w3int)
        b[9] = 2*(w4max - w4int)

    # solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A
    # matrix is not required; not including it should give same result, except
    # Lagrange multipliers will differ by a factor of two)
    import scipy.linalg
    x = scipy.linalg.solve(A,b,overwrite_a=1,overwrite_b=1)

    nasa_low = NASA(Tmin=0, Tmax=0, coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
    nasa_high = NASA(Tmin=0, Tmax=0, coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

    return nasa_low, nasa_high

def Wilhoit2NASA_TintOpt(wilhoit, tmin, tmax, weighting, contCons):
    #input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
    #output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
    #1. vary Tint, bounded by tmin and tmax, to minimize TintOpt_objFun
    #cf. http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html and http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html#scipy.optimize.fminbound)
    import scipy.optimize
    tint = scipy.optimize.fminbound(TintOpt_objFun, tmin, tmax, args=(wilhoit, tmin, tmax, weighting, contCons))
    tint = float(tint) # fminbound returns a numpy.ndarray object
    #note that we have not used any guess when using this minimization routine
    #2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
    (nasa1, nasa2) = Wilhoit2NASA(wilhoit, tmin, tmax, tint, weighting, contCons)
    return nasa1, nasa2, tint

def TintOpt_objFun(tint, wilhoit, tmin, tmax, weighting, contCons):
    #input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
    #output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
    if (weighting == 1):
        result = TintOpt_objFun_W(tint, wilhoit, tmin, tmax, contCons)
    else:
        result = TintOpt_objFun_NW(tint, wilhoit, tmin, tmax, contCons)

    # numerical errors could accumulate to give a slightly negative result
    # this is unphysical (it's the integral of a *squared* error) so we
    # set it to zero to avoid later problems when we try find the square root.
    if result < 0:
        if result<-1E-13:
            logging.error("Greg thought he fixed the numerical problem, but apparently it is still an issue; please e-mail him with the following results:")
            logging.error(tint)
            logging.error(wilhoit)
            logging.error(tmin)
            logging.error(tmax)
            logging.error(weighting)
            logging.error(result)
        logging.info("Negative ISE of %f reset to zero."%(result))
        result = 0

    return result

def TintOpt_objFun_NW(tint, wilhoit, tmin, tmax, contCons):
    """
    Evaluate the objective function - the integral of the square of the error in the fit.

    input: Tint (intermediate temperature, in kiloKelvin)
            Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3,
            Tmin (minimum temperature (in kiloKelvin),
            Tmax (maximum temperature (in kiloKelvin)
    output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
    """
    nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 0, contCons)
    b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
    b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

    q0=Wilhoit_integral_T0(wilhoit, tint)
    q1=Wilhoit_integral_T1(wilhoit, tint)
    q2=Wilhoit_integral_T2(wilhoit, tint)
    q3=Wilhoit_integral_T3(wilhoit, tint)
    q4=Wilhoit_integral_T4(wilhoit, tint)
    result = (Wilhoit_integral2_T0(wilhoit, tmax) - Wilhoit_integral2_T0(wilhoit, tmin) +
                 NASA_integral2_T0(nasa_low, tint) - NASA_integral2_T0(nasa_low, tmin) +
                 NASA_integral2_T0(nasa_high, tmax) - NASA_integral2_T0(nasa_high, tint)
                 - 2* (b6*(Wilhoit_integral_T0(wilhoit, tmax)-q0)+b1*(q0-Wilhoit_integral_T0(wilhoit, tmin))
                 +b7*(Wilhoit_integral_T1(wilhoit, tmax) - q1) +b2*(q1 - Wilhoit_integral_T1(wilhoit, tmin))
                 +b8*(Wilhoit_integral_T2(wilhoit, tmax) - q2) +b3*(q2 - Wilhoit_integral_T2(wilhoit, tmin))
                 +b9*(Wilhoit_integral_T3(wilhoit, tmax) - q3) +b4*(q3 - Wilhoit_integral_T3(wilhoit, tmin))
                 +b10*(Wilhoit_integral_T4(wilhoit, tmax) - q4)+b5*(q4 - Wilhoit_integral_T4(wilhoit, tmin))))

    return result

def TintOpt_objFun_W(tint, wilhoit, tmin, tmax, contCons):
    """
    Evaluate the objective function - the integral of the square of the error in the fit.

    If fit is close to perfect, result may be slightly negative due to numerical errors in evaluating this integral.
    input: Tint (intermediate temperature, in kiloKelvin)
            Wilhoit parameters: Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3,
            Tmin (minimum temperature (in kiloKelvin),
            Tmax (maximum temperature (in kiloKelvin)
    output: the quantity Integrate[1/t*(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
    """
    nasa_low, nasa_high = Wilhoit2NASA(wilhoit,tmin,tmax,tint, 1, contCons)
    b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
    b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

    qM1=Wilhoit_integral_TM1(wilhoit, tint)
    q0=Wilhoit_integral_T0(wilhoit, tint)
    q1=Wilhoit_integral_T1(wilhoit, tint)
    q2=Wilhoit_integral_T2(wilhoit, tint)
    q3=Wilhoit_integral_T3(wilhoit, tint)
    result = (Wilhoit_integral2_TM1(wilhoit, tmax) - Wilhoit_integral2_TM1(wilhoit, tmin) +
                 NASA_integral2_TM1(nasa_low, tint) - NASA_integral2_TM1(nasa_low, tmin) +
                 NASA_integral2_TM1(nasa_high, tmax) - NASA_integral2_TM1(nasa_high, tint)
                 - 2* (b6*(Wilhoit_integral_TM1(wilhoit, tmax)-qM1)+b1*(qM1 - Wilhoit_integral_TM1(wilhoit, tmin))
                 +b7*(Wilhoit_integral_T0(wilhoit, tmax)-q0)+b2*(q0 - Wilhoit_integral_T0(wilhoit, tmin))
                 +b8*(Wilhoit_integral_T1(wilhoit, tmax)-q1)+b3*(q1 - Wilhoit_integral_T1(wilhoit, tmin))
                 +b9*(Wilhoit_integral_T2(wilhoit, tmax)-q2)+b4*(q2 - Wilhoit_integral_T2(wilhoit, tmin))
                 +b10*(Wilhoit_integral_T3(wilhoit, tmax)-q3)+b5*(q3 - Wilhoit_integral_T3(wilhoit, tmin))))

    return result

################################################################################

#a faster version of the integral based on H from Yelvington's thesis; it differs from the original (see above) by a constant (dependent on parameters but independent of t)
def Wilhoit_integral_T0(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(y=cython.double, y2=cython.double, logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    y = t/(t+B)
    y2 = y*y
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = cp0*t - (cpInf-cp0)*t*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplust))
    return result

#a faster version of the integral based on S from Yelvington's thesis; it differs from the original by a constant (dependent on parameters but independent of t)
def Wilhoit_integral_TM1(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^-1, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(y=cython.double, logt=cython.double, logy=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    y = t/(t+B)
    if cython.compiled:
        logy = log(y); logt = log(t)
    else:
        logy = math.log(y); logt = math.log(t)
    result = cpInf*logt-(cpInf-cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
    return result

def Wilhoit_integral_T1(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t + (cpInf*t**2)/2. + (a3*B**7*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 6*a3)*B**6*(cp0 - cpInf))/(4.*(B + t)**4) -
        ((a1 + 5*(a2 + 3*a3))*B**5*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(cp0 - cpInf))/(2.*(B + t)**2) -
        ((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(cp0 - cpInf))/(B + t) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T2(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^2, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*t) + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**2)/2. + (cpInf*t**3)/3. + (a3*B**8*(cp0 - cpInf))/(5.*(B + t)**5) -
        ((a2 + 7*a3)*B**7*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(cp0 - cpInf))/(2.*(B + t)**2) +
        ((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(cp0 - cpInf))/(B + t) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T3(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^3, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**3)/3. +
        (cpInf*t**4)/4. + (a3*B**9*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 8*a3)*B**8*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(cp0 - cpInf))/(3.*(B + t)**3) +
        ((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(cp0 - cpInf))/(B + t) -
        (5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral_T4(wilhoit, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^4, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*t) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t**2)/2. +
        ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**4)/4. + (cpInf*t**5)/5. + (a3*B**10*(cp0 - cpInf))/(5.*(B + t)**5) -
        ((a2 + 9*a3)*B**9*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(cp0 - cpInf))/(2.*(B + t)**2) +
        ((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(cp0 - cpInf))/(B + t) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(cp0 - cpInf)*logBplust)
    return result

def Wilhoit_integral2_T0(wilhoit, t):
    #output: the quantity Integrate[(Cp(Wilhoit)/R)^2, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t)
    else:
        logBplust = math.log(B + t)
    result = (cpInf**2*t - (a3**2*B**12*(cp0 - cpInf)**2)/(11.*(B + t)**11) + (a3*(a2 + 5*a3)*B**11*(cp0 - cpInf)**2)/(5.*(B + t)**10) -
        ((a2**2 + 18*a2*a3 + a3*(2*a1 + 45*a3))*B**10*(cp0 - cpInf)**2)/(9.*(B + t)**9) + ((4*a2**2 + 36*a2*a3 + a1*(a2 + 8*a3) + a3*(a0 + 60*a3))*B**9*(cp0 - cpInf)**2)/(4.*(B + t)**8) -
        ((a1**2 + 14*a1*(a2 + 4*a3) + 2*(14*a2**2 + a3 + 84*a2*a3 + 105*a3**2 + a0*(a2 + 7*a3)))*B**8*(cp0 - cpInf)**2)/(7.*(B + t)**7) +
        ((3*a1**2 + a2 + 28*a2**2 + 7*a3 + 126*a2*a3 + 126*a3**2 + 7*a1*(3*a2 + 8*a3) + a0*(a1 + 6*a2 + 21*a3))*B**7*(cp0 - cpInf)**2)/(3.*(B + t)**6) -
        (B**6*(cp0 - cpInf)*(a0**2*(cp0 - cpInf) + 15*a1**2*(cp0 - cpInf) + 10*a0*(a1 + 3*a2 + 7*a3)*(cp0 - cpInf) + 2*a1*(1 + 35*a2 + 70*a3)*(cp0 - cpInf) +
         2*(35*a2**2*(cp0 - cpInf) + 6*a2*(1 + 21*a3)*(cp0 - cpInf) + a3*(5*(4 + 21*a3)*cp0 - 21*(cpInf + 5*a3*cpInf)))))/(5.*(B + t)**5) +
        (B**5*(cp0 - cpInf)*(14*a2*cp0 + 28*a2**2*cp0 + 30*a3*cp0 + 84*a2*a3*cp0 + 60*a3**2*cp0 + 2*a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) +
         a0*(1 + 10*a1 + 20*a2 + 35*a3)*(cp0 - cpInf) + a1*(5 + 35*a2 + 56*a3)*(cp0 - cpInf) - 15*a2*cpInf - 28*a2**2*cpInf - 35*a3*cpInf - 84*a2*a3*cpInf - 60*a3**2*cpInf))/
         (2.*(B + t)**4) - (B**4*(cp0 - cpInf)*((1 + 6*a0**2 + 15*a1**2 + 32*a2 + 28*a2**2 + 50*a3 + 72*a2*a3 + 45*a3**2 + 2*a1*(9 + 21*a2 + 28*a3) + a0*(8 + 20*a1 + 30*a2 + 42*a3))*cp0 -
         (1 + 6*a0**2 + 15*a1**2 + 40*a2 + 28*a2**2 + 70*a3 + 72*a2*a3 + 45*a3**2 + a0*(8 + 20*a1 + 30*a2 + 42*a3) + a1*(20 + 42*a2 + 56*a3))*cpInf))/(3.*(B + t)**3) +
        (B**3*(cp0 - cpInf)*((2 + 2*a0**2 + 3*a1**2 + 9*a2 + 4*a2**2 + 11*a3 + 9*a2*a3 + 5*a3**2 + a0*(5 + 5*a1 + 6*a2 + 7*a3) + a1*(7 + 7*a2 + 8*a3))*cp0 -
         (2 + 2*a0**2 + 3*a1**2 + 15*a2 + 4*a2**2 + 21*a3 + 9*a2*a3 + 5*a3**2 + a0*(6 + 5*a1 + 6*a2 + 7*a3) + a1*(10 + 7*a2 + 8*a3))*cpInf))/(B + t)**2 -
        (B**2*((2 + a0 + a1 + a2 + a3)**2*cp0**2 - 2*(5 + a0**2 + a1**2 + 8*a2 + a2**2 + 9*a3 + 2*a2*a3 + a3**2 + 2*a0*(3 + a1 + a2 + a3) + a1*(7 + 2*a2 + 2*a3))*cp0*cpInf +
         (6 + a0**2 + a1**2 + 12*a2 + a2**2 + 14*a3 + 2*a2*a3 + a3**2 + 2*a1*(5 + a2 + a3) + 2*a0*(4 + a1 + a2 + a3))*cpInf**2))/(B + t) +
        2*(2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*cpInf*logBplust)
    return result

def Wilhoit_integral2_TM1(wilhoit, t):
    #output: the quantity Integrate[(Cp(Wilhoit)/R)^2*t^-1, t'] evaluated at t'=t
    cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double, a2=cython.double, a3=cython.double)
    cython.declare(logBplust=cython.double, logt=cython.double, result=cython.double)
    cp0, cpInf, B, a0, a1, a2, a3 = wilhoit.cp0.value, wilhoit.cpInf.value, wilhoit.B.value, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3
    if cython.compiled:
        logBplust = log(B + t); logt = log(t)
    else:
        logBplust = math.log(B + t); logt = math.log(t)
    result = ( (a3**2*B**11*(cp0 - cpInf)**2)/(11.*(B + t)**11) - (a3*(2*a2 + 9*a3)*B**10*(cp0 - cpInf)**2)/(10.*(B + t)**10) +
        ((a2**2 + 16*a2*a3 + 2*a3*(a1 + 18*a3))*B**9*(cp0 - cpInf)**2)/(9.*(B + t)**9) -
        ((7*a2**2 + 56*a2*a3 + 2*a1*(a2 + 7*a3) + 2*a3*(a0 + 42*a3))*B**8*(cp0 - cpInf)**2)/(8.*(B + t)**8) +
        ((a1**2 + 21*a2**2 + 2*a3 + 112*a2*a3 + 126*a3**2 + 2*a0*(a2 + 6*a3) + 6*a1*(2*a2 + 7*a3))*B**7*(cp0 - cpInf)**2)/(7.*(B + t)**7) -
        ((5*a1**2 + 2*a2 + 30*a1*a2 + 35*a2**2 + 12*a3 + 70*a1*a3 + 140*a2*a3 + 126*a3**2 + 2*a0*(a1 + 5*(a2 + 3*a3)))*B**6*(cp0 - cpInf)**2)/(6.*(B + t)**6) +
        (B**5*(cp0 - cpInf)*(10*a2*cp0 + 35*a2**2*cp0 + 28*a3*cp0 + 112*a2*a3*cp0 + 84*a3**2*cp0 + a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) + 2*a1*(1 + 20*a2 + 35*a3)*(cp0 - cpInf) +
        4*a0*(2*a1 + 5*(a2 + 2*a3))*(cp0 - cpInf) - 10*a2*cpInf - 35*a2**2*cpInf - 30*a3*cpInf - 112*a2*a3*cpInf - 84*a3**2*cpInf))/(5.*(B + t)**5) -
        (B**4*(cp0 - cpInf)*(18*a2*cp0 + 21*a2**2*cp0 + 32*a3*cp0 + 56*a2*a3*cp0 + 36*a3**2*cp0 + 3*a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) +
        2*a0*(1 + 6*a1 + 10*a2 + 15*a3)*(cp0 - cpInf) + 2*a1*(4 + 15*a2 + 21*a3)*(cp0 - cpInf) - 20*a2*cpInf - 21*a2**2*cpInf - 40*a3*cpInf - 56*a2*a3*cpInf - 36*a3**2*cpInf))/
        (4.*(B + t)**4) + (B**3*(cp0 - cpInf)*((1 + 3*a0**2 + 5*a1**2 + 14*a2 + 7*a2**2 + 18*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(5 + 6*a2 + 7*a3))*cp0 -
        (1 + 3*a0**2 + 5*a1**2 + 20*a2 + 7*a2**2 + 30*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(6 + 6*a2 + 7*a3))*cpInf))/(3.*(B + t)**3) -
        (B**2*((3 + a0**2 + a1**2 + 4*a2 + a2**2 + 4*a3 + 2*a2*a3 + a3**2 + 2*a1*(2 + a2 + a3) + 2*a0*(2 + a1 + a2 + a3))*cp0**2 -
        2*(3 + a0**2 + a1**2 + 7*a2 + a2**2 + 8*a3 + 2*a2*a3 + a3**2 + 2*a1*(3 + a2 + a3) + a0*(5 + 2*a1 + 2*a2 + 2*a3))*cp0*cpInf +
        (3 + a0**2 + a1**2 + 10*a2 + a2**2 + 12*a3 + 2*a2*a3 + a3**2 + 2*a1*(4 + a2 + a3) + 2*a0*(3 + a1 + a2 + a3))*cpInf**2))/(2.*(B + t)**2) +
        (B*(cp0 - cpInf)*(cp0 - (3 + 2*a0 + 2*a1 + 2*a2 + 2*a3)*cpInf))/(B + t) + cp0**2*logt + (-cp0**2 + cpInf**2)*logBplust)
    return result

################################################################################

def NASA_integral2_T0(polynomial, T):
    #output: the quantity Integrate[(Cp(NASA)/R)^2, t'] evaluated at t'=t
    cython.declare(c0=cython.double, c1=cython.double, c2=cython.double, c3=cython.double, c4=cython.double)
    cython.declare(T2=cython.double, T4=cython.double, T8=cython.double)
    c0, c1, c2, c3, c4 = polynomial.c0, polynomial.c1, polynomial.c2, polynomial.c3, polynomial.c4
    T2=T*T; T4=T2*T2; T8=T4*T4
    result = (
        c0*c0*T + c0*c1*T2 + 2./3.*c0*c2*T2*T + 0.5*c0*c3*T4 + 0.4*c0*c4*T4*T +
        c1*c1*T2*T/3. + 0.5*c1*c2*T4 + 0.4*c1*c3*T4*T + c1*c4*T4*T2/3. +
        0.2*c2*c2*T4*T + c2*c3*T4*T2/3. + 2./7.*c2*c4*T4*T2*T +
        c3*c3*T4*T2*T/7. + 0.25*c3*c4*T8 +
        c4*c4*T8*T/9.
    )
    return result

def NASA_integral2_TM1(polynomial, T):
    #output: the quantity Integrate[(Cp(NASA)/R)^2*t^-1, t'] evaluated at t'=t
    cython.declare(c0=cython.double, c1=cython.double, c2=cython.double, c3=cython.double, c4=cython.double)
    cython.declare(T2=cython.double, T4=cython.double, logT=cython.double)
    c0, c1, c2, c3, c4 = polynomial.c0, polynomial.c1, polynomial.c2, polynomial.c3, polynomial.c4
    T2=T*T; T4=T2*T2
    if cython.compiled:
        logT = log(T)
    else:
        logT = math.log(T)
    result = (
        c0*c0*logT + 2*c0*c1*T + c0*c2*T2 + 2./3.*c0*c3*T2*T + 0.5*c0*c4*T4 +
        0.5*c1*c1*T2 + 2./3.*c1*c2*T2*T + 0.5*c1*c3*T4 + 0.4*c1*c4*T4*T +
        0.25*c2*c2*T4 + 0.4*c2*c3*T4*T + c2*c4*T4*T2/3. +
        c3*c3*T4*T2/6. + 2./7.*c3*c4*T4*T2*T +
        c4*c4*T4*T4/8.
    )
    return result
