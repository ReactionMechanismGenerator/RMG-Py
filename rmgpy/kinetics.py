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
This module contains classes and methods for working with kinetics models. All 
such models derive from the :class:`KineticsModel` base class:

* :class:`Arrhenius` - A kinetics model based on the modified Arrhenius
  equation

* :class:`ArrheniusEP` - A kinetics model based on the modified Arrhenius
  equation with Evans-Polanyi correction to the activation energy

* :class:`MultiKinetics` - A kinetics model based on a summation of
  several other kinetics expressions

* :class:`ThirdBody` - A pressure-dependent kinetics model based on the  
  modified Arrhenius equation, but with an additional factor for the third body 
  concentration

* :class:`Lindemann` - A pressure-dependent kinetics model based on the
  Lindemann equation

* :class:`Troe` - A pressure-dependent kinetics model based on the
  Lindemann equation with improved Troe falloff factor

* :class:`PDepArrhenius` - A pressure-dependent kinetics model based on a
  set of modified Arrhenius equations at various pressures, which are then
  interpolated between on a logarithmic pressure scale

* :class:`Chebyshev` - A pressure-dependent kinetics model using an array
  of Chebyshev polynomials in inverse temperature and logarithmic pressure

"""

################################################################################

import math
import numpy
import numpy.linalg
import cython

from quantity import Quantity, constants
from molecule import Molecule

################################################################################

class KineticsError(Exception):
    """
    An exception to be raised when an error occurs while working with 
    kinetics models and data. Pass a string describing the circumstances of the
    exceptional behavior.
    """
    pass

################################################################################

class KineticsModel:
    """
    A base class for kinetics models, containing several attributes common to 
    all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          :class:`Quantity`   The minimum absolute temperature in K at which the model is valid
    `Tmax`          :class:`Quantity`   The maximum absolute temperature in K at which the model is valid
    `Pmin`          :class:`Quantity`   The minimum absolute pressure in Pa at which the model is valid
    `Pmax`          :class:`Quantity`   The maximum absolute pressure in Pa at which the model is valid
    `comment`       :class:`str`        A string containing information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Quantity(Tmin)
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Quantity(Tmax)
        else:
            self.Tmax = None
        if Pmin is not None:
            self.Pmin = Quantity(Pmin)
        else:
            self.Pmin = None
        if Pmax is not None:
            self.Pmax = Quantity(Pmax)
        else:
            self.Pmax = None
        self.comment = comment

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        return 'KineticsModel(Tmin={0!r}, Tmax={1!r}, Pmin={2!r}, Pmax={3!r}, comment="""{4}""")'.format(self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsModel object.
        """
        return (KineticsModel, (self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isTemperatureValid(self, T):
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the kinetic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return self.Tmin is None or self.Tmax is None or (self.Tmin.value <= T and T <= self.Tmax.value)

    def isPressureValid(self, P):
        """
        Return ``True`` if the pressure `P` in Pa is within the valid
        pressure range of the kinetic data, or ``False`` if not. If
        the minimum and maximum pressure are not defined, ``True`` is 
        returned.
        """
        return self.Pmin is None or self.Pmax is None or (self.Pmin.value <= P and P <= self.Pmax.value)

    def isPressureDependent(self):
        """
        Return ``True`` if the kinetics are pressure-dependent or ``False`` if
        they are pressure-independent. This method must be overloaded in the 
        derived class.
        """
        raise KineticsError('Unexpected call to KineticsModel.isPressureDependent(); you should be using a class derived from KineticsModel.')

    def getRateCoefficients(self, Tlist):
        """
        Return the rate coefficients k(T) in SI units at temperatures
        `Tlist` in K.
        """
        return numpy.array([self.getRateCoefficient(T) for T in Tlist], numpy.float64)
    
    def toHTML(self):
        """
        Return an HTML rendering.
        """    
        cython.declare(T=cython.double, k=cython.double)
        cython.declare(string=str)
        cython.declare(kdata1=list, Tdata=list, kdata10=list)
        string = '<table class="KineticsData"><tr class="KineticsData_Tdata"><th>T/[K]</th>\n   '
        kdata1 = []
        kdata10 = []
        Tdata = [500,1000,1500,2000]
        for T in Tdata:
            string += '<td>{0:.0f}</td>'.format(T)
            kdata1.append(self.getRateCoefficient(T,P=1e5))
            if self.isPressureDependent():
                kdata10.append(self.getRateCoefficient(T,P=1e6))
                
        if self.isPressureDependent():
            string += '\n</tr><tr class="KineticsData_kdata"><th>log<sub>10</sub>(k(1 bar)/[mole,m,s])\n    '
            for k in kdata1:
                string += '<td>{0:+.1f}'.format(math.log10(k))
            string += '\n</tr><tr class="KineticsData_kdata"><th>log<sub>10</sub>(k(10 bar)/[mole,m,s])\n    '
            for k in kdata10:
                string += '<td>{0:+.1f}'.format(math.log10(k))
        else:
            string += '\n</tr><tr class="KineticsData_kdata"><th>log<sub>10</sub>(k/[mole,m,s])\n    '
            for k in kdata1:
                string += '<td>{0:+.1f}'.format(math.log10(k))
     
        string += '\n</tr></table>'
        string += "<span class='KineticsData_repr'>{0!r}</span>".format(self)
        return string
        

################################################################################

class KineticsData(KineticsModel):
    """
    A kinetics model based around a set of discrete (high-pressure limit)
    rate coefficients at various temperatures. The attributes are:

    =========== =================== ============================================
    Attribute   Type                Description
    =========== =================== ============================================
    `Tdata`     :class:`Quantity`   The temperatures at which the heat capacity data is provided
    `kdata`     :class:`Quantity`   The rate coefficients in SI units at each temperature in `Tdata`
    =========== =================== ============================================
    
    """

    def __init__(self, Tdata=None, kdata=None, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Tdata = Quantity(Tdata)
        self.kdata = Quantity(kdata)

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsData object.
        """
        string = 'KineticsData(Tdata={0!r}, kdata={1!r}'.format(self.Tdata, self.kdata)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsData object.
        """
        return (KineticsData, (self.Tdata, self.kdata, self.Tmin, self.Tmax, self.comment))
    
    def toHTML(self):
        """
        Return an HTML rendering.
        """
        cython.declare(T=cython.double, k=cython.double)
        cython.declare(string=str)
        string = '<table class="KineticsData"><tr class="KineticsData_Tdata"><th>T/[{0!s}]</th>\n   '.format(self.Tdata.units)
        for T in self.Tdata.values:
            string += '<td>{0:.0f}</td>'.format(T)
        string += '\n</tr><tr class="KineticsData_kdata"><th>log<sub>10</sub>(k/[{0!s}])\n    '.format(self.kdata.units)
        for k in self.kdata.values:
            string += '<td>{0:+.1f}'.format(math.log10(k))
        string += '\n</tr></table>'
        return string

    def isPressureDependent(self):
        """
        Returns ``False`` since KineticsData kinetics are not 
        pressure-dependent.
        """
        return False

    def getRateCoefficient(self, T, P=1e5):
        """
        Return the rate coefficient k(T) in SI units at temperature
        `T` in K.
        """
        cython.declare(Tmin=cython.double, Tmax=cython.double, kmin=cython.double, kmax=cython.double)
        cython.declare(k=cython.double, slope=cython.double)
        k = 0.0
        if not self.isTemperatureValid(T):
            raise KineticsError('Invalid temperature "{0:g} K" for heat capacity estimation.'.format(T))
        if T < numpy.min(self.Tdata.values):
            k = self.kdata.values[0]
        elif T >= numpy.max(self.Tdata.values):
            k = self.kdata.values[-1]
        else:
            for Tmin, Tmax, kmin, kmax in zip(self.Tdata.values[:-1], self.Tdata.values[1:], self.kdata.values[:-1], self.kdata.values[1:]):
                if Tmin <= T and T < Tmax:
                    slope = (1.0/T - 1.0/Tmin) / (1.0/Tmax - 1.0/Tmin)
                    k = kmin * (kmax / kmin)**slope
        return k
    
    def toArrhenius(self):
        """
        Return an :class:`Arrhenius` expression fitted to this data
        """
        return Arrhenius(comment=self.comment, Tmin=self.Tmin, Tmax=self.Tmax).fitToData(Tlist=self.Tdata.values, klist=self.kdata.values, kunits=self.kdata.units)

################################################################################

class Arrhenius(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics. The kinetic expression has
    the form

    .. math:: k(T) = A \\left( \\frac{T}{T_0} \\right)^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

    where :math:`A`, :math:`n`, :math:`E_\\mathrm{a}`, and :math:`T_0` are the
    parameters to be set, :math:`T` is absolute temperature, and :math:`R` is
    the gas law constant. The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `A`             :class:`Quantity`   The preexponential factor in s^-1, m^3/mol*s, etc.
    `T0`            :class:`Quantity`   The reference temperature in K
    `n`             :class:`Quantity`   The temperature exponent
    `Ea`            :class:`Quantity`   The activation energy in J/mol
    =============== =================== ========================================
    
    """
    
    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=1.0, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.A = Quantity(A)
        self.T0 = Quantity(T0)
        self.n = Quantity(n)
        self.Ea = Quantity(Ea)
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string
    
    def __str__(self):
        """
        Return a string representation that is a bit shorter and prettier than __repr__.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r})'.format(self.A, self.n, self.Ea, self.T0)
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``False`` since Arrhenius kinetics are not pressure-dependent.
        """
        return False

    def getRateCoefficient(self, T, P=1e5):
        """
        Return the rate coefficient k(T) in SI units at temperature 
        `T` in K.
        """
        return self.A.value * (T / self.T0.value)** self.n.value * math.exp(-self.Ea.value / constants.R / T)

    def changeT0(self, T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K, and
        adjusts the preexponential accordingly.
        """
        self.A.value /= (self.T0.value / T0)**self.n.value
        self.T0.value = T0

    def fitToData(self, Tlist, klist, kunits, T0=300):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in 
        K. A linear least-squares fit is used, which guarantees that the 
        resulting parameters provide the best possible approximation to the 
        data.
        """
        import numpy.linalg
        A = numpy.zeros((len(Tlist),3), numpy.float64)
        A[:,0] = numpy.ones_like(Tlist)
        A[:,1] = numpy.log(Tlist / T0)
        A[:,2] = -1.0 / constants.R / Tlist
        b = numpy.log(klist)
        x = numpy.linalg.lstsq(A,b)[0]
        
        self.A = Quantity(math.exp(x[0]), kunits)
        self.n = Quantity(x[1])
        self.Ea = Quantity((x[2], "J/mol"))
        self.T0 = Quantity((T0, "K"))
        return self
    
################################################################################

class ArrheniusEP(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics with Evans-Polanyi data. The
    kinetic expression has the form

    .. math:: k(T) = A T^n \\exp \\left( - \\frac{E_0 + \\alpha \\Delta H_\\mathrm{rxn}}{RT} \\right)

    The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `A`             :class:`Quantity`   The preexponential factor in s^-1, m^3/mol*s, etc.
    `n`             :class:`Quantity`   The temperature exponent
    `E0`            :class:`Quantity`   The activation energy at zero enthalpy of reaction in J/mol
    `alpha`         :class:`Quantity`   The linear dependence of activation energy on enthalpy of reaction
    =============== =================== ========================================

    """

    def __init__(self, A=0.0, n=0.0, alpha=0.0, E0=0.0, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.A = Quantity(A)
        self.n = Quantity(n)
        self.alpha = Quantity(alpha)
        self.E0 = Quantity(E0)

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusEP object.
        """
        string = 'ArrheniusEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``False`` since ArrheniusEP kinetics are not pressure-dependent.
        """
        return False

    def getActivationEnergy(self, dHrxn):
        """
        Return the activation energy in J/mol using the enthalpy of reaction 
        `dHrxn` in J/mol.
        """
        return self.E0.value + self.alpha.value * dHrxn
    
    def getRateCoefficient(self, T, dHrxn):
        """
        Return the rate coefficient k(T, P) in SI units at a 
        temperature `T` in K for a reaction having an enthalpy of reaction 
        `dHrxn` in J/mol.
        """
        Ea = cython.declare(cython.double)
        Ea = self.getActivationEnergy(dHrxn)
        return self.A.value * (T ** self.n.value) * math.exp(-Ea / constants.R / T)

    def toArrhenius(self, dHrxn):
        """
        Return an :class:`Arrhenius` object corresponding to this object
        by using the provided enthalpy of reaction `dHrxn` in J/mol to calculate
        the activation energy.
        """
        return Arrhenius(A=self.A, n=self.n, Ea=(self.getActivationEnergy(dHrxn),"J/mol"), T0=(1.0,"K"), Tmin=self.Tmin, Tmax=self.Tmax, comment=self.comment)

################################################################################

class MultiKinetics(KineticsModel):
    """
    Represent a rate coefficient as multiple sets of kinetics expressions

    .. math:: k(T,P) = \\sum_{i=1}^N k_i(T,P)

    The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `kineticsList` ``list``        A list of the :class:`KineticsModel` objects that sum to represent the kinetics
    =============== =============== ============================================

    """

    def __init__(self, kineticsList=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.kineticsList = kineticsList or []

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiKinetics object.
        """
        string = 'MultiKinetics(kineticsList=[{0}]'.format(', '.join([repr(k) for k in self.kineticsList]))
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a MultiKinetics object.
        """
        return (MultiKinetics, (self.kineticsList, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` if any of the multiple kinetics expressions are
        pressure-dependent.
        """
        return any([k.isPressureDependent() for k in self.kineticsList])

    def getRateCoefficient(self, T, P=1e5):
        """
        Return the rate coefficient k(T) in SI units at temperature
        `T` in K.
        """
        cython.declare(k=cython.double, kinetics=KineticsModel)
        k = 0.0
        for kin in self.kineticsList:
            k += kin.getRateCoefficient(T,P)
        return k

################################################################################

class PDepArrhenius(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = A(P) T^{n(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

    where the modified Arrhenius parameters are stored at a variety of pressures
    and interpolated between on a logarithmic scale. The attributes are:

    =============== ================== ============================================
    Attribute       Type               Description
    =============== ================== ============================================
    `pressures`     :class:`list`      The list of pressures in Pa
    `arrhenius`     :class:`list`      The list of :class:`Arrhenius` objects at each pressure
    `highPlimit`    :class:`Arrhenius` The high (infinite) pressure limiting :class:`Arrhenius` expression
    =============== ================== ============================================
    
    Note that `highPlimit` is not used in evaluating k(T,P).
    """

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.pressures = Quantity(pressures)
        self.arrhenius = arrhenius or []
        self.highPlimit = highPlimit or None

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        string = 'PDepArrhenius(\n pressures={0!r},\n arrhenius=[\n  {1}]'.format(self.pressures, ',\n  '.join([repr(arrh) for arrh in self.arrhenius]))
        if self.highPlimit is not None: string += ",\n highPlimit={0!r}".format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ',\n comment="""{0}"""'.format(self.comment)
        string += '\n)'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.highPlimit, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since PDepArrhenius kinetics are pressure-dependent.
        """
        return True

    def __getAdjacentExpressions(self, P):
        """
        Returns the pressures and Arrhenius expressions for the pressures that
        most closely bound the specified pressure `P` in Pa.
        """
        cython.declare(pressures=list)
        cython.declare(arrh=Arrhenius)
        cython.declare(i=cython.int, ilow=cython.int, ihigh=cython.int)

        pressures = [pressure for pressure in self.pressures.values]
        if P in pressures:
            arrh = self.arrhenius[pressures.index(P)]
            return P, P, arrh, arrh
        else:
            ilow = 0; ihigh = -1
            for i in range(1, len(pressures)):
                if pressures[i] <= P:
                    ilow = i
                if pressures[i] > P and ihigh == -1:
                    ihigh = i
            return pressures[ilow], pressures[ihigh], self.arrhenius[ilow], self.arrhenius[ihigh]
    
    def getRateCoefficient(self, T, P):
        """
        Return the rate constant k(T, P) in SI units at a temperature 
        `T` in K and pressure `P` in Pa by evaluating the pressure-
        dependent Arrhenius expression.
        
        If k(P+) and k(P-) (the values either side of the P requested) are zero,
        then zero is returned. If only k(P-)==0 then it is replaced with 
        k(P+)/1e10, and vice versa. This allows the logarithmic interpolation
        to procede without zero-division errors. (The expression is not defined when
        one of them is zero and the other is not, so we have to assume something.)
        """
        cython.declare(Plow=cython.double, Phigh=cython.double)
        cython.declare(alow=Arrhenius, ahigh=Arrhenius)
        cython.declare(j=cython.int, klist=cython.double, klow=cython.double, khigh=cython.double)
        
        k = 0.0
        Plow, Phigh, alow, ahigh = self.__getAdjacentExpressions(P)
        if Plow == Phigh:
            k = alow.getRateCoefficient(T)
        else:
            klow = alow.getRateCoefficient(T)
            khigh = ahigh.getRateCoefficient(T)
            if klow == khigh == 0.0: return 0.0
            if klow == 0: klow = khigh/1e10
            if khigh == 0: khigh = klow/1e10
            k = klow * 10**(math.log10(P/Plow)/math.log10(Phigh/Plow)*math.log10(khigh/klow))
        return k

    def fitToData(self, Tlist, Plist, K, kunits, T0=298.0):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` with units of `kunits` corresponding to a set of 
        temperatures `Tlist` in K and pressures `Plist` in Pa. An Arrhenius 
        model is fit at each pressure.
        """
        cython.declare(i=cython.int)
        self.pressures = Quantity(Plist/1e5,"bar")
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = Arrhenius().fitToData(Tlist, K[:,i], kunits, T0)
            self.arrhenius.append(arrhenius)
        return self

################################################################################

class Chebyshev(KineticsModel):
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
    `kunits`        ``str``         The units of the generated k(T, P) values
    `degreeT`       :class:`int`    The number of terms in the inverse temperature direction
    `degreeP`       :class:`int`    The number of terms in the log pressure direction
    =============== =============== ============================================
    
    """

    def __init__(self, coeffs=None, kunits='', Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        if coeffs is not None:
            self.coeffs = Quantity(numpy.array(coeffs, numpy.float64)).values
            self.degreeT = self.coeffs.shape[0]
            self.degreeP = self.coeffs.shape[1]
        else:
            self.coeffs = None
            self.degreeT = 0
            self.degreeP = 0
        self.kunits = kunits
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Chebyshev object.
        """
        coeffs = '['
        for i in range(self.degreeT):
            if i > 0: coeffs += ', '
            coeffs += '[{0}]'.format(','.join(['{0:g}'.format(self.coeffs[i,j]) for j in range(self.degreeP)]))
        coeffs += ']'
        
        string = 'Chebyshev(coeffs={0}'.format(coeffs)
        if self.kunits != '': string += ', kunits="{0}"'.format(self.kunits)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Chebyshev object.
        """
        return (Chebyshev, (self.coeffs, self.kunits, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since Chebyshev polynomial kinetics are 
        pressure-dependent.
        """
        return True

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
        """
        Return the reduced temperature corresponding to the given temperature
        `T` in K. This maps the inverse of the temperature onto the domain 
        [-1, 1] using the `Tmin` and `Tmax` attributes as the limits.
        """
        return (2.0/T - 1.0/self.Tmin.value - 1.0/self.Tmax.value) / (1.0/self.Tmax.value - 1.0/self.Tmin.value)
    
    def __getReducedPressure(self, P):
        """
        Return the reduced pressure corresponding to the given pressure
        `P` in Pa. This maps the logarithm of the pressure onto the domain 
        [-1, 1] using the `Pmin` and `Pmax` attributes as the limits.
        """
        if cython.compiled:
            return (2.0*log10(P) - log10(self.Pmin.value) - log10(self.Pmax.value)) / (log10(self.Pmax.value) - log10(self.Pmin.value))
        else:
            return (2.0*math.log(P) - math.log(self.Pmin.value) - math.log(self.Pmax.value)) / (math.log(self.Pmax.value) - math.log(self.Pmin.value))
    
    def getRateCoefficient(self, T, P):
        """
        Return the rate constant k(T, P) in SI units at a temperature 
        `T` in K and pressure `P` in Pa by evaluating the Chebyshev 
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

    def fitToData(self, Tlist, Plist, K, kunits, degreeT, degreeP, Tmin, Tmax, Pmin, Pmax):
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

        self.kunits = kunits; self.degreeT = degreeT; self.degreeP = degreeP

        # Set temperature and pressure ranges
        self.Tmin = Quantity(Tmin,"K"); self.Tmax = Quantity(Tmax,"K")
        self.Pmin = Quantity(Pmin/1e5,"bar"); self.Pmax = Quantity(Pmax/1e5,"bar")

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
        
        return self

################################################################################

class ThirdBody(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k(T) [\\ce{M}]

    where :math:`k(T)` is an Arrhenius expression and
    :math:`[\\ce{M}] \\approx P/RT` is the concentration of the third body
    (i.e. the bath gas). A collision efficiency can be used to further correct
    the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusHigh=None, efficiencies=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusHigh = arrheniusHigh
        self.efficiencies = {}
        if efficiencies is not None:
            for mol, eff in efficiencies.iteritems():
                if isinstance(mol, Molecule):
                    self.efficiencies[mol] = eff
                else:
                    self.efficiencies[Molecule().fromSMILES(mol)] = eff

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThirdBody object.
        """
        string = 'ThirdBody(arrheniusHigh={0!r}'.format(self.arrheniusHigh)
        molecules = [(molecule.toSMILES(), molecule) for molecule in self.efficiencies]
        molecules.sort()
        string += ', efficiencies={{{0}}}'.format(', '.join(['"{0}": {1:g}'.format(smiles, self.efficiencies[molecule]) for smiles, molecule in molecules]))
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ThirdBody object.
        """
        return (ThirdBody, (self.arrheniusHigh, self.efficiencies, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def isPressureDependent(self):
        """
        Returns ``True`` since third-body kinetics are pressure-dependent.
        """
        return True

    def getColliderEfficiency(self, collider):
        """
        Return the collider efficiency for the specified `collider`, which can
        take one of two forms:

        * A single collider species. If the collider exists in the in the set
          of efficiencies, its efficiency will be returned. If not, an
          efficiency of unity will be returned.

        * A ``dict`` mapping collider species to mole fractions. The overall
          efficiency will be a weighted sum of the efficiencies of the collider
          species, using the mole fractions as the weights. Collider species not
          present in the set of efficiencies will be assumed to have an
          efficiency of unity.

        If collider is ``None`` or otherwise invalid, an efficiency of unity
        will be returned.
        """
        if isinstance(collider, dict):
            # Assume collider is a dict mapping species to weights
            efficiency = 0.0
            for spec, frac in collider.iteritems():
                try:
                    eff = self.efficiencies[spec]
                except KeyError:
                    eff = 1.0
                efficiency += eff * frac
            efficiency /= sum(collider.values())
        else:
            # Assume collider is a single species
            try:
                efficiency = self.efficiencies[collider]
            except KeyError:
                efficiency = 1.0

        return efficiency

    def getRateCoefficient(self, T, P, collider=None):
        """
        Return the rate constant k(T, P) in SI units at a temperature
        `T` in K and pressure `P` in Pa by evaluating the third-body expression.
        If a `collider` is specified the rate coefficient will be modified
        accordingly.
        """
        cython.declare(C=cython.double, k=cython.double, efficiency=cython.double)
        C = P / constants.R / T # bath gas concentration in mol/m^3
        k = self.arrheniusHigh.getRateCoefficient(T)
        efficiency = self.getColliderEfficiency(collider)
        return efficiency * k * C

################################################################################

class Lindemann(ThirdBody):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_\\infty(T) \\left[ \\frac{P_\\mathrm{r}}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the
    bath gas. The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter. For
    the Lindemann model, :math:`F = 1`. A collision efficiency can be used to
    further correct the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        ThirdBody.__init__(self, arrheniusHigh=arrheniusHigh, efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusLow = arrheniusLow

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Lindemann object.
        """
        string = 'Lindemann(arrheniusHigh={0!r}, arrheniusLow={1!r}'.format(self.arrheniusHigh, self.arrheniusLow)
        molecules = [(molecule.toSMILES(), molecule) for molecule in self.efficiencies]
        molecules.sort()
        string += ', efficiencies={{{0}}}'.format(', '.join(['"{0}": {1:g}'.format(smiles, self.efficiencies[molecule]) for smiles, molecule in molecules]))
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Lindemann object.
        """
        return (Lindemann, (self.arrheniusLow, self.arrheniusHigh, self.efficiencies, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def getRateCoefficient(self, T, P, collider=None):
        """
        Return the rate constant k(T, P) in SI units at a temperature
        `T` in K and pressure `P` in Pa by evaluating the Lindemann expression.
        If a `collider` is specified the rate coefficient will be modified
        accordingly.
        """
        cython.declare(C=cython.double, k0=cython.double, kinf=cython.double,
            Pr=cython.double, F=cython.double, efficiency=cython.double)
        C = P / constants.R / T # bath gas concentration in mol/m^3
        k0 = self.arrheniusLow.getRateCoefficient(T)
        kinf = self.arrheniusHigh.getRateCoefficient(T)
        Pr = k0 * C / kinf
        F = 1.0
        efficiency = self.getColliderEfficiency(collider)
        return efficiency * kinf * (Pr / (1 + Pr)) * F

################################################################################

class Troe(Lindemann):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_\\infty(T) \\left[ \\frac{P_\\mathrm{r}}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the
    bath gas. The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter. A
    collision efficiency can be used to further correct the value of
    :math:`k(T,P)`.

    For the Troe model the parameter :math:`F` is computed via

    .. math::

        \\log F &= \\left\\{1 + \\left[ \\frac{\\log P_\\mathrm{r} + c}{n - d (\\log P_\\mathrm{r} + c)} \\right]^2 \\right\\}^{-1} \\log F_\\mathrm{cent}

        c &= -0.4 - 0.67 \\log F_\\mathrm{cent}

        n &= 0.75 - 1.27 \\log F_\\mathrm{cent}

        d &= 0.14

        F_\\mathrm{cent} &= (1 - \\alpha) \\exp \\left( -T/T_3 \\right) + \\alpha \\exp \\left( -T/T_1 \\right) + \\exp \\left( -T_2/T \\right)

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `alpha`         :class:`Quantity`       The :math:`\\alpha` parameter
    `T1`            :class:`Quantity`       The :math:`T_1` parameter
    `T2`            :class:`Quantity`       The :math:`T_2` parameter
    `T3`            :class:`Quantity`       The :math:`T_3` parameter
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None, alpha=0.0, T3=0.0, T1=0.0, T2=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        Lindemann.__init__(self, arrheniusLow=arrheniusLow, arrheniusHigh=arrheniusHigh, efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.alpha = Quantity(alpha)
        self.T3 = Quantity(T3)
        self.T1 = Quantity(T1)
        if T2 is not None:
            self.T2 = Quantity(T2)
        else:
            self.T2 = None
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Troe object.
        """
        string = 'Troe(arrheniusHigh={0!r}, arrheniusLow={1!r}'.format(self.arrheniusHigh, self.arrheniusLow)
        string += ', alpha={0!r}, T3={1!r}, T1={2!r}'.format(self.alpha, self.T3, self.T1)
        if self.T2 is not None: string += ', T2={0!r}'.format(self.T2)
        molecules = [(molecule.toSMILES(), molecule) for molecule in self.efficiencies]
        molecules.sort()
        string += ', efficiencies={{{0}}}'.format(', '.join(['"{0}": {1:g}'.format(smiles, self.efficiencies[molecule]) for smiles, molecule in molecules]))
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Troe object.
        """
        return (Troe, (self.arrheniusLow, self.arrheniusHigh, self.efficiencies, self.alpha, self.T3, self.T1, self.T2, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    def getRateCoefficient(self, T, P, collider=None):
        """
        Return the rate constant k(T, P) in SI units at a temperature
        `T` in K and pressure `P` in Pa by evaluating the Lindemann expression.
        If a `collider` is specified the rate coefficient will be modified
        accordingly.
        """
        cython.declare(C=cython.double, k0=cython.double, kinf=cython.double,
            Pr=cython.double, F=cython.double, efficiency=cython.double,
            d=cython.double, n=cython.double, c=cython.double, Fcent=cython.double)
        
        C = P / constants.R / T # bath gas concentration in mol/m^3
        k0 = self.arrheniusLow.getRateCoefficient(T)
        kinf = self.arrheniusHigh.getRateCoefficient(T)
        Pr = k0 * C / kinf
        efficiency = self.getColliderEfficiency(collider)
        
        Fcent = (1 - self.alpha.value) * math.exp(-T / self.T3.value) + self.alpha.value * math.exp(-T / self.T1.value)
        if self.T2 is not None: Fcent += math.exp(-self.T2.value / T)
        d = 0.14
        n = 0.75 - 1.27 * math.log10(Fcent)
        c = -0.4 - 0.67 * math.log10(Fcent)
        F = 10.0**(math.log10(Fcent)/(1 + ((math.log10(Pr) + c)/(n - d * (math.log10(Pr))))**2))

        return efficiency * kinf * (Pr / (1 + Pr)) * F
