# cython: embedsignature=True, cdivision=True

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

import numpy
from libc.math cimport exp, log, sqrt, log10

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity

################################################################################

cdef class Arrhenius(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0,"K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K. 
        """
        cdef double A, n, Ea, T0
        A = self._A.value_si
        n = self._n.value_si
        Ea = self._Ea.value_si
        T0 = self._T0.value_si
        return A * (T / T0)**n * exp(-Ea / (constants.R * T))

    cpdef changeT0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K, 
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0)**self._n.value_si
        self._T0.value_si = T0

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray klist, str kunits, double T0=1, numpy.ndarray weights=None, bint threeParams=True):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in 
        K. A linear least-squares fit is used, which guarantees that the 
        resulting parameters provide the best possible approximation to the 
        data.
        """
        import numpy.linalg
        import scipy.stats
        if threeParams:
            A = numpy.zeros((len(Tlist),3), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = numpy.log(Tlist / T0)
            A[:,2] = -1.0 / constants.R / Tlist
        else:
            A = numpy.zeros((len(Tlist),2), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = -1.0 / constants.R / Tlist
        b = numpy.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n,:] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = numpy.linalg.lstsq(A,b)
        
        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * numpy.linalg.inv(numpy.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)
            
        if not threeParams:
            x = numpy.array([x[0], 0, x[1]])
            cov = numpy.array([[cov[0,0], 0, cov[0,1]], [0,0,0], [cov[1,0], 0, cov[1,1]]])
        
        self.A = (exp(x[0]),kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001,"kJ/mol")
        self.T0 = (T0,"K")
        self.Tmin = (numpy.min(Tlist),"K")
        self.Tmax = (numpy.max(Tlist),"K")
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0,0])),
            sqrt(cov[1,1]),
            sqrt(cov[2,2]) * 0.001,
        )
        
        return self

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(otherKinetics,Arrhenius):
            return False
        if not KineticsModel.isIdenticalTo(self, otherKinetics):
            return False
        if (not self.A.equals(otherKinetics.A) or not self.n.equals(otherKinetics.n)
            or not self.Ea.equals(otherKinetics.Ea) or not self.T0.equals(otherKinetics.T0)):
            return False
                
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes A factor in Arrhenius expression by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor
    

    def toCanteraKinetics(self):
        """
        Converts the Arrhenius object to a cantera Arrhenius object 

        Arrhenius(A,b,E) where A is in units of m^3/kmol/s, b is dimensionless, and E is in J/kmol
        """

        import cantera as ct

        rateUnitsDimensionality = {'1/s':0,
                                   's^-1':0, 
                               'm^3/(mol*s)':1,
                               'm^6/(mol^2*s)':2,
                               'cm^3/(mol*s)':1,
                               'cm^6/(mol^2*s)':2,
                               'm^3/(molecule*s)': 1, 
                               'm^6/(molecule^2*s)': 2,
                               'cm^3/(molecule*s)': 1,
                               'cm^6/(molecule^2*s)': 2,
                               }

        if self._T0.value_si != 1:
            A = self._A.value_si/(self._T0.value_si)**self._n.value_si
        else:
            A = self._A.value_si

        try:
            A *= 1000**rateUnitsDimensionality[self._A.units]
        except KeyError:
            raise Exception('Arrhenius A-factor units {0} not found among accepted units for converting to Cantera Arrhenius object.'.format(self._A.units))
        
        b = self._n.value_si
        E = self._Ea.value_si*1000   # convert from J/mol to J/kmol
        return ct.Arrhenius(A,b,E)

    def setCanteraKinetics(self, ctReaction, speciesList):
        """
        Passes in a cantera ElementaryReaction() object and sets its
        rate to a Cantera Arrhenius() object.
        """
        import cantera as ct
        assert isinstance(ctReaction, ct.ElementaryReaction), "Must be a Cantera ElementaryReaction object"

        # Set the rate parameter to a cantera Arrhenius object
        ctReaction.rate = self.toCanteraKinetics()
    
################################################################################

cdef class ArrheniusEP(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `alpha`         The Evans-Polanyi slope
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusEP object.
        """
        string = 'ArrheniusEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property alpha:
        """The Evans-Polanyi slope."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    property E0:
        """The activation energy for a thermoneutral reaction."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    cpdef double getRateCoefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol. 
        """
        cdef double A, n, Ea        
        Ea = self.getActivationEnergy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T**n * exp(-Ea / (constants.R * T))

    cpdef double getActivationEnergy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double Ea
        Ea = self._alpha.value_si * dHrxn + self._E0.value_si
        if self._E0.value_si > 0:
            if dHrxn < 0.0 and Ea < 0.0:
                Ea = 0.0
            elif dHrxn > 0.0 and Ea < dHrxn:
                Ea = dHrxn
        return Ea
    
    cpdef Arrhenius toArrhenius(self, double dHrxn):
        """
        Return an :class:`Arrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        """
        return Arrhenius(
            A = self.A,
            n = self.n,
            Ea = (self.getActivationEnergy(dHrxn)*0.001,"kJ/mol"),
            T0 = (1,"K"),
            Tmin = self.Tmin,
            Tmax = self.Tmax,
            comment = self.comment,
        )

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(otherKinetics,ArrheniusEP):
            return False
        if not KineticsModel.isIdenticalTo(self, otherKinetics):
            return False
        if (not self.A.equals(otherKinetics.A) or not self.n.equals(otherKinetics.n)
            or not self.alpha.equals(otherKinetics.alpha) or not self.E0.equals(otherKinetics.E0)):
            return False
                
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def setCanteraKinetics(self, ctReaction, speciesList):
        """
        Sets a cantera ElementaryReaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('setCanteraKinetics() is not implemented for ArrheniusEP class kinetics.')

################################################################################

cdef class PDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    a set of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `pressures`     The list of pressures
    `arrhenius`     The list of :class:`Arrhenius` objects at each pressure
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure in bar at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure in bar at which the model is valid, or zero if unknown or undefined
    `efficiencies`  A dict associating chemical species with associated efficiencies
    `order`         The reaction order (1 = first, 2 = second, etc.)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================
    
    """

    def __init__(self, pressures=None, arrhenius=[], highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, highPlimit=highPlimit, comment=comment)
        self.pressures = pressures
        self.arrhenius = arrhenius

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        string = 'PDepArrhenius(pressures={0!r}, arrhenius={1!r}'.format(self.pressures, self.arrhenius)
        if self.highPlimit is not None: string += ', highPlimit={0!r}'.format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.highPlimit, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))
    
    property pressures:
        """The list of pressures."""
        def __get__(self):
            return self._pressures
        def __set__(self, value):
            self._pressures = quantity.Pressure(value)

    cdef getAdjacentExpressions(self, double P):
        """
        Returns the pressures and Arrhenius expressions for the pressures that
        most closely bound the specified pressure `P` in Pa.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] pressures
        cdef int i, ilow, ihigh
        
        pressures = self._pressures.value_si
        
        ilow = 0; ihigh = -1
        for i in range(pressures.shape[0]):
            if pressures[i] <= P:
                ilow = i
            if pressures[i] >= P and ihigh == -1:
                ihigh = i
        return pressures[ilow], pressures[ihigh], self.arrhenius[ilow], self.arrhenius[ihigh]
    
    cpdef double getRateCoefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double Plow, Phigh, klow, khigh, k
        cdef KineticsModel alow, ahigh
        cdef int j
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent PDepArrhenius.getRateCoefficient().')
        
        k = 0.0
        Plow, Phigh, alow, ahigh = self.getAdjacentExpressions(P)
        if Plow == Phigh:
            k = alow.getRateCoefficient(T)
        else:
            klow = alow.getRateCoefficient(T)
            khigh = ahigh.getRateCoefficient(T)
            if klow == khigh == 0.0: return 0.0
            k = klow * 10**(log10(P/Plow)/log10(Phigh/Plow)*log10(khigh/klow))
        return k
    
    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, str kunits, double T0=1):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` with units of `kunits` corresponding to a set of 
        temperatures `Tlist` in K and pressures `Plist` in Pa. An Arrhenius 
        model is fit cpdef changeRate(self, double factor)at each pressure.
        """
        cdef int i
        self.pressures = (Plist*1e-5,"bar")
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = Arrhenius().fitToData(Tlist, K[:,i], kunits, T0)
            self.arrhenius.append(arrhenius)
        return self

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other PDepArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(otherKinetics, PDepArrhenius):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if len(self.arrhenius) != len(otherKinetics.arrhenius):
            return False
        if not self.pressures.equals(otherKinetics.pressures):
            return False
        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].isIdenticalTo(otherKinetics.arrhenius[index]):
                return False
        if self.highPlimit and not self.highPlimit.equals(otherKinetics.highPlimit):
            return False
        
        return True
    
    cpdef changeRate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.changeRate(factor)
        if self.highPlimit is not None:
            self.highPlimit.changeRate(factor)

    def setCanteraKinetics(self, ctReaction, speciesList):
        """
        Sets a Cantera PlogReaction()'s `rates` attribute with
        A list of tuples containing [(pressure in Pa, cantera arrhenius object), (..)]
        """
        import cantera as ct
        import copy
        assert isinstance(ctReaction, ct.PlogReaction), "Must be a Cantera PlogReaction object"

        pressures = copy.deepcopy(self._pressures.value_si)
        ctArrhenius = [arr.toCanteraKinetics() for arr in self.arrhenius]

        ctReaction.rates = zip(pressures, ctArrhenius)
            
################################################################################

cdef class MultiArrhenius(KineticsModel):
    """
    A kinetics model based on a set of (modified) Arrhenius equations, which
    are summed to obtain the overall rate. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`Arrhenius` kinetics
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, arrhenius=[], Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiArrhenius object.
        """
        string = 'MultiArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiArrhenius object.
        """
        return (MultiArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K. 
        """
        cdef double k
        cdef Arrhenius arrh
        k = 0.0
        for arrh in self.arrhenius:
            k += arrh.getRateCoefficient(T)
        return k

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(otherKinetics, MultiArrhenius):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if len(self.arrhenius) != len(otherKinetics.arrhenius):
            return False
        
        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].isIdenticalTo(otherKinetics.arrhenius[index]):
                return False
        
        return True
    
    cpdef Arrhenius toArrhenius(self, double Tmin=-1, double Tmax=-1 ):
        """
        Return an :class:`Arrhenius` instance of the kinetics model 

        Fit the Arrhenius parameters to a set of rate coefficient data generated
        from the MultiArrhenius kinetics, over the temperature range
        Tmin to Tmax, in Kelvin. If Tmin or Tmax are unspecified (or -1)
        then the MultiArrhenius's Tmin and Tmax are used.
        A linear least-squares fit is used, which guarantees that the 
        resulting parameters provide the best possible approximation to the 
        data.
        """
        cdef Arrhenius arrh
        cdef numpy.ndarray Tlist, klist
        cdef str kunits
        if Tmin == -1: Tmin = self.Tmin.value_si
        if Tmax == -1: Tmax = self.Tmax.value_si
        kunits = str(quantity.pq.Quantity(1.0, self.arrhenius[0].A.units).simplified).split()[-1] # is this the best way to get the units returned by k??
        Tlist = numpy.logspace(log10(Tmin), log10(Tmax), num=25)
        klist = numpy.array( map(self.getRateCoefficient, Tlist) , numpy.float64)
        arrh = Arrhenius().fitToData(Tlist, klist, kunits)
        arrh.comment = "Fitted to Multiple Arrhenius kinetics over range {Tmin}-{Tmax} K. {comment}".format(Tmin=Tmin, Tmax=Tmax, comment=self.comment)
        return arrh
    
    cpdef changeRate(self, double factor):
        """
        Change kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.changeRate(factor)

    def setCanteraKinetics(self, ctReaction, speciesList):
        """
        Sets the kinetic rates for a list of cantera `Reaction` objects
        Here, ctReaction must be a list rather than a single cantera reaction.
        """
        if len(ctReaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of Arrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.setCanteraKinetics(ctReaction[i], speciesList)
    
    
################################################################################

cdef class MultiPDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    sets of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`PDepArrhenius` kinetics at each temperature
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, arrhenius=[], Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiPDepArrhenius object.
        """
        string = 'MultiPDepArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiPDepArrhenius object.
        """
        return (MultiPDepArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double k, klow, khigh, Plow, Phigh
        cdef PDepArrhenius arrh
        cdef KineticsModel arrh_low, arrh_high
        cdef numpy.ndarray Plist1, Plist2
        cdef int i
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent MultiPDepArrhenius.getRateCoefficient().')
        
        Plist1 = self.arrhenius[0].pressures.value_si
        for arrh in self.arrhenius[1:]:
            Plist2 = arrh.pressures.value_si
            assert Plist1.shape[0] == Plist2.shape[0]
            for i in range(Plist1.shape[0]):
                assert 0.99 < (Plist2[i] / Plist1[i]) < 1.01            
        
        klow = 0.0; khigh = 0.0
        for arrh in self.arrhenius:
            Plow, Phigh, arrh_low, arrh_high = arrh.getAdjacentExpressions(P)
            klow += arrh_low.getRateCoefficient(T)
            khigh += arrh_high.getRateCoefficient(T)
            
        if klow == khigh == 0.0: 
            return 0.0
        elif Plow == Phigh:
            k = klow
        else:
            k = klow * 10**(log10(P/Plow)/log10(Phigh/Plow)*log10(khigh/klow))
        
        return k

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(otherKinetics, MultiPDepArrhenius):
            return False
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False
        if len(self.arrhenius) != len(otherKinetics.arrhenius):
            return False
        
        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].isIdenticalTo(otherKinetics.arrhenius[index]):
                return False
        
        return True
    
    cpdef changeRate(self, double factor):
        """
        Change kinetic rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.changeRate(factor)

    def setCanteraKinetics(self, ctReaction, speciesList):
        """
        Sets the PLOG kinetics for multiple cantera `Reaction` objects, provided in a list.
        ctReaction is a list of cantera reaction objects.
        """
        if len(ctReaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of PdepArrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.setCanteraKinetics(ctReaction[i], speciesList)