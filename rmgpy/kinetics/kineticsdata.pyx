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

cdef class KineticsData(KineticsModel):
    """
    A kinetics model based on an array of rate coefficient data vs. temperature. 
    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `Tdata`         An array of temperatures at which rate coefficient values are known
    `kdata`         An array of rate coefficient values
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, Tdata=None, kdata=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.Tdata = Tdata
        self.kdata = kdata
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsData object.
        """
        string = 'KineticsData(Tdata={0!r}, kdata={1!r}'.format(self.Tdata, self.kdata)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsData object.
        """
        return (KineticsData, (self.Tdata, self.kdata, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    property Tdata:
        """An array of temperatures at which rate coefficient values are known."""
        def __get__(self):
            return self._Tdata
        def __set__(self, value):
            self._Tdata = quantity.Temperature(value)

    property kdata:
        """An array of rate coefficient values."""
        def __get__(self):
            return self._kdata
        def __set__(self, value):
            self._kdata = quantity.RateCoefficient(value)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K. 
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata, kdata
        cdef double Tlow, Thigh, klow, khigh
        cdef double k
        cdef int i, N
        
        Tdata = self._Tdata.value_si
        kdata = self._kdata.value_si
        N = kdata.shape[0]
        k = 0.0
        
        # Make sure we are interpolating and not extrapolating
        if T < Tdata[0]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K using KineticsData model.'.format(T))
        elif T > Tdata[N-1]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K using KineticsData model.'.format(T))
        else:
            for i in range(N-1):
                Tlow = Tdata[i]; Thigh = Tdata[i+1]
                if Tlow <= T and T <= Thigh:
                    klow = kdata[i]; khigh = kdata[i+1]
                    k = klow * (khigh / klow)**((T - Tlow) / (Thigh - Tlow))
                    break
        return k
    
    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if the kdata and Tdata match. Returns ``False`` otherwise.
        """
        if not isinstance(otherKinetics,KineticsData):
            return False        
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False        
        if not self.Tdata.equals(otherKinetics.Tdata) or not self.kdata.equals(otherKinetics.kdata):
            return False
        return True

################################################################################

cdef class PDepKineticsData(PDepKineticsModel):
    """
    A kinetics model based on an array of rate coefficient data vs. temperature
    and pressure. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `Tdata`         An array of temperatures at which rate coefficient values are known
    `Pdata`         An array of pressures at which rate coefficient values are known
    `kdata`         An array of rate coefficient values at each temperature and pressure
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    
    def __init__(self, Tdata=None, Pdata=None, kdata=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.Tdata = Tdata
        self.Pdata = Pdata
        self.kdata = kdata
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepKineticsData object.
        """
        string = 'PDepKineticsData(Tdata={0!r}, Pdata={1!r}, kdata={2!r}'.format(self.Tdata, self.Pdata, self.kdata)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepKineticsData object.
        """
        return (PDepKineticsData, (self.Tdata, self.Pdata, self.kdata, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    property Tdata:
        """An array of temperatures at which rate coefficient values are known."""
        def __get__(self):
            return self._Tdata
        def __set__(self, value):
            self._Tdata = quantity.Temperature(value)

    property Pdata:
        """An array of pressures at which rate coefficient values are known."""
        def __get__(self):
            return self._Pdata
        def __set__(self, value):
            self._Pdata = quantity.Pressure(value)

    property kdata:
        """An array of rate coefficient values at each temperature and pressure."""
        def __get__(self):
            return self._kdata
        def __set__(self, value):
            self._kdata = quantity.RateCoefficient(value)

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3, 
        mol, and s at temperature `T` in K and pressure `P` in Pa. 
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] Tdata, Pdata
        cdef numpy.ndarray[numpy.float64_t,ndim=2] kdata
        cdef double Tlow, Thigh, Plow, Phigh, klow, khigh
        cdef double k
        cdef int i, j, M, N
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent PDepKineticsData.getRateCoefficient().')

        Tdata = self._Tdata.value_si
        Pdata = self._Pdata.value_si
        kdata = self._kdata.value_si
        M = kdata.shape[0]
        N = kdata.shape[1]
        k = 0.0
        
        # Make sure we are interpolating and not extrapolating
        if T < Tdata[0]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K and {1:g} Pa using PDepKineticsData model.'.format(T, P))
        elif T > Tdata[M-1]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K and {1:g} Pa using PDepKineticsData model.'.format(T, P))
        if P < Pdata[0]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K and {1:g} Pa using PDepKineticsData model.'.format(T, P))
        elif P > Pdata[N-1]:
            raise ValueError('Unable to compute rate coefficient at {0:g} K and {1:g} Pa using PDepKineticsData model.'.format(T, P))
        else:
            for i in range(M-1):
                Tlow = Tdata[i]; Thigh = Tdata[i+1]
                if Tlow <= T and T <= Thigh:
                    for j in range(N-1):
                        Plow = Pdata[j]; Phigh = Pdata[j+1]
                        if Plow <= P and P <= Phigh:
                            klow = kdata[i,j] * (kdata[i+1,j] / kdata[i,j])**((T - Tlow) / (Thigh - Tlow))
                            khigh = kdata[i,j+1] * (kdata[i+1,j+1] / kdata[i,j+1])**((T - Tlow) / (Thigh - Tlow))
                            k = klow * (khigh / klow)**(log(P / Plow) / log(Phigh / Plow))
                            break
        
        return k

    cpdef bint isIdenticalTo(self, KineticsModel otherKinetics) except -2:
        """
        Returns ``True`` if the kdata and Tdata match. Returns ``False`` otherwise.
        """
        if not isinstance(otherKinetics,KineticsData):
            return False        
        if not KineticsModel.isIdenticalTo(self,otherKinetics):
            return False        
        if not self.Tdata.equals(otherKinetics.Tdata) or not self.Pdata.equals(otherKinetics.Pdata) or not self.kdata.equals(otherKinetics.kdata):
            return False
        return True
