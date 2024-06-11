# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import rmgpy.quantity as quantity
from rmgpy.rmgobject cimport RMGObject

"""
This module contains base classes that represent various thermodynamic
property calculation methods.
"""


################################################################################

cdef class HeatCapacityModel(RMGObject):
    """
    A base class for heat capacity models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or ``None`` if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or ``None`` if unknown or undefined
    `E0`            The energy at zero Kelvin (including zero point energy)
    `Cp0`           The heat capacity at zero Kelvin
    `CpInf`         The heat capacity at infinity
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, E0=None, Cp0=None, CpInf=None, label='', comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.E0 = E0
        self.Cp0 = Cp0
        self.CpInf = CpInf
        self.label = label
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HeatCapacityModel object.
        """
        return 'HeatCapacityModel(Tmin={0!r}, Tmax={1!r}, E0={2!r}, Cp0={3!r}, Cp0={4!r}, label="""{5}""", comment="""{6}""")'.format(self.Tmin,
                                 self.Tmax, self.E0, self.Cp0, self.CpInf, self.label, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a HeatCapacityModel object.
        """
        return (HeatCapacityModel, (self.Tmin, self.Tmax, self.E0, self.Cp0, self.CpInf, self.label, self.comment))

    property E0:
        """The ground state energy (J/mol) at zero Kelvin, including zero point energy, or ``None`` if not yet specified."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Enthalpy(value)

    property Tmin:
        """The minimum temperature at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Tmin
        def __set__(self, value):
            self._Tmin = quantity.Temperature(value)

    property Tmax:
        """The maximum temperature at which the model is valid, or ``None`` if not defined."""
        def __get__(self):
            return self._Tmax
        def __set__(self, value):
            self._Tmax = quantity.Temperature(value)

    property Cp0:
        """The heat capacity at zero temperature."""
        def __get__(self):
            return self._Cp0
        def __set__(self, value):
            self._Cp0 = quantity.HeatCapacity(value)

    property CpInf:
        """The heat capacity at infinite temperature."""
        def __get__(self):
            return self._CpInf
        def __set__(self, value):
            self._CpInf = quantity.HeatCapacity(value)

    cpdef bint is_temperature_valid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the thermodynamic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self._Tmin is None or self._Tmin.value_si <= T) and (self.Tmax is None or T <= self._Tmax.value_si)

    cpdef double get_heat_capacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the specified
        temperature `T` in K. This method must be overloaded in the derived
        class.
        """
        raise NotImplementedError('Unexpected call to HeatCapacityModel.get_heat_capacity(); you should be using a class derived from HeatCapacityModel.')

    cpdef double get_enthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at the specified temperature `T` in K.
        This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to HeatCapacityModel.get_enthalpy(); you should be using a class derived from HeatCapacityModel.')

    cpdef double get_entropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to HeatCapacityModel.get_entropy(); you should be using a class derived from HeatCapacityModel.')

    cpdef double get_free_energy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at the specified temperature `T`
        in K. This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to HeatCapacityModel.get_free_energy(); you should be using a class derived from HeatCapacityModel.')

    cpdef bint is_similar_to(self, HeatCapacityModel other) except -2:
        """
        Returns ``True`` if `self` and `other` report similar thermo values
        for heat capacity, enthalpy, entropy, and free energy over a wide
        range of temperatures, or ``False`` otherwise.
        """
        cdef double T
        cdef list Tdata
        
        Tdata = [300,400,500,600,800,1000,1500,2000]
        for T in Tdata:
            # Do exact comparison in addition to relative in case both are zero (surface site)
            if self.get_heat_capacity(T) != other.get_heat_capacity(T) and \
                not (0.8 < self.get_heat_capacity(T) / other.get_heat_capacity(T) < 1.25):
                return False
            elif self.get_enthalpy(T) != other.get_enthalpy(T) and \
                not (0.8 < self.get_enthalpy(T) / other.get_enthalpy(T) < 1.25):
                return False
            elif self.get_entropy(T) != other.get_entropy(T) and \
                not (0.8 < self.get_entropy(T) / other.get_entropy(T) < 1.25):
                return False
            elif self.get_free_energy(T) != other.get_free_energy(T) and \
                not (0.8 < self.get_free_energy(T) / other.get_free_energy(T) < 1.25):
                return False

        return True

    cpdef bint is_identical_to(self, HeatCapacityModel other) except -2:
        """
        Returns ``True`` if `self` and `other` report very similar thermo values
        for heat capacity, enthalpy, entropy, and free energy over a wide
        range of temperatures, or ``False`` otherwise.
        """
        cdef double T
        cdef list Tdata
        
        Tdata = [300,400,500,600,800,1000,1500,2000]
        for T in Tdata:
            # Do exact comparison in addition to relative in case both are zero (surface site)
            if self.get_heat_capacity(T) != other.get_heat_capacity(T) and \
                not (0.95 < self.get_heat_capacity(T) / other.get_heat_capacity(T) < 1.05):
                return False
            elif self.get_enthalpy(T) != other.get_enthalpy(T) and \
                not (0.95 < self.get_enthalpy(T) / other.get_enthalpy(T) < 1.05):
                return False
            elif self.get_entropy(T) != other.get_entropy(T) and \
                not (0.95 < self.get_entropy(T) / other.get_entropy(T) < 1.05):
                return False
            elif self.get_free_energy(T) != other.get_free_energy(T) and \
                not (0.95 < self.get_free_energy(T) / other.get_free_energy(T) < 1.05):
                return False

        return True
    
    cpdef double discrepancy(self, HeatCapacityModel other) except -2:
        """
        Return some measure of how dissimilar `self` is from `other`.
        
        The measure is arbitrary, but hopefully useful for sorting purposes.
        Discrepancy of 0 means they are identical
        """
        cdef double T
        cdef double discrepancy
        cdef list Tdata
        
        Tdata = [300,400,500,600,800,1000,1500,2000]
        discrepancy = 0.0
        for T in Tdata:
            #discrepancy += abs(self.get_heat_capacity(T) / other.get_heat_capacity(T) - 1.0)
            #discrepancy += abs(self.get_enthalpy(T) - other.get_enthalpy(T) )/1000
            #discrepancy += abs(self.get_entropy(T) - other.get_entropy(T) )
            discrepancy += abs(self.get_free_energy(T) - other.get_free_energy(T) )/1000
        return discrepancy
