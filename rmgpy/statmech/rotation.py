# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains classes that represent various models of rotational
motion. For most molecular systems, a classical treatment of rotational motion
is sufficient since typical rotational energies are much smaller in magnitude
than :math:`k_\\mathrm{B} T`.
"""

import numpy as np
from math import log, sqrt

import rmgpy.constants as constants
import rmgpy.quantity as quantity
import rmgpy.statmech.schrodinger as schrodinger
from rmgpy.rmgobject import recursive_make_object
from rmgpy.statmech.mode import Mode

################################################################################

class Rotation(Mode):
    """
    A base class for all rotational degrees of freedom. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    In the majority of chemical applications, the rotational energy levels are
    very close together compared to :math:`k_\\mathrm{B} T`, which makes the
    rotational motion well-approximated by the classical limit at all relevant
    temperatures. Therefore, the classical model is used by default. 
    """

    def __init__(self, symmetry=1, quantum=False):
        Mode.__init__(self, quantum)
        self.symmetry = symmetry

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Rotation object.
        """
        result = 'Rotation(symmetry={0:d}'.format(self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a Rotation object.
        """
        return (Rotation, (self.symmetry, self.quantum))

    def make_object(self, data, class_dict):
        kwargs = recursive_make_object(data, class_dict, make_final_object=False)
        if ('inertia' in kwargs) and ('rotationalConstant' in kwargs):  # Only one of these can be specified
            del kwargs['inertia']
        self.__init__(**kwargs)

################################################################################

def get_rotational_constant_energy(inertia):
    """
    Return the value of the rotational constant in energy units (J/mol)
    corresponding to the given moment of `inertia` in kg*m^2.
    """
    return constants.hbar * constants.hbar / (2. * inertia) * constants.Na

################################################################################

class LinearRotor(Rotation):
    """
    A statistical mechanical model of a two-dimensional (linear) rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.

    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """

    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        LinearRotor object.
        """
        result = 'LinearRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a LinearRotor object.
        """
        return (LinearRotor, (self.inertia, self.symmetry, self.quantum, None))

    @property
    def inertia(self):
        """The moment of inertia of the rotor."""
        return self._inertia

    @inertia.setter
    def inertia(self, value):
        self._inertia = quantity.Inertia(value)

    @property
    def rotationalConstant(self):
        """The rotational constant of the rotor."""
        I = self._inertia.value_si
        B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
        return quantity.Quantity(B, "cm^-1")

    @rotationalConstant.setter
    def rotationalConstant(self, value):
        B = quantity.Frequency(value)
        I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
        self._inertia = quantity.ScalarQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    def get_level_energy(self, J):
        """
        Return the energy of level `J` in kJ/mol.
        """
        return get_rotational_constant_energy(self.inertia.value_si) * J * (J + 1)

    def get_level_degeneracy(self, J):
        """
        Return the degeneracy of level `J`.
        """
        return 2 * J + 1

    def get_partition_function(self, T):
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Q = schrodinger.get_partition_function(T, self.get_level_energy, self.get_level_degeneracy, 0) / self.symmetry
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            Q = constants.R * T / B / self.symmetry
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Cv = schrodinger.get_heat_capacity(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            Cv = 1.0
        return Cv * constants.R

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            H = schrodinger.get_enthalpy(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            H = 1.0
        return H * constants.R * T

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            S = schrodinger.get_entropy(T, self.get_level_energy, self.get_level_degeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.get_partition_function(T)) + 1.0
        return S * constants.R

    def get_sum_of_states(self, e_list, sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        if self.quantum:
            sum_states = schrodinger.get_sum_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                    sum_states_0) / self.symmetry
        elif sum_states_0 is not None:
            sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list, dens_states_0=None))
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            sum_states = e_list / B / self.symmetry
        return sum_states

    def get_density_of_states(self, e_list, dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        if self.quantum:
            dens_states = schrodinger.get_density_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                         dens_states_0) / self.symmetry
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            dE = e_list[1] - e_list[0]
            dens_states = np.ones_like(e_list) * dE / B / self.symmetry
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states

################################################################################

class NonlinearRotor(Rotation):
    """
    A statistical mechanical model of an N-dimensional nonlinear rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moments of inertia of the rotor
    `rotationalConstant`     The rotational constants of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================

    Note that the moments of inertia and the rotational constants are simply two
    ways of representing the same quantity; only one set of these can be 
    specified independently.
    
    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. In the
    current implementation, the quantum mechanical model has not been 
    implemented, and a :class:`NotImplementedError` will be raised if you try
    to use it. 
    """

    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        NonlinearRotor object.
        """
        result = 'NonlinearRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a NonlinearRotor object.
        """
        return (NonlinearRotor, (self.inertia, self.symmetry, self.quantum, None))

    @property
    def inertia(self):
        """The moments of inertia of the rotor."""
        return self._inertia

    @inertia.setter
    def inertia(self, value):
        self._inertia = quantity.Inertia(value)

    @property
    def rotationalConstant(self):
        """The rotational constant of the rotor."""
        I = self._inertia.value_si
        B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
        return quantity.Quantity(B, "cm^-1")

    @rotationalConstant.setter
    def rotationalConstant(self, value):
        B = quantity.Frequency(value)
        I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
        self._inertia = quantity.ArrayQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    def get_rotational_constant_energy(self):
        """
        Return the values of the rotational constants in kJ/mol.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia.value_si) * constants.Na

    def get_partition_function(self, T):
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        theta = 1.0
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            rotational_constants = self.get_rotational_constant_energy()
            for i in range(rotational_constants.shape[0]):
                theta *= rotational_constants[i] / constants.R
            Q = sqrt(constants.pi * T ** rotational_constants.shape[0] / theta) / self.symmetry
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            Cv = 0.5 * self._inertia.value_si.shape[0]
        return Cv * constants.R

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            H = 0.5 * self._inertia.value_si.shape[0]
        return H * constants.R * T

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            S = np.log(self.get_partition_function(T)) + 0.5 * self._inertia.value_si.shape[0]
        return S * constants.R

    def get_sum_of_states(self, e_list, sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        theta = 1.0
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            if sum_states_0 is not None:
                sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list, dens_states_0=None))
            else:
                rotational_constants = self.get_rotational_constant_energy()
                assert rotational_constants.shape[0] == 3
                for i in range(rotational_constants.shape[0]):
                    theta *= rotational_constants[i]
                sum_states = 4.0 / 3.0 * e_list * np.sqrt(e_list / theta) / self.symmetry
        return sum_states

    def get_density_of_states(self, e_list, dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        theta = 1.0
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            dE = e_list[1] - e_list[0]
            rotational_constants = self.get_rotational_constant_energy()
            assert rotational_constants.shape[0] == 3
            for i in range(rotational_constants.shape[0]):
                theta *= rotational_constants[i]
            dens_states = 2.0 * np.sqrt(e_list / theta) / self.symmetry * dE
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states

################################################################################

class KRotor(Rotation):
    """
    A statistical mechanical model of an active K-rotor (a one-dimensional
    rigid rotor). The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor in amu*angstrom^2
    `rotationalConstant`     The rotational constant of the rotor in cm^-1
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    
    In the majority of chemical applications, the energies involved in the
    K-rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """

    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KRotor object.
        """
        result = 'KRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a KRotor object.
        """
        return (KRotor, (self.inertia, self.symmetry, self.quantum, None))

    @property
    def inertia(self):
        """The moments of inertia of the rotor."""
        return self._inertia

    @inertia.setter
    def inertia(self, value):
        self._inertia = quantity.Inertia(value)

    @property
    def rotationalConstant(self):
        """The rotational constant of the rotor."""
        I = self._inertia.value_si
        B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
        return quantity.Quantity(B, "cm^-1")

    @rotationalConstant.setter
    def rotationalConstant(self, value):
        B = quantity.Frequency(value)
        I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
        self._inertia = quantity.ArrayQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    def get_level_energy(self, J):
        """
        Return the energy of level `J` in kJ/mol.
        """
        return get_rotational_constant_energy(self.inertia.value_si) * J * J

    def get_level_degeneracy(self, J):
        """
        Return the degeneracy of level `J`.
        """
        return 1 if J == 0 else 2

    def get_partition_function(self, T):
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Q = schrodinger.get_partition_function(T, self.get_level_energy, self.get_level_degeneracy, 0) / self.symmetry
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            Q = sqrt(constants.pi * constants.R * T / B) / self.symmetry
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Cv = schrodinger.get_heat_capacity(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            Cv = 0.5
        return Cv * constants.R

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            H = schrodinger.get_enthalpy(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            H = 0.5
        return H * constants.R * T

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            S = schrodinger.get_entropy(T, self.get_level_energy, self.get_level_degeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.get_partition_function(T)) + 0.5
        return S * constants.R

    def get_sum_of_states(self, e_list, sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        if self.quantum:
            sum_states = schrodinger.get_sum_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                    sum_states_0) / self.symmetry
        elif sum_states_0 is not None:
            sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list, dens_states_0=None))
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            sum_states = 2.0 * np.sqrt(e_list / B) / self.symmetry
        return sum_states

    def get_density_of_states(self, e_list, dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        if self.quantum:
            dens_states = schrodinger.get_density_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                         dens_states_0) / self.symmetry
        else:
            dens_states = np.zeros(e_list.shape[0], dtype=np.float)
            B = get_rotational_constant_energy(self._inertia.value_si)
            dE = e_list[1] - e_list[0]
            for r in range(e_list.shape[0]):
                if e_list[r] == 0: continue
                dens_states[r] = 1.0 / sqrt(B * e_list[r]) * dE / self.symmetry
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states

################################################################################

class SphericalTopRotor(Rotation):
    """
    A statistical mechanical model of a three-dimensional rigid rotor with a
    single rotational constant: a spherical top. The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.

    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """

    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None):
        Rotation.__init__(self, symmetry, quantum)
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SphericalTopRotor object.
        """
        result = 'SphericalTopRotor(inertia={0!r}, symmetry={1:d}'.format(self.inertia, self.symmetry)
        if self.quantum:
            result += ', quantum=True'
        result += ')'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling a SphericalTopRotor object.
        """
        return (SphericalTopRotor, (self.inertia, self.symmetry, self.quantum, None))

    @property
    def inertia(self):
        """The moments of inertia of the rotor."""
        return self._inertia

    @inertia.setter
    def inertia(self, value):
        self._inertia = quantity.Inertia(value)

    @property
    def rotationalConstant(self):
        """The rotational constant of the rotor."""
        I = self._inertia.value_si
        B = constants.h / (8 * constants.pi * constants.pi * I) / (constants.c * 100.)
        return quantity.Quantity(B, "cm^-1")

    @rotationalConstant.setter
    def rotationalConstant(self, value):
        B = quantity.Frequency(value)
        I = constants.h / (8 * constants.pi * constants.pi * (B.value_si * constants.c * 100.))
        self._inertia = quantity.ArrayQuantity(I / (constants.amu * 1e-20), "amu*angstrom^2")

    def get_level_energy(self, J):
        """
        Return the energy of level `J` in kJ/mol.
        """
        return get_rotational_constant_energy(self.inertia.value_si) * J * (J + 1)

    def get_level_degeneracy(self, J):
        """
        Return the degeneracy of level `J`.
        """
        return (2 * J + 1) * (2 * J + 1)

    def get_partition_function(self, T):
        """
        Return the value of the partition function :math:`Q(T)` at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Q = schrodinger.get_partition_function(T, self.get_level_energy, self.get_level_degeneracy, 0) / self.symmetry
        else:
            B = get_rotational_constant_energy(self._inertia.value_si)
            theta = constants.R * T / B
            Q = sqrt(theta * theta * theta * constants.pi) / self.symmetry
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            Cv = schrodinger.get_heat_capacity(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            Cv = 1.5
        return Cv * constants.R

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            H = schrodinger.get_enthalpy(T, self.get_level_energy, self.get_level_degeneracy, 0)
        else:
            H = 1.5
        return H * constants.R * T

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        if self.quantum:
            S = schrodinger.get_entropy(T, self.get_level_energy, self.get_level_degeneracy, 0) - log(self.symmetry)
        else:
            S = log(self.get_partition_function(T)) + 1.5
        return S * constants.R

    def get_sum_of_states(self, e_list, sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        if self.quantum:
            sum_states = schrodinger.get_sum_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                    sum_states_0) / self.symmetry
        elif sum_states_0 is not None:
            sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list, dens_states_0=None))
        else:
            if sum_states_0 is not None:
                sum_states = schrodinger.convolve(sum_states_0, self.get_density_of_states(e_list, dens_states_0=None))
            else:
                B = get_rotational_constant_energy(self._inertia.value_si)
                theta = B * B * B
                sum_states = 4.0 / 3.0 * e_list * np.sqrt(e_list / theta) / self.symmetry
        return sum_states

    def get_density_of_states(self, e_list, dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        if self.quantum:
            dens_states = schrodinger.get_density_of_states(e_list, self.get_level_energy, self.get_level_degeneracy, 0,
                                                         dens_states_0) / self.symmetry
        else:
            dE = e_list[1] - e_list[0]
            B = get_rotational_constant_energy(self._inertia.value_si)
            theta = B * B * B
            dens_states = 2.0 * np.sqrt(e_list / theta) / self.symmetry * dE
            if dens_states_0 is not None:
                dens_states = schrodinger.convolve(dens_states_0, dens_states)
        return dens_states
