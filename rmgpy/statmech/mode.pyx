# cython: embedsignature=True

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
This module contains the :class:`Mode` class, a base class for all molecular
degrees of freedom.
"""

import numpy as np
cimport numpy as np

################################################################################

cdef class Mode(RMGObject):
    """
    A base class for representing molecular degrees of freedom. The attributes
    are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `quantum`       ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    =============== ============================================================
    
    The choice of classical or quantum mechanical treatments generally depends
    on how the energy levels of the quantum model compare with the temperature
    of interest. In general, translational and rotational modes can be treated
    classically, torsional modes can be treated semiclassically, and vibrational
    modes require a quantum mechanical treatment. Note that not all molecular
    degrees of freedom may support both the classical and quantum mechanical
    models.
    """

    def __init__(self, quantum):
        self.quantum = quantum

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Mode object.
        """
        return 'Mode(quantum={0!r})'.format(self.quantum)

    def __reduce__(self):
        """
        A helper function used when pickling a Mode object.
        """
        return (Mode, (self.quantum,))

    cpdef double get_partition_function(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        raise NotImplementedError('Unexpected call to Mode.get_partition_function(); '
                                  'you should be using a class derived from Mode.')

    cpdef double get_heat_capacity(self, double T) except -100000000:
        """
        Return the heat capacity in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        raise NotImplementedError('Unexpected call to Mode.get_heat_capacity(); '
                                  'you should be using a class derived from Mode.')

    cpdef double get_enthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in J/mol for the degree of freedom at the
        specified temperature `T` in K.
        """
        raise NotImplementedError('Unexpected call to Mode.get_enthalpy(); '
                                  'you should be using a class derived from Mode.')

    cpdef double get_entropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K for the degree of freedom at the
        specified temperature `T` in K.
        """
        raise NotImplementedError('Unexpected call to Mode.get_entropy(); '
                                  'you should be using a class derived from Mode.')

    cpdef np.ndarray get_sum_of_states(self, np.ndarray e_list, np.ndarray sum_states_0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol above the ground state. If an initial sum of states 
        `sum_states_0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        raise NotImplementedError('Unexpected call to Mode.get_sum_of_states(); '
                                  'you should be using a class derived from Mode.')

    cpdef np.ndarray get_density_of_states(self, np.ndarray e_list, np.ndarray dens_states_0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state. If an initial density
        of states `dens_states_0` is given, the rotor density of states will be
        convoluted into these states.
        """
        raise NotImplementedError('Unexpected call to Mode.get_density_of_states(); '
                                  'you should be using a class derived from Mode.')
