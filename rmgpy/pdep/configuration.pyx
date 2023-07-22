###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains classes and functions for working with molecular
configurations on a potential energy surface. This includes local minima
(unimolecular and bimolecular) and first-order saddle points (transition
states).
"""

import logging

import cython
import numpy as np
cimport numpy as np
from libc.math cimport log, exp, sqrt

import rmgpy.constants as constants
from rmgpy.statmech import LinearRotor, NonlinearRotor, IdealGasTranslation, HarmonicOscillator
from rmgpy.statmech.conformer import get_density_of_states_forst
from rmgpy.species import Species, TransitionState
from rmgpy.transport import TransportData

################################################################################


cdef class Configuration(object):
    """
    A representation of a molecular configuration on a potential energy surface.
    """

    def __init__(self, *species):
        self.species = list(species)
        self.e_list = None
        self.dens_states = None
        self.sum_states = None
        self.active_j_rotor = False
        self.active_k_rotor = False
        self.energy_correction = 0.0

    def __str__(self):
        return ' + '.join([str(spec) for spec in self.species])

    def __repr__(self):
        string = 'Configuration('
        string += 'species="{0!r}", '.format(self.species)
        if self.e_list is not None: string += 'e_list={0}, '.format(self.e_list)
        if self.dens_states is not None: string += 'dens_states={0}, '.format(self.dens_states)
        if self.sum_states is not None: string += 'sum_states={0}, '.format(self.sum_states)
        string += 'active_k_rotor={0}, '.format(self.active_k_rotor)
        string += 'active_j_rotor={0}, '.format(self.active_j_rotor)
        string += ')'
        return string

    property E0:
        """The ground-state energy of the configuration in J/mol."""
        def __get__(self):
            return sum([float(spec.conformer.E0.value_si) for spec in self.species]) + self.energy_correction 

    cpdef cleanup(self):
        """
        Delete intermediate arrays used in computing k(T,P) values.
        """
        self.e_list = None
        self.dens_states = None
        self.sum_states = None

    cpdef bint is_unimolecular(self) except -2:
        """
        Return ``True`` if the configuration represents a unimolecular isomer,
        or ``False`` otherwise.
        """
        return len(self.species) == 1 and isinstance(self.species[0], Species)

    cpdef bint is_bimolecular(self) except -2:
        """
        Return ``True`` if the configuration represents a bimolecular reactant
        or product channel, or ``False`` otherwise.
        """
        return len(self.species) == 2

    cpdef bint is_termolecular(self) except -2:
        """
        Return ``True`` if the configuration represents a termolecular reactant
        or product channel, or ``False`` otherwise.
        """
        return len(self.species) == 3

    cpdef bint is_transition_state(self) except -2:
        """
        Return ``True`` if the configuration represents a transition state,
        or ``False`` otherwise.
        """
        return len(self.species) == 1 and isinstance(self.species[0], TransitionState)

    cpdef bint has_statmech(self) except -2:
        """
        Return ``True`` if all species in the configuration have statistical
        mechanics parameters, or ``False`` otherwise.
        """
        return all([spec.has_statmech() for spec in self.species])

    cpdef bint has_thermo(self) except -2:
        """
        Return ``True`` if all species in the configuration have thermodynamics
        parameters, or ``False`` otherwise.
        """
        return all([spec.has_thermo() for spec in self.species])
    
    cpdef double get_heat_capacity(self, double T) except -100000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        cdef double cp = 0.0
        for spec in self.species:
            cp += spec.get_heat_capacity(T)
        return cp

    cpdef double get_enthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in kJ/mol at the specified temperature `T` in K.
        """
        cdef double h = 0.0
        for spec in self.species:
            h += spec.get_enthalpy(T)
        h += self.energy_correction
        return h

    cpdef double get_entropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double s = 0.0
        for spec in self.species:
            s += spec.get_entropy(T)
        return s

    cpdef double get_free_energy(self, double T) except 100000000:
        """
        Return the Gibbs free energy in kJ/mol at the specified temperature
        `T` in K.
        """
        cdef double g = 0.0
        for spec in self.species:
            g += spec.get_free_energy(T)
        g += self.energy_correction
        return g

    cpdef double calculate_collision_frequency(self, double T, double P, dict bath_gas) except -1:
        """
        Return the value of the collision frequency in Hz at the given 
        temperature `T` in K and pressure `P` in Pa. If a dictionary `bath_gas`
        of bath gas species and corresponding mole fractions is given, the
        collision parameters of the bas gas species will be averaged with those
        of the species before computing the collision frequency. 
        
        Only the Lennard-Jones collision model is currently supported.
        """
        cdef double bath_gas_sigma, bath_gas_epsilon, bath_gas_mw
        cdef double sigma, epsilon, mu, gas_concentration, frac, tred, omega22

        assert self.is_unimolecular()
        assert isinstance(self.species[0].get_transport_data(), TransportData)
        for spec, frac in bath_gas.items():
            assert isinstance(spec.get_transport_data(), TransportData)

        bath_gas_sigma, bath_gas_epsilon, bath_gas_mw = 0.0, 1.0, 0.0
        for spec, frac in bath_gas.items():
            bath_gas_sigma += spec.get_transport_data().sigma.value_si * frac
            bath_gas_epsilon *= spec.get_transport_data().epsilon.value_si ** frac
            bath_gas_mw += spec.molecular_weight.value_si * frac

        sigma = 0.5 * (self.species[0].get_transport_data().sigma.value_si + bath_gas_sigma)
        epsilon = sqrt((self.species[0].get_transport_data().epsilon.value_si * bath_gas_epsilon))
        mu = 1.0 / (1.0 / self.species[0].molecular_weight.value_si + 1.0 / bath_gas_mw)
        gas_concentration = P / constants.kB / T

        # Evaluate configuration integral
        tred = constants.R * T / epsilon
        omega22 = 1.16145 * tred ** (-0.14874) + 0.52487 * exp(-0.77320 * tred) + 2.16178 * exp(-2.43787 * tred)
        # Evaluate collision frequency
        collision_freq = omega22 * sqrt(8 * constants.kB * T / constants.pi / mu) * constants.pi * sigma*sigma * gas_concentration
        logging.debug("Obtained the following collision parameters: sigma {0}, epsilon {1}, "
                      "omega22 {2}, collision rate {3}".format(sigma, epsilon, omega22, collision_freq))

        return collision_freq

    cpdef np.ndarray generate_collision_matrix(self, double T, np.ndarray dens_states, np.ndarray e_list,
                                               np.ndarray j_list=None):
        """
        Return the collisional energy transfer probabilities matrix for the
        configuration at the given temperature `T` in K using the given
        energies `e_list` in kJ/mol and total angular momentum quantum numbers
        `j_list`. The density of states of the configuration `dens_states` in
        mol/kJ is also required.
        """
        assert self.is_unimolecular()
        assert self.species[0].energy_transfer_model is not None
        return self.species[0].energy_transfer_model.generate_collision_matrix(T, dens_states, e_list, j_list)

    cpdef calculate_density_of_states(self, np.ndarray e_list, bint active_j_rotor=True, bint active_k_rotor=True,
                                      bint rmgmode=False):
        """
        Calculate the density (and sum) of states for the configuration at the
        given energies above the ground state `e_list` in J/mol. The
        `active_j_rotor` and `active_k_rotor` flags control whether the J-rotor
        and/or K-rotor are treated as active (and therefore included in the
        density and sum of states). The computed density and sum of states
        arrays are stored on the object for future use.
        """
        cdef list modes
        cdef int i

        logging.debug('calculating density of states for {}'.format(self.__str__()))

        self.e_list = e_list
        self.active_j_rotor = active_j_rotor
        self.active_k_rotor = active_k_rotor

        # Get the active rovibrational modes for each species in the configuration
        modes = []
        for i, species in enumerate(self.species):
            modes.extend(species.conformer.get_active_modes(active_k_rotor=self.active_k_rotor,
                                                            active_j_rotor=self.active_j_rotor))

        if rmgmode or len(modes) == 0:
            # Include an arbitrary active rigid rotor if needed
            # The moments of inertia cancel in all subsequent calculations
            for mode in modes:
                if isinstance(mode, (LinearRotor,NonlinearRotor)):
                    break
            else:
                linear = False
                for species in self.species:
                    for molecule in species.molecule:
                        if molecule.is_linear():
                            linear = True
                            break
                if linear:
                    modes.insert(0, LinearRotor(inertia=(1.0, "amu*angstrom^2"), symmetry=1))
                else:
                    modes.insert(0, NonlinearRotor(inertia=([1.0, 1.0, 1.0], "amu*angstrom^2"), symmetry=1))

        if len(modes) == 0:
            self.dens_states = None
            self.sum_states = None
        else:        
            # If the configuration is bimolecular, also include the relative
            # translational motion of the two molecules
            if self.is_bimolecular() or self.is_termolecular():
                mass = []
                for species in self.species:
                    for mode in species.conformer.modes:
                        if isinstance(mode, IdealGasTranslation):
                            mass.append(mode.mass.value_si)
                            break
                    else:
                        if species.molecular_weight is not None:
                            mass.append(species.molecular_weight.value_si)
                        else:
                            m = 0
                            for atom in species.molecule[0].atoms:
                                m += atom.element.mass
                            mass.append(m * constants.amu * 1000)
                if self.is_bimolecular():
                    if len(mass) != 2:
                        raise AttributeError('Length of masses should be two for bimolecular reactants. '
                                             'We got {0}.'.format(len(mass)))
                    mu = 1.0 / (1.0 / mass[0] + 1.0 / mass[1])
                    modes.insert(0, IdealGasTranslation(mass=(mu/constants.amu, "amu")))
                else:
                    if len(mass) != 3:
                        raise AttributeError('Length of masses should be three for termolecular reactants. '
                                             'We got {0}.'.format(len(mass)))
                    mu = 1.0 / (1.0 / mass[0] + 1.0 / mass[1])
                    mu2 = 1.0 / (1.0 / mass[1] + 1.0 / mass[2])
                    mu3 = 1.0 / (1.0 / mass[0] + 1.0 / mass[2])
                    mu_avg = (mu + mu2 + mu3) / 3 * 1.1447142
                    modes.insert(0, IdealGasTranslation(mass=(mu_avg / constants.amu, "amu")))
                    modes.insert(0, IdealGasTranslation(mass=(mu_avg / constants.amu, "amu")))
            if rmgmode:
                # Compute the density of states by direct count
                # This is currently faster than the method of steepest descents,
                # but requires classical hindered rotors
                dens_states = None
                for mode in modes:
                    if not isinstance(mode,HarmonicOscillator):
                        dens_states = mode.get_density_of_states(self.e_list, dens_states)
                    # Fix a numerical artifact that occurs when two modes have
                    # density of states expressions that are zero at the
                    # ground state
                    # Convoluting these modes gives a zero at the first excited
                    # energy grain as well, which is unphysical
                    # Instead, fill in an approximate value by extrapolation
                    # This should only occur in systems with IdealGasTranslation
                    # and NonlinearRotor modes
                    if dens_states[1] == 0:
                        dens_states[1] = dens_states[2] * dens_states[2] / dens_states[3]
                for mode in modes:
                    if isinstance(mode,HarmonicOscillator):
                        dens_states = mode.get_density_of_states(self.e_list, dens_states)
                self.dens_states = dens_states
                for spec in self.species:
                    self.dens_states *= spec.conformer.spin_multiplicity * spec.conformer.optical_isomers

            else:
                # Since the evaluation of quantum hindered rotors is slow, it is
                # actually faster (and probably negligibly less accurate) to use
                # interpolation in the steepest descents algorithm
                import scipy.interpolate

                log_t_data = np.linspace(log(10.), log(10000.), 250)
                t_data = np.exp(log_t_data)
                q_data = np.ones_like(t_data)
                for i in range(t_data.shape[0]):
                    t = t_data[i]
                    for mode in modes:
                        q_data[i] = q_data[i] * mode.get_partition_function(t)
                log_q = scipy.interpolate.InterpolatedUnivariateSpline(t_data, np.log(q_data))
                # log_q = LinearInterpolator(t_data, np.log(q_data))

                self.dens_states, self.sum_states = get_density_of_states_forst(self.e_list, log_q)

                for spec in self.species:
                    self.dens_states *= spec.conformer.spin_multiplicity * spec.conformer.optical_isomers
                    self.sum_states *= spec.conformer.spin_multiplicity * spec.conformer.optical_isomers
        if self.dens_states is None:
            raise ValueError("Species {} has no active modes".format(species.label))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def map_density_of_states(self, np.ndarray[float_t, ndim=1] e_list, np.ndarray[np.int_t, ndim=1] j_list=None):
        """
        Return a mapping of the density of states for the configuration to the
        given energies `e_list` in J/mol and, if the J-rotor is not active, the
        total angular momentum quantum numbers `j_list`.
        """
        cdef np.ndarray[float_t,ndim=2] dens_states
        cdef double E0, de0, b1, b2, e, d_j
        cdef int r0, r, s, t, n_grains, n_j, j, j1, j2
        cdef list b_list

        import scipy.interpolate
        for r in range(self.e_list.shape[0]):
            if self.dens_states[r] > 0:
                break
        f = scipy.interpolate.InterpolatedUnivariateSpline(self.e_list[r:], np.log(self.dens_states[r:]))

        E0 = self.E0
        n_grains = e_list.shape[0]
        d_e = e_list[1] - e_list[0]
        d_e0 = self.e_list[1] - self.e_list[0]

        if self.active_j_rotor:
            dens_states = np.zeros((n_grains,1))
            for r0 in range(n_grains):
                if e_list[r0] >= E0: break
            for r in range(r0, n_grains):
                dens_states[r, 0] = f(e_list[r] - E0)
            dens_states[r0:, 0] = np.exp(dens_states[r0:, 0])
        else:
            assert j_list is not None
            n_j = j_list.shape[0]
            d_j = j_list[1] - j_list[0]
            dens_states = np.zeros((n_grains, n_j))

            b_list = []
            for spec in self.species:
                j_rotor, k_rotor = spec.conformer.get_symmetric_top_rotors()
                b_list.append(float(j_rotor.rotationalConstant.value_si))

            for r0 in range(n_grains):
                if e_list[r0] >= E0: break

            if len(b_list) == 1:
                b1 = b_list[0] * 11.962  # cm^-1 to J/mol
                for r in range(r0, n_grains):
                    for s in range(n_j):
                        j1 = j_list[s]
                        e = e_list[r] - E0 - b1 * j1 * (j1 + 1)
                        if e < 0: break
                        dens_states[r,s] = (2 * j1 + 1) * exp(f(e)) * d_j

            elif len(b_list) == 2:
                b1 = b_list[0] * 11.962  # cm^-1 to J/mol
                b2 = b_list[1] * 11.962
                for r in range(r0, n_grains):
                    for s in range(n_j):
                        j = j_list[s]
                        for t in range(s + 1):
                            j1 = j_list[t]
                            j2 = j - j1
                            e = e_list[r] - E0 - b1 * j1 * (j1 + 1) - b2 * j2 * (j2 + 1)
                            if e > 0:
                                dens_states[r, s] += (2 * j1 + 1) * (2 * j2 + 2) * exp(f(e)) * d_j * d_j

        return dens_states * d_e / d_e0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def map_sum_of_states(self, np.ndarray[float_t, ndim=1] e_list, np.ndarray[np.int_t, ndim=1] j_list=None):
        """
        Return a mapping of the density of states for the configuration to the
        given energies `e_list` in J/mol and, if the J-rotor is not active, the
        total angular momentum quantum numbers `j_list`.
        """
        cdef np.ndarray[float_t,ndim=2] sum_states
        cdef double E0, b1, b2, d_j
        cdef int r0, r, s, n_grains, n_j, j1, j2

        import scipy.interpolate
        for r in range(self.e_list.shape[0]):
            if self.sum_states[r] > 0:
                break
        f = scipy.interpolate.InterpolatedUnivariateSpline(self.e_list[r:], np.log(self.sum_states[r:]))

        E0 = self.E0
        n_grains = len(e_list)

        if self.active_j_rotor:
            sum_states = np.zeros((n_grains,1))
            for r0 in range(n_grains):
                if e_list[r0] >= E0:
                    break
            for r in range(r0, n_grains):
                sum_states[r, 0] = f(e_list[r] - E0)
            sum_states[r0:, 0] = np.exp(sum_states[r0:, 0])
        else:
            assert j_list is not None
            n_j = len(j_list)
            d_j = j_list[1] - j_list[0]
            sum_states = np.zeros((n_grains, n_j))

            b_list = []
            for spec in self.species:
                j_rotor, k_rotor = spec.conformer.get_symmetric_top_rotors()
                b_list.append(float(j_rotor.rotationalConstant.value_si))

            for r0 in range(n_grains):
                if e_list[r0] >= E0:
                    break

            if len(b_list) == 1:
                b1 = b_list[0] * 11.962  # cm^-1 to J/mol
                for r in range(r0, n_grains):
                    for s in range(n_j):
                        j1 = j_list[s]
                        e = e_list[r] - E0 - b1 * j1 * (j1 + 1)
                        if e < 0:
                            break
                        sum_states[r,s] = (2 * j1 + 1) * exp(f(e)) * d_j

            elif len(b_list) == 2:
                b1 = b_list[0] * 11.962  # cm^-1 to J/mol
                b2 = b_list[1] * 11.962
                for r in range(r0, n_grains):
                    for s in range(n_j):
                        j = j_list[s]
                        for t in range(s + 1):
                            j1 = j_list[t]
                            j2 = j - j1
                            e = e_list[r] - E0 - b1 * j1 * (j1 + 1) - b2 * j2 * (j2 + 1)
                            if e > 0:
                                sum_states[r, s] += (2 * j1 + 1) * (2 * j2 + 1) * exp(f(e)) * d_j * d_j

        return sum_states
