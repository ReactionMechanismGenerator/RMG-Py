# cython: embedsignature=True, cdivision=True

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
This module contains functions for computing the microcanonical rate 
coefficient :math:`k(E)` using RRKM theory (the microcanonical analogue of
TST) and the inverse Laplace transform method. The former is more accurate,
but requires detailed information about the transition state and reactants.
"""

import logging

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport abs, exp, sqrt, log

cimport rmgpy.constants as constants
from rmgpy.exceptions import PressureDependenceError
from rmgpy.kinetics.arrhenius cimport Arrhenius
from rmgpy.statmech.schrodinger import convolve

################################################################################


@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_microcanonical_rate_coefficient(reaction,
                                              np.ndarray[float_t,ndim=1] e_list,
                                              np.ndarray[np.int_t,ndim=1] j_list,
                                              np.ndarray[float_t,ndim=2] reac_dens_states,
                                              np.ndarray[float_t,ndim=2] prod_dens_states=None,
                                              double T=0.0):
    """
    Calculate the microcanonical rate coefficient :math:`k(E)` for the reaction
    `reaction` at the energies `e_list` in J/mol. `reac_dens_states` and
    `prod_dens_states` are the densities of states of the reactant and product
    configurations for this reaction. If the reaction is irreversible, only the
    reactant density of states is required; if the reaction is reversible, then
    both are required. This function will try to use the best method that it
    can based on the input data available:
    
    * If detailed information has been provided for the transition state (i.e.
      the molecular degrees of freedom), then RRKM theory will be used.
    
    * If the above is not possible but high-pressure limit kinetics
      :math:`k_\\infty(T)` have been provided, then the inverse Laplace 
      transform method will be used.

    The density of states for the product `prod_dens_states` and the temperature
    of interest `T` in K can also be provided. For isomerization and association
    reactions `prod_dens_states` is required; for dissociation reactions it is
    optional. The temperature is used if provided in the detailed balance
    expression to determine the reverse kinetics, and in certain cases in the
    inverse Laplace transform method.
    """        
    cdef int n_grains, n_j, r, s
    cdef np.ndarray[float_t,ndim=2] kf, kr
    cdef double c0_inv
    cdef list modes
    cdef bint reactant_states_known, product_states_known, forward

    n_grains = e_list.shape[0]
    n_j = j_list.shape[0]
    kf = np.zeros((n_grains,n_j))
    kr = np.zeros_like(kf)
    active_j_rotor = j_list is None
    active_k_rotor = False

    reactant_states_known = reac_dens_states is not None and reac_dens_states.any()
    product_states_known = prod_dens_states is not None and prod_dens_states.any()
    forward = True

    c0_inv = constants.R * T / 1.0e5

    if reaction.can_tst():
        modes = reaction.transition_state.conformer.get_active_modes(active_j_rotor=active_j_rotor,
                                                                     active_k_rotor=active_k_rotor)

        # We've been provided with molecular degree of freedom data for the
        # transition state, so let's use the more accurate RRKM theory
        logging.debug('Calculating microcanonical rate coefficient using RRKM theory for %s...', reaction)
        if reactant_states_known and (reaction.is_isomerization() or reaction.is_dissociation()):
            kf = apply_rrkm_theory(reaction.transition_state, e_list, j_list, reac_dens_states)
            kf *= c0_inv ** (len(reaction.reactants) - 1)
            forward = True
        elif product_states_known and reaction.is_association():
            kr = apply_rrkm_theory(reaction.transition_state, e_list, j_list, prod_dens_states)
            kr *= c0_inv ** (len(reaction.products) - 1)        
            forward = False
        else:
            raise PressureDependenceError('Unable to compute k(E) values via RRKM theory for path reaction '
                                          '"{0}".'.format(reaction))

    elif reaction.kinetics is not None:
        # We've been provided with high-pressure-limit rate coefficient data,
        # so let's use the less accurate inverse Laplace transform method
        logging.debug('Calculating microcanonical rate coefficient using ILT method for %s...', reaction)
        if reactant_states_known:
            kinetics = reaction.kinetics if reaction.network_kinetics is None else reaction.network_kinetics
            kf = apply_inverse_laplace_transform_method(reaction.transition_state, kinetics, e_list, j_list, reac_dens_states, T)
            forward = True
        elif product_states_known:
            kinetics = reaction.generate_reverse_rate_coefficient(network_kinetics=True)
            kr = apply_inverse_laplace_transform_method(reaction.transition_state, kinetics, e_list, j_list, prod_dens_states, T)
            forward = False
        else:
            raise PressureDependenceError('Unable to compute k(E) values via ILT method for path reaction '
                                          '"{0}".'.format(reaction))
    else:
        raise PressureDependenceError('Unable to compute k(E) values for path reaction "{0}".'.format(reaction))

    # If the reaction is endothermic and barrierless, it is possible that the
    # forward k(E) will have a nonzero value at an energy where the product
    # density of states is zero (but the reactant density of states is not),
    # which violates detailed balance
    # To fix, we set the forward k(E) to zero wherever this is true
    # (This is correct within the accuracy of discretizing the energy grains)
    if reactant_states_known and product_states_known:
        for r in range(n_grains):
            for s in range(n_j):
                if reac_dens_states[r, s] != 0 and prod_dens_states[r, s] != 0:
                    break
                kf[r, s] = 0
                kr[r, s] = 0

    # Get the reverse microcanonical rate coefficient if possible
    if forward and product_states_known:
        # We computed the forward rate coefficient above
        # Thus we need to compute the reverse rate coefficient here
        kr = np.zeros_like(kf)        
        for s in range(n_j):
            for r in range(n_grains):
                if prod_dens_states[r, s] != 0:
                    kr[r, s] = kf[r, s] * reac_dens_states[r, s] / prod_dens_states[r, s]
        kr *= c0_inv ** (len(reaction.products) - len(reaction.reactants))
    elif not forward and reactant_states_known:
        # We computed the reverse rate coefficient above
        # Thus we need to compute the forward rate coefficient here
        kf = np.zeros_like(kr)        
        for s in range(n_j):
            for r in range(n_grains):
                if reac_dens_states[r, s] != 0:
                    kf[r, s] = kr[r, s] * prod_dens_states[r, s] / reac_dens_states[r, s]
        kf *= c0_inv ** (len(reaction.reactants) - len(reaction.products))
    logging.debug('Finished finding microcanonical rate coefficients for path reaction %s', reaction)
    logging.debug('The forward and reverse rates are found to be %g and %g respectively.', kf, kr)

    return kf, kr

@cython.boundscheck(False)
@cython.wraparound(False)
def apply_rrkm_theory(transition_state,
                      np.ndarray[float_t,ndim=1] e_list,
                      np.ndarray[np.int_t,ndim=1] j_list,
                      np.ndarray[float_t,ndim=2] dens_states):
    """
    Calculate the microcanonical rate coefficient for a reaction using RRKM
    theory, where `transition_state` is the transition state of the reaction,
    `e_list` is the array of energies in J/mol at which to evaluate the
    microcanonial rate, and `dens_states` is the density of states of the
    reactant.
    """
    cdef np.ndarray[float_t,ndim=2] k0, k, sum_states
    cdef int n_grains, n_j
    cdef bint active_j_rotor
    cdef double d_e, e0_ts
    cdef int r, s

    from rmgpy.pdep import Configuration

    n_grains = e_list.shape[0]
    n_j = j_list.shape[0]
    active_j_rotor = (n_j == 1)
    k0 = np.zeros((n_grains,n_j))
    k = np.zeros_like(k0)
    d_e = e_list[1] - e_list[0]
    e0_ts = transition_state.conformer.E0.value_si

    conf = Configuration(transition_state)
    conf.calculate_density_of_states(e_list - e_list[0], active_j_rotor=active_j_rotor)

    # Compute tunneling function
    kappa = transition_state.calculate_tunneling_function(e_list)

    # Convolve with transition state density of states to get new transition
    # state sum of states that includes tunneling
    conf.sum_states = convolve(conf.dens_states, kappa)
    conf.e_list += e_list[0] - e0_ts

    E0 = None
    for r in range(n_grains):
        if conf.sum_states[r] > 0:
            E0 = conf.e_list[r]
            break
    if E0 is None:
        raise ValueError, "Could not find a positive sum of states for {0}".format(conf)
    conf.e_list -= E0

    sum_states = conf.map_sum_of_states(e_list - E0, j_list)

    # Generate k(E) using RRKM formula (with tunneling)
    d_e /= constants.h * constants.Na  # J/mol -> s^-1
    for s in range(n_j):
        for r in range(n_grains):
            if sum_states[r, s] > 0 and dens_states[r, s] > 0:
                k[r, s] = sum_states[r, s] / dens_states[r, s] * d_e
    logging.debug('Finished applying RRKM for path transition state %s', transition_state)
    logging.debug('The rate constant is found to be %s', k)
    return k

@cython.boundscheck(False)
@cython.wraparound(False)
def apply_inverse_laplace_transform_method(transition_state,
                                           Arrhenius kinetics,
                                           np.ndarray[float_t,ndim=1] e_list,
                                           np.ndarray[np.int_t,ndim=1] j_list,
                                           np.ndarray[float_t,ndim=2] dens_states,
                                           double T=0.0):
    """
    Calculate the microcanonical rate coefficient for a reaction using the
    inverse Laplace transform method, where `kinetics` is the high pressure 
    limit rate coefficient, `E0` is the ground-state energy of the transition
    state, `e_list` is the array of energies in kJ/mol at which to evaluate the
    microcanonial rate, and `dens_states` is the density of states of the
    reactant. The temperature `T` in K is not required, and is only used when
    the temperature exponent of the Arrhenius expression is negative (for which
    the inverse transform is undefined).
    """
    cdef np.ndarray[float_t,ndim=2] k
    cdef np.ndarray[float_t,ndim=1] phi0, phi
    cdef int n_grains, n_j
    cdef bint active_j_rotor
    cdef double d_e, gas_constant, freq_factor, n, e_a, m0, rem, E0, num, energy, n_crit
    cdef int r, s, m

    n_grains = e_list.shape[0]
    n_j = j_list.shape[0]
    k = np.zeros((n_grains,n_j))
    d_e = e_list[1] - e_list[0]
    gas_constant = constants.R
    E0 = transition_state.conformer.E0.value_si

    n = kinetics.n.value_si
    freq_factor = kinetics.A.value_si / (kinetics.T0.value_si ** n)
    e_a = kinetics.Ea.value_si

    if isinstance(kinetics, Arrhenius) and (T != 0.0 or (e_a >= 0 and n >= 0)):
        
        # The inverse Laplace transform is not defined for Ea < 0 or n < 0
        # In these cases we move the offending portion into the preexponential
        # at the temperature of interest
        # This is an approximation, but it's not worth a more robust procedure

        # Including the T^n piece explicitly also has numerical difficulties
        # for small positive n; it turns out that using this approximation is
        # actually more accurate than trying to handle the T^n piece "properly"
        # For now the implementation is to use this approximation for all n
        # below some critical value, which is purposely placed a bit above zero
        n_crit = 0.25

        if e_a < 0:
            freq_factor *= exp(-e_a / gas_constant / T)
            e_a = 0.0
        if n < n_crit:
            freq_factor *= T ** n
            n = 0.0

        if n < n_crit:
            # Determine the microcanonical rate directly
            m0, rem = divmod(e_a, d_e)
            m = int(m0)
            if rem == 0:
                for s in range(n_j):
                    for r in range(m, n_grains):
                        if e_list[r] > E0 and dens_states[r, s] != 0:
                            k[r, s] = freq_factor * dens_states[r - m, s] / dens_states[r, s]
            else:
                for s in range(n_j):
                    for r in range(m + 1, n_grains):
                        if e_list[r] > E0 and dens_states[r, s] != 0 \
                                and abs(dens_states[r - m, s]) > 1e-12 \
                                and abs(dens_states[r - m - 1, s]) > 1e-12:
                            num = dens_states[r - m, s] * (dens_states[r - m - 1, s] / dens_states[r - m, s]) \
                                  ** (-rem / (e_list[r - m - 1] - e_list[r - m]))
                            k[r, s] = freq_factor * num / dens_states[r, s]

        elif n >= n_crit:
            import scipy.special
            # Evaluate the inverse Laplace transform of the T**n exp(-Ea/RT) piece, which only exists for n >= 0
            phi0 = np.zeros(n_grains, float)
            for r in range(n_grains):
                energy = e_list[r] - e_list[0] - e_a
                if energy > 1:
                    phi0[r] = (energy / gas_constant) ** (n - 1.0)
            phi0 = phi0 * (d_e / gas_constant) / scipy.special.gamma(n)
            # Evaluate the convolution
            for s in range(n_j):
                phi = convolve(phi0, dens_states[:, s])
                # Apply to determine the microcanonical rate
                for r in range(n_grains):
                    if dens_states[r, s] != 0:
                        k[r, s] = freq_factor * phi[r] / dens_states[r, s]

    else:
        raise PressureDependenceError('Unable to use inverse Laplace transform method for non-Arrhenius kinetics or for n < 0.')
    logging.debug('Finished applying inverse lapace transform for path transition state %s', transition_state)
    logging.debug('The rate constant is found to be %s', k)

    return k

def fit_interpolation_model(reaction, Tlist, Plist, K, model, Tmin, Tmax, Pmin, Pmax, error_check=False):
    """
    For a set of phenomenological rate coefficients `K` computed at a grid of
    temperatures `Tlist` in K and pressures `Plist` in Pa, fit a :math:`k(T,P)`
    interpolation `model`, a tuple where the first item is a string describing
    the type of model - either ``'chebyshev'`` or ``'pdeparrhenius'`` - and the
    remaining elements contain parameters for that model. For Chebyshev
    polynomials, the parameters are the number of terms to use in each of the
    temperature and pressure dimensions. For pressure-dependent Arrhenius models
    there are no additional parameters. `Tmin`, `Tmax`, `Pmin`, and `Pmax`
    specify the temperature and pressure ranges in K and Pa, respectively,
    over which the interpolation model is valid. If `error_check` is ``True``,
    a check will be performed to ensure that the interpolation model does not
    deviate too much from the data; as this is not necessarily a fast process,
    it is optional.
    """

    from rmgpy.quantity import Quantity
    from rmgpy.kinetics import PDepArrhenius, Chebyshev

    # Determine units for k(T,P)
    if len(reaction.reactants) == 1:
        kunits = 's^-1'
    elif len(reaction.reactants) == 2:
        kunits = 'm^3/(mol*s)'
    elif len(reaction.reactants) == 3:
        kunits = 'm^6/(mol^2*s)'
    else:
        kunits = ''

    # Set/update the net reaction kinetics using interpolation model
    if model[0].lower() == 'chebyshev':
        model_type, degree_t, degree_p = model
        chebyshev = Chebyshev()
        chebyshev.fit_to_data(Tlist, Plist, K, kunits, degree_t, degree_p, Tmin, Tmax, Pmin, Pmax)
        kinetics = chebyshev
    elif model[0].lower() == 'pdeparrhenius':
        pdep_arrhenius = PDepArrhenius()
        pdep_arrhenius.fit_to_data(Tlist, Plist, K, kunits=kunits, T0=298.0)
        kinetics = pdep_arrhenius
    else:
        return None

    # Set temperature and pressure ranges explicitly (as they may be different
    # from min(Tlist), max(Tlist), min(Plist), max(Plist))
    kinetics.Tmin = Quantity(Tmin, "K")
    kinetics.Tmax = Quantity(Tmax, "K")
    kinetics.Pmin = Quantity(Pmin / 1e5, "bar")
    kinetics.Pmax = Quantity(Pmax / 1e5, "bar")

    # Compute log RMS error for fit
    if error_check:
        log_rms = 0.0
        # Check that fit is within an order of magnitude at all points
        for t, T in enumerate(Tlist):
            for p, P in enumerate(Plist):
                log_k_model = log(kinetics.get_rate_coefficient(T, P))
                log_k_data = log(K[t,p])
                log_rms += (log_k_model - log_k_data) * (log_k_model - log_k_data)
        log_rms = sqrt(log_rms / len(Tlist) / len(Plist))
        if log_rms > 0.5:
            logging.warning('RMS error for k(T,P) fit = %g for reaction %s.', log_rms, reaction)
    logging.debug('Finished fitting model for path reaction %s', reaction)
    logging.debug('The kinetics fit is %r', kinetics)
    return kinetics
