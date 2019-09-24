#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
from __future__ import division

import logging
import math

import numpy as np

import rmgpy.constants as constants


class DiffusionLimited(object):

    def __init__(self):
        # default is false, enabled if there is a solvent
        self.enabled = False

    def enable(self, solvent_data, solvation_database, comment=''):
        # diffusion_limiter is enabled if a solvent has been added to the RMG object.
        logging.info("Enabling diffusion-limited kinetics...")
        diffusion_limiter.enabled = True
        diffusion_limiter.database = solvation_database
        diffusion_limiter.solvent_data = solvent_data

    def disable(self):
        "Turn it off. Mostly useful for unit testing teardown"
        diffusion_limiter.enabled = False
        del diffusion_limiter.database
        del diffusion_limiter.solvent_data

    def get_solvent_viscosity(self, T):
        return self.solvent_data.get_solvent_viscosity(T)

    def get_effective_rate(self, reaction, T):
        """
        Returns the effective rate of reaction, accounting for diffusion.
        For 1<=>2 and 1<=>3 reactions, the reverse rate is limited.
        For 2<=>1 and 3<=>1 reactions, the forward rate is limited.
        For 2<=>2, 2<=>3, 3<=>2, and 3<=>3 reactions, the faster direction is limited.
        """
        intrinsic_kinetics = reaction.kinetics
        reactants = len(reaction.reactants)
        products = len(reaction.products)
        k_forward = intrinsic_kinetics.get_rate_coefficient(T, P=100e5)
        Keq = reaction.get_equilibrium_constant(T)  # Kc
        k_reverse = k_forward / Keq
        k_eff = k_forward

        if reactants == 1:
            if products == 1:
                k_eff = k_forward
            else:  # 2 or 3 products; reverse rate is limited
                k_diff = self.get_diffusion_limit(T, reaction, forward=False)
                k_eff_reverse = k_reverse * k_diff / (k_reverse + k_diff)
                k_eff = k_eff_reverse * Keq
        else:  # 2 or 3 reactants
            if products == 1:
                k_diff = self.get_diffusion_limit(T, reaction, forward=True)
                k_eff = k_forward * k_diff / (k_forward + k_diff)
            else:  # 2 or 3 products
                kf_diff = self.get_diffusion_limit(T, reaction, forward=True)
                krev_diff = self.get_diffusion_limit(T, reaction, forward=False)
                kff = k_forward * kf_diff / (k_forward + kf_diff)
                krevr = k_reverse * krev_diff / (k_reverse + krev_diff)
                kfr = Keq * krevr
                k_eff = min(kff, kfr)
        return k_eff

    def get_diffusion_factor(self, reaction, T):
        """
        Return the diffusion factor of the specified reaction
        This is the ratio of k_eff to k_intrinsic, which is between 0 and 1.
        It is 1.0 if diffusion has no effect.
        """
        return self.get_effective_rate(reaction, T) / reaction.kinetics.get_rate_coefficient(T, P=0)

    def get_diffusion_limit(self, T, reaction, forward=True):
        """
        Return the diffusive limit on the rate coefficient, k_diff.

        This is the upper limit on the rate, in the specified direction.
        (ie. forward direction if forward=True [default] or reverse if forward=False)
        Returns the rate coefficient k_diff in m3/mol/s.
        """
        if forward:
            reacting = reaction.reactants
        else:
            reacting = reaction.products

        if len(reacting) < 2:
            raise Exception("Cannot calculate diffusion limit for a unimolecular reaction")

        radii = 0.0
        diffusivities = []
        for spec in reacting:
            solute_data = self.database.get_solute_data(spec)
            # calculate radius with the McGowan volume and assuming sphere
            radius = ((75 * solute_data.V / constants.pi / constants.Na) ** (1. / 3)) / 100  # m
            diff = solute_data.get_stokes_diffusivity(T, self.get_solvent_viscosity(T))
            radii += radius  # m
            diffusivities.append(diff)  # m^2/s

        # Calculate Smoluchowski kinetics for any reaction order
        # Flegg, SIAM J. Appl. Math., 76 (4), 2016
        N = len(reacting)

        Dinv = 1.0 / np.array(diffusivities)
        Dhat = np.empty(N - 1)

        alpha = (3.0 * N - 5.0) / 2.0
        denom = 0.0

        for i in range(N - 1):
            Dbar = 1.0 / np.sum(Dinv[:(i + 1)])
            Dhat[i] = diffusivities[i + 1] + Dbar

            for j in range(i + 1, N):
                denom += 1.0 / (diffusivities[i] * diffusivities[j])

        delta = np.sum(Dinv) / denom
        k_diff = (np.prod(Dhat) ** 1.5
                  * 4.0 * constants.pi ** (alpha + 1.0) / math.gamma(alpha)
                  * (radii / np.sqrt(delta)) ** (2.0 * alpha)
                  * constants.Na ** (N - 1.0))  # m^(3*(N-1))/mol^(N-1)/s

        return k_diff


# module level variable. There should only ever be one. It starts off disabled
diffusion_limiter = DiffusionLimited()
