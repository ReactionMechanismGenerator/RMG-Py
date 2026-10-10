#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains classes for sensitivity analysis
of kinetics and pressure-dependent jobs.
"""

import logging
import os
import string
import yaml
from copy import deepcopy
import numpy as np

from rmgpy.pdep import Configuration
import rmgpy.quantity as quantity
from rmgpy.species import TransitionState
from rmgpy.kinetics import Arrhenius
from rmgpy.exceptions import InvalidMicrocanonicalRateError, ModifiedStrongCollisionError, PressureDependenceError

################################################################################


class KineticsSensitivity(object):
    """
    The :class:`KineticsSensitivity` class represents an instance of a sensitivity analysis job
    performed for a KineticsJob. The attributes are:

    =================== ================================================================================================
    Attribute           Description
    =================== ================================================================================================
    `conditions`        A list of the conditions at which the sensitivity coefficients are calculated
    `job`               The KineticsJob object
    `f_rates`           A list of forward rates from `job` at the respective `conditions` in the appropriate units
    `r_rates`           A list of reverse rates from `job` at the respective `conditions` in the appropriate units
    `f_sa_rates`        A dictionary with Species as keys and each value is a list of forward rates from `job` at the
                        respective `conditions` in the appropriate units after perturbing the corresponding Species' E0
    `r_sa_rates`        Same as `f_sa_rates`, only for the reverse direction
    `f_sa_coefficients` A dictionary with Species keys and sensitivity coefficients in the forward direction as values
    `r_sa_coefficients` A dictionary with Species keys and sensitivity coefficients in the reverse direction as values
    =================== ================================================================================================
    """

    def __init__(self, job, output_directory):
        self.job = job
        self.output_directory = output_directory
        self.sensitivity_path = os.path.join(output_directory, 'sensitivity')
        self.conditions = self.job.sensitivity_conditions
        self.f_rates = [self.job.reaction.kinetics.get_rate_coefficient(condition.value_si)
                        for condition in self.conditions]
        kr = self.job.reaction.generate_reverse_rate_coefficient()
        self.r_rates = [kr.get_rate_coefficient(condition.value_si) for condition in self.conditions]
        self.f_sa_rates = {}
        self.r_sa_rates = {}
        self.f_sa_coefficients = {}
        self.r_sa_coefficients = {}
        self.perturbation = quantity.Quantity(1, 'kcal/mol')
        self.execute()

    def execute(self):
        """
        Execute the sensitivity analysis for a :class:KineticsJob: object
        """
        for species in [self.job.reaction.reactants[0], self.job.reaction.products[0],
                        self.job.reaction.transition_state]:
            self.perturb(species)
            self.job.execute(output_file=None, plot=False)  # run the perturbed job
            self.f_sa_rates[species] = [self.job.reaction.kinetics.get_rate_coefficient(condition.value_si)
                                        for condition in self.conditions]
            kr = self.job.reaction.generate_reverse_rate_coefficient()
            self.r_sa_rates[species] = [kr.get_rate_coefficient(condition.value_si)
                                        for condition in self.conditions]
            self.unperturb(species)
            # Calculate the sensitivity coefficients according to dln(r) / dln(E0) = (E0 * dr) / (r * dE0)
            self.f_sa_coefficients[species] = [(self.f_sa_rates[species][i] - self.f_rates[i]) /
                                               (self.perturbation.value_si * self.f_rates[i])
                                               for i in range(len(self.conditions))]
            self.r_sa_coefficients[species] = [(self.r_sa_rates[species][i] - self.r_rates[i]) /
                                               (self.perturbation.value_si * self.r_rates[i])
                                               for i in range(len(self.conditions))]
        self.save()
        self.plot()

    def perturb(self, species):
        """Perturb a species' E0"""
        species.conformer.E0.value_si += self.perturbation.value_si

    def unperturb(self, species):
        """Return the species' E0 to its original value"""
        species.conformer.E0.value_si -= self.perturbation.value_si  # restore E0 to its original value

    def save(self):
        """Save the SA results as tabulated data as well as in YAML format"""
        if not os.path.exists(self.sensitivity_path):
            os.mkdir(self.sensitivity_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.job.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.job.reaction.products]))
        filename = ''.join(c for c in reaction_str if c in valid_chars)
        path = os.path.join(self.sensitivity_path, filename + '.txt')
        sa_data = dict()
        sa_data['structures'] = dict()
        with open(path, 'w') as sa_f:
            sa_f.write("Sensitivity analysis for reaction {0}\n\n"
                       "The semi-normalized sensitivity coefficients are calculated as dln(r)/dE0\n"
                       "by perturbing E0 of each well or TS by {1}, and are given in "
                       "`mol/J` units.\n\n\n".format(reaction_str, self.perturbation))
            for species in self.job.reaction.reactants + self.job.reaction.products:
                if species.label not in sa_data['structures']:
                    sa_data['structures'][species.label] = species.to_adjacency_list()
            reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
            ts_label = self.job.reaction.transition_state.label
            products_label = ' + '.join([reactant.label for reactant in self.job.reaction.products])
            sa_data[reactants_label], sa_data[ts_label], sa_data[products_label] = dict(), dict(), dict()
            max_label = max(len(reactants_label), len(products_label), len(ts_label), 10)
            sa_f.write('========================={0}=============================================\n'
                       '| Direction | Well or TS {1}| Temperature (K) | Sensitivity coefficient |\n'
                       '|-----------+------------{2}+-----------------+-------------------------|\n'
                       .format('=' * (max_label - 10), ' ' * (max_label - 10), '-' * (max_label - 10)))
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    reactants_label, ' ' * (max_label - len(reactants_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]))
                sa_data[reactants_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    products_label, ' ' * (max_label - len(products_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]))
                sa_data[products_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    ts_label, ' ' * (max_label - len(ts_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]))
                sa_data[ts_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]
            sa_f.write('|-----------+------------{0}+-----------------+-------------------------|\n'.format(
                '-' * (max_label - 10)))
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    reactants_label, ' ' * (max_label - len(reactants_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.reactants[0]][i]))
                sa_data[reactants_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    products_label, ' ' * (max_label - len(products_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.products[0]][i]))
                sa_data[products_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    ts_label, ' ' * (max_label - len(ts_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.transition_state][i]))
                sa_data[ts_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]
            sa_f.write('========================={0}=============================================\n'.format(
                '=' * (max_label - 10)))

        with open(os.path.join(self.sensitivity_path, filename + '.yml'), 'w') as f:
            yaml.dump(data=sa_data, stream=f)

    def plot(self):
        """Plot the SA results as horizontal bars"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
        ts_label = self.job.reaction.transition_state.label
        products_label = ' + '.join([reactant.label for reactant in self.job.reaction.products])

        plt.rcdefaults()
        ax = plt.subplots(nrows=len(self.conditions), ncols=2, tight_layout=True)[1]
        labels = [reactants_label, ts_label, products_label]
        min_sa = min(min(min(self.f_sa_coefficients.values())), min(min(self.r_sa_coefficients.values())))
        max_sa = max(max(max(self.f_sa_coefficients.values())), max(max(self.r_sa_coefficients.values())))
        for i, condition in enumerate(self.conditions):
            f_values = [self.f_sa_coefficients[self.job.reaction.reactants[0]][i],
                        self.f_sa_coefficients[self.job.reaction.transition_state][i],
                        self.f_sa_coefficients[self.job.reaction.products[0]][i]]
            r_values = [self.r_sa_coefficients[self.job.reaction.reactants[0]][i],
                        self.r_sa_coefficients[self.job.reaction.transition_state][i],
                        self.r_sa_coefficients[self.job.reaction.products[0]][i]]
            y_pos = np.arange(3)
            ax[i][0].barh(y_pos, f_values, align='center', color='green')
            ax[i][0].set_yticks(y_pos)
            ax[i][0].set_yticklabels(labels)
            ax[i][0].invert_yaxis()  # labels read top-to-bottom
            ax[i][0].set_xlabel(r'Sensitivity: $\frac{\partial\:\ln{k}}{\partial\:E0}$, ($\frac{mol}{J}$)')
            ax[i][0].set_title('Forward, {0}'.format(condition))
            ax[i][0].set_xlim([min_sa, max_sa])
            ax[i][1].barh(y_pos, r_values, align='center', color='blue')
            ax[i][1].set_yticks(y_pos)
            ax[i][1].set_yticklabels(labels)
            ax[i][1].invert_yaxis()  # labels read top-to-bottom
            ax[i][1].set_xlabel(r'Sensitivity: $\frac{\partial\:\ln{k}}{\partial\:E0}$, ($\frac{mol}{J}$)')
            ax[i][1].set_title('Reverse, {0}'.format(condition))
            ax[i][1].set_xlim([min_sa, max_sa])
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

        if not os.path.exists(self.sensitivity_path):
            os.mkdir(self.sensitivity_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
        products_label = ' + '.join([product.label for product in self.job.reaction.products])
        reaction_str = f'{reactants_label} <=> {products_label}'
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        path = os.path.join(self.sensitivity_path, filename)
        plt.savefig(path)
        plt.close()


class PDepSensitivity(object):
    """
    The :class:`Sensitivity` class represents an instance of a sensitivity analysis job
    performed for a PressureDependenceJob. The attributes are:

    =================== ================================================================================================
    Attribute           Description
    =================== ================================================================================================
    `conditions`        A list of the conditions (each entry is a list of one T and one P quantities) at which the
                        sensitivity coefficients are calculated
    `job`               The PressureDependenceJob object
    `rates`             A dictionary with net_reactions as keys. Values are lists of forward rates from `job` for the
                        respective path reaction at the respective `conditions` in the appropriate units
    `sa_rates`          A dictionary with string representations of net_reactions as keys. Values are dictionaries with
                        Wells or TransitionStates as keys and each value is a list of forward rates from `job` at the
                        respective `conditions` after perturbing the corresponding well's E0, or, for a TS, its E0
                        and (if the owning path reaction is ILT-based rather than RRKM-based) the Ea of the
                        Arrhenius kinetics used to derive its microcanonical rate, both by the same amount. For
                        such TS rows the reported coefficient is therefore NOT a plain dln(r)/dE0(TS): it is a
                        derivative along the coordinate that raises the barrier (TS E0 and Ea together)
    `sa_coefficients`   A dictionary with similar structure as `sa_rates`, containing the sensitivity coefficients
                        in the forward direction
    =================== ================================================================================================
    """

    def __init__(self, job, output_directory, perturbation, max_iters=5):
        self.job = job
        self.output_directory = output_directory
        self.sensitivity_path = os.path.join(output_directory, 'sensitivity')
        self.conditions = self.job.sensitivity_conditions
        self.perturbation = perturbation
        
        base_wells = []
        base_wells.extend(self.job.network.reactants)
        base_wells.extend(self.job.network.isomers)
        base_wells.extend(self.job.network.products)
        Emin = min([x.E0 for x in base_wells])
        
        base_transition_states = []
        Emax = Emin
        for rxn in self.job.network.path_reactions:
            base_transition_states.append(rxn.transition_state)
            if rxn.transition_state.conformer.E0.value_si > Emax:
                Emax = rxn.transition_state.conformer.E0.value_si
        if self.perturbation.value_si > 0:
            self.job.network.Emin = Emin
            self.job.network.Emax = Emax + self.perturbation.value_si
        else:
            self.job.network.Emin = Emin + self.perturbation.value_si
            self.job.network.Emax = Emax
        
        self.job.execute(output_file=None, plot=False, print_summary=False)
        
        self.rates = {}
        for rxn in self.job.network.net_reactions:
            self.rates[str(rxn)] = []
            for condition in self.conditions:
                self.rates[str(rxn)].append(rxn.kinetics.get_rate_coefficient(condition[0].value_si,
                                                                              condition[1].value_si))
        self.sa_rates = {}
        self.sa_coefficients = {}
        for rxn in self.job.network.net_reactions:
            self.sa_rates[str(rxn)] = {}
            self.sa_coefficients[str(rxn)] = {}
        
        self.max_iters = max_iters
        self.execute()

    def execute(self):
        """
        Execute the sensitivity analysis for a :class:PressureDependenceJob: object
        """
        base_wells = []
        base_wells.extend(self.job.network.reactants)
        base_wells.extend(self.job.network.isomers)
        base_wells.extend(self.job.network.products)
        base_transition_states = []
        for rxn in self.job.network.path_reactions:
            # if rxn.transition_state is not None:
            base_transition_states.append(rxn.transition_state)
        base_perturbation = deepcopy(self.perturbation)
        base_job = self.job
        for j in range(len(base_wells + base_transition_states)):
            self.perturbation = base_perturbation
            base_entry = (base_wells + base_transition_states)[j]
            failed = False
            c = 0
            wells = []
            transition_states = []
            while c < self.max_iters:
                # Each attempt works on a fresh deep copy of the base job, so the in-place E0/Ea
                # shifts that perturb() applies never touch base_job. This is what makes the
                # perturb/unperturb pair exception-safe: if self.job.execute() below raises, the
                # corrupted copy is simply discarded on the next iteration -- the invariant does
                # not rely on unperturb() running on every code path.
                self.job = deepcopy(base_job)
                wells = []
                wells.extend(self.job.network.reactants)
                wells.extend(self.job.network.isomers)
                wells.extend(self.job.network.products)
                transition_states = []
                for rxn in self.job.network.path_reactions:
                    # if rxn.transition_state is not None:
                    transition_states.append(rxn.transition_state)
                entry = (wells+transition_states)[j]
                if entry in wells:
                    logging.info("\n\nPerturbing well '{0}' by {1}:".format(entry, self.perturbation))
                else:
                    logging.info("\n\nPerturbing TS '{0}' by {1}:".format(entry.label, self.perturbation))
                self.perturb(entry)
                try:
                    self.job.execute(output_file=None, plot=False, print_summary=False)  # run the perturbed job
                    self.unperturb(entry)
                    break
                except (InvalidMicrocanonicalRateError, ModifiedStrongCollisionError) as e:
                    self.unperturb(entry)
                    c += 1
                    self.perturbation = quantity.Quantity(self.perturbation.value/2.0, self.perturbation.units)
                    logging.error("Decreasing perturbation to {}".format(self.perturbation))
                    
            if c == self.max_iters:
                if entry in wells:
                    logging.error("Perturbation of well '{0}' has failed".format(entry))
                else:
                    logging.error("Perturbation of TS '{0}' has failed".format(entry.label))
                failed = True
            for rxn in self.job.network.net_reactions:
                if failed:
                    self.sa_rates[str(rxn)][base_entry] = [np.NaN for condition in self.conditions]
                    self.sa_coefficients[str(rxn)][base_entry] = [np.NaN for i in range(len(self.conditions))]
                else:
                    self.sa_rates[str(rxn)][base_entry] = [rxn.kinetics.get_rate_coefficient(
                        condition[0].value_si, condition[1].value_si) for condition in self.conditions]
                    self.sa_coefficients[str(rxn)][base_entry] = [((self.sa_rates[str(rxn)][base_entry][i]
                                                               - self.rates[str(rxn)][i])) /
                                                             (self.perturbation.value_si * self.rates[str(rxn)][i])
                                                             for i in range(len(self.conditions))]
        self.perturbation = base_perturbation
        self.job = base_job
        self.save(base_wells, base_transition_states)
        self.plot(base_wells, base_transition_states)

    def perturb(self, entry, unperturb=False):
        """
        Perturb E0 of `entry` which could be either a :class:TransitionState or a :class:Configuration
        In the latter case, only the first species in the Configuration.species list is perturbed.
        The perturbation is done by addition of the energy amount in self.perturbation.
        If unperturb is `False`, the perturbation is addition of the energy amount in self.perturbation.
        If unperturb is `False`, this is done by subtracting.

        For a :class:TransitionState belonging to one or more path reactions whose microcanonical
        rate is computed via the inverse Laplace transform (ILT) method rather than RRKM theory
        (i.e. ``not rxn.can_tst()``), perturbing the TS's E0 alone is a structural no-op for
        n >= 0.25 and very nearly one otherwise. In ``apply_inverse_laplace_transform_method``
        (rmgpy/pdep/reaction.pyx) the TS E0 is never referenced when the Arrhenius temperature
        exponent n >= 0.25, so E0 is genuinely dead there; for n < 0.25 it enters only as a
        grain-threshold cutoff, and in every case the resulting k(E) is renormalized
        (rmgpy/pdep/network.py, ``calculate_microcanonical_rates``) by a single scalar so that its
        canonical average reproduces ``kf_expected``, the high-pressure-limit Arrhenius k(T), which
        depends on `Ea` and not on the TS E0. (For n < 0.25 the threshold shift does survive that
        scalar renormalization into the pressure-dependent rate, so E0 is not strictly a no-op
        there -- but that is precisely why perturbing E0 alongside Ea is correct rather than a
        double-count: for an ILT-based reaction the barrier physically lives in `Ea`, with
        E0(TS) = sum(E0(reactants)) + Ea, so raising the barrier necessarily raises both together.)
        The `Ea` perturbation is therefore applied here in addition to the E0 perturbation for every
        ILT-based path reaction with Arrhenius kinetics that shares this TS object (Arkane input
        files may point more than one `reaction` block at the same `transitionState` label, so a
        single TS entry can own several path reactions; since the E0 perturbation above already
        affects all of them -- the object is shared -- the Ea perturbation must match that same
        blast radius or the two representations of the barrier would desynchronize). This assumes
        each such path reaction owns its own Arrhenius kinetics object: a single Arrhenius instance
        shared across a TS's reactions is one barrier and is perturbed only once, while one shared
        across reactions of *different* transition states is an unsupported input that would
        cross-contaminate their sensitivities -- a warning is emitted if that is detected. RRKM-based
        (`can_tst()` is `True`) path reactions are unaffected, since RRKM genuinely uses the TS's
        E0 and sum of states, and perturbing both would double-count.
        """
        perturbation = self.perturbation.value_si
        if unperturb:
            perturbation *= -1
        if isinstance(entry, TransitionState):
            entry.conformer.E0 = quantity.Energy(entry.conformer.E0.value_si + perturbation, 'J/mol')
            # Shift Ea by the same amount for every ILT-based path reaction owning this TS. The gate
            # below MUST return the same verdict on a perturb call and its matching unperturb call,
            # or Ea would be left shifted (or shifted twice) and every later iteration would run on a
            # corrupted network. It trivially does: it depends only on the presence of statmech modes
            # on the TS (`can_tst`) and on the type of the kinetics object, neither of which is
            # affected by the E0/Ea shifts applied here.
            perturbed_kinetics = set()  # id()s of Arrhenius objects already shifted this call
            for rxn in self.job.network.path_reactions:
                if rxn.transition_state is not None and rxn.transition_state is entry and not rxn.can_tst():
                    kinetics = rxn.kinetics if rxn.network_kinetics is None else rxn.network_kinetics
                    if isinstance(kinetics, Arrhenius):
                        if id(kinetics) in perturbed_kinetics:
                            # Several of this TS's path reactions share one Arrhenius object: it is a
                            # single barrier, already shifted once to match the single E0 shift above,
                            # so do not shift it again (that would move Ea by 2x the E0 shift).
                            continue
                        perturbed_kinetics.add(id(kinetics))
                        kinetics.Ea = quantity.Energy(kinetics.Ea.value_si + perturbation, 'J/mol')
                    elif not unperturb:
                        # Defensive only: a non-Arrhenius ILT path reaction cannot actually reach
                        # sensitivity analysis, because apply_inverse_laplace_transform_method
                        # (rmgpy/pdep/reaction.pyx) declares a Cython-typed `Arrhenius kinetics`
                        # parameter, so the unperturbed base run in __init__ would already raise a
                        # TypeError before we ever get here. The warning is guarded by `not unperturb`
                        # so it is logged once per perturbation rather than twice (perturb + unperturb).
                        logging.warning(
                            "Not perturbing Ea for ILT-based path reaction '{0}' with kinetics of type "
                            "'{1}': the resulting TS sensitivity coefficient will be meaningless.".format(
                                rxn, kinetics.__class__.__name__))
            if not unperturb:
                # Aliased-kinetics guard: if an Arrhenius object just shifted for this TS is also the
                # kinetics of an ILT path reaction owned by a *different* TS, that reaction's rate has
                # been collaterally perturbed and its sensitivity row is contaminated. This is an
                # unsupported input (distinct reactions should not share one mutable kinetics object),
                # so surface it loudly rather than silently reporting a corrupted coefficient.
                for rxn in self.job.network.path_reactions:
                    if (rxn.transition_state is not None and rxn.transition_state is not entry
                            and not rxn.can_tst()):
                        other_kinetics = rxn.kinetics if rxn.network_kinetics is None else rxn.network_kinetics
                        if id(other_kinetics) in perturbed_kinetics:
                            logging.warning(
                                "Path reaction '{0}' (transition state '{1}') shares its Arrhenius kinetics "
                                "object with a path reaction of the perturbed transition state '{2}'; its "
                                "sensitivity coefficients will be contaminated. Give each path reaction its "
                                "own kinetics object.".format(rxn, rxn.transition_state.label, entry.label))
        elif isinstance(entry, Configuration):
            entry.species[0].conformer.E0 = quantity.Energy(entry.species[0].conformer.E0.value_si + perturbation,
                                                            'J/mol')

    def unperturb(self, entry):
        """A helper function for calling self.perturb cleanly when unperturbing"""
        self.perturb(entry, unperturb=True)

    def _classify_ilt_transition_states(self):
        """
        Return ``(ilt_labels, contaminated_labels)``, two ordered lists of transition-state labels.

        ``ilt_labels`` are the TSes whose sensitivity row carries the combined-E0+Ea (ILT) semantics
        rather than a plain dln(r)/dE0(TS): a TS qualifies when at least one of its path reactions is
        both ILT-based (``not can_tst()``) AND carries Arrhenius kinetics, which is exactly when
        perturb() shifts Ea alongside E0.

        ``contaminated_labels`` (a subset of ``ilt_labels``) are TSes whose Arrhenius kinetics object
        is shared with a path reaction of a *different* TS: perturbing one collaterally shifts the
        other, so their coefficients are contaminated. Such sharing is unsupported input; perturb()
        also warns about it at run time.

        Both save() and plot() classify rows from this single source so the table, the YAML metadata
        and the figure never disagree about which rows are ILT-based or contaminated.
        """
        ilt_labels = []
        kinetics_to_ts = {}  # id(Arrhenius) -> set of TS labels among ILT path reactions
        for rxn in self.job.network.path_reactions:
            if rxn.transition_state is not None and not rxn.can_tst():
                kinetics = rxn.kinetics if rxn.network_kinetics is None else rxn.network_kinetics
                if not isinstance(kinetics, Arrhenius):
                    continue
                if rxn.transition_state.label not in ilt_labels:
                    ilt_labels.append(rxn.transition_state.label)
                kinetics_to_ts.setdefault(id(kinetics), set()).add(rxn.transition_state.label)
        contaminated_labels = [label for label in ilt_labels
                               if any(label in ts_set and len(ts_set) > 1
                                      for ts_set in kinetics_to_ts.values())]
        return ilt_labels, contaminated_labels

    def save(self, wells, transition_states):
        """Save the SA output as tabulated data as well as in YAML format"""
        if not os.path.exists(os.path.join(self.output_directory, 'sensitivity')):
            os.mkdir(os.path.join(self.output_directory, 'sensitivity'))
        valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
        network_str = self.job.network.label
        filename = os.path.join('sensitivity', ''.join(c for c in network_str if c in valid_chars) + '.txt')
        path = os.path.join(self.output_directory, filename)
        sa_data = dict()
        sa_data['structures'] = dict()
        # Classify TS rows once (shared with plot()): which carry the combined-E0+Ea (ILT) semantics
        # rather than a plain dln(r)/dE0(TS), and which are contaminated by an Arrhenius kinetics
        # object shared across transition states. The '(TS) ' prefix matches the entry keys used below.
        ilt_raw_labels, contaminated_raw_labels = self._classify_ilt_transition_states()
        ilt_ts_labels = ['(TS) ' + label for label in ilt_raw_labels]
        contaminated_ts_labels = ['(TS) ' + label for label in contaminated_raw_labels]
        # Only emit the top-level `metadata` key when it actually carries something, so an all-RRKM
        # network's YAML keeps the pre-existing shape (`structures` + reaction strings) byte-for-byte
        # and consumers that never dealt with ILT markers see no new key.
        metadata = {}
        if ilt_ts_labels:
            metadata['ilt_transition_states'] = ilt_ts_labels
        if contaminated_ts_labels:
            metadata['contaminated_transition_states'] = contaminated_ts_labels
        if metadata:
            sa_data['metadata'] = metadata
        with open(path, 'w') as sa_f:
            sa_f.write("Sensitivity analysis for network {0}\n\n"
                       "The semi-normalized sensitivity coefficients are calculated as dln(r)/dE0\n"
                       "by perturbing E0 of each well or TS by {1}, and are given in `mol/J` units.\n"
                       "For a TS owned by an ILT-based path reaction (no statmech modes, so its\n"
                       "microcanonical rate comes from the high-pressure-limit Arrhenius kinetics),\n"
                       "the Arrhenius Ea is perturbed by the same amount alongside E0, so that row's\n"
                       "coefficient is a derivative along the coordinate that raises the barrier\n"
                       "(TS E0 and Ea together), not a plain dln(r)/dE0(TS).\n\n\n".format(network_str, self.perturbation))
            if contaminated_ts_labels:
                sa_f.write("WARNING: rows marked '(!)' belong to a transition state whose Arrhenius\n"
                           "kinetics object is shared with a path reaction of a DIFFERENT transition\n"
                           "state. Perturbing one collaterally shifts the other, so these coefficients\n"
                           "are contaminated and unreliable. Give each path reaction its own kinetics\n"
                           "object. (The same labels are listed under metadata.contaminated_transition_"
                           "states in sa_coefficients.yml.)\n\n\n")
            for rxn in self.job.network.net_reactions:
                for species in rxn.reactants + rxn.products:
                    if species.label not in sa_data['structures']:
                        sa_data['structures'][species.label] = species.to_adjacency_list()
                reactants_label = ' + '.join([reactant.label for reactant in rxn.reactants])
                products_label = ' + '.join([product.label for product in rxn.products])
                reaction_str = f'{reactants_label} <=> {products_label}'
                sa_data[reaction_str] = dict()
                sa_f.write('  Sensitivity of network reaction ' + reaction_str + ' :' + '\n')
                max_label = 40
                sa_f.write('========================={0}==================================================\n'
                           '| Well or TS {1}| Temperature (K) | Pressure (bar) | Sensitivity coefficient |\n'
                           '|------------{2}+-----------------+----------------+-------------------------|\n'
                           .format('=' * (max_label - 10), ' ' * (max_label - 10), '-' * (max_label - 10)))
                for entry in wells + transition_states:
                    if isinstance(entry, TransitionState):
                        entry_label = '(TS) ' + entry.label
                    elif isinstance(entry, Configuration):
                        entry_label = ' + '.join([species.label for species in entry.species])
                    # Flag contaminated TS rows in the table with a '(!)' suffix. The marker goes on
                    # the display label only; the YAML data key (entry_label) stays clean so consumers
                    # keying on '(TS) <label>' are unaffected -- the marker is carried in the YAML via
                    # metadata.contaminated_transition_states instead.
                    display_label = entry_label
                    if isinstance(entry, TransitionState) and entry_label in contaminated_ts_labels:
                        display_label = entry_label + ' (!)'
                    mod_entry_label = display_label + ' ' * max(1, max_label - len(display_label))
                    for i, condition in enumerate(self.conditions):
                        sa_f.write('| {0} | {1:6.1f}          | {2:8.2f}       | {3:+1.2e}               |\n'.format(
                            mod_entry_label, condition[0].value_si, condition[1].value_si * 1e-5,
                            self.sa_coefficients[str(rxn)][entry][i]))
                        condition_tuple = (condition[0].value_si, 'K', condition[1].value_si * 1e-5, 'bar')
                        if condition_tuple not in sa_data[reaction_str]:
                            sa_data[reaction_str][condition_tuple] = dict()
                        sa_data[reaction_str][condition_tuple][entry_label] = self.sa_coefficients[str(rxn)][entry][i]
                sa_f.write('========================={0}=================================================='
                           '\n\n\n'.format('=' * (max_label - 10)))

        with open(os.path.join(self.output_directory, 'sensitivity', 'sa_coefficients.yml'), 'w') as f:
            yaml.dump(data=sa_data, stream=f)

    def plot(self, wells, transition_states):
        """Draw the SA results as horizontal bars"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        for rxn in self.job.network.net_reactions:
            plt.rcdefaults()
            ax = plt.subplots(nrows=len(self.conditions), ncols=1, tight_layout=True)[1]
            # Annotate the TS tick labels so the figure distinguishes the two derivative types it
            # plots on one axis: a plain dln(k)/dE0 for wells and RRKM-based TSes, versus the
            # combined-E0+Ea derivative for ILT-based TSes (see save()). Contaminated rows (shared
            # kinetics across TSes) are flagged too, matching the '(!)' note in the text table.
            ilt_raw_labels, contaminated_raw_labels = self._classify_ilt_transition_states()
            labels = [str(entry) for entry in wells]
            for ts in transition_states:
                if ts.label in contaminated_raw_labels:
                    labels.append(ts.label + ' (E0+Ea, contaminated)')
                elif ts.label in ilt_raw_labels:
                    labels.append(ts.label + ' (E0+Ea)')
                else:
                    labels.append(ts.label)
            max_sa = min_sa = self.sa_coefficients[str(rxn)][wells[0]][0]
            for conformer_sa in self.sa_coefficients[str(rxn)].values():
                for sa_condition in conformer_sa:
                    if min_sa > sa_condition:
                        min_sa = sa_condition
                    if max_sa < sa_condition:
                        max_sa = sa_condition
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
            for i, condition in enumerate(self.conditions):
                values = [self.sa_coefficients[str(rxn)][conf][i] for conf in wells + transition_states]
                y_pos = np.arange(len(labels))
                if len(self.conditions) > 1:
                    axis = ax[i]
                else:
                    axis = ax
                axis.barh(y_pos, values, align='center', color=colors[i % len(colors)])
                axis.set_yticks(y_pos)
                axis.set_yticklabels(labels)
                axis.invert_yaxis()  # labels read top-to-bottom
                # A generic label: for wells and RRKM-based TS rows this is dln(k)/dE0, but for a
                # TS row belonging to an ILT-based path reaction it is a derivative along the
                # coordinate that raises E0 and Ea together, not a plain dln(k)/dE0(TS) (see save()).
                # The coefficient has units of mol/J (the reciprocal of the perturbed energy).
                axis.set_xlabel(r'Sensitivity coefficient ($\frac{mol}{J}$)')
                # axis.ticklabel_format('sci')
                axis.set_title('{0}, {1}'.format(condition[0], condition[1]))
                try:
                    axis.set_xlim([min_sa, max_sa])
                except:
                    logging.error("Could not set sensitivity plotting limits may be NaNs.")
                    pass
                axis.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

            if not os.path.exists(self.sensitivity_path):
                os.mkdir(self.sensitivity_path)
            valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
            reaction_str = str(rxn)
            filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
            path = os.path.join(self.sensitivity_path, filename)
            plt.savefig(path)
            plt.close()
