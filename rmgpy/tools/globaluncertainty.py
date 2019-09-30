#!/usr/bin/env python3

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

import copy
import logging
from time import time

import numpy as np

import muq.Modeling as muqm
import muq.Approximation as muqa
import muq.Utilities as muqu


class ReactorModPiece(muqm.PyModPiece):
    def __init__(self, cantera, output_species_list, k_params, k_uncertainty, g_params, g_uncertainty,
                 correlated=False, logx=True):
        """
        ======================= ====================================================
        Attribute               Description
        ======================= ====================================================
        `cantera`               A Cantera() object containing CanteraConditions and initialized species and reactions
        `output_species_list`   A list of Species() objects corresponding to the desired observables for uncertainty analysis
        `k_params`              Uncorrelated: A list of indices of the Reaction() objects in cantera.reaction_list corresponding to the uncertain input rate coefficients
                                Correlated: this is a list of strings of the uncertain kinetic parameter sources to be propagated in the model. i.e. 'H_Abstraction CHO/Oa'
        `k_uncertainty`         Uncorrelated: A list of uncertainties dlnk corresponding to the reactions in kReactions
                                Correlated: A list of dictionaries corresponding to each reaction's partial uncertainties with respect to various kinetic sources.
                                            This is the output object from the Uncertainty.kineticInputUncertainties 
        `g_params `             Uncorrelated: A list of indices for the Species() objects in cantera.species_list corresponding to the uncertain input free energies of individual species
                                Correlated: A list of strings corresponding to the uncertain thermo sources to be propagated.  i.e. 'Group(group) C=O'
        `g_uncertainty`         Uncorrelated: A list of uncertainties dG corresponding to the species in gSpecies in units of kcal/mol
                                Correlated: A list of dictionaries corresponding to each specie's partial uncertainties with respect to various thermo sources.
                                            This is the output object from the Uncertainty.thermoInputUncertainties 
        `correlated`            A flag set to either `True` or `False` depending on whether the uncertainties are correlated or not
                                If correlated, the correlated uncertainties will propagate and affect multiple species thermo and reaction kinetics parameters
        
        `affected_reactions`     Used only in the correlated case, this is a list of the indices of the affected rxns for k_params
        `affected_species`       Used only in the correlated case, this is a list of the indices of the affected species for g_params
        `logx`                  A flag to set whether to use mole fraction or ln(mole fraction) as the output variable
        ============================================================================
        """
        self.cantera = cantera
        self.output_species_list = output_species_list
        self.output_species_indices = [cantera.species_list.index(output_species) for output_species in output_species_list]
        self.correlated = correlated
        self.logx = logx

        # The size of the uncertain inputs: [parameters affecting k, parameters affecting free energy G]
        input_size = [len(k_params) + len(g_params)]

        if not self.correlated:
            # Convert input indices (RMG indices) into species list indices
            new_k_params = []
            for ind in k_params:
                for i, rxn in enumerate(cantera.reaction_list):
                    if rxn.index == ind:
                        new_k_params.append(i)
                        logging.debug('Replacing reaction index {0} with {1} for reaction {2!s}'.format(ind, i, rxn))
                        break
                else:
                    raise ValueError('Could not find requested index {0} in reaction list.'.format(ind))
            new_g_params = []
            for ind in g_params:
                for i, spc in enumerate(cantera.species_list):
                    if spc.index == ind:
                        new_g_params.append(i)
                        logging.debug('Replacing species index {0} with {1} for species {2!s}'.format(ind, i, spc))
                        break
                else:
                    raise ValueError('Could not find requested index {0} in species list.'.format(ind))
            k_params = new_k_params
            g_params = new_g_params

        # for uncorrelated case, these are indices of reactions
        # for correlated case, this is the list of labels for the uncertain parameters to be perturbed, i.e. 'H_Abstraction CHO/Oa'
        self.k_params = k_params
        # for uncorrelated case, these are indices of the species
        # for correlated case, this is the list of labels for the uncertain thermo parameters to be perturbed, i.e. 'Group(ring) cyclohexane'
        self.g_params = g_params

        if not self.correlated:
            k_uncertainty_factors = [val * np.sqrt(3) / np.log(10) for val in k_uncertainty]
            self.k_uncertainty_factors = {}
            for i, rxn_index in enumerate(k_params):
                self.k_uncertainty_factors[rxn_index] = k_uncertainty_factors[rxn_index]
                logging.debug('For {0}, set uncertainty factor to {1}'.format(cantera.reaction_list[rxn_index], k_uncertainty_factors[rxn_index]))
            
            g_uncertainty_factors = [val * np.sqrt(3) for val in g_uncertainty]
            self.g_uncertainty_factors = {}
            for i, spc_index in enumerate(g_params):
                self.g_uncertainty_factors[spc_index] = g_uncertainty_factors[spc_index]
                logging.debug('For {0}, set uncertainty factor to {1}'.format(cantera.species_list[spc_index], g_uncertainty_factors[spc_index]))
            
        else:
            # In the correlated case, keep track of which reactions and species each 
            # uncertain parameter affects 
            self.k_uncertainty_factors = {}
            affected_reactions = set()
            for k_param in k_params:
                rxn_partial_uncertainty = []
                for rxn_index, rxn_input_dict in enumerate(k_uncertainty):
                    if k_param in rxn_input_dict:
                        affected_reactions.add(rxn_index)
                        uncertainty_factor = rxn_input_dict[k_param] * np.sqrt(3) / np.log(10)
                        # If this parameter string contributes to the reaction, add this reaction index
                        rxn_partial_uncertainty.append((rxn_index, uncertainty_factor))
                self.k_uncertainty_factors[k_param] = rxn_partial_uncertainty  # list of indices of the reactions affected
            self.affected_reactions = list(affected_reactions)

            self.g_uncertainty_factors = {}
            affected_species = set()
            for g_param in g_params:
                spc_partial_uncertainty = []
                for spc_index, spc_input_dict in enumerate(g_uncertainty):
                    if g_param in spc_input_dict:
                        affected_species.add(spc_index)
                        uncertainty_factor = spc_input_dict[g_param] * np.sqrt(3)
                        # If this parameter contributes to the species thermo, add this species index
                        spc_partial_uncertainty.append((spc_index, uncertainty_factor))
                self.g_uncertainty_factors[g_param] = spc_partial_uncertainty
            self.affected_species = list(affected_species)

        # The size of the vector corresponding to the outputs to be analyzed for uncertainty analysis
        # is equal to the number of cantera conditions involved multiplied by the number of desired observables
        output_size = [len(cantera.conditions) * len(output_species_list)]

        self.num_output_species = len(output_species_list)
        self.num_conditions = len(cantera.conditions)

        # Initialize the ModPiece with some input and output size specifications
        muqm.PyModPiece.__init__(self, input_size, output_size)

    def EvaluateImpl(self, inputs):
        """
        Evaluate the desired output mole fractions based on a set of inputs
            inputs = [[k_rv], [g_rv]]
        which contains the random normal variables attributed to the uncertain
        kinetics and free energy parameters, respectively.

        Overrides methods in parent class, PyModPiece.
        
        The output returned contains
            [Condition1_outputMoleFraction1,
             Condition1_outputMoleFraction2,
             Condition2_output...,
             ConditionN_output...,]
        """
        assert len(inputs[0]) == self.inputSizes[0], "Number of inputs matches number of uncertain parameters"

        k_rv = inputs[0][0:len(self.k_params)]
        g_rv = inputs[0][len(self.k_params):]

        if not self.correlated:
            # Make deepcopies of the thermo and kinetics so as to not modify the originals in the species_list and reaction_list
            original_thermo = [copy.deepcopy(self.cantera.species_list[index].thermo) for index in self.g_params]
            original_kinetics = [copy.deepcopy(self.cantera.reaction_list[index].kinetics) for index in self.k_params]

            # Scale the thermo and kinetics of the current objects
            for i, rv in enumerate(k_rv):
                rxn_index = self.k_params[i]
                self.scale_to_kinetics(rv, self.k_uncertainty_factors[rxn_index], rxn_index)
            for i, rv in enumerate(g_rv):
                spc_index = self.g_params[i]
                self.scale_to_thermo(rv, self.g_uncertainty_factors[spc_index], spc_index)

        else:
            # Make deepcopies of the thermo and kinetics so as to not modify the originals in the species_list and reaction_list
            original_kinetics = [copy.deepcopy(self.cantera.reaction_list[index].kinetics) 
                                 for index in self.affected_reactions]
            original_thermo = [copy.deepcopy(self.cantera.species_list[index].thermo) 
                               for index in self.affected_species]

            mapped_reaction_scaling = {index: 0.0 for index in self.affected_reactions}
            mapped_species_scaling = {index: 0.0 for index in self.affected_species}

            for i, k_param in enumerate(self.k_params):
                for rxn_index, uncertainty_factor in self.k_uncertainty_factors[k_param]:
                    mapped_reaction_scaling[rxn_index] += k_rv[i] * uncertainty_factor

            for i, g_param in enumerate(self.g_params):
                for spc_index, uncertainty_factor in self.g_uncertainty_factors[g_param]:
                    mapped_species_scaling[spc_index] += g_rv[i] * uncertainty_factor

            for rxn_index, uncertainty_factor in mapped_reaction_scaling.items():
                self.scale_to_kinetics(1.0, uncertainty_factor, rxn_index)
            for spc_index, uncertainty_factor in mapped_species_scaling.items():
                self.scale_to_thermo(1.0, uncertainty_factor, spc_index)

        # The model must be refreshed when there are any thermo changes
        # kinetics can be refreshed automatically so we don't need to recreate the Solution() object.
        if np.any(g_rv):
            self.cantera.refresh_model()

        # Run the cantera simulation
        all_data = self.cantera.simulate()

        # Create a vector to hold the ModPiece output, which will be the mole fraction of the output species of interest
        output = np.zeros(self.outputSizes[0])

        # Extract the final time point for each of the mole fractions within the output_species_list

        for i in range(self.num_conditions):
            for j in range(self.num_output_species):
                species_index = self.output_species_indices[j]
                species_generic_data = all_data[i][1][2:]
                if self.logx:
                    output[i * self.num_output_species + j] = np.log(species_generic_data[species_index].data[-1])
                else:
                    output[i * self.num_output_species + j] = species_generic_data[species_index].data[-1]

        if not self.correlated:
            # Now reset the cantera object's species_list and reaction_list back to original thermo and kinetics 
            for i, thermo in enumerate(original_thermo):
                index = self.g_params[i]
                self.cantera.species_list[index].thermo = thermo

            for i, kinetics in enumerate(original_kinetics):
                index = self.k_params[i]
                self.cantera.reaction_list[index].kinetics = kinetics
        else:
            for i, thermo in enumerate(original_thermo):
                index = self.affected_species[i]
                self.cantera.species_list[index].thermo = thermo
            for i, kinetics in enumerate(original_kinetics):
                index = self.affected_reactions[i]
                self.cantera.reaction_list[index].kinetics = kinetics

        self.outputs = [output]

    def scale_to_kinetics(self, random_input, uncertainty_factor, reaction_index):
        """
        This function takes a random uniform input X = Unif(-1,1) and scales 
        the kinetics within a reaction to that value, given that the kinetics 
        has a loguniform distribution where ln(k) = Unif[ln(k_min), ln(k_max)]
        
        k_sampled = 10^(random_input * uncertainty_factor) * k0
        
        The kinetics is permanently altered in the cantera model and must be 
        reset to its original value after the evaluation is finished.
        """

        rxn = self.cantera.reaction_list[reaction_index]
        factor = random_input * uncertainty_factor

        # The rate is loguniform in k
        rxn.kinetics.change_rate(10 ** factor)
        self.cantera.modify_reaction_kinetics(reaction_index, rxn)

    def scale_to_thermo(self, random_input, uncertainty_factor, species_index):
        """
        This function takes a random normal input X = Unif(-1,1) and scales 
        the thermodynamics for a species to that value, given that the thermo 
        has a uniform distribution G = Unif(-Gmin,Gmax)
        
        G_sampled = random_input * uncertainty_factor + G0
        
        The thermo is permanently altered in the cantera model and must be 
        reset to its original value after the evaluation is finished.
        """

        species = self.cantera.species_list[species_index]
        delta_h = random_input * uncertainty_factor * 4184.0  # Convert kcal/mol to J/mol

        species.thermo.change_base_enthalpy(delta_h)
        self.cantera.modify_species_thermo(species_index, species, use_chemkin_identifier=True)


class ReactorPCEFactory(object):
    """
    This class uses MUQ to generate adaptive Polynomial Chaos Expansions for 
    global uncertainty analysis in chemical reaction systems. 
    It uses RMG, Cantera, and MUQ dependencies.
    
    Methodology
        1. Set up reactor conditions and load chemical kinetic mechanism. Select desired outputs
        2. Run local uncertainty analysis
        3. Extract top N most uncertain parameters
        4. Input these set of parameters into a MUQ model class (which is a child of the ModPiece class)
        5. Create EvaluateImpl function within model class that runs simulation based on reactor conditions through Cantera
        6. Perform PCE analysis of desired outputs
    """

    def __init__(self, cantera, output_species_list, k_params, k_uncertainty, g_params, g_uncertainty, 
                 correlated=False, logx=True):

        self.reactor_mod = ReactorModPiece(
            cantera=cantera,
            output_species_list=output_species_list,
            k_params=k_params,
            k_uncertainty=k_uncertainty,
            g_params=g_params,
            g_uncertainty=g_uncertainty,
            correlated=correlated,
            logx=logx,
        )

        # Select the polynomial and quadrature families
        quad_family = muqa.GaussPattersonQuadrature()  # We select the Gauss-Patterson quadrature as it is recommended
                                                       # as the fastest in the Patrick Conrad, Youssef Marzouk paper
        poly_family = muqa.Legendre()                  # Uniform random variables used chemical kinetics
                                                       # uncertainty propagation uses Legendre polynomials

        # Initialize the PCE Factory
        self.dim = self.reactor_mod.inputSizes[0]
        self.factory = muqa.AdaptiveSmolyakPCE(self.reactor_mod, [quad_family] * self.dim, [poly_family] * self.dim)

        self.pce = None
        self.logx = logx

    def generate_pce(self, run_time=None, start_order=2, tolerance=None, fixed_terms=False):
        """
        Generate the PCEs adaptively. There are three methods for doing so. 
        `run_time` should be given in seconds
        Option 1: Adaptive for a pre-specified amount of time
        Option 2: Adaptively construct PCE to error tolerance
        Option 3: Use a fixed order
        """

        if run_time is None and tolerance is None and not fixed_terms:
            raise ValueError('Must define at least one termination criteria')

        multis = muqu.MultiIndexFactory.CreateTotalOrder(self.dim, start_order)

        options = dict()
        options['ShouldAdapt'] = 1
        if fixed_terms:
            options['ShouldAdapt'] = 0
        if run_time is not None:
            options['MaximumAdaptTime'] = run_time
        if tolerance is not None:
            options['ErrorTol'] = tolerance

        # Also monitor the amount of time it takes
        start_time = time()
        self.pce = self.factory.Compute(multis, options)
        end_time = time()

        time_taken = end_time - start_time
        logging.info('Polynomial Chaos Expansion construction took {0:2f} seconds.'.format(time_taken))
        logging.info('Number of Model Evaluations: {0}'.format(self.factory.NumEvals()))
        logging.info('Estimated L2 Error: {0:.4e}'.format(self.factory.Error()))

    def compare_output(self, test_point, log=True):
        """
        Evaluate the PCEs against what the real output might give for a test point.
        test_point is an array of all the values in terms of factor of f
        
        Returns a tuple containing the 
        (true output mole fractions, pce output mole fractions) evaluated at the test point.
        """

        true_output = self.reactor_mod.Evaluate([test_point])[0]
        pce_output = self.pce.Evaluate([test_point])[0]

        output = ''
        for i in range(self.reactor_mod.num_conditions):
            output += """============================================================
Condition {0}
------------------------------------------------------------
{1!s}
============================================================
Condition {0} {2}Mole Fractions Evaluated at Test Point
------------------------------------------------------------
Species                      True Output          PCE Output
------------------------------------------------------------
""".format(i + 1, self.reactor_mod.cantera.conditions[i], 'Log ' if self.logx else '')

            for j, outputSpecies in enumerate(self.reactor_mod.output_species_list):
                output_index = i * self.reactor_mod.num_output_species + j
                output += '{0:<20}{1:>20.3f}{2:>20.3f}\n'.format(outputSpecies.to_chemkin(),
                                                                 true_output[output_index],
                                                                 pce_output[output_index])
            output += '============================================================\n'

        if log:
            logging.info(output)
        else:
            print(output)

        return true_output, pce_output

    def analyze_results(self, log=True):
        """
        Obtain the results: the prediction mean and variance, as well as the global sensitivity indices
        Returns a tuple containing the following statistics
        
        (mean species mole fractions, variance, covariance, main sensitivity indices, total sensitivity indices)
        """
        # Compute the mean and variance for each of the uncertain parameters
        mean = np.array(self.pce.Mean())

        var = np.array(self.pce.Variance())
        stddev = np.sqrt(var)
        stddev_percent = stddev / mean * 100.0

        cov = self.pce.Covariance()

        # Extract the global sensitivity indices
        main_sens = np.array(self.pce.MainSensitivity())
        total_sens = np.array(self.pce.TotalSensitivity())

        output = ''
        for i in range(self.reactor_mod.num_conditions):
            output += """============================================================
Condition {0}
------------------------------------------------------------
{1!s}
============================================================
Condition {0} {2}Mole Fractions
------------------------------------------------------------
Species                   Mean         Stddev     Stddev (%)
------------------------------------------------------------
""".format(i + 1, self.reactor_mod.cantera.conditions[i], 'Log ' if self.logx else '')

            for j, output_species in enumerate(self.reactor_mod.output_species_list):
                output_index = i * self.reactor_mod.num_output_species + j
                output += '{0:<15}{1:>15.3e}{2:>15.3e}{3:>15.3f}\n'.format(output_species.to_chemkin(),
                                                                           mean[output_index],
                                                                           stddev[output_index],
                                                                           stddev_percent[output_index])
            output += '============================================================\n\n'

            if self.reactor_mod.k_params:
                output += """====================================================================================================
Condition {0} Reaction Rate Sensitivity Indices
----------------------------------------------------------------------------------------------------
Description                                                                 sens_main     sens_total
""".format(i + 1)

                for j, output_species in enumerate(self.reactor_mod.output_species_list):
                    output += '----------------------------------------------------------------------------------------------------\n'
                    output_index = i * self.reactor_mod.num_output_species + j
                    for k, descriptor in enumerate(self.reactor_mod.k_params):
                        parameter_index = k
                        if not self.reactor_mod.correlated:
                            description = 'dln[{0}]/dln[{1}]'.format(output_species.to_chemkin(),
                                                                     self.reactor_mod.cantera.reaction_list[descriptor].to_chemkin(kinetics=False))
                        else:
                            description = 'dln[{0}]/dln[{1}]'.format(output_species.to_chemkin(), descriptor)

                        output += '{0:<70}{1:>14.3f}%{2:>14.3f}%\n'.format(description,
                                                                           100 * main_sens[output_index][parameter_index],
                                                                           100 * total_sens[output_index][parameter_index])
                output += '====================================================================================================\n\n'

            if self.reactor_mod.g_params:
                output += """====================================================================================================
Condition {0} Thermochemistry Sensitivity Indices
----------------------------------------------------------------------------------------------------
Description                                                                 sens_main     sens_total
""".format(i + 1)
                for j, output_species in enumerate(self.reactor_mod.output_species_list):
                    output += '----------------------------------------------------------------------------------------------------\n'
                    output_index = i * self.reactor_mod.num_output_species + j
                    for g, descriptor in enumerate(self.reactor_mod.g_params):
                        parameter_index = len(self.reactor_mod.k_params) + g
                        if not self.reactor_mod.correlated:
                            description = 'dln[{0}]/dG[{1}]'.format(output_species.to_chemkin(),
                                                                    self.reactor_mod.cantera.species_list[descriptor].to_chemkin())
                        else:
                            description = 'dln[{0}]/dG[{1}]'.format(output_species.to_chemkin(), descriptor)

                        output += '{0:<70}{1:>14.3f}%{2:>14.3f}%\n'.format(description,
                                                                           100 * main_sens[output_index][parameter_index],
                                                                           100 * total_sens[output_index][parameter_index])
                output += '====================================================================================================\n\n'

        if log:
            logging.info(output)
        else:
            print(output)

        return mean, var, cov, main_sens, total_sens
