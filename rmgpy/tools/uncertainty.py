#!/usr/bin/env python3

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

import os
import pickle

import numpy as np
import warnings

import rmgpy.data.thermo
import rmgpy.util as util
from rmgpy.species import Species
from rmgpy.tools.data import GenericData
from rmgpy.tools.plot import parse_csv_data, plot_sensitivity, ReactionSensitivityPlot, ThermoSensitivityPlot


# This represents the ratio k_max/k_0 for a reaction if using a uniform distribution for uncertainty
# or the the ratio k_95/k_0 where k_95 is the 95% confidence bound for a lognormal distribution.
KINETIC_RANK_MULTIPLIERS = {
    1: 1.2,     # Experiment/FCI
    2: 1.5,     # W4/HEAT with very good (2-d if necessary) rotors
    3: 2.0,     # CCSD(T)-F12/cc-PVnZ with n>2 or CCSD(T)-F12/CBS with good (2-d if necessary) rotors
    4: 3.0,     # CCSD(T)-F12/DZ, with good (2-d if necessary) rotors
    5: 5.0,     # CBS-QB3 with 1-d rotors
    6: 7.5,     # Double-hybrid DFT with 1-D rotors
    7: 10.0,    # Hybrid DFT (w/ dispersion) (rotors if necessary)
    8: 30,      # B3LYP & lower DFT (rotors if necessary)
    9: 100,     # Direct Estimate/Guess
    10: 1000,   # Average of Rates
    11: 10000,  # General Estimate (Never used in generation)
    0: 10000,   # General Estimate (Never used in generation)
}

# Use this for uniform distribution
KINETIC_RANK_UNCERTAINTIES = {k: np.log(v) / np.sqrt(3.0) for k, v in KINETIC_RANK_MULTIPLIERS.items()}

# # Use this for lognormal distribution
# KINETIC_RANK_UNCERTAINTIES = {k: np.log(v) / 1.96 for k, v in KINETIC_RANK_MULTIPLIERS.items()}

THERMO_RANK_UNCERTAINTIES = {  # THESE ARE FILLER. PLEASE UPDATE WITH BETTER UNCERTAINTIES BASED ON DATA ANALYSIS
    1: 0.1,     # Experiment/FCI
    2: 0.2,     # W4/HEAT with very good (2-d if necessary) rotors
    3: 0.3,     # CCSD(T)-F12/cc-PVnZ with n>2 or CCSD(T)-F12/CBS with good (2-d if necessary) rotors
    4: 0.5,     # CCSD(T)-F12/DZ, with good (2-d if necessary) rotors
    5: 0.8,     # CBS-QB3 with 1-d rotors
    6: 1.0,     # Double-hybrid DFT with 1-D rotors
    7: 1.5,     # Hybrid DFT (w/ dispersion) (rotors if necessary)
    8: 2.0,     # B3LYP & lower DFT (rotors if necessary)
    9: 2.5,     # Direct Estimate/Guess
    10: 3.0,    # Average of Rates
    11: 3.0,    # General Estimate (Never used in generation)
    0: 3.0,     # General Estimate (Never used in generation)
}


class ThermoParameterUncertainty(object):
    """
    This class is an engine that generates the species uncertainty based on its thermo sources.
    """

    def __init__(self, dG_library=1.5, dG_QM=3.0, dG_GAV=1.5, dG_group=0.7159, dG_ADS_correction=6.918, dG_surf_lib=6.918, other_covariances=None, rank_dictionary=True):
        """
        Initialize the different uncertainties dG_library, dG_QM, dG_GAV, and dG_other with set values
        in units of kcal/mol.

        We expect a uniform distribution for some species free energy G in [Gmin, Gmax].
        dG = (Gmax-Gmin)/2

        rank_dictionary is set to True or False. If True, uses the default RMG ranks in THERMO_RANK_UNCERTAINTIES to assign rank-based uncertainties.
        """
        self.dG_library = dG_library
        self.dG_QM = dG_QM
        self.dG_GAV = dG_GAV
        self.dG_group = dG_group
        self.dG_ADS_correction = dG_ADS_correction
        self.dG_surf_lib = dG_surf_lib
        self.other_covariances = other_covariances  # storage of covariances as a dict. Keys are sorted tuples of parameter labels and values are covariances
        if rank_dictionary:  # default is to use default rank dictionary
            self.rank_dictionary = THERMO_RANK_UNCERTAINTIES
        else:  # ignore rank
            self.rank_dictionary = None

    def _get_rank_uncertainty(self, entry, default_uncertainty):
        """Helper function for returning the rank-based uncertainty for a given item,
        which may be a training reaction or a rate rule entry.
        If the item does not have a rank attribute or if the rank is not in the rank dictionary, return the default uncertainty."""
        if not self.rank_dictionary:
            return default_uncertainty
        rank = getattr(entry, 'rank', None)
        if rank is None or rank not in self.rank_dictionary:
            return default_uncertainty
        return self.rank_dictionary[rank]

    def get_uncertainty_value(self, source):
        """
        Retrieve the uncertainty value in kcal/mol when the source of the thermo of a species is given.
        """
        varG = 0.0
        if 'Library' in source:
            library_entry = source['Library'][1]
            default_uncertainty = self.dG_library
            dG_library = self._get_rank_uncertainty(library_entry, default_uncertainty)
            varG += dG_library * dG_library
        if 'Surface_Library' in source:
            surf_lib_entry = source['Surface_Library'][1]
            default_uncertainty = self.dG_surf_lib
            dG_surf_lib = self._get_rank_uncertainty(surf_lib_entry, default_uncertainty)
            surf_lib_varG = dG_surf_lib * dG_surf_lib

            # covariance libraries should overrule the default uncertainties when available,
            if self.other_covariances is not None:
                label = source["Surface_Library"][2]  # match the covariance dict label format
                cov_label = (label, label)
                if cov_label in self.other_covariances:
                    surf_lib_varG = self.other_covariances[cov_label]
            varG += surf_lib_varG  # Add the variance of the surface library parameter if covariance is not specified in the covariance libraries
        if 'QM' in source:
            varG += self.dG_QM * self.dG_QM
        if 'GAV' in source:
            varG += self.dG_GAV * self.dG_GAV  # Add a fixed uncertainty for the GAV method
            for group_type, group_entries in source['GAV'].items():
                group_weights = [groupTuple[-1] for groupTuple in group_entries]
                varG += np.sum([weight * weight * self.dG_group * self.dG_group for weight in group_weights])
        if 'ADS' in source:
            varG += self.dG_ADS_correction * self.dG_ADS_correction  # Add adsorption correction uncertainty

        return np.sqrt(varG)

    def get_partial_uncertainty_value(self, source, corr_source_type, corr_param=None, corr_group_type=None):
        """
        Obtain the partial uncertainty dG/dG_corr*dG_corr, where dG_corr is the correlated parameter

        `corr_param` is the parameter identifier itself, which is
            (method_name, species) for QM
            (library_name, library_entry, entry_label) for library parameters
            or a string for group values
        `corr_source_type` is a string, being either 'Library', 'QM', 'GAV', 'ADS', or 'Estimation'
        `corr_group_type` is a string used only when the source type is 'GAV' and indicates grouptype
        """

        if corr_source_type == 'Library':
            if 'Library' in source:  # check if same library and index match, could do isomorphism check on the entries, but this is faster
                if source['Library'][0] == corr_param[0] and source['Library'][1].index == corr_param[1].index:
                    # Correlated parameter is a source of the overall parameter
                    library_entry = source['Library'][1]
                    dG_library = self._get_rank_uncertainty(library_entry, self.dG_library)
                    return dG_library

        elif corr_source_type == 'Surface_Library':
            if 'Surface_Library' in source:
                if source['Surface_Library'][0] == corr_param[0] and source['Surface_Library'][1].index == corr_param[1].index:
                    # Correlated parameter is a source of the overall parameter
                    surf_lib_entry = source['Surface_Library'][1]
                    dG_surf_lib = self._get_rank_uncertainty(surf_lib_entry, self.dG_surf_lib)
                    return dG_surf_lib

        elif corr_source_type == 'QM':
            if 'QM' in source:  # here we do an isomorphism check since no index is saved for QM
                if source['QM'][0] == corr_param[0] and source['QM'][1].is_isomorphic(corr_param[1]):
                    # Correlated parameter is a source of the overall parameter
                    return self.dG_QM

        elif corr_source_type == 'ADS':
            if 'ADS' in source:
                if corr_group_type in source['ADS']:
                    group_list = source['ADS'][corr_group_type]
                    for group, weight in group_list:  # there should be only one group entry for adsorption corrections
                        if group == corr_param:
                            return weight * self.dG_ADS_correction

        elif corr_source_type == 'GAV':
            if 'GAV' in source:
                if corr_group_type in source['GAV']:
                    group_list = source['GAV'][corr_group_type]
                    for group, weight in group_list:
                        if group == corr_param:
                            return weight * self.dG_group

        elif corr_source_type == 'Estimation':
            if 'GAV' in source:
                return self.dG_GAV

        else:
            raise Exception('Thermo correlated source must be GAV, QM, Library, Surface_Library, ADS, or Estimation')

        # If we get here, it means the correlated parameter was not found
        return None

    def get_intermediate_uncertainty_value(self, source_type, source):
        """
        Get the uncertainty of an intermediate parameter, which is an underlying parameter that contributes to the overall uncertainty but is not directly a source parameter. For example, this could be the uncertainty of a specific group contribution value that contributes to the overall GAV uncertainty for a species.

        `source_type` is a string indicating the type of intermediate parameter, such as 'Library', 'Surface_Library', 'QM', 'GAV', or 'ADS'
        `source` is the identifier for the intermediate parameter, such as (method_name, species) for QM or (library_name, library_entry, entry_label) for library parameters
        """
        if source_type == 'Library':
            library_entry = source[1]
            return self._get_rank_uncertainty(library_entry, self.dG_library)
        elif source_type == 'Surface_Library':
            surf_lib_entry = source[1]
            return self._get_rank_uncertainty(surf_lib_entry, self.dG_surf_lib)
        elif source_type == 'QM':
            return self.dG_QM
        elif source_type == 'GAV':
            return self.dG_group
        elif source_type == 'ADS':
            return self.dG_ADS_correction
        elif source_type == 'Estimation':
            return self.dG_GAV
        else:
            raise Exception('Thermo intermediate parameter source must be Library, Surface_Library, QM, GAV, or ADS')

    def get_uncertainty_factor(self, source):
        """
        Retrieve the uncertainty factor f in kcal/mol when the source of the thermo of a species is given.

        This is equivalent to sqrt(3)*dG in a uniform uncertainty interval
        """
        dG = self.get_uncertainty_value(source)
        f = np.sqrt(3) * dG
        return f
    
    def _get_covariance_qq(self, q_label1, q_label2, dG_source1=None):
        """
        Gets the covariance between two intermediate sources q1 and q2 in units of (kcal/mol)^2
        Where q_label1 and q_label2 are both labels for a given thermo source and dG_source1 is the partial uncertainty of source 1

        The possible intermediate parameter types are:
        Library, Surface_Library, QM, Estimation, AdsorptionCorrection, or Group

        Sources are assumed to be independent unless the labels match exactly or an explicit covariance exists for the two
        """
        # Explicit covariance value overrides all other possibilities
        if self.other_covariances is not None:
            # check if covariance is specified in other_covariances dict
            sorted_labels = tuple(sorted([q_label1, q_label2]))
            if sorted_labels in self.other_covariances:
                return self.other_covariances[sorted_labels]

        if q_label1 == q_label2:
            # If the two correlated parameters are exactly the same, return the variance of that parameter
            return dG_source1 * dG_source1
        return 0


class KineticParameterUncertainty(object):
    """
    This class is an engine that generates the reaction uncertainty based on its kinetic sources.
    """

    def __init__(self, dlnk_library=0.5, dlnk_training=0.5, dlnk_pdep=2.0, dlnk_family=1.0, dlnk_nonexact=3.5,
                 dlnk_rule=0.5, dlnk_surf_library=2.659, dlnk_surf_training=2.659, dlnk_surf_rule=2.659, rank_dictionary=True, other_covariances=None):
        """
        Initialize the different uncertainties dlnk

        We expect a uniform distribution for some reaction kinetics  about ln(k0) in [ln(kmin), ln(kmax)].
        dlnk = (ln(kmax)-ln(kmin))/2

        rank_dictionary is set to True or False. If True, uses the default RMG ranks in KINETIC_RANK_UNCERTAINTIES to assign rank-based uncertainties."""
        self.dlnk_library = dlnk_library
        self.dlnk_training = dlnk_training
        self.dlnk_pdep = dlnk_pdep
        self.dlnk_family = dlnk_family
        self.dlnk_nonexact = dlnk_nonexact
        self.dlnk_rule = dlnk_rule
        self.dlnk_surf_library = dlnk_surf_library
        self.dlnk_surf_training = dlnk_surf_training
        self.dlnk_surf_rule = dlnk_surf_rule
        self.other_covariances = other_covariances  # storage of covariances as a dict. Keys are sorted tuples of parameter labels and values are covariances
        if rank_dictionary:
            self.rank_dictionary = KINETIC_RANK_UNCERTAINTIES  # use default rank dictionary
        else:  # ignore rank
            self.rank_dictionary = None

    def _get_rank_uncertainty(self, entry, default_uncertainty):
        """Helper function for returning the rank-based uncertainty for a given item,
        which may be a training reaction or a rate rule entry.
        If the item does not have a rank attribute or if the rank is not in the rank dictionary, return the default uncertainty."""
        if not self.rank_dictionary:
            return default_uncertainty
        rank = getattr(entry, 'rank', None)
        if rank is None or rank not in self.rank_dictionary:
            return default_uncertainty
        return self.rank_dictionary[rank]

    def get_uncertainty_value(self, source):
        """
        Retrieve the dlnk uncertainty when the source of the reaction kinetics are given
        """
        varlnk = 0.0
        if 'Library' in source:
            # Should be a single library reaction source
            library_entry = source['Library'][1]
            default_uncertainty = self.dlnk_library
            dlnk_library = self._get_rank_uncertainty(library_entry, default_uncertainty)
            varlnk += dlnk_library * dlnk_library
        elif 'Surface_Library' in source:
            # Should be a single library reaction source
            surface_library_entry = source['Surface_Library'][1]
            default_uncertainty = self.dlnk_surf_library
            dlnk_surf_library = self._get_rank_uncertainty(surface_library_entry, default_uncertainty)
            varlnk += dlnk_surf_library * dlnk_surf_library
        elif 'PDep' in source:
            # Should be a single pdep reaction source
            varlnk += self.dlnk_pdep * self.dlnk_pdep
        elif 'Training' in source:
            # Should be a single training reaction
            # Although some training entries may be used in reverse,
            # We still consider the kinetics to be directly dependent
            training_entry = source['Training'][1]
            if training_entry.item.is_surface_reaction():
                default_uncertainty = self.dlnk_surf_training
            else:
                default_uncertainty = self.dlnk_training
            dlnk_training = self._get_rank_uncertainty(training_entry, default_uncertainty)
            varlnk += dlnk_training * dlnk_training
        elif 'Rate Rules' in source:
            family_label = source['Rate Rules'][0]
            source_dict = source['Rate Rules'][1]
            exact = source_dict['exact']
            rule_weights = [ruleTuple[-1] for ruleTuple in source_dict['rules']]
            training_weights = [trainingTuple[-1] for trainingTuple in source_dict['training']]

            default_rule_uncertainty = self.dlnk_surf_rule if 'surface' in family_label.lower() else self.dlnk_rule
            # Note that even if these source from training reactions, we actually
            # use the uncertainty for rate rules, since these are now approximations
            # of the original reaction. We consider these to be independent of original the training
            # parameters because the rate rules may be reversing the training reactions,
            # which leads to more complicated dependence
            rule_uncertainties = [self._get_rank_uncertainty(ruleEntry, default_rule_uncertainty) for ruleEntry, weight in source_dict['rules']]
            training_uncertainties = [self._get_rank_uncertainty(trainingEntry, default_rule_uncertainty) for ruleEntry, trainingEntry, weight in source_dict['training']]

            varlnk += self.dlnk_family * self.dlnk_family

            N = len(rule_weights) + len(training_weights)
            if not exact:
                # nonexactness contribution increases as N increases
                varlnk += (np.log10(N + 1) * self.dlnk_nonexact) * (np.log10(N + 1) * self.dlnk_nonexact)

            # Add the contributions from rules
            varlnk += np.sum([weight * weight * rule_uncertainty * rule_uncertainty for weight, rule_uncertainty in zip(rule_weights, rule_uncertainties)])
            
            # Add the contributions from training
            varlnk += np.sum([weight * weight * training_uncertainty * training_uncertainty for weight, training_uncertainty in zip(training_weights, training_uncertainties)])

        return np.sqrt(varlnk)

    def get_partial_uncertainty_value(self, source, corr_source_type, corr_param=None, corr_family=None):
        """
        Obtain the partial uncertainty dlnk/dlnk_corr*dlnk_corr, where dlnk_corr is the correlated parameter

        `corr_param` is the parameter identifier itself, which is the string identifier of the rate rule
        `corr_source_type` is a string, being either 'Rate Rules', 'Library', 'PDep', 'Training', 'Estimation Nonexact', or 'Estimation Family'
        `corr_family` is a string used only when the source type is 'Rate Rules' and indicates the family
        """

        if corr_source_type == 'Rate Rules':
            if 'Rate Rules' in source:
                family_label = source['Rate Rules'][0]
                if corr_family == family_label:
                    surface_family = 'surface' in family_label.lower()
                    source_dict = source['Rate Rules'][1]
                    rules = source_dict['rules']
                    training = source_dict['training']
                    if rules:
                        for ruleEntry, weight in rules:
                            if corr_param == ruleEntry:
                                if surface_family:
                                    dlnk_rule = self.dlnk_surf_rule
                                else:
                                    dlnk_rule = self.dlnk_rule

                                # get rank-based uncertainty if rank_dictionary is being used
                                dlnk_rule = self._get_rank_uncertainty(ruleEntry, dlnk_rule)

                                return weight * dlnk_rule
                    if training:
                        for ruleEntry, trainingEntry, weight in training:
                            if corr_param == ruleEntry:
                                # We use the uncertainty for rate rules, since these are now approximations
                                # of the original reaction. We consider these to be independent of original the training
                                # parameters because the rate rules may be reversing the training reactions,
                                # which leads to more complicated dependence
                                if surface_family:
                                    dlnk_rule = self.dlnk_surf_rule
                                else:
                                    dlnk_rule = self.dlnk_rule

                                # get rank-based uncertainty if rank_dictionary is being used
                                dlnk_rule = self._get_rank_uncertainty(ruleEntry, dlnk_rule)

                                return weight * dlnk_rule

        # Writing it this way in the function is not the most efficient, but makes it easy to use, and
        # testing a few if statements is not too costly
        elif corr_source_type == 'Library':
            if 'Library' in source:
                if corr_param == source['Library']:
                    # Should be a single library reaction source
                    library_entry = source['Library'][1]
                    return self._get_rank_uncertainty(library_entry, self.dlnk_library)
        elif corr_source_type == 'Surface_Library':
            if 'Surface_Library' in source:
                if corr_param == source['Surface_Library']:
                    # Should be a single library reaction source
                    surface_library_entry = source['Surface_Library'][1]
                    return self._get_rank_uncertainty(surface_library_entry, self.dlnk_surf_library)
        elif corr_source_type == 'PDep':
            if 'PDep' in source:
                if corr_param == source['PDep']:
                    return self.dlnk_pdep
        elif corr_source_type == 'Training':
            if 'Training' in source:
                # Should be a unique single training reaction
                if corr_param == source['Training']:
                    training_entry = source['Training'][1]
                    if training_entry.item.is_surface_reaction():
                        return self._get_rank_uncertainty(training_entry, self.dlnk_surf_training)
                    else:
                        return self._get_rank_uncertainty(training_entry, self.dlnk_training)
        elif corr_source_type == 'Estimation Nonexact':
            # Return the uncorrelated uncertainty associated with using a non-exact rate rule
            if 'Rate Rules' in source:
                source_dict = source['Rate Rules'][1]
                exact = source_dict['exact']
                N = len(source_dict['rules']) + len(source_dict['training'])
                if not exact:
                    # nonexactness contribution increases as N increases
                    return np.log10(N + 1) * self.dlnk_nonexact
        elif corr_source_type == 'Estimation Family':
            # Return the uncertainty associated with using a family decision tree
            if 'Rate Rules' in source:
                return self.dlnk_family
        else:
            raise Exception('Kinetics correlated source must be Rate Rules, Library, PDep, Training, Estimation Nonexact, or Estimation Family')

        # If we get here, it means that we did not find the correlated parameter in the source
        return None

    def get_intermediate_uncertainty_value(self, source_type, source):
        """
        Get the uncertainty of an intermediate parameter, which is an underlying parameter that contributes to the overall uncertainty but is not directly a source parameter. For example, this could be the uncertainty of a specific rate rule that contributes to the overall uncertainty for a reaction.

        `source_type` is a string indicating the type of intermediate parameter, such as 'Rate Rules', 'Library', 'PDep', or 'Training'
        `source` is the identifier for the intermediate parameter, such as the string identifier of the rate rule
        """
        if source_type == 'Rate Rules':
            return self._get_rank_uncertainty(source, self.dlnk_rule)
        elif source_type == 'Library':
            return self._get_rank_uncertainty(source, self.dlnk_library)
        elif source_type == 'Surface_Library':
            return self._get_rank_uncertainty(source, self.dlnk_surf_library)
        elif source_type == 'PDep':
            return self.dlnk_pdep
        elif source_type == 'Training':
            return self._get_rank_uncertainty(source, self.dlnk_training)
        elif source_type == 'Estimation Nonexact':
            return self.dlnk_nonexact
        elif source_type == 'Estimation Family':
            return self.dlnk_family
        else:
            raise Exception('Kinetic intermediate parameter source must be Rate Rules, Library, PDep, or Training')


    def get_uncertainty_factor(self, source):
        """
        Retrieve the uncertainty factor f when the source of the reaction kinetics are given.

        This is equivalent to sqrt(3)/ln(10) * dlnk  in a uniform uncertainty interval
        """
        dlnk = self.get_uncertainty_value(source)
        f = np.sqrt(3) / np.log(10) * dlnk
        return f
    
    def _get_covariance_qq(self, q_label1, q_label2, dlnk_source1=None):
        """
        Gets the covariance between two intermediate sources q1 and q2
        Where q_label1 and q_label2 are both labels for a given kinetics source and dlnk_source1 is the partial uncertainty of source 1

        The possible intermediate parameter types are:
        Library, Surface_Library, PDEP, Training, Rate Rule,
        Estimation Family, Estimation Nonexact, Surface Training, and Surface Rate Rule
        """
        # covariance libraries should overrule the default uncertainties when available
        if self.other_covariances is not None:
            # check if covariance is specified in other_covariances dict
            sorted_labels = tuple(sorted([q_label1, q_label2]))
            if sorted_labels in self.other_covariances:
                return self.other_covariances[sorted_labels]

        if q_label1 == q_label2:
            # If the two correlated parameters are exactly the same, return the variance of that parameter
            return dlnk_source1 * dlnk_source1
        
        return 0


class Uncertainty(object):
    """
    This class contains functions associated with running uncertainty analyses
    for a single RMG-generated mechanism.
    """

    def __init__(self, species_list=None, reaction_list=None, output_directory='', thermo_covariance_libraries=None,
                 thermo_covariance_groups=None, kinetic_covariance_rules=None):
        """
        `species_list`: list of RMG species objects
        `reaction_list`: list of RMG reaction objects
        `outputDirectoy`: directory path for saving output files from the analyses
        `thermo_covariance_libraries`: list of library paths to pull additional thermo covariances from
        `thermo_covariance_groups`: list of groups to get additional thermo covariances from
        """
        self.database = None
        self.species_list = species_list
        self.reaction_list = reaction_list
        self.species_sources_dict = None
        self.reaction_sources_dict = None
        self.all_thermo_sources = None
        self.all_kinetic_sources = None
        self.thermo_input_uncertainties = None          # thermo parameter uncertainties
        self.kinetic_input_uncertainties = None         # kinetic parameter uncertainties
        self.thermo_intermediate_uncertainties = None   # list of dictionaries of intermediate parameter uncertainties for each species
        self.kinetic_intermediate_uncertainties = None  # list of dictionaries of intermediate parameter uncertainties for each reaction
        self.dG_dqs = None                              # derivatives of Gibbs free energy with respect to underlying thermo parameters
        self.dlnk_dqs = None                            # derivatives of ln(k) with respect to underlying kinetic parameters
        self.thermo_covariance_matrix = None            # covariance matrix of all species thermo uncertainties
        self.kinetic_covariance_matrix = None           # covariance matrix of all reaction kinetic uncertainties
        self.Sigma_ww_thermo = None                     # covariance matrix of all underlying thermo parameter uncertainties
        self.Sigma_ww_kinetics = None                   # covariance matrix of all underlying kinetics parameter uncertainties
        self.thermo_intermediates = None                # list of labels of underlying thermo parameters (may be truncated to relevant parameters based on sensitivity)
        self.kinetic_intermediates = None               # list of labels of underlying kinetic parameters (may be truncated to relevant parameters based on sensitivity)
        self.output_directory = output_directory if output_directory else os.getcwd()
        self.thermo_covariance_libraries = thermo_covariance_libraries
        self.thermo_covariance_groups = thermo_covariance_groups
        self.thermo_covariances_dict = {}  # dictionary to store covariances from covariance libraries
        self.used_rank = False  # keep track of whether we apply rank-based uncertainties
        self.kinetic_covariance_rules = kinetic_covariance_rules  # autogenerated tree covariances
        self.kinetic_covariances_dict = {}

        # Make output directory if it does not yet exist:
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except OSError:
                raise Exception('Uncertainty output directory could not be created.')

    def load_database(self, kinetics_families='all', kinetics_depositories=None, thermo_libraries=None, reaction_libraries=None):
        """
        This function loads a single copy of the RMGDatabase with full verbose averaging
        of the rate rule to trace kinetics sources.

        By default, this function loads all the kinetics families, only the training kinetics depository,
        the primaryThermoLibrary, and no reaction libraries.
        """
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy import settings

        if not kinetics_depositories:
            kinetics_depositories = ['training']
        if not thermo_libraries:
            thermo_libraries = ['primaryThermoLibrary']
        if not reaction_libraries:
            reaction_libraries = []

        self.database = RMGDatabase()
        self.database.load(
            settings['database.directory'],
            kinetics_families=kinetics_families,
            kinetics_depositories=kinetics_depositories,
            thermo_libraries=thermo_libraries,
            reaction_libraries=reaction_libraries,
        )

        # Prepare the database by loading training reactions but not averaging the rate rules
        for familyLabel, family in self.database.kinetics.families.items():
            if not family.auto_generated:
                family.add_rules_from_training(thermo_database=self.database.thermo)
                family.fill_rules_by_averaging_up(verbose=True)

    def load_model(self, chemkin_path, dictionary_path, transport_path=None, surface_path=None):
        """
        Load a RMG-generated model into the Uncertainty class
        `chemkin_path`: path to the chem_annotated.inp CHEMKIN mechanism
        `dictionary_path`: path to the species_dictionary.txt file
        `transport_path`: path to the tran.dat file (optional)

        Then create dictionaries stored in self.thermoGroups and self.rateRules
        containing information about the source of the thermodynamic and kinetic
        parameters
        """
        from rmgpy.chemkin import load_chemkin_file

        self.species_list, self.reaction_list = load_chemkin_file(chemkin_path,
                                                                  dictionary_path=dictionary_path,
                                                                  transport_path=transport_path,
                                                                  surface_path=surface_path)

    def load_thermo_covariances_from_libraries(self):
        """
        This function populates the self.thermo_covariances_dict with covariance data (in units of (kcal/mol)^2) from the given covariance libraries

        For each library, it expects:
        1. a covariance.npy file containing the thermo covariance matrix and 
        2. a species_dictionary.txt file containing the species corresponding to the covariance data. In the same order.

        See the RMG-database/scripts/compile_BEEF_cov.ipynb Jupyter notebook for more details on how to generate these covariance libraries.

        This function only adds covariance data for species that are actually in the model, (which may still include a subset of library species that are not in the model as in the case of the radical/HBI correction)
        and only for the thermo source associated with that library. The goal is to keep the dictionary as small as possible because the lookups scale badly.

        Note: the covariance.npy matrix is in units of (kJ/mol)^2, but gets converted to (kcal/mol)^2 in this function to match the rest of the analysis
        """
        from rmgpy.chemkin import load_species_dictionary
        if self.database is None:
            raise RuntimeError('Must load database before loading covariance libraries, since we need the path to the covariance libraries from the database')
        if self.thermo_covariance_libraries is not None:
            for cov_lib in self.thermo_covariance_libraries:
                library_name = os.path.basename(cov_lib)
                if library_name in self.database.thermo.libraries:
                    library = self.database.thermo.libraries[library_name]
                else:
                    raise ValueError(f'Thermo covariance library {library_name} not found in the loaded database')
                covariance_file = os.path.join(cov_lib, 'covariance.npy')
                covariance_species = os.path.join(cov_lib, 'species_dictionary.txt')

                if not os.path.isfile(covariance_file):
                    raise ValueError(f'Thermo covariance file {covariance_file} not found in library {cov_lib}')
                if not os.path.isfile(covariance_species):
                    raise ValueError(f'Thermo species file {covariance_species} not found in library {cov_lib}')

                # warn the user if the covariance library is older than the thermo library, since this likely means the covariance data is out of date
                covariance_data_time = os.path.getmtime(covariance_file)
                library_file_path = os.path.join(rmgpy.settings['database.directory'], 'thermo', 'libraries', f'{library_name}.py')
                if not os.path.isfile(library_file_path):
                    # just warn the user instead of raising an error, since the covariance data may still be relevant even if the library file is not available
                    warnings.warn(f'Could not find the file for library {library_name} to compare modification times with the covariance library {cov_lib}. It may or may not be up to date.')
                else:
                    library_file_time = os.path.getmtime(library_file_path)
                    if covariance_data_time < library_file_time:
                        warnings.warn(f'Thermo covariance library {cov_lib} is older than the thermo library {library_name}, which may mean the covariance data is out of date and not consistent with the current library data')

                # Load covariance data and species
                cov_data = np.load(covariance_file) / 4.184 / 4.184  # convert from (kJ/mol)^2 to (kcal/mol)^2
                cov_species_dict = load_species_dictionary(covariance_species)
                cov_specs = [item for _, item in cov_species_dict.items()]
                
                # quick check to make sure the covariance data and molecule data are consistent with each other
                if cov_data.shape[0] != len(cov_specs):
                    raise ValueError(f'Covariance data and molecule data in library {cov_lib} are inconsistent: covariance data has shape {cov_data.shape} but molecule data has length {len(cov_specs)}')

                # load the labels, but only include species in the model
                subset_indices = []  # keep track of indices relevant to the model
                for i_lib, lib_species in enumerate(cov_specs):
                    i_sp = get_i_thing(lib_species, self.species_list)
                    if i_sp < 0:
                        continue

                    # make sure the species actually comes from this library, otherwise skip
                    result = self.database.thermo.get_thermo_data_from_library(lib_species, library)
                    if result is not None:
                        surface_prefix = 'Surface_' if lib_species.contains_surface_site() else ''
                        # match the label as constructed in assign_intermediate_uncertainties,
                        # where the number corresponds to the index of the species in species_list
                        try:
                            label = f'{surface_prefix}Library {self.species_list[i_sp].to_chemkin()}'
                        except IndexError:
                            # uses InChI now
                            label = f'{surface_prefix}Library {self.species_list[i_sp].molecule[0].to_inchi()}'
                        lib_species.label = label
                        subset_indices.append(i_lib)

                # fill in the dictionary of covariances from the covariance libraries,
                # with keys being sorted tuples of the labels of the correlated parameters
                # and values being the covariance between those parameters
                # only go through upper triangle of the covariance matrix, since the other half is just the transpose
                tolerance = 1e-12  # consider anything with covariance less than this to be uncorrelated
                for i, index_i in enumerate(subset_indices):
                    for j in range(i, len(subset_indices)):
                        index_j = subset_indices[j]
                        if abs(cov_data[index_i, index_j]) > tolerance:
                            label1 = cov_specs[index_i].label
                            label2 = cov_specs[index_j].label
                            covariance = cov_data[index_i, index_j]
                            self.thermo_covariances_dict[tuple(sorted([label1, label2]))] = covariance

    def load_thermo_covariances_from_groups(self):
        """
        This function populates the self.thermo_covariances_dict with covariance data (in units of (kcal/mol)^2) from the given covariance group trees

        For each group tree, it expects:
        1. a covariance.npy file containing the thermo covariance matrix, 
        2. a groups.py file containing the group definitions for the covariance data, in the same order as covariance.npy, and
        3. (optional) a species_dictionary.txt file containing the species in an associated library containing the training data.
           For example, the adsorptionPt111 correction tree uses the same species as surfaceThermoPt111, so it includes a species_dictionary.txt file to be able to get correlations with that library.

        See the RMG-database/scripts/compile_BEEF_cov.ipynb Jupyter notebook for more details on how to generate these covariance libraries.

        This function only adds covariance data for groups and species that are actually in the model, (which may still include a subset of library entries for species not in the model as in the case of the radical/HBI correction)
        and only for the thermo source associated with that group/library. The goal is to keep the dictionary as small as possible because the lookups scale badly.

        Note: the covariance.npy matrix is in units of (kJ/mol)^2, but gets converted to (kcal/mol)^2 in this function to match the rest of the analysis
        """
        from rmgpy.chemkin import load_species_dictionary
        # assumes there might also be covariances associated with library entries

        # associated library is hardcoded for now
        associated_libraries = {
            'adsorptionPt111': 'surfaceThermoPt111',
        }

        if self.database is None:
            raise RuntimeError('Must load database before loading covariance groups, since we need the database to find the associated libraries for the covariance groups')
        if self.thermo_covariance_groups is not None:
            self.compile_all_sources()
            for cov_group_tree in self.thermo_covariance_groups:
                associated_library = None
                cov_group_tree_name = os.path.basename(cov_group_tree)
                if cov_group_tree_name in self.database.thermo.groups:
                    grouptree = self.database.thermo.groups[cov_group_tree_name]
                    rmg_group_tree_file = os.path.join(rmgpy.settings['database.directory'], 'thermo', 'groups', f'{cov_group_tree_name}.py')
                else:
                    raise ValueError(f'Thermo covariance library {cov_group_tree_name} not found in the loaded database')
                if cov_group_tree_name in associated_libraries:
                    library_name = associated_libraries[cov_group_tree_name]
                    if library_name in self.database.thermo.libraries:
                        associated_library = self.database.thermo.libraries[library_name]
                    else:
                        raise ValueError(f'Associated library {library_name} for covariance group {cov_group_tree_name} not found in the loaded database')

                covariance_file = os.path.join(cov_group_tree, 'covariance.npy')
                if not os.path.isfile(covariance_file):
                    raise ValueError(f'Thermo covariance file {covariance_file} not found in {cov_group_tree}')
                group_database_file = os.path.join(cov_group_tree, 'groups.py')
                covariance_molecules = None
                if associated_library is not None:
                    covariance_molecules = os.path.join(cov_group_tree, 'species_dictionary.txt')
                    if not os.path.isfile(covariance_molecules):
                        raise ValueError(f'Thermo molecules file {covariance_molecules} not found in {cov_group_tree}')

                # warn the user if the covariance group tree is older than the associated library, since this likely means the covariance data is out of date
                covariance_data_time = os.path.getmtime(covariance_file)
                rmg_group_file_time = os.path.getmtime(rmg_group_tree_file)
                if covariance_data_time < rmg_group_file_time:
                    warnings.warn(f'Thermo covariance group tree {cov_group_tree} is older than the RMG group tree file {rmg_group_tree_file}, which may mean the covariance data is out of date and not consistent with the current group definitions in RMG')
                if associated_library is not None:
                    library_file_path = os.path.join(rmgpy.settings['database.directory'], 'thermo', 'libraries', f'{associated_library.label}.py')
                    if not os.path.isfile(library_file_path):
                        # just warn the user instead of raising an error, since the covariance data may still be relevant even if the library file is not available
                        warnings.warn(f'Could not find the file for associated library {associated_library.label} to compare modification times with the covariance group tree {cov_group_tree}. It may or may not be up to date.')
                    else:
                        library_file_time = os.path.getmtime(library_file_path)
                        if covariance_data_time < library_file_time:
                            warnings.warn(f'Thermo covariance group tree {cov_group_tree} is older than the associated library {associated_library.label}, which may mean the covariance data is out of date and not consistent with the current library data')

                # load data
                cov_data = np.load(covariance_file) / 4.184 / 4.184  # convert from (kJ/mol)^2 to (kcal/mol)^2
                group_database = rmgpy.data.thermo.ThermoGroups()
                group_database.load(group_database_file)

                # reconstruct the groups and molecules stored in the molecules.pickle file
                cov_specs = []
                if covariance_molecules is not None:
                    cov_species_dict = load_species_dictionary(covariance_molecules)
                    cov_specs = [item for _, item in cov_species_dict.items()]

                group_items = [x.item for _, x in group_database.entries.items()]
                group_labels = [x.label for _, x in group_database.entries.items()]
                n_groups = len(group_items)
                n_mols = len(cov_specs)
                if cov_data.shape[0] != n_groups + n_mols:
                    raise ValueError(f'Covariance data size {cov_data.shape[0]} does not match the number of groups and molecules {n_mols} + {n_groups} in covariance group tree {cov_group_tree}')


                # add all groups in the model
                tolerance = 1e-12  # consider anything with covariance less than this to be uncorrelated

                # --------------------------------- add group-group correlations ---------------------------------
                # only consider groups actually in the model to keep this a reasonable size
                groups_in_model = []
                if cov_group_tree_name == 'adsorptionPt111' and 'ADS' in self.all_thermo_sources:
                    groups_in_model = self.all_thermo_sources['ADS'][cov_group_tree_name]
                else:
                    groups_in_model = self.all_thermo_sources['GAV'][cov_group_tree_name]
                group_items_in_model = [x.item for x in groups_in_model]
                group_labels_in_model = [x.label for x in groups_in_model]

                # get labels through isomorphism with the groups in the model, so it doesn't matter if names change, it matches structure
                valid_group_indices = [i for i in range(n_groups) if get_i_thing(group_items[i], group_items_in_model) >= 0]
                valid_group_labels = [group_labels_in_model[get_i_thing(group_items[i], group_items_in_model)] for i in valid_group_indices]

                groupname_prefix = f'Group({cov_group_tree_name})'
                if cov_group_tree_name == 'adsorptionPt111':
                    groupname_prefix = f'AdsorptionCorrection({cov_group_tree_name})'

                # only iterate through groups in the model
                # go through upper triangle of the covariance matrix, since the other half is just the transpose
                for i, index_i in enumerate(valid_group_indices):
                    label_i = f'{groupname_prefix} {valid_group_labels[i]}'
                    for j in range(i, len(valid_group_indices)):
                        index_j = valid_group_indices[j]
                        label_j = f'{groupname_prefix} {valid_group_labels[j]}'
                        if abs(cov_data[index_i, index_j]) > tolerance:
                            covariance = cov_data[index_i, index_j]
                            self.thermo_covariances_dict[tuple(sorted([label_i, label_j]))] = covariance
                
                # ---------------------------- add group-molecule correlations from associated library (if applicable) ----------------------------
                if associated_library is not None:
                    # Figure out the associated library labels, but only include species in the model
                    subset_species_indices = []  # keep track of indices relevant to the model
                    for i_lib, lib_species in enumerate(cov_specs):
                        i_sp = get_i_thing(lib_species, self.species_list)
                        if i_sp < 0:
                            continue

                        # make sure the species actually comes from this library, otherwise skip
                        result = self.database.thermo.get_thermo_data_from_library(lib_species, associated_library)
                        if result is not None:
                            surface_prefix = 'Surface_' if lib_species.contains_surface_site() else ''
                            # match the label as constructed in assign_intermediate_uncertainties,
                            # where the number corresponds to the index of the species in species_list
                            try:
                                label = f'{surface_prefix}Library {self.species_list[i_sp].to_chemkin()}'
                            except IndexError:
                                label = f'{surface_prefix}Library {self.species_list[i_sp].molecule[0].to_inchi()}'
                            lib_species.label = label
                            subset_species_indices.append(i_lib)

                    # fill in the dictionary of covariances from the covariance libraries
                    for i, index_i in enumerate(valid_group_indices):
                        label_i = f'{groupname_prefix} {valid_group_labels[i]}'
                        for j, index_j in enumerate(subset_species_indices):
                            label_j = cov_specs[index_j].label
                            if abs(cov_data[index_i, n_groups + index_j]) > tolerance:
                                covariance = cov_data[index_i, n_groups + index_j]
                                self.thermo_covariances_dict[tuple(sorted([label_i, label_j]))] = covariance


    def load_kinetic_covariances_from_rules(self):
        """
        This function populates the self.kinetic_covariances_dict with covariance data (in units of (ln k)^2) from the given covariance rate rule trees

        For each kinetics tree, it expects:
        1. a covariance.npy file containing the kinetic covariance matrix
        2. a reactions.py file containing the training reactions used to train the tree
        3. a species_dictionary.txt file containing the species structures associated with the training reactions
        (4). it expects the kinetics tree in the database to have the same nodes as the covariance.npy, but not necessarily
             the same training reactions, because people add training reactions a lot more frequently than they retrain the trees


        # TODO, get the database PR together so that this notebook exists
        See the RMG-database/scripts/compile_BEEF_cov.ipynb Jupyter notebook for more details on how to generate these covariance libraries.

        This function only adds covariance data for rule nodes and training reactions that are actually in the model,
        The goal is to keep the dictionary as small as possible because the lookups scale badly.

        TODO: figure out if this function needs SIDT in the name, is specific to SIDTs
        """

        if self.database is None:
            raise RuntimeError('Must load database before loading covariance groups, since we need the database to find the associated libraries for the covariance groups')
        if self.kinetic_covariance_rules is not None:
            self.compile_all_sources()

            for correlated_kinetics_tree_path in self.kinetic_covariance_rules:
                if (correlated_kinetics_tree_path.endswith(os.sep)):  # help out the poor user who puts a / at the end of the path
                    correlated_kinetics_tree_path = correlated_kinetics_tree_path[:-1]
                family_name = os.path.basename(correlated_kinetics_tree_path)

                if family_name in self.database.kinetics.families:
                    kineticstree = self.database.kinetics.families[family_name]
                    rmg_kinetics_tree_file = os.path.join(rmgpy.settings['database.directory'], 'kinetics', 'families', family_name, f'groups.py')
                else:
                    raise ValueError(f'Kinetic covariance family {family_name} not found in the loaded database')

                covariance_file = os.path.join(correlated_kinetics_tree_path, 'covariance.npy')
                if not os.path.isfile(covariance_file):
                    raise ValueError(f'Kinetic covariance file {covariance_file} not found')
                training_reactions_file = os.path.join(correlated_kinetics_tree_path, 'reactions.pkl')

                with open(training_reactions_file, 'rb') as f:
                    training_reactions = pickle.load(f)

                # warn the user if the covariance group tree is older than the associated library, since this likely means the covariance data is out of date
                covariance_data_time = os.path.getmtime(covariance_file)
                rmg_kinetics_tree_file_time = os.path.getmtime(rmg_kinetics_tree_file)
                if covariance_data_time < rmg_kinetics_tree_file_time:
                    warnings.warn(f'Kinetic covariance group tree {family_name} is older than the RMG kinetics tree file {rmg_kinetics_tree_file}, which may mean the covariance data is out of date and not consistent with the current group definitions in RMG')

                # load data
                cov_data = np.load(covariance_file)

                # load the ordered list of nodes in the tree
                tree_nodes = list(kineticstree.rules.entries.keys())

                assert len(tree_nodes) + len(training_reactions) == cov_data.shape[0], f'Covariance data size {cov_data.shape[0]} does not match the number of nodes {len(tree_nodes)} in covariance kinetics tree {family_name}'

                # ---------------------------- add rule-rule correlations ---------------------------------
                # only consider rules actually in the model to keep this a reasonable size
                if 'Rate Rules' not in self.all_kinetic_sources or family_name not in self.all_kinetic_sources['Rate Rules']:
                    continue
                rules_in_model = self.all_kinetic_sources['Rate Rules'][family_name]
                rule_labels_in_model = [x.label for x in rules_in_model]

                # get labels through isomorphism with the rules in the model, so it doesn't matter if names change, it matches structure
                valid_rule_indices = [i for i in range(len(tree_nodes)) if tree_nodes[i] in rule_labels_in_model]

                # go through upper triangle of the covariance matrix, since the other half is just the transpose
                tolerance = 1e-12  # consider anything with covariance less than this to be uncorrelated
                for index_i in valid_rule_indices:
                    label_i = f'Rate Rule {family_name} {tree_nodes[index_i]}'
                    for index_j in (valid_rule_indices):
                        label_j = f'Rate Rule {family_name} {tree_nodes[index_j]}'
                        if abs(cov_data[index_i, index_j]) > tolerance:
                            covariance = cov_data[index_i, index_j]
                            self.kinetic_covariances_dict[tuple(sorted([label_i, label_j]))] = covariance


    def extract_sources_from_model(self):
        """
        Extract the source data from the model using its comments.
        Must be done after loading model and database to work.
        """
        self.species_sources_dict = {}
        allowed_source_keys = {'Library', 'QM', 'GAV', 'ADS'}
        for species in self.species_list:
            source = self.database.thermo.extract_source_from_comments(species)
            unexpected_source_keys = set(source.keys()) - allowed_source_keys
            if unexpected_source_keys:
                raise ValueError(
                    f'Source of thermo must be either Library, QM, GAV, or ADS; '
                    f'got unexpected source keys {unexpected_source_keys} for species {species.label}'
                )

            # Prep the source data:
            # extend the library/QM sources to include species labels specific to the model
            # Specify the source as a Surface Library (if it has surface sites and comes from a library), for better differentiation when assigning uncertainties
            if 'Library' in source:
                lib_species = rmgpy.species.Species(molecule=[source["Library"][1].item])
                label = self.get_species_label(lib_species)
                if lib_species.contains_surface_site():
                    extended_source = source['Library'] + (f'Surface_Library {label}',)
                    source['Surface_Library'] = extended_source
                    source.pop('Library')
                else:
                    extended_source = source['Library'] + (f'Library {label}',)
                    source['Library'] = extended_source

            if 'QM' in source:
                qm_species = source["QM"][1]
                label = self.get_species_label(qm_species)
                extended_source = source['QM'] + (f'QM {label}',)
                source['QM'] = extended_source

            self.species_sources_dict[species] = source

        self.reaction_sources_dict = {}
        for reaction in self.reaction_list:
            source = self.database.kinetics.extract_source_from_comments(reaction)
            # Prep the source data
            # Consider any library or PDep reaction to be an independent parameter for now
            if 'Library' in source:
                # add a label to the end of the tuple with the official name to be used to match covariances
                if reaction.is_surface_reaction():
                    extended_source = source['Library'] + (f'Surface_Library {reaction.to_chemkin(self.species_list, kinetics=False)}',)
                    source['Surface_Library'] = extended_source
                    source.pop('Library')
                else:
                    extended_source = source['Library'] + (f'Library {reaction.to_chemkin(self.species_list, kinetics=False)}',)
                    source['Library'] = extended_source
            elif 'PDep' in source:
                source['PDep'] = self.reaction_list.index(reaction)
            elif 'Training' in source:
                # Do nothing here because training source already saves the entry from the training reaction
                pass
            elif 'Rate Rules' in source:
                pass
            else:
                raise Exception('Source of kinetics must be either Library, PDep, Training, or Rate Rules')
            self.reaction_sources_dict[reaction] = source

        # -------------------- load covariance libraries ------------------------#
        self.load_thermo_covariances_from_libraries()
        self.load_thermo_covariances_from_groups()
        self.load_kinetic_covariances_from_rules()

    def compile_all_sources(self):
        """
        Compile two dictionaries composed of all the thermo and kinetic sources.  Must
        be performed after extract_sources_from_model function
        """
        # Account for all the thermo sources
        all_thermo_sources = {'GAV': {}, 'Library': set(), 'QM': set(), 'ADS': {}, 'Surface_Library': set()}
        for source in self.species_sources_dict.values():
            if 'GAV' in source:
                for groupType in source['GAV'].keys():
                    group_entries = [groupTuple[0] for groupTuple in source['GAV'][groupType]]
                    if groupType not in all_thermo_sources['GAV']:
                        all_thermo_sources['GAV'][groupType] = set(group_entries)
                    else:
                        all_thermo_sources['GAV'][groupType].update(group_entries)
            if 'Library' in source:
                all_thermo_sources['Library'].add(source['Library'])
            if 'QM' in source:
                all_thermo_sources['QM'].add(source['QM'])
            if 'ADS' in source:
                for ads_group in source['ADS'].keys():
                    ads_group_entries = [groupTuple[0] for groupTuple in source['ADS'][ads_group]]
                    if ads_group not in all_thermo_sources['ADS']:
                        all_thermo_sources['ADS'][ads_group] = set(ads_group_entries)
                    else:
                        all_thermo_sources['ADS'][ads_group].update(ads_group_entries)
            if 'Surface_Library' in source:
                all_thermo_sources['Surface_Library'].add(source['Surface_Library'])

                # Convert to lists
        self.all_thermo_sources = {}
        self.all_thermo_sources['Library'] = list(all_thermo_sources['Library'])
        self.all_thermo_sources['QM'] = list(all_thermo_sources['QM'])
        self.all_thermo_sources['GAV'] = {}
        for groupType in all_thermo_sources['GAV'].keys():
            self.all_thermo_sources['GAV'][groupType] = list(all_thermo_sources['GAV'][groupType])
        self.all_thermo_sources['ADS'] = {}
        for ads_group in all_thermo_sources['ADS'].keys():
            self.all_thermo_sources['ADS'][ads_group] = list(all_thermo_sources['ADS'][ads_group])
        self.all_thermo_sources['Surface_Library'] = list(all_thermo_sources['Surface_Library'])

        # Account for all the kinetics sources
        all_kinetic_sources = {'Rate Rules': {}, 'Training': {}, 'Library': [], 'PDep': [], 'Surface_Library': []}
        for source in self.reaction_sources_dict.values():
            if 'Training' in source:
                family_label = source['Training'][0]
                training_entry = source['Training'][1]
                if family_label not in all_kinetic_sources['Training']:
                    all_kinetic_sources['Training'][family_label] = set([training_entry])
                else:
                    all_kinetic_sources['Training'][family_label].add(training_entry)
            elif 'Library' in source:
                all_kinetic_sources['Library'].append(source['Library'])
            elif 'Surface_Library' in source:
                all_kinetic_sources['Surface_Library'].append(source['Surface_Library'])
            elif 'PDep' in source:
                all_kinetic_sources['PDep'].append(source['PDep'])
            elif 'Rate Rules' in source:
                family_label = source['Rate Rules'][0]
                source_dict = source['Rate Rules'][1]
                rules = source_dict['rules']
                training = source_dict['training']
                if rules:
                    rule_entries = [ruleTuple[0] for ruleTuple in rules]
                    if family_label not in all_kinetic_sources['Rate Rules']:
                        all_kinetic_sources['Rate Rules'][family_label] = set(rule_entries)
                    else:
                        all_kinetic_sources['Rate Rules'][family_label].update(rule_entries)
                if training:
                    # Even though they are from training reactions, we consider the rate rules derived from the training
                    # reactions to be noncorrelated, due to the fact that some may be reversed.
                    training_rules = [trainingTuple[0] for trainingTuple in training]  # Pick the rate rule entries
                    if family_label not in all_kinetic_sources['Rate Rules']:
                        all_kinetic_sources['Rate Rules'][family_label] = set(training_rules)
                    else:
                        all_kinetic_sources['Rate Rules'][family_label].update(training_rules)

        self.all_kinetic_sources = {}
        self.all_kinetic_sources['Library'] = all_kinetic_sources['Library']
        self.all_kinetic_sources['Surface_Library'] = all_kinetic_sources['Surface_Library']
        self.all_kinetic_sources['PDep'] = all_kinetic_sources['PDep']
        # Convert to lists
        self.all_kinetic_sources['Rate Rules'] = {}
        for family_label in all_kinetic_sources['Rate Rules'].keys():
            self.all_kinetic_sources['Rate Rules'][family_label] = list(all_kinetic_sources['Rate Rules'][family_label])

        self.all_kinetic_sources['Training'] = {}
        for family_label in all_kinetic_sources['Training'].keys():
            self.all_kinetic_sources['Training'][family_label] = list(all_kinetic_sources['Training'][family_label])

    def get_species_label(self, species):
        i_sp = get_i_thing(species, self.species_list)
        if i_sp >= 0:
            return self.species_list[i_sp].to_chemkin()
        else:
            return species.to_chemkin()

    def assign_parameter_uncertainties(self, g_param_engine=None, k_param_engine=None, correlated=False, use_rank=False):
        """
        Assign uncertainties based on the sources of the species thermo and reaction kinetics.

        use_rank: whether to apply the rank-based uncertainties for library entries
        """
        if g_param_engine is None:
            g_param_engine = ThermoParameterUncertainty(other_covariances=self.thermo_covariances_dict, rank_dictionary=use_rank)
        if k_param_engine is None:
            k_param_engine = KineticParameterUncertainty(other_covariances=self.kinetic_covariances_dict, rank_dictionary=use_rank)
        self.used_rank = use_rank  # keep track of whether we applied rank-based uncertainties in case it's needed for later analysis or plotting

        # reset variables
        self.thermo_input_uncertainties = []
        self.kinetic_input_uncertainties = []
        self.thermo_intermediate_uncertainties = []
        self.kinetic_intermediate_uncertainties = []
        self.thermo_intermediates = None
        self.kinetic_intermediates = None

        self.dG_dqs = []  # a list of dictionaries to store the intermediate derivatives dG_i/dq for each parameter q that contributes to the uncertainty of species i's G
        self.dlnk_dqs = []  # a list of dictionaries to store the intermediate derivatives dlnk_i/dq for each parameter q that contributes to the uncertainty of reaction i's k

        for species in self.species_list:
            if not correlated:
                dG = g_param_engine.get_uncertainty_value(self.species_sources_dict[species])
                self.thermo_input_uncertainties.append(dG)
            else:
                source = self.species_sources_dict[species]
                intermediate_source = {}
                dG = {}
                dG_dq = {}
                if 'Library' in source:
                    pdG = g_param_engine.get_partial_uncertainty_value(source, 'Library', corr_param=source['Library'])
                    label = source['Library'][2]  # the label we added to the end of the library source tuple in extract_sources_from_model
                    dG[label] = pdG
                    dG_dq[label] = 1  # derivative of G with respect to the library parameter is 1, since it's a direct contribution
                    intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('Library', source['Library'])
                if 'Surface_Library' in source:
                    pdG = g_param_engine.get_partial_uncertainty_value(source, 'Surface_Library', corr_param=source['Surface_Library'])
                    label = source['Surface_Library'][2]
                    dG[label] = pdG
                    dG_dq[label] = 1  # derivative of G with respect to the surface library parameter is 1, since it's a direct contribution
                    intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('Surface_Library', source['Surface_Library'])
                if 'QM' in source:
                    pdG = g_param_engine.get_partial_uncertainty_value(source, 'QM', corr_param=source['QM'])
                    label = source['QM'][2]
                    dG[label] = pdG
                    dG_dq[label] = 1  # derivative of G with respect to the QM parameter is 1, since it's a direct contribution
                    intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('QM', source['QM'])
                if 'ADS' in source:
                    for adsGroupType, groupList in source['ADS'].items():
                        for group, weight in groupList:
                            pdG = g_param_engine.get_partial_uncertainty_value(source, 'ADS', group, adsGroupType)
                            label = 'AdsorptionCorrection({}) {}'.format(adsGroupType, group.label)
                            dG[label] = pdG
                            if weight != 1:
                                raise ValueError('Weight for adsorption group contribution to thermo should be 1, but got weight={weight} for {adsGroupType} in species {species}'.format(weight=weight, adsGroupType=adsGroupType, species=species))
                            dG_dq[label] = weight  # This should be 1 because there's only one group contribution per adsorption correction
                            intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('ADS', source['ADS'])
                if 'GAV' in source:
                    for groupType, groupList in source['GAV'].items():
                        for group, weight in groupList:
                            pdG = g_param_engine.get_partial_uncertainty_value(source, 'GAV', group, groupType)
                            label = 'Group({}) {}'.format(groupType, group.label)
                            dG[label] = pdG
                            dG_dq[label] = weight  # the derivative of G with respect to the group contribution is the weight of that group contribution to the overall G
                            intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('GAV', source['GAV'])
                    # We also know if there is group additivity used, there will be uncorrelated estimation error
                    est_pdG = g_param_engine.get_partial_uncertainty_value(source, 'Estimation')
                    if est_pdG:
                        label = 'Estimation {}'.format(species.to_chemkin())
                        dG[label] = est_pdG
                        dG_dq[label] = 1  # the derivative of G with respect to the estimation error is 1, since we add the term directly
                        intermediate_source[label] = g_param_engine.get_intermediate_uncertainty_value('Estimation', source['GAV'])
                self.thermo_input_uncertainties.append(dG)
                self.dG_dqs.append(dG_dq)
                self.thermo_intermediate_uncertainties.append(intermediate_source)

        for reaction in self.reaction_list:
            if not correlated:
                dlnk = k_param_engine.get_uncertainty_value(self.reaction_sources_dict[reaction])
                self.kinetic_input_uncertainties.append(dlnk)
            else:
                source = self.reaction_sources_dict[reaction]
                dlnk = {}
                dlnk_dq = {}
                intermediate_source = {}
                if 'Rate Rules' in source:
                    family = source['Rate Rules'][0]
                    source_dict = source['Rate Rules'][1]
                    rules = source_dict['rules']
                    training = source_dict['training']
                    surface_prefix = ''
                    if reaction.is_surface_reaction():
                        surface_prefix = 'Surface '

                    for ruleEntry, weight in rules:
                        dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Rate Rules', corr_param=ruleEntry,
                                                                             corr_family=family)
                        label = '{}Rate Rule {} {}'.format(surface_prefix, family, ruleEntry)
                        dlnk[label] = dplnk
                        dlnk_dq[label] = weight  # the derivative of ln(k) with respect to the rate rule contribution
                        intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Rate Rules', source['Rate Rules'])

                        # TODO this is a bug: adding this should not be changing things, but it is and I don't know why
                        # should probably overwrite with covariance matrix value if exists
                        # just to make kinetic_intermediate_uncertainties more consistent
                        # but it won't change the computation, which calls the _get_covariance_qq function later
                        # if self.kinetic_covariances_dict and (label, label) in self.kinetic_covariances_dict:
                        #     intermediate_source[label] = np.sqrt(self.kinetic_covariances_dict[(label, label)])
                    for ruleEntry, trainingEntry, weight in training:
                        dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Rate Rules', corr_param=ruleEntry,
                                                                             corr_family=family)
                        label = '{}Rate Rule {} {}'.format(surface_prefix, family, ruleEntry)
                        dlnk[label] = dplnk
                        dlnk_dq[label] = weight  # the derivative of ln(k) with respect to the training rule contribution
                        intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Rate Rules', source['Rate Rules'])

                    N = len(source_dict['rules']) + len(source_dict['training'])

                    # There is also estimation error if rate rules are used (nonexact and family contribute to this)
                    nonexact_dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Estimation Nonexact', corr_family=family)
                    if nonexact_dplnk:
                        label = 'Estimation Nonexact {}'.format(reaction.to_chemkin(self.species_list, kinetics=False))
                        dlnk[label] = nonexact_dplnk
                        dlnk_dq[label] = np.log10(N + 1)  # the derivative of ln(k) with respect to the nonexact estimation error
                        intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Estimation Nonexact', source['Rate Rules'])

                    family_dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Estimation Family', corr_family=family)
                    if family_dplnk:
                        label = 'Estimation Family {}'.format(reaction.to_chemkin(self.species_list, kinetics=False))
                        dlnk[label] = family_dplnk
                        dlnk_dq[label] = 1  # the derivative of ln(k) with respect to the family estimation error
                        intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Estimation Family', source['Rate Rules'])
                elif 'PDep' in source:
                    dplnk = k_param_engine.get_partial_uncertainty_value(source, 'PDep', source['PDep'])
                    label = 'PDep {}'.format(reaction.to_chemkin(self.species_list, kinetics=False))
                    dlnk[label] = dplnk
                    dlnk_dq[label] = 1
                    intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('PDep', source['PDep'])

                elif 'Library' in source:
                    dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Library', source['Library'])
                    label = source['Library'][2]
                    dlnk[label] = dplnk
                    dlnk_dq[label] = 1
                    intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Library', source['Library'])

                elif 'Surface_Library' in source:
                    dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Surface_Library', source['Surface_Library'])
                    label = source['Surface_Library'][2]
                    dlnk[label] = dplnk
                    dlnk_dq[label] = 1
                    intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Surface_Library', source['Surface_Library'])

                elif 'Training' in source:
                    dplnk = k_param_engine.get_partial_uncertainty_value(source, 'Training', source['Training'])
                    family = source['Training'][0]
                    label = 'Training {} {}'.format(family, reaction.to_chemkin(self.species_list, kinetics=False))
                    dlnk[label] = dplnk
                    dlnk_dq[label] = 1
                    intermediate_source[label] = k_param_engine.get_intermediate_uncertainty_value('Training', source['Training'])

                self.kinetic_input_uncertainties.append(dlnk)
                self.dlnk_dqs.append(dlnk_dq)
                self.kinetic_intermediate_uncertainties.append(intermediate_source)

    def sensitivity_analysis(self, initial_mole_fractions, sensitive_species, T, P, termination_time,
                             sensitivity_threshold=1e-3, number=10, fileformat='.png', initial_surface_coverages=None,
                             surface_volume_ratio=None, surface_site_density=2.72e-5):
        """
        Run sensitivity analysis using the RMG solver in a single ReactionSystem object

        initial_mole_fractions is a dictionary with Species objects as keys and mole fraction initial conditions
        sensitive_species is a list of sensitive Species objects
        number is the number of top species thermo or reaction kinetics desired to be plotted
        """

        from rmgpy.solver import SimpleReactor, SurfaceReactor, TerminationTime
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
        from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
        T = Quantity(T)
        P = Quantity(P)
        termination = [TerminationTime(Quantity(termination_time))]

        surface_mech = any([x.contains_surface_site() for x in self.species_list])
        if surface_mech:
            assert surface_volume_ratio is not None, 'Must provide surface_volume_ratio for sensitivity analysis of surface mechanisms'
            surface_volume_ratio = Quantity(surface_volume_ratio)

        if not surface_mech:
            reaction_system = SimpleReactor(T=T,
                                            P=P,
                                            initial_mole_fractions=initial_mole_fractions,
                                            termination=termination,
                                            sensitive_species=sensitive_species,
                                            sensitivity_threshold=sensitivity_threshold)
        else:
            reaction_system = SurfaceReactor(T=T,
                                             P_initial=P,
                                             initial_gas_mole_fractions=initial_mole_fractions,
                                             initial_surface_coverages=initial_surface_coverages,
                                             surface_volume_ratio=surface_volume_ratio,
                                             surface_site_density=surface_site_density,
                                             n_sims=1,
                                             termination=termination,
                                             sensitive_species=sensitive_species,
                                             sensitivity_threshold=sensitivity_threshold)

        # Create the csv worksheets for logging sensitivity
        util.make_output_subdirectory(self.output_directory, 'solver')
        sens_worksheet = []
        reaction_system_index = 0
        for spec in reaction_system.sensitive_species:
            csvfile_path = os.path.join(self.output_directory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1, spec.index))
            sens_worksheet.append(csvfile_path)

        reaction_system.attach(SimulationProfileWriter(
            self.output_directory, reaction_system_index, self.species_list))
        reaction_system.attach(SimulationProfilePlotter(
            self.output_directory, reaction_system_index, self.species_list))

        simulator_settings = SimulatorSettings()  # defaults

        model_settings = ModelSettings()  # defaults
        model_settings.tol_move_to_core = 0.1
        model_settings.tol_interrupt_simulation = 1.0
        model_settings.tol_keep_in_edge = 0.0

        reaction_system.simulate(
            core_species=self.species_list,
            core_reactions=self.reaction_list,
            edge_species=[],
            edge_reactions=[],
            surface_species=[],
            surface_reactions=[],
            model_settings=model_settings,
            simulator_settings=simulator_settings,
            sensitivity=True,
            sens_worksheet=sens_worksheet,
        )

        plot_sensitivity(self.output_directory, reaction_system_index, reaction_system.sensitive_species,
                         number=number, fileformat=fileformat)

    def local_analysis(self, sensitive_species, reaction_system_index=0, correlated=False, number=10,
                       fileformat='.png', t=None):
        """
        Conduct local uncertainty analysis on the reaction model.
        sensitive_species is a list of sensitive Species objects
        number is the number of highest contributing uncertain parameters desired to be plotted
        fileformat can be either .png, .pdf, or .svg
        t is the time at which to analyze the sensitivity (if None, uses the final time)
        """
        output = {}
        for sens_species in sensitive_species:
            csvfile_path = os.path.join(self.output_directory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1,
                                                                             sens_species.index))
            time, data_list = parse_csv_data(csvfile_path)

            if t is not None:  # truncate to the time of interest if specified, otherwise use the full time range
                t_index = int(np.argmin(np.abs(time.data - t)))
                if t_index < len(time.data) - 1:
                    time.data = time.data[:t_index + 1]
                    for data in data_list:
                        data.data = data.data[:t_index + 1]

            # Assign uncertainties
            thermo_data_list = []
            reaction_data_list = []
            for data in data_list:
                if data.species:
                    for species in self.species_list:
                        if species.to_chemkin() == data.species:
                            index = self.species_list.index(species)
                            break
                    else:
                        raise Exception('Chemkin name {} of species in the CSV file does not match anything in the '
                                        'species list.'.format(data.species))

                    data.uncertainty = self.thermo_input_uncertainties[index]
                    thermo_data_list.append(data)

                if data.reaction:
                    rxn_index = int(data.index) - 1
                    data.uncertainty = self.kinetic_input_uncertainties[rxn_index]
                    reaction_data_list.append(data)

            if correlated:
                correlated_thermo_data = {}
                correlated_reaction_data = {}
                for data in thermo_data_list:
                    for label, dpG in data.uncertainty.items():
                        if label in correlated_thermo_data:
                            # Unpack the labels and partial uncertainties
                            correlated_thermo_data[label].data[-1] += data.data[-1] * dpG  # Multiply the sensitivity with the partial uncertainty
                        else:
                            correlated_thermo_data[label] = GenericData(data=[data.data[-1] * dpG],
                                                                        uncertainty=1, label=label, species='dummy')
                for data in reaction_data_list:
                    for label, dplnk in data.uncertainty.items():
                        if label in correlated_reaction_data:
                            correlated_reaction_data[label].data[-1] += data.data[-1] * dplnk
                        else:
                            correlated_reaction_data[label] = GenericData(data=[data.data[-1] * dplnk],
                                                                          uncertainty=1, label=label, reaction='dummy')

                thermo_data_list = list(correlated_thermo_data.values())
                reaction_data_list = list(correlated_reaction_data.values())

            # Compute total variance
            total_variance = 0.0
            for data in thermo_data_list:
                total_variance += (data.data[-1] * data.uncertainty) ** 2
            for data in reaction_data_list:
                total_variance += (data.data[-1] * data.uncertainty) ** 2

            if not correlated:
                # Add the reaction index to the data label of the reaction uncertainties
                # data.index stores the physical index of the reaction + 1, so we convert it to the RMG index here
                for data in reaction_data_list:
                    data.label = 'k' + str(self.reaction_list[data.index - 1].index) + ': ' + data.label.split()[-1]

            if correlated:
                folder = os.path.join(self.output_directory, 'correlated')
            else:
                folder = os.path.join(self.output_directory, 'uncorrelated')
            if not os.path.exists(folder):
                try:
                    os.makedirs(folder)
                except OSError as e:
                    raise OSError('Uncertainty output directory could not be created: {0!s}'.format(e))

            r_path = os.path.join(folder, 'kineticsLocalUncertainty_{0}'.format(sens_species.to_chemkin()) + fileformat)
            t_path = os.path.join(folder, 'thermoLocalUncertainty_{0}'.format(sens_species.to_chemkin()) + fileformat)
            reaction_uncertainty = ReactionSensitivityPlot(x_var=time, y_var=reaction_data_list, num_reactions=number).uncertainty_plot(total_variance, filename=r_path)
            thermo_uncertainty = ThermoSensitivityPlot(x_var=time, y_var=thermo_data_list, num_species=number).uncertainty_plot(total_variance, filename=t_path)

            output[sens_species] = (total_variance, reaction_uncertainty, thermo_uncertainty)

        return output

    def local_analysis_intermediate(self, sensitive_species, reaction_system_index=0, correlated=False, number=10,
                                    fileformat='.png', t=None):
        """
        local uncertainty analysis using new formulation where parameters might not be fully independent

        sensitive_species is a list of sensitive Species objects
        number is the number of highest contributing uncertain parameters desired to be plotted
        fileformat can be either .png, .pdf, or .svg

        t is the time in seconds at which to perform the analysis.
        The default (None) uses the final timestep in the sensitivity analysis

        returns a dictionary of tuples for every sensitive species with:
        - the total variance in concentration of that species
        - the kinetic contributions to variance of that species's concentration
        - the thermo contributions to variance of that species's concentration
        """

        output = {}
        for sens_species in sensitive_species:
            # 1. ------------------------- Get sensitivities --------------------------
            csvfile_path = os.path.join(self.output_directory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1,
                                                                             sens_species.index))
            time, data_list = parse_csv_data(csvfile_path)

            if t is None:
                t_index = -1
            else:
                t_index = int(np.argmin(np.abs(time.data - t)))

            # get the sensitivities and compile a list of the species and reaction indices actually used
            # record sensitivities for all time, so sensitivity array is #time steps x #species
            species_sensitivity_full_array = np.zeros(len(self.species_list))
            reaction_sensitivity_full_array = np.zeros(len(self.reaction_list))

            # keeping track of which species/reactions are sensitive gives us a speedup in computing covariance matrices later on
            species_used = []
            reactions_used = []
            for data in data_list:
                if data.species:
                    for species in self.species_list:
                        if species.to_chemkin() == data.species:
                            index = self.species_list.index(species)
                            break
                    else:
                        raise ValueError(f'Chemkin name {data.species} of species in the CSV file does not match anything in the species list.')
                    species_sensitivity_full_array[index] = data.data[t_index]
                    species_used.append(index)

                if data.reaction:
                    rxn_index = int(data.index) - 1
                    reaction_sensitivity_full_array[rxn_index] = data.data[t_index]
                    reactions_used.append(rxn_index)
            species_used = sorted(species_used)
            reactions_used = sorted(reactions_used)
            
            # shorten the sensitivity vectors to only the nonzero elements
            species_sensitivity = np.array([species_sensitivity_full_array[i] for i in species_used])
            reaction_sensitivity = np.array([reaction_sensitivity_full_array[i] for i in reactions_used])

            # 2. ------------------------- Compute intermediate covariance --------------------------
            # Now get the covariance matrix of the intermediate parameters (if uncorrelated, these are just covariance(G_i, G_j) or covariance(lnk_i, lnk_j)
            Sigma_qq_thermo = self._get_intermediate_thermo_covariance_matrix(subset_indices=species_used)
            Sigma_qq_kinetics = self._get_intermediate_kinetics_covariance_matrix(subset_indices=reactions_used)

            # 3. ------------------------- Compute partial derivatives dG/dq or dlnk/dq --------------------------
            # get the parital derivative of G or lnk with respect to intermediates (if uncorrelated, these are identity matrices)
            dG_dq = self._get_dG_dq_matrix(subset_indices=species_used)
            dlnkdq = self._get_dlnk_dq_matrix(subset_indices=reactions_used)
            N_q_thermo = dG_dq.shape[1]  # the number of thermo contributions
            N_q_kinetics = dlnkdq.shape[1]  # the number of kinetics contributions

            # 4. ------------------------- Multiply all matrices together to get total variance --------------------------
            # total variance =
            #   species_sensitivity * dG_dq * Sigma_qq_thermo * dG_dq' * species_sensitivity' +
            #   reaction_sensitivity * dlnk_dq * Sigma_qq_kinetics * dlnk_dq' * reaction_sensitivity' +

            # we split this into:
            # species_sensitivity' * dG_dq' and Sigma_qq_thermo * dG_dq * species_sensitivity
            # because this gives us the contributions of each intermediate parameter
            thermo_contributions = np.multiply(np.dot(species_sensitivity, dG_dq), np.dot(Sigma_qq_thermo, np.dot(dG_dq.T, species_sensitivity.T)).T)
            kinetic_contributions = np.multiply(np.dot(reaction_sensitivity, dlnkdq), np.dot(Sigma_qq_kinetics, np.dot(dlnkdq.T, reaction_sensitivity.T)).T) 

            missed_negative_thermo = 0
            for i in range(len(thermo_contributions)):
                if thermo_contributions[i] < 0:
                    print(f'Warning: negative contribution to variance from {self.thermo_intermediates[i]} of {thermo_contributions[i]}. Setting contribution to 0 for plotting purposes.')
                    missed_negative_thermo += thermo_contributions[i]
                    thermo_contributions[i] = 0
            if missed_negative_thermo < 0:
                print(f'Warning: total missed negative contribution to variance from thermo intermediates is {missed_negative_thermo}.')

            # missed_negative_kinetic = 0
            # for i in range(len(kinetic_contributions)):
            #     if kinetic_contributions[i] < 0:
            #         print(f'Warning: negative contribution to variance from {self.kinetic_intermediates[i]} of {kinetic_contributions[i]}. Setting contribution to 0 for plotting purposes.')
            #         missed_negative_kinetic += kinetic_contributions[i]
            #         kinetic_contributions[i] = 0
            # if missed_negative_kinetic < 0:
            #     print(f'Warning: total missed negative contribution to variance from kinetic intermediates is {missed_negative_kinetic}.')

            total_variance = np.sum(thermo_contributions) + np.sum(kinetic_contributions)

            # 5. ------------------------- Make plots --------------------------
            # define labels if they're not already listed in all_thermo/kinetic_intermediates
            if self.thermo_intermediates is None or len(self.thermo_intermediates) != N_q_thermo:
                self.thermo_intermediates = [f'dln[{sens_species.to_chemkin()}]/dG[{self.species_list[sp_idx].to_chemkin()}]' for sp_idx in species_used]
            if self.kinetic_intermediates is None or len(self.kinetic_intermediates) != N_q_kinetics:
                self.kinetic_intermediates = ['k' + str(self.reaction_list[rxn_idx].index) + ': ' + self.reaction_list[rxn_idx].to_chemkin(kinetics=False) for rxn_idx in reactions_used]

            # append all data points
            thermo_plotting_data = []
            kinetics_plotting_data = []
            for i in range(N_q_thermo):
                label = self.thermo_intermediates[i]
                thermo_plotting_data.append(GenericData(data=[np.sqrt(thermo_contributions[i])], uncertainty=1.0, label=label, species='dummy'))
            for i in range(N_q_kinetics):
                label = self.kinetic_intermediates[i]
                kinetics_plotting_data.append(GenericData(data=[np.sqrt(kinetic_contributions[i])], uncertainty=1.0, label=label, reaction='dummy'))

            # set up the folders and filenames for plotting
            folder = os.path.join(self.output_directory, 'uncorrelated')
            if correlated:
                folder = os.path.join(self.output_directory, 'correlated')
            os.makedirs(folder, exist_ok=True)

            r_path = os.path.join(folder, f'kineticsLocalUncertainty_{sens_species.to_chemkin()}{fileformat}')
            t_path = os.path.join(folder, f'thermoLocalUncertainty_{sens_species.to_chemkin()}{fileformat}')
            reaction_uncertainty = ReactionSensitivityPlot(x_var=time, y_var=kinetics_plotting_data, num_reactions=number).uncertainty_plot(total_variance, filename=r_path)
            thermo_uncertainty = ThermoSensitivityPlot(x_var=time, y_var=thermo_plotting_data, num_species=number).uncertainty_plot(total_variance, filename=t_path)

            output[sens_species] = (total_variance, reaction_uncertainty, thermo_uncertainty)

        return output

    def get_variance_across_x(self, sensitive_species, reaction_system_index=0, correlated=False):
        """
        Like local_analysis_intermediate, but a more computationally efficient way to compute variance across a range of time/temperature points
        by computing all values at once and not keeping track of individual contributions for plotting/ranking purposes

        sensitive_species is a single sensitive Species object

        returns an array of variances for every time/temperature point in the sensitivity analysis
        """

        if isinstance(sensitive_species, list):
            print('Warning: get_variance_across_x is designed to take in a single sensitive species, but got a list. Using only the first species in the list for analysis.')
            sensitive_species = sensitive_species[0]

        # 1. ------------------------- Get sensitivities --------------------------
        csvfile_path = os.path.join(self.output_directory, 'solver', f'sensitivity_{reaction_system_index + 1}_SPC_{sensitive_species.index}.csv')
        x, data_list = parse_csv_data(csvfile_path)  # x is usually time or temperature, depending on how the sensitivity analysis was set up

        # get the sensitivities and compile a list of the species and reaction indices actually used
        # record sensitivities for all time, so sensitivity array is #time steps x #species
        species_sensitivity_full_array = np.zeros((len(x.data), len(self.species_list)))
        reaction_sensitivity_full_array = np.zeros((len(x.data), len(self.reaction_list)))

        # keeping track of which species/reactions are sensitive gives us a speedup in computing covariance matrices later on
        species_used = []
        reactions_used = []
        for data in data_list:
            if data.species:
                for species in self.species_list:
                    if species.to_chemkin() == data.species:
                        index = self.species_list.index(species)
                        break
                else:
                    raise ValueError(f'Chemkin name {data.species} of species in the CSV file does not match anything in the species list.')
                species_sensitivity_full_array[:, index] = data.data
                species_used.append(index)

            if data.reaction:
                rxn_index = int(data.index) - 1
                reaction_sensitivity_full_array[:, rxn_index] = data.data
                reactions_used.append(rxn_index)
        species_used = sorted(species_used)
        reactions_used = sorted(reactions_used)
        
        # shorten the sensitivity vectors to only the nonzero elements
        species_sensitivity = np.zeros((len(x.data), len(species_used)))
        reaction_sensitivity = np.zeros((len(x.data), len(reactions_used)))
        for i in range(len(species_used)):
            species_sensitivity[:, i] = species_sensitivity_full_array[:, species_used[i]]
        for i in range(len(reactions_used)):
            reaction_sensitivity[:, i] = reaction_sensitivity_full_array[:, reactions_used[i]]

        # 2. ------------------------- Compute intermediate covariance --------------------------
        # Now get the covariance matrix of the intermediate parameters (if uncorrelated, these are just covariance(G_i, G_j) or covariance(lnk_i, lnk_j)
        Sigma_qq_thermo = self._get_intermediate_thermo_covariance_matrix(subset_indices=species_used)
        Sigma_qq_kinetics = self._get_intermediate_kinetics_covariance_matrix(subset_indices=reactions_used)

        # 3. ------------------------- Compute partial derivatives dG/dq or dlnk/dq --------------------------
        # get the parital derivative of G or lnk with respect to intermediates (if uncorrelated, these are identity matrices)
        dG_dq = self._get_dG_dq_matrix(subset_indices=species_used)
        dlnkdq = self._get_dlnk_dq_matrix(subset_indices=reactions_used)

        # 4. ------------------------- Multiply all matrices together to get total variance --------------------------
        # total variance =
        #   species_sensitivity * dG_dq * Sigma_qq_thermo * dG_dq' * species_sensitivity' +
        #   reaction_sensitivity * dlnk_dq * Sigma_qq_kinetics * dlnk_dq' * reaction_sensitivity' +
        total_variance = np.dot(species_sensitivity, np.dot(dG_dq, np.dot(Sigma_qq_thermo, np.dot(dG_dq.T, species_sensitivity.T)))) + \
                         np.dot(reaction_sensitivity, np.dot(dlnkdq, np.dot(Sigma_qq_kinetics, np.dot(dlnkdq.T, reaction_sensitivity.T))))
        total_variance = np.diag(total_variance)  # we only care about the diagonal elements, which give us the variance at each time/temperature point

        return total_variance


    def get_thermo_covariance_matrix(self, g_param_engine=None):
        """
        Return the thermo covariance matrix as a numpy array.
        NxN square matrix where N is the number of species in the model,
        with the covariance between species i and j in the ith row and jth column. 
        Units are in (kcal/mol)^2.
        Must call assign_parameter_uncertainties first to populate the source dictionaries.

        TODO speed this up with sparse matrix multiplication?
        """
        assert self.thermo_input_uncertainties is not None, 'Must call assign_parameter_uncertainties first'
        assert len(self.thermo_input_uncertainties) > 0, 'No thermodynamic parameters found'
        if isinstance(self.thermo_input_uncertainties[0], np.float64):  # uncorrelated case
            self.thermo_covariance_matrix = np.float_power(np.diag(self.thermo_input_uncertainties), 2.0)
            return self.thermo_covariance_matrix

        self.thermo_covariance_matrix = np.zeros((len(self.species_list), len(self.species_list)))

        if g_param_engine is None:
            g_param_engine = ThermoParameterUncertainty(other_covariances=self.thermo_covariances_dict, rank_dictionary=self.used_rank)
        
        for i in range(len(self.species_list)):
            for j in range((len(self.species_list))):
                for q in self.dG_dqs[i].keys():
                    dG_i_dq = self.dG_dqs[i][q]
                    q_uncertainty = self.thermo_intermediate_uncertainties[i][q]
                    for r in self.dG_dqs[j].keys():
                        dG_j_dr = self.dG_dqs[j][r]
                        self.thermo_covariance_matrix[i, j] += dG_i_dq * g_param_engine._get_covariance_qq(q, r, q_uncertainty) * dG_j_dr

        return self.thermo_covariance_matrix

    def get_kinetic_covariance_matrix(self, k_param_engine=None):
        """
        Return the kinetic covariance matrix as a numpy array.
        MxM square matrix where M is the number of reactions in the model,
        with the covariance between reaction i and j in the ith row and jth column.
        Units are in (ln(k))^2.
        Must call assign_parameter_uncertainties first to populate the source dictionaries.

        TODO speed this up with sparse matrix multiplication?
        """
        assert self.kinetic_input_uncertainties is not None, 'Must call assign_parameter_uncertainties first'
        assert len(self.kinetic_input_uncertainties) > 0, 'No kinetic parameters found'
        if isinstance(self.kinetic_input_uncertainties[0], np.float64):  # uncorrelated case
            self.kinetic_covariance_matrix = np.float_power(np.diag(self.kinetic_input_uncertainties), 2.0)
            return self.kinetic_covariance_matrix

        if k_param_engine is None:
            k_param_engine = KineticParameterUncertainty(other_covariances=self.kinetic_covariances_dict, rank_dictionary=self.used_rank)

        self.kinetic_covariance_matrix = np.zeros((len(self.reaction_list), len(self.reaction_list)))
        
        for i in range(len(self.reaction_list)):
            for j in range(len(self.reaction_list)):
                for q in self.dlnk_dqs[i].keys():
                    dlnk_i_dq = self.dlnk_dqs[i][q]
                    q_uncertainty = self.kinetic_intermediate_uncertainties[i][q]
                    for r in self.dlnk_dqs[j].keys():
                        dlnk_j_dr = self.dlnk_dqs[j][r]
                        self.kinetic_covariance_matrix[i, j] += dlnk_i_dq * k_param_engine._get_covariance_qq(q, r, q_uncertainty) * dlnk_j_dr

        return self.kinetic_covariance_matrix
    
    def _get_intermediate_thermo_covariance_matrix(self, g_param_engine=None, subset_indices=None):
        """
        Make an explicit covariance matrix of all the qs (intermediate thermo parameters, like specific groups or library entries)

        Requires calling assign_parameter_uncertainties first

        if subset_indices is None, computes the full matrix.
        Otherwise, only computes the matrices relevant to the species indicated by subset_indices
        """
        if subset_indices is None:
            subset_indices = np.arange(len(self.species_list))

        if isinstance(self.thermo_input_uncertainties[0], np.float64):  # uncorrelated case
            self.thermo_intermediates = None  # in the uncorrelated case, we don't have specific intermediate parameters to track, so we can set this to None
            self.Sigma_ww_thermo = np.diag(np.float_power([self.thermo_input_uncertainties[i] for i in subset_indices], 2.0))
            return self.Sigma_ww_thermo

        if g_param_engine is None:
            g_param_engine = ThermoParameterUncertainty(other_covariances=self.thermo_covariances_dict, rank_dictionary=self.used_rank)
        

        # the goal is to construct a list of unique tuples (q_label, q_partial_uncertainty)
        # but first, find the unique list of intermediate thermo sources used
        self.thermo_intermediates = {}
        for sp_idx in subset_indices:
            for q_label in self.thermo_intermediate_uncertainties[sp_idx].keys():
                self.thermo_intermediates[q_label] = self.thermo_intermediate_uncertainties[sp_idx][q_label]
        q_uncertainties = [uncertainty for q_label, uncertainty in self.thermo_intermediates.items()]
        self.thermo_intermediates = [q_label for q_label in self.thermo_intermediates.keys()]
        W = len(self.thermo_intermediates)

        self.Sigma_ww_thermo = np.zeros((W, W))
        for i in range(W):
            q_i_label = self.thermo_intermediates[i]
            q_i_uncertainty = q_uncertainties[i]
            for j in range(i + 1):
                q_j_label = self.thermo_intermediates[j]
                self.Sigma_ww_thermo[i, j] = g_param_engine._get_covariance_qq(q_i_label, q_j_label, q_i_uncertainty)
                self.Sigma_ww_thermo[j, i] = self.Sigma_ww_thermo[i, j]  # symmetric matrix
        return self.Sigma_ww_thermo

    def _get_intermediate_kinetics_covariance_matrix(self, k_param_engine=None, subset_indices=None):
        """
        Make an explicit covariance matrix of all the qs (intermediate kinetic parameters, like specific rate rules or libraries entries)

        Requires calling assign_parameter_uncertainties first

        if subset_indices is None, computes the full matrix.
        Otherwise, only computes the matrices relevant to the reactions indicated by subset_indices
        """
        if subset_indices is None:
            subset_indices = np.arange(len(self.reaction_list))

        if isinstance(self.kinetic_input_uncertainties[0], np.float64):  # uncorrelated case
            self.kinetic_intermediates = None  # in the uncorrelated case, we don't have specific intermediate parameters to track, so we can set this to None
            self.Sigma_ww_kinetics = np.diag(np.float_power([self.kinetic_input_uncertainties[i] for i in subset_indices], 2.0))
            return self.Sigma_ww_kinetics

        if k_param_engine is None:
            k_param_engine = KineticParameterUncertainty(other_covariances=self.kinetic_covariances_dict, rank_dictionary=self.used_rank)

        # the goal is to construct a list of unique tuples (q_label, q_partial_uncertainty)
        # but first, find the unique list of intermediate kinetic sources used
        self.kinetic_intermediates = {}
        for rxn_idx in subset_indices:
            for q in self.kinetic_intermediate_uncertainties[rxn_idx].keys():
                self.kinetic_intermediates[q] = self.kinetic_intermediate_uncertainties[rxn_idx][q]
        q_uncertainties = [uncertainty for q_label, uncertainty in self.kinetic_intermediates.items()]
        self.kinetic_intermediates = [q_label for q_label in self.kinetic_intermediates.keys()]

        W = len(self.kinetic_intermediates)

        self.Sigma_ww_kinetics = np.zeros((W, W))
        for i in range(W):
            q_i_label = self.kinetic_intermediates[i]
            q_i_uncertainty = q_uncertainties[i]
            for j in range(i + 1):
                q_j_label = self.kinetic_intermediates[j]
                self.Sigma_ww_kinetics[i, j] = k_param_engine._get_covariance_qq(q_i_label, q_j_label, q_i_uncertainty)
                self.Sigma_ww_kinetics[j, i] = self.Sigma_ww_kinetics[i, j]  # symmetric matrix
        return self.Sigma_ww_kinetics

    def _get_dG_dq_matrix(self, subset_indices=None):
        """
        Returns an nxW matrix of partial derivatives where
        n is number of species (or subset species) and W is number of relevant intermediate paramaters

        subset_indices is the set of species indices that matter
        if subset_indices is None, it computes the full matrix
        use at your own risk!

        assumes that get_intermediate_thermo_covariance_matrix was called with matching subset_indices
        if not, then matrix dimensions will probably not match up for later multiplication,
        which will signal the user that there's a problem
        """
        if subset_indices is None:
            subset_indices = np.arange(len(self.species_list))

        # return a square identity matrix if uncorrelated
        if isinstance(self.thermo_input_uncertainties[0], np.float64):  # uncorrelated case
            return np.eye(len(subset_indices))

        dGdq = np.zeros((len(subset_indices), len(self.thermo_intermediates)))

        for i, sp_idx in enumerate(subset_indices):
            for key in self.dG_dqs[sp_idx].keys():
                q_index = self.thermo_intermediates.index(key)
                dGdq[i, q_index] = self.dG_dqs[sp_idx][key]

        return dGdq
    
    def _get_dlnk_dq_matrix(self, subset_indices=None):
        """
        Returns an mxW matrix of partial derivatives where
        m is number of reactions and W is number of intermediate paramaters

        assumes that get_intermediate_kinetic_covariance_matrix was called with matching subset_indices
        if not, then matrix dimensions will probably not match up for later multiplication,
        which will signal the user that there's a problem
        """
        if subset_indices is None:
            subset_indices = np.arange(len(self.reaction_list))

        # return a square identity matrix if uncorrelated
        if isinstance(self.kinetic_input_uncertainties[0], np.float64):  # uncorrelated case
            return np.eye(len(subset_indices))

        dlnkdq = np.zeros((len(subset_indices), len(self.kinetic_intermediates)))

        for i, rxn_idx in enumerate(subset_indices):
            for key in self.dlnk_dqs[rxn_idx].keys():
                q_index = self.kinetic_intermediates.index(key)
                dlnkdq[i, q_index] = self.dlnk_dqs[rxn_idx][key]

        return dlnkdq


def process_local_results(results, sensitive_species, number=10):
    """
    Return a dictionary of processed results along with a formatted string
    given results from local uncertainty analysis.
    """
    processed_results = {}
    for spc in sensitive_species:
        total_var, reaction_u, thermo_u = results[spc]
        reaction_c = []
        for label, reaction, u in reaction_u:
            reaction_c.append((label, reaction, u / total_var * 100))
        reaction_c.sort(key=lambda x: abs(x[2]), reverse=True)

        thermo_c = []
        for label, species, u in thermo_u:
            thermo_c.append((label, species, u / total_var * 100))
        thermo_c.sort(key=lambda x: abs(x[2]), reverse=True)

        processed_results[spc] = (total_var, reaction_c, thermo_c)

    output = ''
    for spc in sensitive_species:
        output += '================================================================================\n'
        total_var, reaction_c, thermo_c = processed_results[spc]
        output += 'Total variance [(d ln(c))^2] for species {0} is {1:.6f}\n'.format(spc.label, total_var)
        output += '--------------------------------------------------------------------------------\n'
        output += 'Top {0:2} reaction rate contributors                              Sensitivity Index\n'.format(number)
        output += '--------------------------------------------------------------------------------\n'
        for label, reaction, c in reaction_c[:number]:
            output += '{0:<65}{1:>14.4f}%\n'.format(label, c)
        output += '--------------------------------------------------------------------------------\n'
        output += 'Top {0:2} thermochemistry contributors                            Sensitivity Index\n'.format(number)
        output += '--------------------------------------------------------------------------------\n'
        for label, species, c in thermo_c[:number]:
            output += '{0:<65}{1:>14.4f}%\n'.format(label, c)
        output += '================================================================================\n\n'

    return processed_results, output

def get_i_thing(thing, thing_list):
    # get index of a species/molecule/group/reaction in a list of those things,
    # where the thing might not be exactly the same object as the one in the list but is isomorphic to it
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1
