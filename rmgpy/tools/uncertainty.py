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

import os

import numpy as np

import rmgpy.util as util
from rmgpy.species import Species
from rmgpy.tools.data import GenericData
from rmgpy.tools.plot import parseCSVData, plot_sensitivity, ReactionSensitivityPlot, ThermoSensitivityPlot


class ThermoParameterUncertainty(object):
    """
    This class is an engine that generates the species uncertainty based on its thermo sources.
    """

    def __init__(self, dG_library=1.5, dG_QM=3.0, dG_GAV=1.5, dG_group=0.10):
        """
        Initialize the different uncertainties dG_library, dG_QM, dG_GAV, and dG_other with set values
        in units of kcal/mol.
        
        We expect a uniform distribution for some species free energy G in [Gmin, Gmax].
        dG = (Gmax-Gmin)/2
        """
        self.dG_library = dG_library
        self.dG_QM = dG_QM
        self.dG_GAV = dG_GAV
        self.dG_group = dG_group

    def getUncertaintyValue(self, source):
        """
        Retrieve the uncertainty value in kcal/mol when the source of the thermo of a species is given.
        """
        dG = 0.0
        if 'Library' in source:
            dG += self.dG_library
        if 'QM' in source:
            dG += self.dG_QM
        if 'GAV' in source:
            dG += self.dG_GAV  # Add a fixed uncertainty for the GAV method
            for group_type, group_entries in source['GAV'].items():
                group_weights = [groupTuple[-1] for groupTuple in group_entries]
                dG += np.sum([weight * self.dG_group for weight in group_weights])

        return dG

    def getPartialUncertaintyValue(self, source, corrSourceType, corrParam=None, corrGroupType=None):
        """
        Obtain the partial uncertainty dG/dG_corr*dG_corr, where dG_corr is the correlated parameter
        
        `corrParam` is the parameter identifier itself, which is a integer for QM and library parameters, or a string for group values
        `corrSourceType` is a string, being either 'Library', 'QM', 'GAV', or 'Estimation'
        `corrGroupType` is a string used only when the source type is 'GAV' and indicates grouptype
        """

        if corrSourceType == 'Library':
            if 'Library' in source:
                if source['Library'] == corrParam:
                    # Correlated parameter is a source of the overall parameter
                    return self.dG_library

        elif corrSourceType == 'QM':
            if 'QM' in source:
                if source['QM'] == corrParam:
                    # Correlated parameter is a source of the overall parameter
                    return self.dG_QM

        elif corrSourceType == 'GAV':
            if 'GAV' in source:
                if corrGroupType in source['GAV']:
                    groupList = source['GAV'][corrGroupType]
                    for group, weight in groupList:
                        if group == corrParam:
                            return weight * self.dG_group

        elif corrSourceType == 'Estimation':
            if 'GAV' in source:
                return self.dG_GAV
        else:
            raise Exception('Thermo correlated source must be GAV, QM, Library, or Estimation')

        # If we get here, it means the correlated parameter was not found
        return None

    def getUncertaintyFactor(self, source):
        """
        Retrieve the uncertainty factor f in kcal/mol when the source of the thermo of a species is given.
        
        This is equivalent to sqrt(3)*dG in a uniform uncertainty interval
        """
        dG = self.getUncertaintyValue(source)
        f = np.sqrt(3) * dG


class KineticParameterUncertainty(object):
    """
    This class is an engine that generates the reaction uncertainty based on its kinetic sources.
    """

    def __init__(self, dlnk_library=0.5, dlnk_training=0.5, dlnk_pdep=2.0, dlnk_family=1.0, dlnk_nonexact=3.5,
                 dlnk_rule=0.5):
        """
        Initialize the different uncertainties dlnk
        
        We expect a uniform distribution for some reaction kinetics  about ln(k0) in [ln(kmin), ln(kmax)].
        dlnk = (ln(kmax)-ln(kmin))/2
        """
        self.dlnk_library = dlnk_library
        self.dlnk_training = dlnk_training
        self.dlnk_pdep = dlnk_pdep
        self.dlnk_family = dlnk_family
        self.dlnk_nonexact = dlnk_nonexact
        self.dlnk_rule = dlnk_rule

    def getUncertaintyValue(self, source):
        """
        Retrieve the dlnk uncertainty when the source of the reaction kinetics are given
        """
        dlnk = 0.0
        if 'Library' in source:
            # Should be a single library reaction source
            dlnk += self.dlnk_library
        elif 'PDep' in source:
            # Should be a single pdep reaction source
            dlnk += self.dlnk_pdep
        elif 'Training' in source:
            # Should be a single training reaction
            # Although some training entries may be used in reverse,
            # We still consider the kinetics to be directly dependent 
            dlnk += self.dlnk_training
        elif 'Rate Rules' in source:
            family_label = source['Rate Rules'][0]
            source_dict = source['Rate Rules'][1]
            exact = source_dict['exact']
            rule_weights = [ruleTuple[-1] for ruleTuple in source_dict['rules']]
            training_weights = [trainingTuple[-1] for trainingTuple in source_dict['training']]

            dlnk += self.dlnk_family ** 2
            N = len(rule_weights) + len(training_weights)
            if not exact:
                # nonexactness contribution increases as N increases
                dlnk += np.log10(N + 1) * self.dlnk_nonexact

            # Add the contributions from rules
            dlnk += np.sum([weight * self.dlnk_rule for weight in rule_weights])
            # Add the contributions from training
            # Even though these source from training reactions, we actually
            # use the uncertainty for rate rules, since these are now approximations
            # of the original reaction.  We consider these to be independent of original the training
            # parameters because the rate rules may be reversing the training reactions,
            # which leads to more complicated dependence
            dlnk += np.sum([weight * self.dlnk_rule for weight in training_weights])

        return dlnk

    def getPartialUncertaintyValue(self, source, corrSourceType, corrParam=None, corrFamily=None):
        """
        Obtain the partial uncertainty dlnk/dlnk_corr*dlnk_corr, where dlnk_corr is the correlated parameter
        
        `corrParam` is the parameter identifier itself, which is the string identifier of the rate rule
        `corrSourceType` is a string, being either 'Rate Rules', 'Library', 'PDep', 'Training' or 'Estimation'
        `corrFamily` is a string used only when the source type is 'Rate Rules' and indicates the family
        """

        if corrSourceType == 'Rate Rules':
            if 'Rate Rules' in source:
                familyLabel = source['Rate Rules'][0]
                if corrFamily == familyLabel:
                    sourceDict = source['Rate Rules'][1]
                    rules = sourceDict['rules']
                    training = sourceDict['training']
                    if rules:
                        for ruleEntry, weight in rules:
                            if corrParam == ruleEntry:
                                return weight * self.dlnk_rule
                    if training:
                        for ruleEntry, trainingEntry, weight in training:
                            if corrParam == ruleEntry:
                                return weight * self.dlnk_rule

        # Writing it this way in the function is not the most efficient, but makes it easy to use, and
        # testing a few if statements is not too costly
        elif corrSourceType == 'Library':
            if 'Library' in source:
                if corrParam == source['Library']:
                    # Should be a single library reaction source
                    return self.dlnk_library
        elif corrSourceType == 'PDep':
            if 'PDep' in source:
                if corrParam == source['PDep']:
                    return self.dlnk_pdep
        elif corrSourceType == 'Training':
            if 'Training' in source:
                # Should be a unique single training reaction
                if corrParam == source['Training']:
                    return self.dlnk_training

        elif corrSourceType == 'Estimation':
            # Return all the uncorrelated uncertainty associated with using an estimation scheme

            if 'Rate Rules' in source:
                sourceDict = source['Rate Rules'][1]
                exact = sourceDict['exact']

                dlnk = self.dlnk_family  # Base uncorrelated uncertainty just from using rate rule estimation
                # Additional uncertainty from using non-exact rate rule
                N = len(sourceDict['rules']) + len(sourceDict['training'])
                if not exact:
                    # nonexactness contribution increases as N increases
                    dlnk += np.log10(N + 1) * self.dlnk_nonexact
                return dlnk
        else:
            raise Exception('Kinetics correlated source must be Rate Rules, Library, PDep, Training, or Estimation')

        # If we get here, it means that we did not find the correlated parameter in the source
        return None

    def getUncertaintyFactor(self, source):
        """
        Retrieve the uncertainty factor f when the source of the reaction kinetics are given.
        
        This is equivalent to sqrt(3)/ln(10) * dlnk  in a uniform uncertainty interval
        """
        dlnk = self.getUncertaintyValue(source)
        f = np.sqrt(3) / np.log(10) * dlnk


class Uncertainty(object):
    """
    This class contains functions associated with running uncertainty analyses
    for a single RMG-generated mechanism.
    """

    def __init__(self, speciesList=None, reactionList=None, outputDirectory=''):
        """
        `speciesList`: list of RMG species objects
        `reactionList`: list of RMG reaction objects
        `outputDirectoy`: directory path for saving output files from the analyses
        """
        self.database = None
        self.speciesList = speciesList
        self.reactionList = reactionList
        self.speciesSourcesDict = None
        self.reactionSourcesDict = None
        self.allThermoSources = None
        self.allKineticSources = None
        self.thermoInputUncertainties = None
        self.kineticInputUncertainties = None
        self.outputDirectory = outputDirectory if outputDirectory else os.getcwd()

        # For extra species needed for correlated analysis but not in model
        self.extraSpecies = []

        # Make output directory if it does not yet exist:
        if not os.path.exists(self.outputDirectory):
            try:
                os.makedirs(self.outputDirectory)
            except:
                raise Exception('Uncertainty output directory could not be created.')

    def loadDatabase(self, kineticsFamilies='all', kineticsDepositories=None, thermoLibraries=None, reactionLibraries=None):
        """
        This function loads a single copy of the RMGDatabase with full verbose averaging
        of the rate rule to trace kinetics sources.  
        
        By default, this function loads all the kinetics families, only the training kinetics depository,
        the primaryThermoLibrary, and no reaction libraries.  
        """
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy import settings

        if not kineticsDepositories:
            kineticsDepositories = ['training']
        if not thermoLibraries:
            thermoLibraries = ['primaryThermoLibrary']
        if not reactionLibraries:
            reactionLibraries = []

        self.database = RMGDatabase()
        self.database.load(
            settings['database.directory'],
            kineticsFamilies=kineticsFamilies,
            kineticsDepositories=kineticsDepositories,
            thermoLibraries=thermoLibraries,
            reactionLibraries=reactionLibraries,
        )

        # Prepare the database by loading training reactions but not averaging the rate rules
        for familyLabel, family in self.database.kinetics.families.items():
            family.add_rules_from_training(thermoDatabase=self.database.thermo)
            family.fill_rules_by_averaging_up(verbose=True)

    def loadModel(self, chemkinPath, dictionaryPath, transportPath=None):
        """
        Load a RMG-generated model into the Uncertainty class
        `chemkinPath`: path to the chem_annotated.inp CHEMKIN mechanism 
        `dictionaryPath`: path to the species_dictionary.txt file 
        `transportPath`: path to the tran.dat file (optional)

        Then create dictionaries stored in self.thermoGroups and self.rateRules
        containing information about the source of the thermodynamic and kinetic
        parameters
        """
        from rmgpy.chemkin import loadChemkinFile

        self.speciesList, self.reactionList = loadChemkinFile(chemkinPath,
                                                              dictionaryPath=dictionaryPath,
                                                              transportPath=transportPath)

    def retrieveSaturatedSpeciesFromList(self, species):
        """
        Given a radical `species`, this function retrieves the saturated species objects from a list of species objects
        and returns the saturated species object along with a boolean that indicates if the species is not part of the model
        (True->not in the model, False->in the model)
        """

        molecule = species.molecule[0]
        assert molecule.is_radical(), "Method only valid for radicals."
        saturated_struct = molecule.copy(deep=True)
        saturated_struct.saturate_radicals()
        for otherSpecies in self.speciesList:
            if otherSpecies.is_isomorphic(saturated_struct):
                return otherSpecies, False

        # couldn't find saturated species in the model, try libraries
        new_spc = Species(molecule=[saturated_struct])
        thermo = self.database.thermo.get_thermo_data_from_libraries(new_spc)

        if thermo is not None:
            new_spc.thermo = thermo
            self.speciesList.append(new_spc)
            return new_spc, True
        else:
            raise Exception('Could not retrieve saturated species form of {0} from the species list'.format(species))

    def extractSourcesFromModel(self):
        """
        Extract the source data from the model using its comments.
        Must be done after loading model and database to work.
        """
        self.speciesSourcesDict = {}
        for species in self.speciesList:
            if species not in self.extraSpecies:
                source = self.database.thermo.extract_source_from_comments(species)

                # Now prep the source data
                # Do not alter the GAV information, but reassign QM and Library sources to the species indices that they came from
                if len(source) == 1:
                    # The thermo came from a single source, so we know it comes from a value describing the exact species
                    if 'Library' in source:
                        # Use just the species index in self.speciesList, for better shorter printouts when debugging
                        source['Library'] = self.speciesList.index(species)
                    if 'QM' in source:
                        source['QM'] = self.speciesList.index(species)

                elif len(source) == 2:
                    # The thermo has two sources, which indicates it's an HBI correction on top of a library or QM value.
                    # We must retrieve the original saturated molecule's thermo instead of using the radical species as the source of thermo
                    saturated_species, ignore_spc = self.retrieveSaturatedSpeciesFromList(species)

                    if ignore_spc:  # this is saturated species that isn't in the actual model
                        self.extraSpecies.append(saturated_species)

                    if 'Library' in source:
                        source['Library'] = self.speciesList.index(saturated_species)
                    if 'QM' in source:
                        source['QM'] = self.speciesList.index(saturated_species)
                else:
                    raise Exception('Source of thermo should not use more than two sources out of QM, Library, or GAV.')

                self.speciesSourcesDict[species] = source

        self.reactionSourcesDict = {}
        for reaction in self.reactionList:
            source = self.database.kinetics.extract_source_from_comments(reaction)
            # Prep the source data 
            # Consider any library or PDep reaction to be an independent parameter for now
            # and assign the source to the index of the reaction within self.reactionList
            if 'Library' in source:
                source['Library'] = self.reactionList.index(reaction)
            elif 'PDep' in source:
                source['PDep'] = self.reactionList.index(reaction)
            elif 'Training' in source:
                # Do nothing here because training source already saves the entry from the training reaction
                pass
            elif 'Rate Rules' in source:
                # Do nothing
                pass
            else:
                raise Exception('Source of kinetics must be either Library, PDep, Training, or Rate Rules')
            self.reactionSourcesDict[reaction] = source

        for spc in self.extraSpecies:
            self.speciesList.remove(spc)

    def compileAllSources(self):
        """
        Compile two dictionaries composed of all the thermo and kinetic sources.  Must
        be performed after extractSourcesFromModel function
        """
        # Account for all the thermo sources
        all_thermo_sources = {'GAV': {}, 'Library': set(), 'QM': set()}
        for source in self.speciesSourcesDict.values():
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

                # Convert to lists
        self.allThermoSources = {}
        self.allThermoSources['Library'] = list(all_thermo_sources['Library'])
        self.allThermoSources['QM'] = list(all_thermo_sources['QM'])
        self.allThermoSources['GAV'] = {}
        for groupType in all_thermo_sources['GAV'].keys():
            self.allThermoSources['GAV'][groupType] = list(all_thermo_sources['GAV'][groupType])

        # Account for all the kinetics sources
        all_kinetic_sources = {'Rate Rules': {}, 'Training': {}, 'Library': [], 'PDep': []}
        for source in self.reactionSourcesDict.values():
            if 'Training' in source:
                family_label = source['Training'][0]
                training_entry = source['Training'][1]
                if family_label not in all_kinetic_sources['Training']:
                    all_kinetic_sources['Training'][family_label] = set([training_entry])
                else:
                    all_kinetic_sources['Training'][family_label].add(training_entry)
            elif 'Library' in source:
                all_kinetic_sources['Library'].append(source['Library'])
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

        self.allKineticSources = {}
        self.allKineticSources['Library'] = all_kinetic_sources['Library']
        self.allKineticSources['PDep'] = all_kinetic_sources['PDep']
        # Convert to lists
        self.allKineticSources['Rate Rules'] = {}
        for family_label in all_kinetic_sources['Rate Rules'].keys():
            self.allKineticSources['Rate Rules'][family_label] = list(all_kinetic_sources['Rate Rules'][family_label])

        self.allKineticSources['Training'] = {}
        for family_label in all_kinetic_sources['Training'].keys():
            self.allKineticSources['Training'][family_label] = list(all_kinetic_sources['Training'][family_label])

    def assignParameterUncertainties(self, gParamEngine=None, kParamEngine=None, correlated=False):
        """
        Assign uncertainties based on the sources of the species thermo and reaction kinetics.
        """
        if gParamEngine is None:
            gParamEngine = ThermoParameterUncertainty()
        if kParamEngine is None:
            kParamEngine = KineticParameterUncertainty()

        self.thermoInputUncertainties = []
        self.kineticInputUncertainties = []

        for species in self.speciesList:
            if not correlated:
                dG = gParamEngine.getUncertaintyValue(self.speciesSourcesDict[species])
                self.thermoInputUncertainties.append(dG)
            else:
                source = self.speciesSourcesDict[species]
                dG = {}
                if 'Library' in source:
                    pdG = gParamEngine.getPartialUncertaintyValue(source, 'Library', corrParam=source['Library'])
                    try:
                        label = 'Library {}'.format(self.speciesList[source['Library']].to_chemkin())
                    except IndexError:
                        label = 'Library {}'.format(self.extraSpecies[source['Library'] - len(self.speciesList)].to_chemkin())
                    dG[label] = pdG
                if 'QM' in source:
                    pdG = gParamEngine.getPartialUncertaintyValue(source, 'QM', corrParam=source['QM'])
                    label = 'QM {}'.format(self.speciesList[source['QM']].to_chemkin())
                    dG[label] = pdG
                if 'GAV' in source:
                    for groupType, groupList in source['GAV'].items():
                        for group, weight in groupList:
                            pdG = gParamEngine.getPartialUncertaintyValue(source, 'GAV', group, groupType)
                            label = 'Group({}) {}'.format(groupType, group.label)
                            dG[label] = pdG
                    # We also know if there is group additivity used, there will be uncorrelated estimation error
                    est_pdG = gParamEngine.getPartialUncertaintyValue(source, 'Estimation')
                    if est_pdG:
                        label = 'Estimation {}'.format(species.to_chemkin())
                        dG[label] = est_pdG
                self.thermoInputUncertainties.append(dG)

        for reaction in self.reactionList:
            if not correlated:
                dlnk = kParamEngine.getUncertaintyValue(self.reactionSourcesDict[reaction])
                self.kineticInputUncertainties.append(dlnk)
            else:
                source = self.reactionSourcesDict[reaction]
                dlnk = {}
                if 'Rate Rules' in source:
                    family = source['Rate Rules'][0]
                    source_dict = source['Rate Rules'][1]
                    rules = source_dict['rules']
                    training = source_dict['training']
                    for ruleEntry, weight in rules:
                        dplnk = kParamEngine.getPartialUncertaintyValue(source, 'Rate Rules', corrParam=ruleEntry,
                                                                        corrFamily=family)
                        label = '{} {}'.format(family, ruleEntry)
                        dlnk[label] = dplnk

                    for ruleEntry, trainingEntry, weight in training:
                        dplnk = kParamEngine.getPartialUncertaintyValue(source, 'Rate Rules', corrParam=ruleEntry,
                                                                        corrFamily=family)
                        label = '{} {}'.format(family, ruleEntry)
                        dlnk[label] = dplnk

                    # There is also estimation error if rate rules are used
                    est_dplnk = kParamEngine.getPartialUncertaintyValue(source, 'Estimation')
                    if est_dplnk:
                        label = 'Estimation {}'.format(reaction.to_chemkin(self.speciesList, kinetics=False))
                        dlnk[label] = est_dplnk

                elif 'PDep' in source:
                    dplnk = kParamEngine.getPartialUncertaintyValue(source, 'PDep', source['PDep'])
                    label = 'PDep {}'.format(reaction.to_chemkin(self.speciesList, kinetics=False))
                    dlnk[label] = dplnk

                elif 'Library' in source:
                    dplnk = kParamEngine.getPartialUncertaintyValue(source, 'Library', source['Library'])
                    label = 'Library {}'.format(reaction.to_chemkin(self.speciesList, kinetics=False))
                    dlnk[label] = dplnk

                elif 'Training' in source:
                    dplnk = kParamEngine.getPartialUncertaintyValue(source, 'Training', source['Training'])
                    family = source['Training'][0]
                    label = 'Training {} {}'.format(family, reaction.to_chemkin(self.speciesList, kinetics=False))
                    dlnk[label] = dplnk

                self.kineticInputUncertainties.append(dlnk)

    def sensitivityAnalysis(self, initialMoleFractions, sensitiveSpecies, T, P, terminationTime,
                            sensitivityThreshold=1e-3, number=10, fileformat='.png'):
        """
        Run sensitivity analysis using the RMG solver in a single ReactionSystem object
        
        initialMoleFractions is a dictionary with Species objects as keys and mole fraction initial conditions
        sensitiveSpecies is a list of sensitive Species objects
        number is the number of top species thermo or reaction kinetics desired to be plotted
        """

        from rmgpy.solver import SimpleReactor, TerminationTime
        from rmgpy.quantity import Quantity
        from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
        from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
        T = Quantity(T)
        P = Quantity(P)
        termination = [TerminationTime(Quantity(terminationTime))]

        reaction_system = SimpleReactor(T=T,
                                        P=P,
                                        initialMoleFractions=initialMoleFractions,
                                        termination=termination,
                                        sensitiveSpecies=sensitiveSpecies,
                                        sensitivityThreshold=sensitivityThreshold)

        # Create the csv worksheets for logging sensitivity
        util.makeOutputSubdirectory(self.outputDirectory, 'solver')
        sens_worksheet = []
        reaction_system_index = 0
        for spec in reaction_system.sensitiveSpecies:
            csvfile_path = os.path.join(self.outputDirectory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1, spec.index))
            sens_worksheet.append(csvfile_path)

        reaction_system.attach(SimulationProfileWriter(
            self.outputDirectory, reaction_system_index, self.speciesList))
        reaction_system.attach(SimulationProfilePlotter(
            self.outputDirectory, reaction_system_index, self.speciesList))

        simulator_settings = SimulatorSettings()  # defaults

        model_settings = ModelSettings()  # defaults
        model_settings.fluxToleranceMoveToCore = 0.1
        model_settings.fluxToleranceInterrupt = 1.0
        model_settings.fluxToleranceKeepInEdge = 0.0

        reaction_system.simulate(
            coreSpecies=self.speciesList,
            coreReactions=self.reactionList,
            edgeSpecies=[],
            edgeReactions=[],
            surfaceSpecies=[],
            surfaceReactions=[],
            modelSettings=model_settings,
            simulatorSettings=simulator_settings,
            sensitivity=True,
            sensWorksheet=sens_worksheet,
        )

        plot_sensitivity(self.outputDirectory, reaction_system_index, reaction_system.sensitiveSpecies, number=number, fileformat=fileformat)

    def localAnalysis(self, sensitiveSpecies, reactionSystemIndex=0, correlated=False, number=10, fileformat='.png'):
        """
        Conduct local uncertainty analysis on the reaction model.
        sensitiveSpecies is a list of sensitive Species objects
        number is the number of highest contributing uncertain parameters desired to be plotted
        fileformat can be either .png, .pdf, or .svg
        """
        output = {}
        for sensSpecies in sensitiveSpecies:
            csvfile_path = os.path.join(self.outputDirectory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reactionSystemIndex+1, sensSpecies.index))
            time, data_list = parseCSVData(csvfile_path)
            # Assign uncertainties
            thermo_data_list = []
            reaction_data_list = []
            for data in data_list:
                if data.species:
                    for species in self.speciesList:
                        if species.to_chemkin() == data.species:
                            index = self.speciesList.index(species)
                            break
                    else:
                        raise Exception('Chemkin name {} of species in the CSV file does not match anything in the '
                                        'species list.'.format(data.species))

                    data.uncertainty = self.thermoInputUncertainties[index]
                    thermo_data_list.append(data)

                if data.reaction:
                    rxn_index = int(data.index) - 1
                    data.uncertainty = self.kineticInputUncertainties[rxn_index]
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
                            correlated_thermo_data[label] = GenericData(data=[data.data[-1] * dpG], uncertainty=1, label=label, species='dummy')
                for data in reaction_data_list:
                    for label, dplnk in data.uncertainty.items():
                        if label in correlated_reaction_data:
                            correlated_reaction_data[label].data[-1] += data.data[-1] * dplnk
                        else:
                            correlated_reaction_data[label] = GenericData(data=[data.data[-1] * dplnk], uncertainty=1, label=label, reaction='dummy')

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
                    data.label = 'k' + str(self.reactionList[data.index-1].index) + ': ' + data.label.split()[-1]

            if correlated:
                folder = os.path.join(self.outputDirectory, 'correlated')
            else:
                folder = os.path.join(self.outputDirectory, 'uncorrelated')
            if not os.path.exists(folder):
                try:
                    os.makedirs(folder)
                except OSError as e:
                    raise OSError('Uncertainty output directory could not be created: {0}'.format(e.message))

            r_path = os.path.join(folder, 'kineticsLocalUncertainty_{0}'.format(sensSpecies.to_chemkin()) + fileformat)
            t_path = os.path.join(folder, 'thermoLocalUncertainty_{0}'.format(sensSpecies.to_chemkin()) + fileformat)
            reaction_uncertainty = ReactionSensitivityPlot(xVar=time, yVar=reaction_data_list, numReactions=number).uncertaintyPlot(total_variance, filename=r_path)
            thermo_uncertainty = ThermoSensitivityPlot(xVar=time, yVar=thermo_data_list, numSpecies=number).uncertaintyPlot(total_variance, filename=t_path)

            output[sensSpecies] = (total_variance, reaction_uncertainty, thermo_uncertainty)

        return output


def process_local_results(results, sensitiveSpecies, number=10):
    """
    Return a dictionary of processed results along with a formatted string
    given results from local uncertainty analysis.
    """
    processed_results = {}
    for spc in sensitiveSpecies:
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
    for spc in sensitiveSpecies:
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
