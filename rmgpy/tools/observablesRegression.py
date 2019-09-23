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

from __future__ import division, print_function

import os.path

import numpy as np

from rmgpy.chemkin import load_chemkin_file
from rmgpy.tools.canteraModel import Cantera, generateCanteraConditions, getRMGSpeciesFromUserSpecies
from rmgpy.tools.plot import GenericPlot, SimulationPlot, findNearest


def curvesSimilar(t1, y1, t2, y2, tol):
    """
    This function returns True if the two given curves are similar enough within tol. Otherwise returns False.

    t1: time/domain of standard curve we assume to be correct
    y1: values of standard curve, usually either temperature in (K) or log of a mol fraction
    t2: time/domain of test curve
    y2: values of test curve, usually either temperature in (K) or log of a mol fraction

    The test curve is first synchronized to the standard curve using geatNearestTime function. We then calculate the
    value of abs((y1-y2')/y1), giving us a normalized difference for every point. If the average value of these
    differences is less than tol, we say the curves are similar.

    We choose this criteria because it is compatible with step functions we expect to see in ignition systems.
    """
    # Make synchornized version of t2,y2 called t2sync,y2sync.
    t2sync = np.zeros_like(t1)
    y2sync = np.zeros_like(t1)
    for i, timepoint1 in enumerate(t1):
        time_index = findNearest(t2, timepoint1)
        t2sync[i] = t2[time_index]
        y2sync[i] = y2[time_index]

    # Get R^2 value equivalent:
    normalized_error = abs((y1 - y2sync) / y1)
    normalized_error = sum(normalized_error) / len(y1)

    if normalized_error > tol:
        return False
    else:
        return True


class ObservablesTestCase(object):
    """
    We use this class to run regressive tests

    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    `title`                 A string describing the test. For regressive tests, should be same as example's name.
    `oldDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the old model
    `newDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the new model
    `conditions`            A list of the :class: 'CanteraCondition' objects describing reaction conditions
    `observables`           A dictionary of observables
                            key: 'species', value: a list of the "class" 'Species' that correspond with species mole fraction observables
                            key: 'variable', value: a list of state variable observables, i.e. ['Temperature'] or ['Temperature','Pressure']
                            key: 'ignitionDelay', value: a tuple containing (ignition metric, yVar)
                                                         for example: ('maxDerivative','P')
                                                                      ('maxHalfConcentration', '[OH]')
                                                                      ('maxSpeciesConcentrations',['[CH2]','[O]'])
                                                        see findIgnitionDelay function for more details
    'exptData'              An array of GenericData objects
    'ck2tci'                Indicates whether to convert chemkin to cti mechanism.  If set to False, RMG will convert the species
                            and reaction objects to Cantera objects internally
    ======================= ==============================================================================================


    """

    def __init__(self, title='', oldDir='', newDir='', observables=None, exptData=None, ck2cti=True):
        self.title = title
        self.newDir = newDir
        self.oldDir = oldDir
        self.conditions = None
        self.exptData = exptData if exptData else []
        self.observables = observables if observables else {}

        # Detect if the transport file exists
        old_transport_path = None
        if os.path.exists(os.path.join(oldDir, 'tran.dat')):
            old_transport_path = os.path.join(oldDir, 'tran.dat')
        new_transport_path = None
        if os.path.exists(os.path.join(newDir, 'tran.dat')):
            new_transport_path = os.path.join(newDir, 'tran.dat')

        # load the species and reactions from each model
        old_species_list, old_reaction_list = load_chemkin_file(os.path.join(oldDir, 'chem_annotated.inp'),
                                                                os.path.join(oldDir, 'species_dictionary.txt'),
                                                                old_transport_path)

        new_species_list, new_reaction_list = load_chemkin_file(os.path.join(newDir, 'chem_annotated.inp'),
                                                                os.path.join(newDir, 'species_dictionary.txt'),
                                                                new_transport_path)

        self.oldSim = Cantera(speciesList=old_species_list,
                              reactionList=old_reaction_list,
                              outputDirectory=oldDir)
        self.newSim = Cantera(speciesList=new_species_list,
                              reactionList=new_reaction_list,
                              outputDirectory=newDir)

        # load each chemkin file into the cantera model
        if not ck2cti:
            self.oldSim.loadModel()
            self.newSim.loadModel()
        else:
            self.oldSim.loadChemkinModel(os.path.join(oldDir, 'chem_annotated.inp'), transportFile=old_transport_path,
                                         quiet=True)
            self.newSim.loadChemkinModel(os.path.join(newDir, 'chem_annotated.inp'), transportFile=new_transport_path,
                                         quiet=True)

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    def generateConditions(self, reactorTypeList, reactionTimeList, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of conditions from from the lists provided. 
        
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactorTypeList`        A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures

        `reactionTimeList`      A tuple object giving the ([list of reaction times], units)
        `molFracList`           A list of molfrac dictionaries with species object keys
                               and mole fraction values
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
        
        This saves all the reaction conditions into both the old and new cantera jobs.
        """
        # Store the conditions in the observables test case, for bookkeeping
        self.conditions = generateCanteraConditions(reactorTypeList, reactionTimeList, molFracList, Tlist=Tlist,
                                                    Plist=Plist, Vlist=Vlist)

        # Map the mole fractions dictionaries to species objects from the old and new models
        old_mol_frac_list = []
        new_mol_frac_list = []

        for molFracCondition in molFracList:
            old_condition = {}
            new_condition = {}
            old_species_dict = getRMGSpeciesFromUserSpecies(list(molFracCondition.keys()), self.oldSim.speciesList)
            new_species_dict = getRMGSpeciesFromUserSpecies(list(molFracCondition.keys()), self.newSim.speciesList)
            for smiles, molfrac in molFracCondition.items():
                if old_species_dict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the old model!'.format(smiles))
                if new_species_dict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the new model!'.format(smiles))

                old_condition[old_species_dict[smiles]] = molfrac
                new_condition[new_species_dict[smiles]] = molfrac
            old_mol_frac_list.append(old_condition)
            new_mol_frac_list.append(new_condition)

        # Generate the conditions in each simulation
        self.oldSim.generateConditions(reactorTypeList, reactionTimeList, old_mol_frac_list,
                                       Tlist=Tlist, Plist=Plist, Vlist=Vlist)
        self.newSim.generateConditions(reactorTypeList, reactionTimeList, new_mol_frac_list,
                                       Tlist=Tlist, Plist=Plist, Vlist=Vlist)

    def compare(self, tol, plot=False):
        """
        Compare the old and new model
        'tol':  average error acceptable between old and new model for variables
        `plot`: if set to True, it will comparison plots of the two models comparing their species.

        Returns a list of variables failed in a list of tuples in the format:
        
        (CanteraCondition, variable label, variable_old, variable_new)

        """
        # Ignore Inerts
        inert_list = ['[Ar]', '[He]', '[N#N]', '[Ne]']

        old_condition_data, new_condition_data = self.runSimulations()

        conditions_broken = []
        variables_failed = []

        print('')
        print('{0} Comparison'.format(self))
        print('================')
        # Check the species profile observables
        if 'species' in self.observables:
            old_species_dict = getRMGSpeciesFromUserSpecies(self.observables['species'], self.oldSim.speciesList)
            new_species_dict = getRMGSpeciesFromUserSpecies(self.observables['species'], self.newSim.speciesList)

        # Check state variable observables 
        implemented_variables = ['temperature', 'pressure']
        if 'variable' in self.observables:
            for item in self.observables['variable']:
                if item.lower() not in implemented_variables:
                    print('Observable variable {0} not yet implemented'.format(item))

        fail_header = '\nThe following observables did not match:\n'
        fail_header_printed = False
        for i in range(len(old_condition_data)):
            time_old, data_list_old, reaction_sensitivity_data_old = old_condition_data[i]
            time_new, data_list_new, reaction_sensitivity_data_new = new_condition_data[i]

            # Compare species observables
            if 'species' in self.observables:
                smiles_list = []  # This is to make sure we don't have species with duplicate smiles
                multiplicity_list = ['', '(S)', '(D)', '(T)', '(Q)']  # list ot add multiplcity
                for species in self.observables['species']:

                    smiles = species.molecule[0].to_smiles()  # For purpose of naming the plot only
                    if smiles in smiles_list: smiles = smiles + multiplicity_list[species.molecule[0].multiplicity]
                    smiles_list.append(smiles)

                    fail = False
                    old_rmg_species = old_species_dict[species]
                    new_rmg_species = new_species_dict[species]

                    if old_rmg_species:
                        variable_old = next((data for data in data_list_old if data.species == old_rmg_species), None)
                    else:
                        print('No RMG species found for observable species {0} in old model.'.format(smiles))
                        fail = True
                    if new_rmg_species:
                        variable_new = next((data for data in data_list_new if data.species == new_rmg_species), None)
                    else:
                        print('No RMG species found for observable species {0} in new model.'.format(smiles))
                        fail = True

                    if fail is False:
                        if not curvesSimilar(time_old.data, variable_old.data, time_new.data, variable_new.data, tol):
                            fail = True

                        # Try plotting only when species are found in both models
                        if plot:
                            old_species_plot = SimulationPlot(xVar=time_old, yVar=variable_old)
                            new_species_plot = SimulationPlot(xVar=time_new, yVar=variable_new)
                            old_species_plot.comparePlot(new_species_plot,
                                                         title='Observable Species {0} Comparison'.format(smiles),
                                                         ylabel='Mole Fraction',
                                                         filename='condition_{0}_species_{1}.png'.format(i + 1, smiles))

                    # Append to failed variables or conditions if this test failed
                    if fail:
                        if not fail_header_printed:
                            print(fail_header)
                            fail_header_printed = True
                        if i not in conditions_broken: conditions_broken.append(i)
                        print("Observable species {0} varied by more than {1:.3f} on average between old model {2} and "
                              "new model {3} in condition {4:d}.".format(smiles, tol, variable_old.label,
                                                                         variable_new.label, i + 1))
                        variables_failed.append((self.conditions[i], smiles, variable_old, variable_new))

            # Compare state variable observables
            if 'variable' in self.observables:
                for varName in self.observables['variable']:
                    variable_old = next((data for data in data_list_old if data.label == varName), None)
                    variable_new = next((data for data in data_list_new if data.label == varName), None)
                    if not curvesSimilar(time_old.data, variable_old.data, time_new.data, variable_new.data, 0.05):
                        if i not in conditions_broken:
                            conditions_broken.append(i)
                        if not fail_header_printed:
                            fail_header_printed = True
                            print(fail_header)

                        print("Observable variable {0} varied by more than {1:.3f} on average between old model and "
                              "new model in condition {2:d}.".format(variable_old.label, tol, i + 1))
                        variables_failed.append((self.conditions[i], varName, variable_old, variable_new))

                    if plot:
                        old_var_plot = GenericPlot(xVar=time_old, yVar=variable_old)
                        new_var_plot = GenericPlot(xVar=time_new, yVar=variable_new)
                        old_var_plot.comparePlot(new_var_plot,
                                                 title='Observable Variable {0} Comparison'.format(varName),
                                                 filename='condition_{0}_variable_{1}.png'.format(i + 1, varName))

            # Compare ignition delay observables
            if 'ignitionDelay' in self.observables:
                print('Ignition delay observable comparison not implemented yet.')

        if fail_header_printed:
            print('')
            print('The following reaction conditions were had some discrepancies:')
            print('')
            for index in conditions_broken:
                print("Condition {0:d}:".format(index + 1))
                print(str(self.conditions[index]))
                print('')

            return variables_failed
        else:
            print('')
            print('All Observables varied by less than {0:.3f} on average between old model and '
                  'new model in all conditions!'.format(tol))
            print('')

    def runSimulations(self):
        """
        Run a selection of conditions in Cantera and return
        generic data objects containing the time, pressure, temperature,
        and mole fractions from the simulations.

        Returns (old_condition_data, new_condition_data)
        where conditionData is a list of of tuples: (time, dataList) for each condition in the same order as conditions
        time is a GenericData object which gives the time domain for each profile
        dataList is a list of GenericData objects for the temperature, profile, and mole fraction of major species
        """
        old_condition_data = self.oldSim.simulate()
        new_condition_data = self.newSim.simulate()
        return (old_condition_data, new_condition_data)
