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

import os.path

import numpy as np

from rmgpy.chemkin import load_chemkin_file
from rmgpy.tools.canteramodel import Cantera, generate_cantera_conditions, get_rmg_species_from_user_species
from rmgpy.tools.plot import GenericPlot, SimulationPlot, find_nearest


def curves_similar(t1, y1, t2, y2, tol):
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
        time_index = find_nearest(t2, timepoint1)
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
    `old_dir`               A directory path containing the chem_annotated.inp and species_dictionary.txt of the old model
    `new_dir`               A directory path containing the chem_annotated.inp and species_dictionary.txt of the new model
    `conditions`            A list of the :class: 'CanteraCondition' objects describing reaction conditions
    `observables`           A dictionary of observables
                            key: 'species', value: a list of the "class" 'Species' that correspond with species mole fraction observables
                            key: 'variable', value: a list of state variable observables, i.e. ['Temperature'] or ['Temperature','Pressure']
                            key: 'ignitionDelay', value: a tuple containing (ignition metric, y_var)
                                                         for example: ('maxDerivative','P')
                                                                      ('maxHalfConcentration', '[OH]')
                                                                      ('maxSpeciesConcentrations',['[CH2]','[O]'])
                                                        see find_ignition_delay function for more details
    'expt_data'             An array of GenericData objects
    'ck2tci'                Indicates whether to convert chemkin to cti mechanism.  If set to False, RMG will convert the species
                            and reaction objects to Cantera objects internally
    ======================= ==============================================================================================


    """

    def __init__(self, title='', old_dir='', new_dir='', observables=None, expt_data=None, ck2cti=True):
        self.title = title
        self.new_dir = new_dir
        self.old_dir = old_dir
        self.conditions = None
        self.expt_data = expt_data if expt_data else []
        self.observables = observables if observables else {}

        # Detect if the transport file exists
        old_transport_path = None
        if os.path.exists(os.path.join(old_dir, 'tran.dat')):
            old_transport_path = os.path.join(old_dir, 'tran.dat')
        new_transport_path = None
        if os.path.exists(os.path.join(new_dir, 'tran.dat')):
            new_transport_path = os.path.join(new_dir, 'tran.dat')

        # define chemkin file path for old and new models
        old_chemkin_path = os.path.join(old_dir, 'chem_annotated.inp')
        new_chemkin_path = os.path.join(new_dir, 'chem_annotated.inp')
        old_species_dict_path = os.path.join(old_dir, 'species_dictionary.txt')
        new_species_dict_path = os.path.join(new_dir, 'species_dictionary.txt')

        surface = False
        if os.path.exists(os.path.join(old_dir, 'chem_annotated-surface.inp')) and \
            os.path.exists(os.path.join(old_dir, 'chem_annotated-gas.inp')):
            surface = True
            old_chemkin_path = os.path.join(old_dir, 'chem_annotated-gas.inp')
            new_chemkin_path = os.path.join(new_dir, 'chem_annotated-gas.inp')
            old_surface_chemkin_path = os.path.join(old_dir, 'chem_annotated-surface.inp')
            new_surface_chemkin_path = os.path.join(new_dir, 'chem_annotated-surface.inp')
            ck2cti = True

        old_species_list, old_reaction_list = load_chemkin_file(
            old_chemkin_path,
            old_species_dict_path,
            old_transport_path
        )
        new_species_list, new_reaction_list = load_chemkin_file(
            new_chemkin_path,
            new_species_dict_path,
            new_transport_path
        )

        old_surface_species_list = None
        new_surface_species_list = None
        new_surface_reaction_list = None
        old_surface_reaction_list = None
        if surface:
            old_surface_species_list, old_surface_reaction_list = load_chemkin_file(
                old_surface_chemkin_path,
                old_species_dict_path,
                old_transport_path
            )
            new_surface_species_list, new_surface_reaction_list = load_chemkin_file(
                new_surface_chemkin_path,
                new_species_dict_path,
                new_transport_path
            )

        self.old_sim = Cantera(species_list=old_species_list,
                            reaction_list=old_reaction_list,
                            output_directory=old_dir,
                            surface_species_list=old_surface_species_list,
                            surface_reaction_list=old_surface_reaction_list,
                            )
        self.new_sim = Cantera(species_list=new_species_list,
                            reaction_list=new_reaction_list,
                            output_directory=new_dir,
                            surface_species_list=new_surface_species_list,
                            surface_reaction_list=new_surface_reaction_list,
                            )

        # load each chemkin file into the cantera model
        if not ck2cti:
            self.old_sim.load_model()
            self.new_sim.load_model()
        else:
            self.old_sim.load_chemkin_model(old_chemkin_path,
                                            transport_file=old_transport_path,
                                            surface_file=old_surface_chemkin_path,
                                            quiet=True)
            self.new_sim.load_chemkin_model(new_chemkin_path,
                                            transport_file=new_transport_path,
                                            surface_file=new_surface_chemkin_path,
                                            quiet=True)

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    def generate_conditions(self, reactor_type_list, reaction_time_list, mol_frac_list, surface_mol_frac_list=None,
                            Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of conditions from from the lists provided. 
        
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactor_type_list`     A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures

        `reaction_time_list`    A tuple object giving the ([list of reaction times], units)
        `mol_frac_list`         A list of molfrac dictionaries with species object keys
                                and mole fraction values
        `surface_mol_frac_list` A list of molfrac dictionaries with species object keys for the surface
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
        
        This saves all the reaction conditions into both the old and new cantera jobs.
        """
        # Store the conditions in the observables test case, for bookkeeping
        self.conditions = generate_cantera_conditions(reactor_type_list, reaction_time_list, mol_frac_list,
                                                      surface_mol_frac_list=surface_mol_frac_list,
                                                      Tlist=Tlist, Plist=Plist, Vlist=Vlist)

        if surface_mol_frac_list is None:
            surface_mol_frac_list = []  # initialize here as list to avoid mutable default argument in function definition

        # Map the mole fractions dictionaries to species objects from the old and new models
        old_mol_frac_list = []
        new_mol_frac_list = []

        for mol_frac in mol_frac_list:
            old_condition = {}
            new_condition = {}
            old_species_dict = get_rmg_species_from_user_species(list(mol_frac.keys()), self.old_sim.species_list)
            new_species_dict = get_rmg_species_from_user_species(list(mol_frac.keys()), self.new_sim.species_list)
            for smiles, molfrac in mol_frac.items():
                if old_species_dict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the old model!'.format(smiles))
                if new_species_dict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the new model!'.format(smiles))

                old_condition[old_species_dict[smiles]] = molfrac
                new_condition[new_species_dict[smiles]] = molfrac
            old_mol_frac_list.append(old_condition)
            new_mol_frac_list.append(new_condition)

        if self.old_sim.surface:
            old_surf_mol_frac_list = []
            new_surf_mol_frac_list = []
            for surf_mol_frac in surface_mol_frac_list:
                old_surface_condition = {}
                new_surface_condition = {}
                old_surface_species_dict = get_rmg_species_from_user_species(list(surf_mol_frac.keys()), self.old_sim.surface_species_list)
                new_surface_species_dict = get_rmg_species_from_user_species(list(surf_mol_frac.keys()), self.new_sim.surface_species_list)
                for smiles, molfrac in surf_mol_frac.items():
                    if old_surface_species_dict[smiles] is None:
                        raise Exception('SMILES {0} was not found in the old model!'.format(smiles))
                    if new_surface_species_dict[smiles] is None:
                        raise Exception('SMILES {0} was not found in the new model!'.format(smiles))

                    old_surface_condition[old_surface_species_dict[smiles]] = molfrac
                    new_surface_condition[new_surface_species_dict[smiles]] = molfrac
                old_surf_mol_frac_list.append(old_surface_condition)
                new_surf_mol_frac_list.append(new_surface_condition)
        else:
            old_surf_mol_frac_list = [None]
            new_surf_mol_frac_list = [None]

        # Generate the conditions in each simulation
        self.old_sim.generate_conditions(reactor_type_list, reaction_time_list, old_mol_frac_list,
                                         surface_mol_frac_list=old_surf_mol_frac_list,
                                         Tlist=Tlist, Plist=Plist, Vlist=Vlist)
        self.new_sim.generate_conditions(reactor_type_list, reaction_time_list, new_mol_frac_list,
                                         surface_mol_frac_list=new_surf_mol_frac_list,
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

        old_condition_data, new_condition_data = self.run_simulations()

        conditions_broken = []
        variables_failed = []

        print('{0} Comparison'.format(self))
        # Check the species profile observables
        if 'species' in self.observables:
            if self.old_sim.surface_species_list:
                old_species_dict = get_rmg_species_from_user_species(self.observables['species'], self.old_sim.species_list + self.old_sim.surface_species_list)
                new_species_dict = get_rmg_species_from_user_species(self.observables['species'], self.new_sim.species_list + self.new_sim.surface_species_list)
            else:
                old_species_dict = get_rmg_species_from_user_species(self.observables['species'], self.old_sim.species_list)
                new_species_dict = get_rmg_species_from_user_species(self.observables['species'], self.new_sim.species_list)

        # Check state variable observables
        implemented_variables = ['temperature', 'pressure']
        if 'variable' in self.observables:
            for item in self.observables['variable']:
                if item.lower() not in implemented_variables:
                    print('Observable variable {0} not yet implemented'.format(item))

        fail_header = '\nThe following observables did not match:\n'
        fail_header_printed = False
        for i in range(len(old_condition_data)):
            time_old, data_list_old, reaction_sensitivity_data_old, thermodynamic_sensitivity_data_old = old_condition_data[i]
            time_new, data_list_new, reaction_sensitivity_data_new, thermodynamic_sensitivity_data_new = new_condition_data[i]

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
                        print('⁉️ No RMG species found for observable species {0} in old model.'.format(smiles))
                        fail = True
                    if new_rmg_species:
                        variable_new = next((data for data in data_list_new if data.species == new_rmg_species), None)
                    else:
                        print('⁉️ No RMG species found for observable species {0} in new model.'.format(smiles))
                        fail = True

                    if fail is False:
                        if not curves_similar(time_old.data, variable_old.data, time_new.data, variable_new.data, tol):
                            fail = True

                        # Try plotting only when species are found in both models
                        if plot:
                            old_species_plot = SimulationPlot(x_var=time_old, y_var=variable_old)
                            new_species_plot = SimulationPlot(x_var=time_new, y_var=variable_new)
                            old_species_plot.compare_plot(new_species_plot,
                                                          title='Observable Species {0} Comparison'.format(smiles),
                                                          ylabel='Mole Fraction',
                                                          filename='condition_{0}_species_{1}.png'.format(i + 1, smiles))

                    # Append to failed variables or conditions if this test failed
                    if fail:
                        if not fail_header_printed:
                            print(fail_header)
                            fail_header_printed = True
                        if i not in conditions_broken:
                            conditions_broken.append(i)
                        print("❌ Observable species {0} varied by more than {1:.3f} on average between old model {2} and "
                              "new model {3} in condition {4:d}.".format(smiles, tol, variable_old.label,
                                                                         variable_new.label, i + 1))
                        variables_failed.append((self.conditions[i], smiles, variable_old, variable_new))

            # Compare state variable observables
            if 'variable' in self.observables:
                for varName in self.observables['variable']:
                    variable_old = next((data for data in data_list_old if data.label == varName), None)
                    variable_new = next((data for data in data_list_new if data.label == varName), None)
                    if not curves_similar(time_old.data, variable_old.data, time_new.data, variable_new.data, 0.05):
                        if i not in conditions_broken:
                            conditions_broken.append(i)
                        if not fail_header_printed:
                            fail_header_printed = True
                            print(fail_header)

                        print("❌ Observable variable {0} varied by more than {1:.3f} on average between old model and "
                              "new model in condition {2:d}.".format(variable_old.label, tol, i + 1))
                        variables_failed.append((self.conditions[i], varName, variable_old, variable_new))

                    if plot:
                        old_var_plot = GenericPlot(x_var=time_old, y_var=variable_old)
                        new_var_plot = GenericPlot(x_var=time_new, y_var=variable_new)
                        old_var_plot.compare_plot(new_var_plot,
                                                  title='Observable Variable {0} Comparison'.format(varName),
                                                  filename='condition_{0}_variable_{1}.png'.format(i + 1, varName))

            # Compare ignition delay observables
            if 'ignitionDelay' in self.observables:
                print('Ignition delay observable comparison not implemented yet.')

        if fail_header_printed:
            print('')
            print('⚠️ The following reaction conditions had some discrepancies:')
            for index in conditions_broken:
                print("Condition {0:d}:".format(index + 1))
                print(str(self.conditions[index]))

            return variables_failed
        else:
            print('')
            print('✅ All Observables varied by less than {0:.3f} on average between old model and '
                  'new model in all conditions!'.format(tol))
            print('')

    def run_simulations(self):
        """
        Run a selection of conditions in Cantera and return
        generic data objects containing the time, pressure, temperature,
        and mole fractions from the simulations.

        Returns (old_condition_data, new_condition_data)
        where conditionData is a list of of tuples: (time, dataList) for each condition in the same order as conditions
        time is a GenericData object which gives the time domain for each profile
        dataList is a list of GenericData objects for the temperature, profile, and mole fraction of major species
        """
        old_condition_data = self.old_sim.simulate()
        new_condition_data = self.new_sim.simulate()
        return (old_condition_data, new_condition_data)
