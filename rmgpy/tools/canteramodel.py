#!/usr/bin/env python3

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

import os.path

import cantera as ct
import numpy as np

from rmgpy.chemkin import get_species_identifier
from rmgpy.quantity import Quantity
from rmgpy.tools.data import GenericData
from rmgpy.tools.plot import GenericPlot, SimulationPlot, ReactionSensitivityPlot, ThermoSensitivityPlot


class CanteraCondition(object):
    """
    This class organizes the inputs needed for a cantera simulation

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `reactor_type`          A string of the cantera reactor type. List of supported types below:
        IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor 
                            for ideal gas mixtures (the same as RMG's SimpleReactor)

    `reaction_time`         A tuple object giving the (reaction time, units)
    `mol_frac`              A dictionary giving the initial mol Fractions. Keys are species objects and the values are floats

    To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
    `T0`                    A tuple giving the (initial temperature, units) which reconstructs a Quantity object
    'P0'                    A tuple giving the (initial pressure, units) which reconstructs a Quantity object
    'V0'                    A tuple giving the (initial specific volume, units) which reconstructs a Quantity object
    ======================= ====================================================


    """

    def __init__(self, reactor_type, reaction_time, mol_frac, T0=None, P0=None, V0=None):
        self.reactor_type = reactor_type
        self.reaction_time = Quantity(reaction_time)

        # Normalize initialMolFrac if not already done:
        if sum(mol_frac.values()) != 1.00:
            total = sum(mol_frac.values())
            for species, value in mol_frac.items():
                mol_frac[species] = value / total

        self.mol_frac = mol_frac

        # Check to see that one of the three attributes T0, P0, and V0 is less unspecified
        props = [T0, P0, V0]
        total = 0
        for prop in props:
            if prop is None: total += 1

        if not total == 1:
            raise Exception("Cantera conditions must leave one of T0, P0, and V0 state variables unspecified")

        self.T0 = Quantity(T0) if T0 else None
        self.P0 = Quantity(P0) if P0 else None
        self.V0 = Quantity(V0) if V0 else None

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = "CanteraCondition("
        string += 'reactor_type="{0}", '.format(self.reactor_type)
        string += 'reaction_time={}, '.format(self.reaction_time.__repr__())
        string += 'mol_frac={0}, '.format(self.mol_frac.__repr__())
        if self.T0: string += 'T0={}, '.format(self.T0.__repr__())
        if self.P0: string += 'P0={}, '.format(self.P0.__repr__())
        if self.V0: string += 'V0={}, '.format(self.V0__repr__())
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the condition.
        """
        string = ""
        string += 'Reactor Type: {0}\n'.format(self.reactor_type)
        string += 'Reaction Time: {}\n'.format(self.reaction_time)
        if self.T0: string += 'T0: {}\n'.format(self.T0)
        if self.P0: string += 'P0: {}\n'.format(self.P0)
        if self.V0: string += 'V0: {}\n'.format(self.V0)
        # ConvertMolFrac to SMILES for keys for display
        pretty_mol_frac = {}
        for key, value in self.mol_frac.items():
            pretty_mol_frac[key.molecule[0].to_smiles()] = value
        string += 'Initial Mole Fractions: {0}'.format(pretty_mol_frac.__repr__())
        return string


def generate_cantera_conditions(reactor_type_list, reaction_time_list, mol_frac_list, Tlist=None, Plist=None, Vlist=None):
    """
    Creates a list of cantera conditions from from the arguments provided.

    ======================= ====================================================
    Argument                Description
    ======================= ====================================================
    `reactor_type_list`        A list of strings of the cantera reactor type. List of supported types below:
        IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor
                            for ideal gas mixtures (the same as RMG's SimpleReactor)

    `reaction_time_list`      A tuple object giving the ([list of reaction times], units)
    `mol_frac_list`           A list of molfrac dictionaries with species object keys
                           and mole fraction values
    To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
    `T0List`                A tuple giving the ([list of initial temperatures], units)
    'P0List'                A tuple giving the ([list of initial pressures], units)
    'V0List'                A tuple giving the ([list of initial specific volumes], units)


    This saves all the reaction conditions into the Cantera class.
    """

    # Create individual ScalarQuantity objects for Tlist, Plist, Vlist, and reaction_time_list
    if Tlist:
        Tlist = Quantity(Tlist)  # Be able to create a Quantity object from it first
        Tlist = [(Tlist.value[i], Tlist.units) for i in range(len(Tlist.value))]
    if Plist:
        Plist = Quantity(Plist)
        Plist = [(Plist.value[i], Plist.units) for i in range(len(Plist.value))]
    if Vlist:
        Vlist = Quantity(Vlist)
        Vlist = [(Vlist.value[i], Vlist.units) for i in range(len(Vlist.value))]
    if reaction_time_list:
        reaction_time_list = Quantity(reaction_time_list)
        reaction_time_list = [(reaction_time_list.value[i], reaction_time_list.units)
                            for i in range(len(reaction_time_list.value))]

    conditions = []

    if Tlist is None:
        for reactor_type in reactor_type_list:
            for reaction_time in reaction_time_list:
                for mol_frac in mol_frac_list:
                    for P in Plist:
                        for V in Vlist:
                            conditions.append(CanteraCondition(reactor_type, reaction_time, mol_frac, P0=P, V0=V))

    elif Plist is None:
        for reactor_type in reactor_type_list:
            for reaction_time in reaction_time_list:
                for mol_frac in mol_frac_list:
                    for T in Tlist:
                        for V in Vlist:
                            conditions.append(CanteraCondition(reactor_type, reaction_time, mol_frac, T0=T, V0=V))

    elif Vlist is None:
        for reactor_type in reactor_type_list:
            for reaction_time in reaction_time_list:
                for mol_frac in mol_frac_list:
                    for T in Tlist:
                        for P in Plist:
                            conditions.append(CanteraCondition(reactor_type, reaction_time, mol_frac, T0=T, P0=P))

    else:
        raise Exception("Cantera conditions must leave one of T0, P0, and V0 state variables unspecified")
    return conditions


class Cantera(object):
    """
    This class contains functions associated with an entire Cantera job
    """

    def __init__(self, species_list=None, reaction_list=None, canteraFile='', output_directory='', conditions=None,
                 sensitive_species=None, thermo_SA=False):
        """
        `species_list`: list of RMG species objects
        `reaction_list`: list of RMG reaction objects
        `reaction_map`: dict mapping the RMG reaction index within the `reaction_list` to cantera model reaction(s) indices
        `canteraFile` path of the chem.yaml file associated with this job
        `conditions`: a list of `CanteraCondition` objects
        `sensitive_species`: a list of RMG species objects for conductng sensitivity analysis on
        `thermo_SA`: a boolean indicating whether or not to run thermo SA. By default, if sensitive_species is given,
                      only kinetic_SA will be calculated and it must be additionally specified to perform thermo SA.
        """
        self.species_list = species_list
        self.reaction_list = reaction_list
        self.reaction_map = {}
        self.model = ct.Solution(canteraFile) if canteraFile else None
        self.output_directory = output_directory if output_directory else os.getcwd()
        self.conditions = conditions if conditions else []
        self.sensitive_species = sensitive_species if sensitive_species else []
        self.thermo_SA = thermo_SA

        # Make output directory if it does not yet exist:
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except:
                raise Exception('Cantera output directory could not be created.')

    def generate_conditions(self, reactor_type_list, reaction_time_list, mol_frac_list, Tlist=None, Plist=None, Vlist=None):
        """
        This saves all the reaction conditions into the Cantera class.
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactor_type_list`        A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor 
                                for ideal gas mixtures (the same as RMG's SimpleReactor)

        `reaction_time_list`      A tuple object giving the ([list of reaction times], units)
        `mol_frac_list`           A list of molfrac dictionaries with species object keys
                               and mole fraction values
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
        """

        self.conditions = generate_cantera_conditions(reactor_type_list, reaction_time_list, mol_frac_list, Tlist, Plist, Vlist)

    def load_model(self):
        """
        Load a cantera Solution model from the job's own species_list and reaction_list attributes
        """

        ct_species = [spec.to_cantera(use_chemkin_identifier=True) for spec in self.species_list]

        self.reaction_map = {}
        ct_reactions = []
        for rxn in self.reaction_list:
            index = len(ct_reactions)

            converted_reactions = rxn.to_cantera(self.species_list, use_chemkin_identifier=True)

            if isinstance(converted_reactions, list):
                indices = list(range(index, index + len(converted_reactions)))
                ct_reactions.extend(converted_reactions)
            else:
                indices = [index]
                ct_reactions.append(converted_reactions)

            self.reaction_map[self.reaction_list.index(rxn)] = indices

        self.model = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                                 species=ct_species, reactions=ct_reactions)

    def refresh_model(self):
        """
        Modification to thermo requires that the cantera model be refreshed to 
        recalculate reverse rate coefficients and equilibrium constants... 
        As soon as cantera has its own Kinetics().modify_thermo function in place,
        this function may be deprecated.
        """
        ct_reactions = self.model.reactions()
        ct_species = self.model.species()

        self.model = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                                 species=ct_species, reactions=ct_reactions)

    def load_chemkin_model(self, chemkin_file, transport_file=None, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.yaml 
        and save it in the output_directory
        Then load it into self.model
        """
        from cantera import ck2yaml

        base = os.path.basename(chemkin_file)
        base_name = os.path.splitext(base)[0]
        out_name = os.path.join(self.output_directory, base_name + ".yaml")
        if os.path.exists(out_name):
            os.remove(out_name)
        parser = ck2yaml.Parser()
        parser.convert_mech(chemkin_file, transport_file=transport_file, out_name=out_name, **kwargs)
        self.model = ct.Solution(out_name)

    def modify_reaction_kinetics(self, rmg_reaction_index, rmg_reaction):
        """
        Modify the corresponding cantera reaction's kinetics to match 
        the reaction kinetics of an `rmg_reaction`, using the `rmg_reaction_index` to
        map to the corresponding reaction in the cantera model. Note that
        this method only works if there is a reaction_map available (therefore only when the cantera model
        is generated directly from rmg objects and not from a chemkin file)
        """
        indices = self.reaction_map[rmg_reaction_index]
        modified_ct_reactions = rmg_reaction.to_cantera(self.species_list, use_chemkin_identifier=True)
        if not isinstance(modified_ct_reactions, list):
            modified_ct_reactions = [modified_ct_reactions]

        for i in range(len(indices)):
            self.model.modify_reaction(indices[i], modified_ct_reactions[i])

    def modify_species_thermo(self, rmg_species_index, rmg_species, use_chemkin_identifier=False):
        """
        Modify the corresponding cantera species thermo to match that of a
        `rmg_species` object, given the `rmg_species_index` which indicates the
        index at which this species appears in the `species_list`
        """
        modified_ct_species = rmg_species.to_cantera(use_chemkin_identifier=use_chemkin_identifier)
        ct_species = self.model.species(rmg_species_index)
        ct_species.thermo = modified_ct_species.thermo

    def plot(self, data, top_species=10, top_sensitive_reactions=10, top_sensitive_species=10):
        """
        Plots data from the simulations from this cantera job.
        Takes data in the format of a list of tuples containing (time, [list of temperature, pressure, and species data]) 
        
        3 plots will be created for each condition:
        - T vs. time
        - P vs. time
        - Maximum species mole fractions (the number of species plotted is based on the `top_species` argument)
        
        Reaction sensitivity plots will also be plotted automatically if there were sensitivities evaluated.
        The number of reactions to be plotted is defined by the `top_sensitive_reactions` argument.
        
        """
        num_ct_reactions = len(self.model.reactions())
        num_ct_species = len(self.model.species())
        for i, condition_data in enumerate(data):
            time, data_list, reaction_sensitivity_data, thermodynamic_sensitivity_data = condition_data
            # In RMG, any species with an index of -1 is an inert and should not be plotted
            inert_list = [species for species in self.species_list if species.index == -1]

            T_data = data_list[0]
            P_data = data_list[1]
            species_data = [data for data in data_list if data.species not in inert_list]

            # plot
            GenericPlot(x_var=time, y_var=T_data).plot(
                os.path.join(self.output_directory, '{0}_temperature.png'.format(i + 1)))
            GenericPlot(x_var=time, y_var=P_data).plot(
                os.path.join(self.output_directory, '{0}_pressure.png'.format(i + 1)))
            SimulationPlot(x_var=time, y_var=species_data, num_species=top_species, ylabel='Mole Fraction').plot(
                os.path.join(self.output_directory, '{0}_mole_fractions.png'.format(i + 1)))

            for j, species in enumerate(self.sensitive_species):
                ReactionSensitivityPlot(x_var=time,
                                        y_var=reaction_sensitivity_data[j * num_ct_reactions:(j + 1) * num_ct_reactions],
                                        num_reactions=top_sensitive_reactions).barplot(
                    os.path.join(self.output_directory, '{0}_{1}_reaction_sensitivity.png'.format(i + 1, species.to_chemkin())))
                if self.thermo_SA:
                    ThermoSensitivityPlot(x_var=time,
                                          y_var=thermodynamic_sensitivity_data[j * num_ct_species:(j + 1) * num_ct_species]*4184000,
                                          xlabel='dln(c)/d(H_i) [(kcal/mol)^-1]',
                                          num_species=top_sensitive_species).barplot(
                        os.path.join(self.output_directory, '{0}_{1}_thermo_sensitivity.png'.format(i + 1, species.to_chemkin())))

    def simulate(self):
        """
        Run all the conditions as a cantera simulation.
        Returns the data as a list of tuples containing: (time, [list of temperature, pressure, and species data]) 
            for each reactor condition
        """
        # Get all the cantera names for the species
        species_names_list = [get_species_identifier(species) for species in self.species_list]
        inert_index_list = [self.species_list.index(species) for species in self.species_list if species.index == -1]

        all_data = []
        for condition in self.conditions:

            # First translate the mol_frac from species objects to species names
            new_mol_frac = {}
            for key, value in condition.mol_frac.items():
                newkey = get_species_identifier(key)
                new_mol_frac[newkey] = value

            # Set Cantera simulation conditions
            if condition.V0 is None:
                self.model.TPX = condition.T0.value_si, condition.P0.value_si, new_mol_frac
            elif condition.P0 is None:
                self.model.TDX = condition.T0.value_si, 1.0 / condition.V0.value_si, new_mol_frac
            else:
                raise Exception(
                    "Cantera conditions in which T0 and P0 or T0 and V0 are not the specified state variables are not yet implemented.")

            # Choose reactor
            if condition.reactor_type == 'IdealGasReactor':
                cantera_reactor = ct.IdealGasReactor(self.model)
            elif condition.reactor_type == 'IdealGasConstPressureReactor':
                cantera_reactor = ct.IdealGasConstPressureReactor(contents=self.model)
            elif condition.reactor_type == 'IdealGasConstPressureTemperatureReactor':
                cantera_reactor = ct.IdealGasConstPressureReactor(contents=self.model, energy='off')
            else:
                raise Exception('Other types of reactor conditions are currently not supported')

            # Run this individual condition as a simulation
            cantera_simulation = ct.ReactorNet([cantera_reactor])

            num_ct_reactions = len(self.model.reactions())
            num_ct_species = len(self.model.species())
            if self.sensitive_species:
                if ct.__version__ == '2.2.1':
                    print('Warning: Cantera version 2.2.1 may not support sensitivity analysis unless SUNDIALS was used during compilation.')
                    print('Warning: Upgrade to newer of Cantera in anaconda using the command "conda update -c rmg cantera"')
                # Add all the reactions as part of the analysis
                for i in range(num_ct_reactions):
                    cantera_reactor.add_sensitivity_reaction(i)
                # if thermo SA is requested, add all species enthalpies
                if self.thermo_SA:
                    for i in range(num_ct_species):
                        cantera_reactor.add_sensitivity_species_enthalpy(i)
                # Set the tolerances for the sensitivity coefficients
                cantera_simulation.rtol_sensitivity = 1e-4
                cantera_simulation.atol_sensitivity = 1e-6

            # Initialize the variables to be saved
            times = []
            temperature = []
            pressure = []
            species_data = []
            kinetic_sensitivity_data = []
            thermo_sensitivity_data = []

            # Begin integration
            time = 0.0
            # Run the simulation over 100 time points
            while cantera_simulation.time < condition.reaction_time.value_si:

                # Advance the state of the reactor network in time from the current time to time t [s], taking as many integrator timesteps as necessary.
                cantera_simulation.step()
                times.append(cantera_simulation.time)
                temperature.append(cantera_reactor.T)
                pressure.append(cantera_reactor.thermo.P)
                species_data.append(cantera_reactor.thermo[species_names_list].X)

                if self.sensitive_species:
                    # Cantera returns mass-based sensitivities rather than molar concentration or mole fraction based sensitivities.
                    # The equation for converting between them is:
                    # 
                    # d ln xi = d ln wi - sum_(species i) (dln wi) (xi)
                    # 
                    # where xi is the mole fraction of species i and wi is the mass fraction of species i

                    mass_frac_sensitivity_array = cantera_simulation.sensitivities()
                    if condition.reactor_type == 'IdealGasReactor':
                        # Row 0: mass, Row 1: volume, Row 2: internal energy or temperature, Row 3+: mass fractions of species
                        mass_frac_sensitivity_array = mass_frac_sensitivity_array[3:, :]
                    elif condition.reactor_type == 'IdealGasConstPressureReactor' or condition.reactor_type == 'IdealGasConstPressureTemperatureReactor':
                        # Row 0: mass, Row 1: enthalpy or temperature, Row 2+: mass fractions of the species
                        mass_frac_sensitivity_array = mass_frac_sensitivity_array[2:, :]
                    else:
                        raise Exception('Other types of reactor conditions are currently not supported')

                    for i in range(len(mass_frac_sensitivity_array)):
                        mass_frac_sensitivity_array[i] *= species_data[-1][i]

                    # extract kinetics SA
                    kinetics_mass_frac_sa = mass_frac_sensitivity_array[:, 0:num_ct_reactions]
                    sensitivity_array = np.zeros(len(self.sensitive_species) * len(self.model.reactions()))
                    for index, species in enumerate(self.sensitive_species):
                        for j in range(num_ct_reactions):
                            sensitivity_array[num_ct_reactions * index + j] = cantera_simulation.sensitivity(
                                species.to_chemkin(), j)

                            for i in range(len(kinetics_mass_frac_sa)):
                                if i not in inert_index_list:
                                    # massFracSensitivity for inerts are returned as 0.0 in Cantera, so we do not include them here
                                    sensitivity_array[num_ct_reactions * index + j] -= kinetics_mass_frac_sa[i][j]
                    kinetic_sensitivity_data.append(sensitivity_array)

                    # extract thermo SA if requested
                    if self.thermo_SA:
                        # extract thermo SA
                        thermo_mass_frac_sa = mass_frac_sensitivity_array[:, num_ct_reactions:]
                        sensitivity_array = np.zeros(len(self.sensitive_species) * num_ct_species)
                        for index, species in enumerate(self.sensitive_species):
                            for j in range(num_ct_species):
                                sensitivity_array[num_ct_species * index + j] = cantera_simulation.sensitivity(
                                    species.to_chemkin(), j + num_ct_reactions)

                                for i in range(len(mass_frac_sensitivity_array)):
                                    if i not in inert_index_list:
                                        # massFracSensitivity for inerts are returned as 0.0 in Cantera, so we must not include them here
                                        sensitivity_array[num_ct_species * index + j] -= thermo_mass_frac_sa[i][j]
                        thermo_sensitivity_data.append(sensitivity_array)

            # Convert species_data and sensitivity data to numpy arrays
            species_data = np.array(species_data)
            kinetic_sensitivity_data = np.array(kinetic_sensitivity_data)
            thermo_sensitivity_data = np.array(thermo_sensitivity_data)

            # Resave data into generic data objects
            time = GenericData(label='Time',
                               data=times,
                               units='s')
            temperature = GenericData(label='Temperature',
                                      data=temperature,
                                      units='K')
            pressure = GenericData(label='Pressure',
                                   data=pressure,
                                   units='Pa')
            condition_data = []
            condition_data.append(temperature)
            condition_data.append(pressure)

            for index, species in enumerate(self.species_list):
                # Create generic data object that saves the species object into the species object.  To allow easier manipulate later.
                species_generic_data = GenericData(label=species_names_list[index],
                                                   species=species,
                                                   data=species_data[:, index],
                                                   index=species.index
                                                   )
                condition_data.append(species_generic_data)

            # save kinetic data as generic data object
            reaction_sensitivity_data = []
            for index, species in enumerate(self.sensitive_species):
                for j in range(num_ct_reactions):
                    reaction_sensitivity_generic_data = GenericData(
                        label='dln[{0}]/dln[k{1}]: {2}'.format(species.to_chemkin(), j + 1, self.model.reactions()[j]),
                        species=species,
                        reaction=self.model.reactions()[j],
                        data=kinetic_sensitivity_data[:, num_ct_reactions * index + j],
                        index=j + 1,
                        )
                    reaction_sensitivity_data.append(reaction_sensitivity_generic_data)

            # save thermo data as generic data object
            thermodynamic_sensitivity_data = []
            if self.thermo_SA:
                for index, species in enumerate(self.sensitive_species):
                    for j in range(num_ct_species):
                        thermo_sensitivity_generic_data = GenericData(
                            label='dln[{0}]/dH[{1}]'.format(species, self.model.species()[j].name),
                            species=species,
                            data=thermo_sensitivity_data[:, num_ct_species * index + j],
                            index=j + 1,
                            )
                        thermodynamic_sensitivity_data.append(thermo_sensitivity_generic_data)

            all_data.append((time, condition_data, reaction_sensitivity_data, thermodynamic_sensitivity_data))

        return all_data


def get_rmg_species_from_user_species(user_list, rmg_list):
    """
    Args:
        user_list: list of generic Class Species Objects created by user
        species_list: a list of RMG species objects

    This function takes a list of generic species objects and returns the species object generated from a loaded RMG
    dictionary, thereby gaining the correct label for a given mechanism.

    Returns: A dict containing the Species Object from user_list and RMG Species objects as their values
    If the species is not found, the value will be returned as None
    """
    mapping = {}
    for user_species in user_list:
        user_species.generate_resonance_structures()

        for rmg_species in rmg_list:
            if user_species.is_isomorphic(rmg_species):
                if user_species in mapping:
                    raise KeyError("The Species with SMIlES {0} has appeared twice in the species list!".format(
                        user_species.molecule[0].to_smiles()))
                mapping[user_species] = rmg_species
                break
        else:
            mapping[user_species] = None

    return mapping


def find_ignition_delay(time, y_var=None, metric='maxDerivative'):
    """
    Identify the ignition delay point based on the following parameters:

    `time`: an array containing different times
    `yVars`: either a single y array or a list of arrays, typically containing only a single array such as pressure,
             but can contain multiple arrays such as species
    `metric`: can be set to 
        'maxDerivative': This is selected by default for y(t_ign) = max(dY/dt), and is typically used for a y_var containing T or P data
        'maxHalfConcentration': This is selected for the case where a metric like [OH](t_ign) = [OH]_max/2 is desired
        'maxSpeciesConcentrations': This is selected for the case where the metric for ignition
            is y1*y2*...*yn(t_ign) = max(y1*y2*...*yn) such as when the time desired if for max([CH][O]).  This is
            the only metric that requires a list of arrays

    Note that numpy array must be used.
    """

    if not isinstance(y_var, list):
        y_var = [y_var]

    for y in y_var:
        if len(y) != len(time):
            raise Exception('Mismatch of array length for time and y variable.')

    if metric == 'maxDerivative':
        if len(y_var) != 1:
            raise Exception('Maximum derivative metric for ignition delay must be used with a single y variable.')

        y = y_var[0]
        dydt = (y[1:] - y[:-1]) / (time[1:] - time[:-1])
        index = next(i for i, d in enumerate(dydt) if d == max(dydt))

        return 0.5 * (time[index] + time[index + 1])
    elif metric == 'maxHalfConcentration':
        if len(y_var) != 1:
            raise Exception('Max([OH]/2) metric for ignition delay must be used with a single y variable.')

        y = y_var[0]
        max_index = y.argmax()
        oh_metric = max(y) / 2
        min_data = oh_metric - y[0:max_index]
        index = min_data.argmin()
        return time[index]

    elif metric == 'maxSpeciesConcentrations':
        mult_data = np.ones(len(y_var[0]))
        for spec in y_var:
            mult_data *= spec
        index = mult_data.argmax()
        return time[index]


def check_nearly_equal(value1, value2, dE=1e-5):
    """
    Check that two values are nearly equivalent by abs(val1-val2) < abs(dE*val1)
    """

    if (abs(value1 - value2) <= abs(dE * value1) or
            abs(value1 - value2) <= abs(dE * value2) or
            abs(value1 - value2) <= dE):
        return True
    else:
        return False


def check_equivalent_cantera_species(ct_spec1, ct_spec2, dE=1e-5):
    """
    Checks that the two cantera species are nearly equivalent
    """
    try:
        assert ct_spec1.name == ct_spec2.name, "Identical name"
        assert ct_spec1.composition == ct_spec2.composition, "Identical composition"
        assert ct_spec1.size == ct_spec2.size, "Identical species size"
        assert ct_spec1.charge == ct_spec2.charge, "Identical charge"

        if ct_spec1.transport or ct_spec2.transport:
            trans1 = ct_spec1.transport
            trans2 = ct_spec2.transport

            assert check_nearly_equal(trans1.acentric_factor, trans2.acentric_factor, dE), "Identical acentric factor"
            assert check_nearly_equal(trans1.diameter, trans2.diameter, dE), "Identical diameter"
            assert check_nearly_equal(trans1.dipole, trans2.dipole, dE), "Identical dipole moment"
            assert trans1.geometry == trans2.geometry, "Identical geometry"
            assert check_nearly_equal(trans1.polarizability, trans2.polarizability, dE), "Identical polarizibility"
            assert check_nearly_equal(trans1.rotational_relaxation, trans2.rotational_relaxation,
                                      dE), "Identical rotational relaxation number"
            assert check_nearly_equal(trans1.well_depth, trans2.well_depth, dE), "Identical well depth"

        if ct_spec1.thermo or ct_spec2.thermo:
            thermo1 = ct_spec1.thermo
            thermo2 = ct_spec2.thermo

            Tlist = [300, 500, 1000, 1500, 2000]
            for T in Tlist:
                assert check_nearly_equal(thermo1.cp(T), thermo2.cp(T), dE), "Similar heat capacity"
                assert check_nearly_equal(thermo1.h(T), thermo2.h(T), dE), "Similar enthalpy"
                assert check_nearly_equal(thermo1.s(T), thermo2.s(T), dE), "Similar entropy"
    except Exception as e:
        print("Cantera species {0} failed equivalency check on: {1}".format(ct_spec1, e))
        return False

    return True


def check_equivalent_cantera_reaction(ct_rxn1, ct_rxn2, check_id=False, dE=1e-5):
    """
    Checks that the two cantera species are nearly equivalent
    if check_id is True, then ID's for the reactions will also be checked
    """

    def check_equivalent_arrhenius(arr1, arr2):
        assert check_nearly_equal(arr1.activation_energy, arr2.activation_energy, dE), "Similar Arrhenius Ea"
        assert check_nearly_equal(arr1.pre_exponential_factor, arr2.pre_exponential_factor,
                                  dE), "Similar Arrhenius A-factor"
        assert check_nearly_equal(arr1.temperature_exponent, arr2.temperature_exponent,
                                  dE), "Similar Arrhenius temperature exponent"

    def check_equivalent_falloff(fall1, fall2):
        assert len(fall1.parameters) == len(fall2.parameters), "Same number of falloff parameters"
        for i in range(len(fall1.parameters)):
            assert check_nearly_equal(fall1.parameters[i], fall2.parameters[i], dE), "Similar falloff parameters"
        assert fall1.type == fall2.type, "Same falloff parameterization type"

    try:
        assert type(ct_rxn1) == type(ct_rxn2), "Same Cantera reaction type"

        if isinstance(ct_rxn1, list):
            assert len(ct_rxn1) == len(ct_rxn2), "Same number of reactions"
            for i in range(len(ct_rxn1)):
                check_equivalent_cantera_reaction(ct_rxn1[i], ct_rxn2[i], check_id=check_id)

        if check_id:
            assert ct_rxn1.ID == ct_rxn2.ID, "Same reaction ID"

        assert ct_rxn1.duplicate == ct_rxn2.duplicate, "Same duplicate attribute"
        assert ct_rxn1.reversible == ct_rxn2.reversible, "Same reversible attribute"
        assert ct_rxn1.orders == ct_rxn2.orders, "Same orders attribute"
        assert ct_rxn1.allow_negative_orders == ct_rxn2.allow_negative_orders, "Same allow_negative_orders attribute"
        assert ct_rxn1.allow_nonreactant_orders == ct_rxn2.allow_nonreactant_orders, "Same allow_nonreactant_orders attribute"
        assert ct_rxn1.reactants == ct_rxn2.reactants, "Same reactants"
        assert ct_rxn1.products == ct_rxn2.products, "Same products"

        if isinstance(ct_rxn1, ct.ElementaryReaction):
            assert ct_rxn1.allow_negative_pre_exponential_factor == ct_rxn2.allow_negative_pre_exponential_factor, \
                "Same allow_negative_pre_exponential_factor attribute"
            if ct_rxn1.rate or ct_rxn2.rate:
                check_equivalent_arrhenius(ct_rxn1.rate, ct_rxn2.rate)

        elif isinstance(ct_rxn1, ct.PlogReaction):
            if ct_rxn1.rates or ct_rxn2.rates:
                assert len(ct_rxn1.rates) == len(ct_rxn2.rates), "Same number of rates in PLOG reaction"

                for i in range(len(ct_rxn1.rates)):
                    P1, arr1 = ct_rxn1.rates[i]
                    P2, arr2 = ct_rxn2.rates[i]
                    assert check_nearly_equal(P1, P2, dE), "Similar pressures for PLOG rates"
                    check_equivalent_arrhenius(arr1, arr2)

        elif isinstance(ct_rxn1, ct.ChebyshevReaction):
            assert ct_rxn1.Pmax == ct_rxn2.Pmax, "Same Pmax for Chebyshev reaction"
            assert ct_rxn1.Pmin == ct_rxn2.Pmin, "Same Pmin for Chebyshev reaction"
            assert ct_rxn1.Tmax == ct_rxn2.Tmax, "Same Tmax for Chebyshev reaction"
            assert ct_rxn1.Tmin == ct_rxn2.Tmin, "Same Tmin for Chebyshev reaction"
            assert ct_rxn1.nPressure == ct_rxn2.nPressure, "Same number of pressure interpolations"
            assert ct_rxn1.nTemperature == ct_rxn2.nTemperature, "Same number of temperature interpolations"
            for i in range(ct_rxn1.coeffs.shape[0]):
                for j in range(ct_rxn1.coeffs.shape[1]):
                    assert check_nearly_equal(ct_rxn1.coeffs[i, j], ct_rxn2.coeffs[i, j], dE), \
                        "Similar Chebyshev coefficients"

        elif isinstance(ct_rxn1, ct.ThreeBodyReaction):
            assert ct_rxn1.default_efficiency == ct_rxn2.default_efficiency, "Same default efficiency"
            assert ct_rxn1.efficiencies == ct_rxn2.efficiencies, "Same efficienciess"

        elif isinstance(ct_rxn1, ct.FalloffReaction):
            assert ct_rxn1.default_efficiency == ct_rxn2.default_efficiency, "Same default efficiency"
            assert ct_rxn1.efficiencies == ct_rxn2.efficiencies, "Same efficienciess"
            if ct_rxn1.falloff or ct_rxn2.falloff:
                check_equivalent_falloff(ct_rxn1.falloff, ct_rxn2.falloff)
            if ct_rxn1.high_rate or ct_rxn2.high_rate:
                check_equivalent_arrhenius(ct_rxn1.high_rate, ct_rxn2.high_rate)
            if ct_rxn1.low_rate or ct_rxn2.low_rate:
                check_equivalent_arrhenius(ct_rxn1.low_rate, ct_rxn2.low_rate)

    except Exception as e:
        print("Cantera reaction {0} failed equivalency check on: {1}".format(ct_rxn1, e))
        return False

    return True
