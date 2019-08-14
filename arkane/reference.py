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

"""
This module defines the ReferenceSpecies class, which are used in isodesmic reaction calculations

"""

from __future__ import print_function, division

import logging
import os

import yaml

from arkane.common import ArkaneSpecies, ARKANE_CLASS_DICT
from arkane.isodesmic import SpeciesConstraints,ErrorCancelingSpecies,ErrorCancelingReaction,ErrorCancelingScheme
from rmgpy.molecule import Molecule
from rmgpy.rmgobject import RMGObject
from rmgpy.species import Species
from rmgpy.statmech import Conformer
from rmgpy.thermo import ThermoData
from rmgpy.thermo.model import HeatCapacityModel


class ReferenceSpecies(ArkaneSpecies):
    """
    A class for storing high level reference data and quantum chemistry calculations for a variety of model chemistry
    selections for use in isodesmic reaction calculations
    """

    def __init__(self, species=None, smiles=None, adjacency_list=None, inchi=None, reference_data=None,
                 calculated_data=None, preferred_reference=None, index=None, label=None, cas_number=None,
                 symmetry_number=None, **kwargs):
        """
        One of the following must be provided: species, smiles, adjacency_list, inchi.

        Args:
            species (rmgpy.molecule.Species): Molecule object representing the reference species
            smiles (str): SMILES string representing the reference species
            adjacency_list (str): An RMG adjacency list representation of the reference species
            inchi (str): InChI string representing the reference species
            reference_data (dict): Formatted as {'source_string': ReferenceDataEntry, ...}
            calculated_data (dict): Formatted as {'model_chemistry': CalculatedDataEntry, ...}
            preferred_reference (str): The source string key for the reference data to use for isodesmic reactions
            index (int): Index of this species in the database of reference species located at
                `RMG-database/input/reference_sets/`
            label (str): A user defined label for easily identifying the species
            cas_number (str): CAS number associated with the reference species
            symmetry_number (int): The true symmetry number of the species (if not provided will default to the number
                calculated by RMG)
            **kwargs: Arguments passed to the parent ArkaneSpecies class when loading from a YAML file. Not intended for
                user input
        """

        if species is None:
            if smiles:
                species = Species(SMILES=smiles)
            elif inchi:
                species = Species(InChI=inchi)
            elif adjacency_list:
                species = Species().fromAdjacencyList(adjacency_list)
            else:
                raise ValueError('Either an rmgpy species object, smiles string, InChI string, or an adjacency list '
                                 'must be given to create a ReferenceSpecies object')

        super(ReferenceSpecies, self).__init__(species=species, label=label, **kwargs)

        self._reference_data = None
        self._calculated_data = None
        self.reference_data = reference_data
        self.calculated_data = calculated_data
        self.index = index
        self.cas_number = cas_number
        self.preferred_reference = preferred_reference

        # Alter the symmetry number calculated by RMG to the one provided by the user
        if symmetry_number:
            self.symmetry_number = symmetry_number

    def __repr__(self):
        if self.index:
            label = '{0}({1})'.format(self.smiles, self.index)
        else:
            label = '{0}'.format(self.smiles)

        return '<ReferenceSpecies {0}>'.format(label)

    @property
    def reference_data(self):
        return self._reference_data

    @reference_data.setter
    def reference_data(self, value):
        if not value:
            self._reference_data = {}
        elif isinstance(value, dict):
            if all(isinstance(source, str) for source in value.keys()):
                if all(isinstance(data_entry, ReferenceDataEntry) for data_entry in value.values()):
                    self._reference_data = value
        else:
            raise ValueError('Reference data must be given as a dictionary of the data source (string) and associated '
                             'ReferenceDataEntry object')

    @property
    def calculated_data(self):
        return self._calculated_data

    @calculated_data.setter
    def calculated_data(self, value):
        if not value:
            self._calculated_data = {}
        elif isinstance(value, dict):
            if all(isinstance(source, str) for source in value.keys()):
                if all(isinstance(data_entry, CalculatedDataEntry) for data_entry in value.values()):
                    self._calculated_data = value
        else:
            raise ValueError('Calculated data must be given as a dictionary of the model chemistry (string) and '
                             'associated CalculatedDataEntry object')

    def load_yaml(self, path, label=None, pdep=False):
        """
        Load a ReferenceSpecies object from a YAML file.

        Args:
            path (str): Location on disk of the YAML file
            label: Unused argument from parent class ArkaneSpecies
            pdep: Unused argument from parent class ArkaneSpecies
        """
        with open(path, 'r') as f:
            data = yaml.safe_load(stream=f)

        if data['class'] != 'ReferenceSpecies':
            raise ValueError('Cannot create ReferenceSpecies object from yaml file {0}: object defined by this file is'
                             'not a ReferenceSpecies object'.format(path))

        data = {key: data[key] for key in data.keys() if key != 'class'}
        class_dict = ARKANE_CLASS_DICT
        class_dict['ReferenceDataEntry'] = ReferenceDataEntry
        class_dict['CalculatedDataEntry'] = CalculatedDataEntry

        self.make_object(data, class_dict)

    def update_from_arkane_spcs(self, arkane_species):
        """
        Add in calculated data from an existing ArkaneSpecies object.

        Notes:
            If the model chemistry already exists then this calculated data will be overwritten by the data contained
            in arkane_species

        Args:
            arkane_species (ArkaneSpecies):  Matching Arkane species that was run at the desired model chemistry
        """
        conformer = arkane_species.conformer
        thermo = arkane_species.thermo
        thermo_data = arkane_species.thermo_data
        calc_data = CalculatedDataEntry(conformer, thermo, thermo_data,)
        self.calculated_data[arkane_species.model_chemistry] = calc_data

    def to_error_canceling_spcs(self, model_chemistry, source=None):
        """
        Extract calculated and reference data from a specified model chemistry and source and return as a new
        ErrorCancelingSpecies object

        Notes:
            If no source is given, the preferred source for this species. If the `preferred_source` attribute is not set
            then the preferred source is taken as the source with the lowest non-zero uncertainty

        Args:
            model_chemistry (str): Model chemistry (level of theory) to use as the low level data
            source (str): Reference data source to take the high level data from

        Raises:
            KeyError: If `model_chemistry` is not available for this reference species
            ValueError: If there is no reference data for this reference species

        Returns:
            ErrorCancelingSpecies
        """
        if model_chemistry not in self.calculated_data:
            raise KeyError('Model chemistry `{0}` not available for species {1}'.format(model_chemistry, self))
        if not self.reference_data:
            raise ValueError('No reference data is included for species {0}'.format(self))

        molecule = Molecule(SMILES=self.smiles)
        preferred_source = source

        if not preferred_source:
            # Find the preferred source
            if self.preferred_reference is not None:
                preferred_source = self.preferred_reference
            else:  # Choose the source that has the smallest uncertainty
                sources = self.reference_data.keys()
                data = self.reference_data.values()
                preferred_source = sources[0]  # If all else fails, use the first source as the preferred one
                uncertainty = data[0].thermo_data.H298.uncertainty_si
                for i, entry in enumerate(data):
                    if (entry.thermo_data.H298.uncertainty_si > 0) and \
                            (entry.thermo_data.H298.uncertainty_si < uncertainty):
                        uncertainty = entry.thermo_data.H298.uncertainty_si
                        preferred_source = sources[i]
        high_level_h298 = self.reference_data[preferred_source].thermo_data.H298.__reduce__()[1]
        low_level_h298 = self.calculated_data[model_chemistry].thermo_data.H298.__reduce__()[1]
        if 'tpss/def2-tzvp' in self.calculated_data.keys():
            try:
                fod = float(self.calculated_data['tpss/def2-tzvp'].fod)
            except:
                fod = None
        else:
            fod = None


        return ErrorCancelingSpecies(molecule, low_level_h298, model_chemistry, high_level_h298, fod, preferred_source)


class ReferenceDataEntry(RMGObject):
    """
    A class for storing reference data for a specific species from a single source
    """
    def __init__(self, thermo_data, atct_id=None):
        """

        Args:
            thermo_data (rmgpy.thermo.ThermoData): Thermochemistry (Hf298, Cp, ...) from the reference for a species
            atct_id (str): ID number in the Active Thermochemical Tables if the source is ATcT
        """
        super(ReferenceDataEntry, self).__init__()
        self._thermo_data = None
        self.thermo_data = thermo_data
        self.atct_id = atct_id

    def __repr__(self):
        return str(self.as_dict())

    @property
    def thermo_data(self):
        return self._thermo_data

    @thermo_data.setter
    def thermo_data(self, value):
        if value:
            if isinstance(value, ThermoData):
                self._thermo_data = value
            else:
                raise ValueError('thermo_data for a ReferenceDataEntry object must be an rmgpy ThermoData instance')
        else:
            self._thermo_data = None


class CalculatedDataEntry(RMGObject):
    """
    A class for storing a single entry of statistical mechanical and thermochemistry information calculated at a single
    model chemistry or level of theory
    """
    def __init__(self, conformer=None, thermo=None, thermo_data=None, t1_diagnostic=None, fod=None):
        """

        Args:
            conformer (rmgpy.statmech.Conformer): Conformer object generated from an Arkane job. Stores many peices of
                information gained from quantum chemistry calculations, including coordinates, frequencies etc.
            thermo (HeatCapacityModel): `NASA` or `Wilhoit` thermo object to store the fitted polynomials
            thermo_data (rmgpy.thermo.ThermoData): Actual thermochemistry values calculated using statistical mechanics
                at select points. Arkane fits a heat capacity model to this data
            t1_diagnostic (float): T1 diagnostic for coupled cluster calculations to check if single reference methods
                are suitable
            fod (float): Fractional Occupation number weighted electron Density
        """
        super(CalculatedDataEntry, self).__init__()
        self._conformer = None
        self._thermo = None
        self._thermo_data = None
        self.conformer = conformer
        self.thermo = thermo
        self.thermo_data = thermo_data
        self.t1_diagnostic = t1_diagnostic
        self.fod = fod

    def __repr__(self):
        return str(self.as_dict())

    @property
    def conformer(self):
        return self._conformer

    @conformer.setter
    def conformer(self, value):
        if value:
            if isinstance(value, Conformer):
                self._conformer = value
            else:
                raise ValueError('conformer for a CalculatedDataEntry object must be an rmgpy Conformer instance')
        else:
            self._conformer = None

    @property
    def thermo(self):
        return self._thermo

    @thermo.setter
    def thermo(self, value):
        if value:
            if issubclass(type(value), HeatCapacityModel):
                self._thermo = value
            else:
                raise ValueError('thermo for a CalculatedDataEntry object must be an object of a subclass of'
                                 'an rmgpy HeatCapacityModel class')
        else:
            self._thermo = None

    @property
    def thermo_data(self):
        return self._thermo_data

    @thermo_data.setter
    def thermo_data(self, value):
        if value:
            if isinstance(value, ThermoData):
                self._thermo_data = value
            else:
                raise ValueError('thermo_data for a CalculatedDataEntry object must be an rmgpy ThermoData object')


class ReferenceDatabase(object):
    """
    A class for loading and working with database of reference species, located at RMG-database/input/reference_sets/
    """
    def __init__(self,paths = None, parse_model_chemistries=True):
        """
        Attributes:
            self.reference_sets (Dict[str, ReferenceSpecies]): {'set name': [ReferenceSpecies, ...], ...}
        """
        
        self.SpeciesConstraints = {}
        self.reference_sets = {}
        self.paths = {}
        self.model_chemistries = {}
        self.errorCancellingSets = {}
        self.descriptors = {}

        if paths:
            self.load(paths=paths,parse_model_chemistries=parse_model_chemistries)

    def load(self, paths='', parse_model_chemistries=False):
        """
        Load one or more set of reference species and append it on to the database

        Args:
            paths (Union[list, str]): A single path string, or a list of path strings pointing to a set of reference
                species to be loaded into the database. The string should point to the folder that has the name of the
                reference set. The name of sub-folders in a reference set directory should be indices starting from 0
                and should contain a YAML file that defines the ReferenceSpecies object of that index, named {index}.yml
        """
        if not paths:  # Default to the main reference set in RMG-database
            file_dir = os.path.dirname(os.path.abspath(__file__))
            rmg_database_dir = os.path.join(os.path.dirname(os.path.dirname(file_dir)), 'RMG-database')
            paths = [os.path.join(rmg_database_dir, 'input/reference_sets/main')]

        if isinstance(paths, str):  # Convert to a list with one element
            paths = [paths]

        molecule_list = []
        for path in paths:
            model_chemistries = []
            set_name = os.path.basename(path)
            logging.info('Loading in reference set `{0}` from {1} ...'.format(set_name, path))
            spcs_dirs = os.listdir(path)
            reference_set = []
            for spcs in spcs_dirs:
                ref_spcs = ReferenceSpecies.__new__(ReferenceSpecies)
                ref_spcs.load_yaml(os.path.join(path, spcs, '{0}.yml'.format(spcs)))
                for model_chem in ref_spcs.calculated_data.keys():
                    if model_chem not in model_chemistries:
                        model_chemistries.append(model_chem)
                molecule = Molecule(SMILES=ref_spcs.smiles)
                if (len(ref_spcs.calculated_data) == 0) or (len(ref_spcs.reference_data) == 0):
                    logging.warning('Molecule {0} from reference set `{1}` does not have any reference data and/or '
                                    'calculated data. This entry will not be added'.format(ref_spcs.smiles, set_name))
                    continue
                # perform isomorphism checks to prevent duplicate species
                for mol in molecule_list:
                    if molecule.isIsomorphic(mol):
                        logging.warning('Molecule {0} from reference set `{1}` already exists in the reference '
                                        'database. The entry from this reference set will not '
                                        'be added'.format(ref_spcs.smiles, set_name))
                        break
                else:
                    molecule_list.append(molecule)
                    reference_set.append(ref_spcs)
            
            self.model_chemistries[set_name] = model_chemistries
            self.paths[set_name] = path
            self.reference_sets[set_name] = reference_set

        if parse_model_chemistries:
            for set_name,chemistries in self.model_chemistries.items():
                self.errorCancellingSets[set_name] = {}
                for chem in chemistries:
                    errorCancellingSet = self.extract_model_chemistry(chem, [set_name])
                    self.errorCancellingSets[set_name][chem] = errorCancellingSet

    def save(self,path,set_name):

        ref_spcs = self.reference_sets[set_name]

        for ref in ref_spcs:
            ref.save_yaml(path)

    def extract_model_chemistry(self, model_chemistry, sets=None):
        """
        Return a list of ErrorCancelingSpecies objects from the reference species in the database that have entries for
        the requested model chemistry

        Args:
            model_chemistry (str): String that describes the level of chemistry used to calculate the low level data
            sets (list): A list of the names of the reference sets to include (all sets in the database will be used if
                not specified or `None`)

        Returns:
            List[ErrorCancelingSpecies]
        """
        reference_list = []

        if sets is None:  # Load in all of the sets
            sets = self.reference_sets.keys()

        for set_name in sets:
            current_set = self.reference_sets[set_name]
            for ref_spcs in current_set:
                if model_chemistry not in ref_spcs.calculated_data:  # Move on to the next reference species
                    continue
                if not ref_spcs.calculated_data[model_chemistry].thermo: # Make sure refernce species has thermo
                    continue
                if not ref_spcs.reference_data:  # This reference species does not have any sources, continue on
                    continue
                reference_list.append(ref_spcs.to_error_canceling_spcs(model_chemistry))

        return reference_list

    def get_constraint_map(self, model_chemistry, sets=None):
  
        if not sets:
            sets = self.reference_sets.keys()
        else:
            if isinstance(sets,str):
                sets = [sets]
            assert(isinstance(sets,list))
        
        reference_list = []

        for s in sets: 
            if model_chemistry in self.errorCancellingSets[s].keys():
                for spcs in self.errorCancellingSets[s][model_chemistry]:
                    reference_list.append(spcs)
            else:
                for spcs in self.extract_model_chemistry(model_chemistry):
                    reference_list.append(spcs)
  
        constraint = SpeciesConstraints(target=None, reference_list=reference_list, 
        constraint_class='all',conserve_bonds=True, conserve_ring_size=True)
        
        self.SpeciesConstraints[model_chemistry] = constraint

        return constraint.constraint_map

    def get_descriptors(self,species,sets=None):

        constraint = SpeciesConstraints()
        descriptors = constraint.get_descriptors(species)
        self.descriptors[species] = descriptors

        return descriptors

    def test(self, model_chemistry, constraint_class = None, iterate_constraint_classes = True, number_of_reactions=5, sets=None):

        import pandas as pd
        
        reference_list = []

        if not sets:
            sets = self.errorCancellingSets.keys()
        
        for s in sets:
            if s in self.errorCancellingSets.keys():
                if model_chemistry in self.errorCancellingSets[s].keys():
                    for spcs in self.errorCancellingSets[s][model_chemistry]:
                        reference_list.append(spcs)
            else:
                for spcs in self.extract_model_chemistry(model_chemistry, sets=[s]):
                    reference_list.append(spcs)
        
        data = []
        for i,error_canceling_spcs in enumerate(reference_list):
            print('calculating thermo for {}, {} of {}'.format(error_canceling_spcs.molecule.toSMILES(),i+1,len(reference_list)))
            reference_set = reference_list[:]
            target_spcs = error_canceling_spcs
            reference_set.remove(target_spcs)

            high_level_hf298 = target_spcs.high_level_hf298
            ref_h298 = high_level_hf298.value_si/high_level_hf298.conversionFactors['kcal/mol']
            ref_h298_uncertainty = high_level_hf298.uncertainty_si/high_level_hf298.conversionFactors['kcal/mol']
            fod = target_spcs.fod
            
            for constraint in ['class_0','class_1','class_2','class_3','class_4','class_5']:
                if not iterate_constraint_classes:
                    if constraint != constraint_class:
                        continue
                isodesmic_scheme = ErrorCancelingScheme(target=target_spcs,reference_set=reference_set,constraint_class=constraint,
                conserve_bonds=True,conserve_ring_size=True)
                h298_mean, reaction_list, fod_dict = isodesmic_scheme.calculate_target_enthalpy(n_reactions_max=5, milp_software='lpsolve')
                h298 = h298_mean.value_si/h298_mean.conversionFactors['kcal/mol']
                species_data = [target_spcs.molecule.toSMILES(),constraint,reaction_list,fod_dict,fod,h298,ref_h298,ref_h298_uncertainty,h298-ref_h298]
                data.append(species_data)

        columns = ['SMILES','constraint_class','Reactions','fod_dict','fod','H298_mean(kcal/mol)','H298_ref(kcal/mol)','uncertainty','calculated-ref']
        df = pd.DataFrame(data,columns=columns)

        return df



if __name__ == '__main__':
    pass
