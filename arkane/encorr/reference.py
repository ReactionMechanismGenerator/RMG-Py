#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import string

import yaml

from arkane.common import ArkaneSpecies, ARKANE_CLASS_DICT
from rmgpy.rmgobject import RMGObject
from rmgpy.species import Species


class ReferenceSpecies(ArkaneSpecies):
    """
    A class for storing high level reference data and quantum chemistry calculations for a variety of model chemistry
    selections for use in isodesmic reaction calculations
    """

    def __init__(self, species=None, smiles=None, adjacency_list=None, inchi=None, reference_data=None,
                 calculated_data=None, preferred_reference=None, index=None, label=None, cas_number=None,
                 symmetry_number=None, default_xyz_chemistry=None, **kwargs):
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
                ``RMG-database/input/reference_sets/``
            label (str): A user defined label for easily identifying the species
            cas_number (str): CAS number associated with the reference species
            symmetry_number (int): The true symmetry number of the species (if not provided will default to the number
                calculated by RMG)
            default_xyz_chemistry (str): The model chemistry that should be used to get the default XYZ coordinates for
                this species.
            **kwargs: Arguments passed to the parent ArkaneSpecies class when loading from a YAML file. Not intended for
                user input
        """

        if species is None:
            if adjacency_list:
                species = Species().from_adjacency_list(adjacency_list, raise_atomtype_exception=False,
                                                        raise_charge_exception=False)
            elif smiles:
                species = Species(smiles=smiles)
            elif inchi:
                species = Species(inchi=inchi)
            else:
                raise ValueError('Either an rmgpy species object, smiles string, InChI string, or an adjacency list '
                                 'must be given to create a ReferenceSpecies object')

        super().__init__(species=species, label=label, **kwargs)

        self.reference_data = reference_data
        self.calculated_data = calculated_data
        self.index = index
        self.cas_number = cas_number
        self.preferred_reference = preferred_reference
        self.default_xyz_chemistry = default_xyz_chemistry

        # Alter the symmetry number calculated by RMG to the one provided by the user
        if symmetry_number:
            self.symmetry_number = symmetry_number

    def __repr__(self):
        if self.index:
            label = f'{self.smiles}({self.index})'
        else:
            label = f'{self.smiles}'

        return f'<ReferenceSpecies {label}>'

    @property
    def reference_data(self):
        return self._reference_data

    @reference_data.setter
    def reference_data(self, value):
        if not value:
            self._reference_data = {}
        elif isinstance(value, dict) and _is_valid_reference_data(value):
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
        elif isinstance(value, dict) and _is_valid_calculated_data(value):
            self._calculated_data = value
        else:
            raise ValueError('Calculated data must be given as a dictionary of the model chemistry (string) and '
                             'associated CalculatedDataEntry object')

    def as_dict(self):
        dictionary = super().as_dict()

        # Don't include "is_ts" attribute, as no reference species are transition states
        del dictionary['is_ts']
        return dictionary

    def save_yaml(self, path=MAIN_REFERENCE_PATH):
        """
        Save the reference species to a .yml file
        """
        if not os.path.exists(os.path.join(os.path.abspath(path), '')):
            os.mkdir(os.path.join(os.path.abspath(path), ''))
        valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
        filename = os.path.join(''.join(c for c in self.label if c in valid_chars) + '.yml')
        full_path = os.path.join(path, filename)
        with open(full_path, 'w') as f:
            yaml.dump(data=self.as_dict(), stream=f)
        logging.debug(f'Dumping species {self.label} data as {filename}')

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
            raise ValueError(f'Cannot create ReferenceSpecies object from yaml file {path}: object defined by this '
                             f'file is not a ReferenceSpecies object')

        data = {key: data[key] for key in data.keys() if key != 'class'}
        class_dict = ARKANE_CLASS_DICT
        class_dict['ReferenceDataEntry'] = ReferenceDataEntry
        class_dict['CalculatedDataEntry'] = CalculatedDataEntry

        self.make_object(data, class_dict)


def _is_valid_reference_data(data_dictionary):
    """
    Determine if the given reference_data dictionary is supplied in a valid format
    Args:
        data_dictionary (dict): reference_data dictionary

    Returns:
        bool
    """
    if all(isinstance(source, str) for source in data_dictionary.keys()):
        if all(isinstance(data_entry, ReferenceDataEntry) for data_entry in data_dictionary.values()):
            return True
    return False


def _is_valid_calculated_data(data_dictionary):
    """
    Determine if the given calculated_data dictionary is supplied in a valid format
    Args:
        data_dictionary (dict): calculated_data dictionary

    Returns:
        bool
    """
    if all(isinstance(source, str) for source in data_dictionary.keys()):
        if all(isinstance(data_entry, CalculatedDataEntry) for data_entry in data_dictionary.values()):
            return True
    return False


if __name__ == '__main__':
    pass
