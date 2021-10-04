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

"""
Arkane common module
"""

import logging
import os.path
import shutil
import string
import time
from typing import List, Union

import numpy as np
import yaml

import rmgpy.constants as constants
from rmgpy import __version__
from rmgpy.exceptions import InputError
from rmgpy.molecule.element import get_element
from rmgpy.molecule.translator import to_inchi, to_inchi_key
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.quantity import ScalarQuantity, ArrayQuantity
from rmgpy.rmgobject import RMGObject
from rmgpy.species import Species, TransitionState
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.thermo import NASA, Wilhoit, ThermoData, NASAPolynomial
from rmgpy.transport import TransportData

from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory
from arkane.pdep import PressureDependenceJob

################################################################################

# Class dictionary for recreating objects from YAML files. This is needed elsewhere, so store as a module level variable
ARKANE_CLASS_DICT = {'ScalarQuantity': ScalarQuantity,
                     'ArrayQuantity': ArrayQuantity,
                     'Conformer': Conformer,
                     'LinearRotor': LinearRotor,
                     'NonlinearRotor': NonlinearRotor,
                     'KRotor': KRotor,
                     'SphericalTopRotor': SphericalTopRotor,
                     'HinderedRotor': HinderedRotor,
                     'FreeRotor': FreeRotor,
                     'IdealGasTranslation': IdealGasTranslation,
                     'HarmonicOscillator': HarmonicOscillator,
                     'TransportData': TransportData,
                     'SingleExponentialDown': SingleExponentialDown,
                     'Wilhoit': Wilhoit,
                     'NASA': NASA,
                     'NASAPolynomial': NASAPolynomial,
                     'ThermoData': ThermoData,
                     'np_array': np.array,
                     'LevelOfTheory': LevelOfTheory,
                     'CompositeLevelOfTheory': CompositeLevelOfTheory
                     }


# Add a custom string representer to use block literals for multiline strings
def str_repr(dumper, data):
    """
    Repair YAML string representation
    """
    if len(data.splitlines()) > 1:
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)


yaml.add_representer(str, str_repr)


class ArkaneSpecies(RMGObject):
    """
    A class for archiving an Arkane species including its statmech data into .yml files
    """

    def __init__(self, species=None, conformer=None, author='', level_of_theory='', model_chemistry='',
                 frequency_scale_factor=None, use_hindered_rotors=None, use_bond_corrections=None, atom_energies='',
                 chemkin_thermo_string='', smiles=None, adjacency_list=None, inchi=None, inchi_key=None, xyz=None,
                 molecular_weight=None, symmetry_number=None, transport_data=None, energy_transfer_model=None,
                 thermo=None, thermo_data=None, label=None, datetime=None, RMG_version=None, reactants=None,
                 products=None, reaction_label=None, is_ts=None, charge=None, formula=None, multiplicity=None):
        # reactants/products/reaction_label need to be in the init() to avoid error when loading a TS YAML file,
        # but we don't use them
        super(ArkaneSpecies, self).__init__()
        if species is None and conformer is None:
            # Expecting to get a species or a TS when generating the object within Arkane,
            # or a conformer when parsing from YAML.
            raise ValueError('No species (or TS) or conformer was passed to the ArkaneSpecies object')
        if conformer is not None:
            self.conformer = conformer
        if label is None and species is not None:
            self.label = species.label
        else:
            self.label = label
        self.author = author
        self.level_of_theory = level_of_theory
        self.model_chemistry = model_chemistry
        self.frequency_scale_factor = frequency_scale_factor
        self.use_hindered_rotors = use_hindered_rotors
        self.use_bond_corrections = use_bond_corrections
        self.atom_energies = atom_energies
        self.xyz = xyz
        self.molecular_weight = molecular_weight
        self.symmetry_number = symmetry_number
        self.charge = charge
        self.multiplicity = multiplicity
        self.is_ts = is_ts if is_ts is not None else isinstance(species, TransitionState)
        if not self.is_ts:
            self.chemkin_thermo_string = chemkin_thermo_string
            self.smiles = smiles
            self.adjacency_list = adjacency_list
            self.inchi = inchi
            self.inchi_key = inchi_key
            self.transport_data = transport_data
            self.energy_transfer_model = energy_transfer_model
            self.thermo = thermo
            self.thermo_data = thermo_data
            self.formula = formula
        else:
            # initialize TS-related attributes
            self.imaginary_frequency = None
            self.reaction_label = ''
            self.reactants = list()
            self.products = list()
        if species is not None:
            self.update_species_attributes(species)
        self.RMG_version = RMG_version if RMG_version is not None else __version__
        self.datetime = datetime if datetime is not None else time.strftime("%Y-%m-%d %H:%M")

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the object
        """
        result = '{0!r}'.format(self.__class__.__name__)
        result += '{'
        for key, value in self.as_dict().items():
            if key != 'class':
                result += '{0!r}: {1!r}'.format(str(key), str(value))
        result += '}'
        return result

    def update_species_attributes(self, species=None):
        """
        Update the object with a new species/TS (while keeping non-species-dependent attributes unchanged)
        """
        if species is None:
            raise ValueError('No species was passed to ArkaneSpecies')
        # Don't overwrite the label if it already exists
        self.label = self.label or species.label
        if isinstance(species, TransitionState):
            self.imaginary_frequency = species.frequency
            if species.conformer is not None:
                self.conformer = species.conformer
                self.xyz = self.update_xyz_string()
        elif species.molecule is not None and len(species.molecule) > 0:
            self.smiles = species.molecule[0].to_smiles()
            self.adjacency_list = species.molecule[0].to_adjacency_list()
            self.charge = species.molecule[0].get_net_charge()
            self.multiplicity = species.molecule[0].multiplicity
            self.formula = species.molecule[0].get_formula()
            try:
                inchi = to_inchi(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi = ''
            try:
                inchi_key = to_inchi_key(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi_key = ''
            self.inchi = inchi
            self.inchi_key = inchi_key
            if species.conformer is not None:
                self.conformer = species.conformer
                self.xyz = self.update_xyz_string()
            self.molecular_weight = species.molecular_weight
            if species.symmetry_number != -1:
                self.symmetry_number = species.symmetry_number
            if species.transport_data is not None:
                self.transport_data = species.transport_data  # called `collisionModel` in Arkane
            if species.energy_transfer_model is not None:
                self.energy_transfer_model = species.energy_transfer_model
            if species.thermo is not None:
                self.thermo = species.thermo.as_dict()
                data = species.get_thermo_data()
                h298 = data.get_enthalpy(298) / 4184.
                s298 = data.get_entropy(298) / 4.184
                temperatures = np.array([300, 400, 500, 600, 800, 1000, 1500, 2000, 2400])
                cp = []
                for t in temperatures:
                    cp.append(data.get_heat_capacity(t) / 4.184)

                self.thermo_data = ThermoData(H298=(h298, 'kcal/mol'),
                                              S298=(s298, 'cal/(mol*K)'),
                                              Tdata=(temperatures, 'K'),
                                              Cpdata=(cp, 'cal/(mol*K)'),
                                              )

    def update_xyz_string(self):
        """
        Generate an xyz string built from self.conformer, and standardize the result

        Returns:
            str: 3D coordinates in an XYZ format.
        """
        xyz_list = list()
        if self.conformer is not None and self.conformer.number is not None:
            # generate the xyz-format string from self.conformer.coordinates and self.conformer.number
            xyz_list.append(str(len(self.conformer.number.value_si)))
            xyz_list.append(self.label)
            for number, coordinate in zip(self.conformer.number.value_si, self.conformer.coordinates.value_si):
                element_symbol = get_element(int(number)).symbol
                row = '{0:4}'.format(element_symbol)
                row += '{0:14.8f}{1:14.8f}{2:14.8f}'.format(*(coordinate * 1e10).tolist())  # convert m to Angstrom
                xyz_list.append(row)
        return '\n'.join(xyz_list)

    def save_yaml(self, path):
        """
        Save the species with all statMech data to a .yml file
        """
        if not os.path.exists(os.path.join(os.path.abspath(path), 'species', '')):
            os.mkdir(os.path.join(os.path.abspath(path), 'species', ''))
        valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
        filename = os.path.join('species',
                                ''.join(c for c in self.label if c in valid_chars) + '.yml')
        full_path = os.path.join(path, filename)
        with open(full_path, 'w') as f:
            yaml.dump(data=self.as_dict(), stream=f)
        logging.debug('Dumping species {0} data as {1}'.format(self.label, filename))

    def load_yaml(self, path, label=None, pdep=False):
        """
        Load the all statMech data from the .yml file in `path` into `species`
        `pdep` is a boolean specifying whether or not job_list includes a pressureDependentJob.
        """
        yml_file = os.path.basename(path)
        if label:
            logging.info('Loading statistical mechanics parameters for {0} from {1} file...'.format(label, yml_file))
        else:
            logging.info('Loading statistical mechanics parameters from {0} file...'.format(yml_file))
        with open(path, 'r') as f:
            content = f.read()
        content = replace_yaml_syntax(content, label)
        data = yaml.safe_load(stream=content)
        if label:
            # First, warn the user if the label doesn't match
            try:
                if label != data['label']:
                    logging.debug('Found different labels for species: {0} in input file, and {1} in the .yml file. '
                                  'Using the label "{0}" for this species.'.format(label, data['label']))
            except KeyError:
                # Lacking label in the YAML file is strange, but accepted
                logging.debug('Did not find label for species {0} in .yml file.'.format(label))

            # Then, set the ArkaneSpecies label to the user supplied label
            data['label'] = label
        try:
            class_name = data['class']
        except KeyError:
            raise KeyError("Can only make objects if the `class` attribute in the dictionary is known")
        if class_name != 'ArkaneSpecies':
            raise KeyError("Expected a ArkaneSpecies object, but got {0}".format(class_name))
        del data['class']
        freq_data = None
        if 'imaginary_frequency' in data:
            freq_data = data['imaginary_frequency']
            del data['imaginary_frequency']
        if not data['is_ts']:
            if 'smiles' in data:
                data['species'] = Species(smiles=data['smiles'])
            elif 'adjacency_list' in data:
                data['species'] = Species().from_adjacency_list(data['adjacency_list'])
            elif 'inchi' in data:
                data['species'] = Species(inchi=data['inchi'])
            else:
                raise ValueError('Cannot load ArkaneSpecies from YAML file {0}. Either `smiles`, `adjacency_list`, or '
                                 'InChI must be specified'.format(path))
            # Finally, set the species label so that the special attributes are updated properly
            data['species'].label = data['label']

        self.make_object(data=data, class_dict=ARKANE_CLASS_DICT)
        if freq_data is not None:
            self.imaginary_frequency = ScalarQuantity()
            self.imaginary_frequency.make_object(data=freq_data, class_dict=ARKANE_CLASS_DICT)

        if pdep and not self.is_ts and self.smiles is None and self.adjacency_list is None \
                and self.inchi is None and self.molecular_weight is None:
            raise ValueError('The molecular weight was not specified, and a structure was not given so it could '
                             'not be calculated. Specify either the molecular weight or structure if '
                             'pressure-dependent calculations are requested. Check file {0}'.format(path))
        logging.debug("Parsed all YAML objects")


def replace_yaml_syntax(content, label=None):
    """
    PEP8 compliant changes to RMG objects could be backward incompatible with Arkane's YAML files.
    Search for knows phrases which were replace, and fix the format on the fly.

    Args:
        content (str): The content of an Arkane YAML file.

    Returns:
        str: The modified content to be processed via yaml.safe_load().
    """
    syntax_correction_dict = {'spinMultiplicity': 'spin_multiplicity',
                              'opticalIsomers': 'optical_isomers',
                              }
    replaced_keys = list()
    for key, value in syntax_correction_dict.items():
        if key in content:
            content = content.replace(key, value)
            replaced_keys.append(key)
    label = ' for species {0}'.format(label) if label is not None else ''
    if replaced_keys:
        logging.info('\nThe loaded YAML file{0} seems to be from an older version of RMG/Arkane.\n'
                     'Some keywords will be automatically replaced before loading objects from this file.'.format(label))
    for key in replaced_keys:
        logging.info('Replacing keyword "{key}" with "{value}" in the Arkane YAML file.'.format(
            key=key, value=syntax_correction_dict[key]))
    if replaced_keys:
        logging.info('\n')
    return content


def is_pdep(job_list):
    """A helper function to determine whether a job is PressureDependenceJob or not"""
    for job in job_list:
        if isinstance(job, PressureDependenceJob):
            return True
    return False


def check_conformer_energy(energies, path):
    """
    Check to see that the starting energy of the species in the potential energy scan calculation
    is not 0.5 kcal/mol (or more) higher than any other energies in the scan. If so, print and 
    log a warning message.  
    """
    energies = np.array(energies, np.float64)
    e_diff = (energies[0] - np.min(energies)) * constants.E_h * constants.Na / 1000
    if e_diff >= 2:  # we choose 2 kJ/mol to be the critical energy
        logging.warning(f'The species corresponding to {os.path.basename(path)} is different in energy from the '
                        f'lowest energy conformer by {e_diff:.2f} kJ/mol. This can cause significant errors in '
                        f'your computed thermodynamic properties and rate coefficients.')


def get_element_mass(input_element, isotope=None):
    """
    Returns the mass and z number of the requested isotope for a given element.
    'input_element' can be wither the atomic number (integer) or an element symbol.
    'isotope' is an integer of the atomic z number. If 'isotope' is None, returns the most common isotope.
    Data taken from NIST, https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl (accessed October 2018)
    """
    symbol = None
    number = None

    if isinstance(input_element, int):
        symbol = symbol_by_number[input_element]
        number = input_element
    elif isinstance(input_element, str):
        symbol = input_element
        try:
            number = next(key for key, value in symbol_by_number.items() if value == input_element)
        except:
            symbol = input_element[0] + input_element[1].lower()
            number = [key for key, value in symbol_by_number.items() if value == symbol][0]

    if symbol is None or number is None:
        raise ValueError('Could not identify element {0}'.format(input_element))

    mass_list = mass_by_symbol[symbol]

    if isotope is not None:
        # a specific isotope is required
        for iso_mass in mass_list:
            if iso_mass[0] == isotope:
                mass = iso_mass[1]
                break
        else:
            raise ValueError("Could not find requested isotope {0} for element {1}".format(isotope, symbol))
    else:
        # no specific isotope is required
        if len(mass_list[0]) == 2:
            # isotope weight is unavailable, use the first entry
            mass = mass_list[0][1]
            logging.warning('Assuming isotope {0} is representative of element {1}'.format(mass_list[0][0], symbol))
        else:
            # use the most common isotope
            max_weight = mass_list[0][2]
            mass = mass_list[0][1]
            for iso_mass in mass_list:
                if iso_mass[2] > max_weight:
                    max_weight = iso_mass[2]
                    mass = iso_mass[1]
    return mass, number


symbol_by_number = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na',
                    12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
                    22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga',
                    32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb',
                    42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
                    52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm',
                    62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
                    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl',
                    82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
                    92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md',
                    102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
                    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}

# Structure of mass_by_symbol items: list(list(isotope1, mass1, weight1), list(isotope2, mass2, weight2), ...)
mass_by_symbol = {
    'H': [[1, 1.00782503224, 0.999885], [2, 2.01410177812, 0.000115], [3, 3.0160492779, 0]],
    'He': [[3, 3.0160293201, 0.00000134], [4, 4.00260325414, 0.99999866]],
    'Li': [[6, 6.0151228874, 0.0759], [7, 7.0160034366, 0.9241]],
    'Be': [[9, 9.012183066, 1]],
    'B': [[10, 10.01293695, 0.199], [11, 11.00930536, 0.801]],
    'C': [[12, 12.0000000, 0.9893], [13, 13.00335483507, 0.0107], [14, 14.0032419884, 0]],
    'N': [[14, 14.00307400443, 0.99636], [15, 15.00010889889, 0.00364]],
    'O': [[16, 15.99491461957, 0.99757], [17, 16.99913175651, 0.00038], [18, 17.99915961287, 0.00205]],
    'F': [[19, 18.99840316274, 1]],
    'Ne': [[20, 19.9924401762, 0.9048], [21, 20.993846685, 0.0027], [22, 21.991385114, 0.0925]],
    'Na': [[23, 22.9897692820, 1]],
    'Mg': [[24, 23.985041697, 0.7899], [25, 24.985836976, 0.1000], [26, 25.982592968, 0.1101]],
    'Al': [[27, 26.98153853, 1]],
    'Si': [[28, 27.97692653465, 0.92223], [29, 28.97649466491, 0.04685], [30, 29.973770136, 0.03092]],
    'P': [[31, 30.97376199843, 1]],
    'S': [[32, 31.9720711744, 0.9499], [33, 32.9714589098, 0.0075], [34, 33.967867004, 0.0425],
          [36, 35.96708071, 0.0001]],
    'Cl': [[35, 34.968852682, 0.7576], [37, 36.965902603, 0.2424]],
    'Ar': [[36, 35.967545105, 0.003336], [38, 37.96273211, 0.000629], [40, 39.9623831237, 0.996035]],
    'K': [[39, 38.9637064864, 0.932581], [40, 39.963998167, 0.000117], [41, 40.9618252579, 0.067302]],
    'Ca': [[40, 39.962590863, 0.96941], [42, 41.95861783, 0.00647], [43, 42.95876644, 0.00135],
           [44, 43.95548156, 0.02086], [46, 45.9536890, 0.00004], [48, 47.95252276, 0.00187]],
    'Sc': [[45, 44.95590829, 1]],
    'Ti': [[46, 45.95262772, 0.0825], [47, 46.95175879, 0.0744], [48, 47.94794198, 0.7372], [49, 48.94786568, 0.0541],
           [50, 49.94478689, 0.0518]],
    'V': [[50, 49.94715602, 0.00250], [51, 50.94395705, 0.99750]],
    'Cr': [[50, 49.94604184, 0.04345], [52, 51.94050624, 0.83789], [53, 52.94064816, 0.09501],
           [54, 53.93887917, 0.02365]],
    'Mn': [[55, 54.93804391, 1]],
    'Fe': [[54, 53.93960900, 0.05845], [56, 55.93493633, 0.91754], [57, 56.93539284, 0.02119],
           [58, 57.93327444, 0.00282]],
    'Co': [[59, 58.93319430, 1]],
    'Ni': [[58, 57.93534242, 0.68077], [60, 59.93078589, 0.26223], [61, 60.93105558, 0.011399],
           [62, 61.92834538, 0.036346], [64, 63.92796682, 0.009255]],
    'Cu': [[63, 62.92959773, 0.6915], [65, 64.92778971, 0.3085]],
    'Zn': [[64, 63.92914202, 0.4917], [66, 65.92603382, 0.2773], [67, 66.92712776, 0.0404], [68, 67.92484456, 0.1845],
           [70, 69.9253192, 0.0061]],
    'Ga': [[69, 68.9255735, 0.60108], [71, 70.92470259, 0.39892]],
    'Ge': [[70, 69.92424876, 0.2057], [72, 71.922075827, 0.2745], [73, 72.923458957, 0.0775],
           [74, 73.921177761, 0.3650], [76, 75.921402726, 0.0773]],
    'As': [[75, 74.92159458, 1]],
    'Se': [[74, 73.922475934, 0.0089], [76, 75.919213704, 0.0937], [77, 76.919914155, 0.0763],
           [78, 77.91730928, 0.2377], [80, 79.9165218, 0.4961], [82, 81.9166995, 0.0873]],
    'Br': [[79, 78.9183376, 0.5069], [81, 80.9162897, 0.4931]],
    'Kr': [[78, 77.92036495, 0.00355], [80, 79.91637809, 0.02286], [82, 81.91348274, 0.11593],
           [83, 82.91412716, 0.11500], [84, 83.9114977282, 0.56987], [86, 85.9106106269, 0.17279]],
    'Rb': [[85, 84.9117897380, 0.7217], [87, 86.9091805311, 0.2783]],
    'Sr': [[84, 83.9134191, 0.0056], [86, 85.9092606, 0.0986], [87, 86.9088775, 0.0700],
           [88, 87.9056125, 0.8258]],
    'Y': [[89, 88.9058403, 1]],
    'Zr': [[90, 89.9046977, 0.5145], [91, 90.9056396, 0.1122], [92, 91.9050347, 0.1715],
           [94, 93.9063108, 0.1738], [96, 95.9082714, 0.0280]],
    'Nb': [[93, 92.9063730, 1]],
    'Mo': [[92, 91.90680797, 0.1453], [94, 93.90508490, 0.0915], [95, 94.90583877, 0.1584], [96, 95.90467612, 0.1667],
           [97, 96.90601812, 0.0960], [98, 97.90540482, 0.2439], [100, 99.9074718, 0.0982]],
    'Tc': [[97, 96.9063667, ], [98, 97.9072124], [99, 98.9062508]],
    'Ru': [[96, 95.90759025, 0.0554], [98, 97.9052869, 0.0187], [99, 98.9059341, 0.1276], [100, 99.9042143, 0.1260],
           [101, 100.9055769, 0.1706], [102, 101.9043441, 0.3155], [104, 103.9054275, 0.1862]],
    'Rh': [[103, 102.905498, 1]],
    'Pd': [[102, 101.9056022, 0.0102], [104, 103.9040305, 0.1114], [105, 104.9050796, 0.2233],
           [106, 105.9034804, 0.2733], [108, 107.9038916, 0.2646], [110, 109.90517221, 0.1172]],
    'Ag': [[107, 106.9050916, 0.51839], [109, 108.9047553, 0.48161]],
    'Cd': [[106, 105.9064599, 0.0125], [108, 107.9041834, 0.0089], [110, 109.90300662, 0.1249],
           [111, 110.90418288, 0.1280], [112, 111.9027629, 0.2413], [113, 112.9044081, 0.1222],
           [114, 113.9033651, 0.2873], [116, 115.9047632, 0.0749]],
    'In': [[113, 112.9040618, 0.0429], [115, 114.9038788, 0.9571]],
    'Sn': [[112, 111.9048239, 0.0097], [114, 113.9027827, 0.0066], [115, 114.9033447, 0.0034],
           [116, 115.9017428, 0.1454], [117, 116.902954, 0.0768], [118, 117.9016066, 0.2422],
           [119, 118.9033112, 0.0859], [120, 119.9022016, 0.3258], [122, 121.9034438, 0.0463],
           [124, 123.9052766, 0.0579]],
    'Sb': [[121, 120.903812, 0.5721], [123, 122.9042132, 0.4279]],
    'Te': [[120, 119.9040593, 0.0009], [122, 121.9030435, 0.0255], [123, 122.9042698, 0.0089],
           [124, 123.9028171, 0.0474], [125, 124.9044299, 0.0707], [126, 125.9033109, 0.1884],
           [128, 127.9044613, 0.3174], [130, 129.9062227, 0.3408]],
    'I': [[127, 126.9044719, 1]],
    'Xe': [[124, 123.905892, 0.000952], [126, 125.9042983, 0.000890], [128, 127.903531, 0.019102],
           [129, 128.9047809, 0.264006], [130, 129.9035093, 0.040710], [131, 130.9050841, 0.212324],
           [132, 131.9041551, 0.269086], [134, 133.9053947, 0.104357], [136, 135.9072145, 0.088573]],
    'Cs': [[133, 132.905452, 1]],
    'Ba': [[130, 129.9063207, 0.00106], [132, 131.9050611, 0.00101], [134, 133.9045082, 0.02417],
           [135, 134.9056884, 0.06592], [136, 135.9045757, 0.07854], [137, 136.9058271, 0.11232],
           [138, 137.905247, 0.71698]],
    'La': [[138, 137.9071149, 0.0008881], [139, 138.9063563, 0.9991119]],
    'Ce': [[136, 135.9071292, 0.00185], [138, 137.905991, 0.00251], [140, 139.9054431, 0.88450],
           [142, 141.9092504, 0.11114]],
    'Pr': [[141, 140.9076576, 1]],
    'Nd': [[142, 141.907729, 0.27152], [143, 142.90982, 0.12174], [144, 143.910093, 0.23798],
           [145, 144.9125793, 0.08293], [146, 145.9131226, 0.17189], [148, 147.9168993, 0.05756],
           [150, 149.9209022, 0.05638]],
    'Pm': [[145, 144.9127559], [147, 146.915145]],
    'Sm': [[144, 143.9120065, 0.0307], [147, 146.9149044, 0.1499], [148, 147.9148292, 0.1124],
           [149, 148.9171921, 0.1382], [150, 149.9172829, 0.0738], [152, 151.9197397, 0.2675],
           [154, 153.9222169, 0.2275]],
    'Eu': [[151, 150.9198578, 0.4781], [153, 152.921238, 0.5219]],
    'Gd': [[152, 151.9197995, 0.0020], [154, 153.9208741, 0.0218], [155, 154.9226305, 0.1480],
           [156, 155.9221312, 0.2047], [157, 156.9239686, 0.1565], [158, 157.9241123, 0.2484],
           [160, 159.9270624, 0.2186]],
    'Tb': [[159, 158.9253547, 1]],
    'Dy': [[156, 155.9242847, 0.00056], [158, 157.9244159, 0.00095], [160, 159.9252046, 0.02329],
           [161, 160.9269405, 0.18889], [162, 161.9268056, 0.25475], [163, 162.9287383, 0.24896],
           [164, 163.9291819, 0.28260]],
    'Ho': [[165, 164.9303288, 1]],
    'Er': [[162, 161.9287884, 0.00139], [164, 163.9292088, 0.01601], [166, 165.9302995, 0.33503],
           [167, 166.9320546, 0.22869], [168, 167.9323767, 0.26978], [170, 169.9354702, 0.14910]],
    'Tm': [[169, 168.9342179, 1]],
    'Yb': [[168, 167.9338896, 0.00123], [170, 169.9347664, 0.02982], [171, 170.9363302, 0.1409],
           [172, 171.9363859, 0.2168], [173, 172.9382151, 0.16103], [174, 173.9388664, 0.32026],
           [176, 175.9425764, 0.12996]],
    'Lu': [[175, 174.9407752, 0.97401], [176, 175.9426897, 0.02599]],
    'Hf': [[174, 173.9400461, 0.0016], [176, 175.9414076, 0.0526], [177, 176.9432277, 0.1860],
           [178, 177.9437058, 0.2728], [179, 178.9458232, 0.1362], [180, 179.946557, 0.3508]],
    'Ta': [[180, 179.9474648, 0.0001201], [181, 180.9479958, 0.9998799]],
    'W': [[180, 179.9467108, 0.0012], [182, 181.9482039, 0.2650], [183, 182.9502228, 0.1431],
          [184, 183.9509309, 0.3064], [186, 185.9543628, 0.2843]],
    'Re': [[185, 184.9529545, 0.3740], [187, 186.9557501, 0.6260]],
    'Os': [[184, 183.9524885, 0.0002], [186, 185.953835, 0.0159], [187, 186.9557474, 0.0196],
           [188, 187.9558352, 0.1324], [189, 188.9581442, 0.1615], [190, 189.9584437, 0.2626],
           [192, 191.961477, 0.4078]],
    'Ir': [[191, 190.9605893, 0.373], [193, 192.9629216, 0.627]],
    'Pt': [[190, 189.9599297, 0.00012], [192, 191.9610387, 0.00782], [194, 193.9626809, 0.3286],
           [195, 194.9647917, 0.3378], [196, 195.9649521, 0.2521], [198, 197.9678949, 0.07356]],
    'Au': [[197, 196.9665688, 1]],
    'Hg': [[196, 195.9658326, 0.0015], [198, 197.9667686, 0.0997], [199, 198.9682806, 0.1687],
           [200, 199.9683266, 0.2310], [201, 200.9703028, 0.1318], [202, 201.9706434, 0.2986],
           [204, 203.973494, 0.0687]],
    'Tl': [[203, 202.9723446, 0.2952], [205, 204.9744278, 0.7048]],
    'Pb': [[204, 203.973044, 0.014], [206, 205.9744657, 0.241], [207, 206.9758973, 0.221],
           [208, 207.9766525, 0.524]],
    'Bi': [[209, 208.9803991, 1]],
    'Po': [[209, 208.9824308], [210, 209.9828741]],
    'At': [[210, 209.9871479], [211, 210.9874966]],
    'Rn': [[211, 210.9906011], [220, 220.0113941], [222, 222.0175782]],
    'Fr': [[223, 223.019736]],
    'Ra': [[223, 223.0185023], [224, 224.020212], [226, 226.0254103], [228, 228.0310707]],
    'Ac': [[227, 227.0277523]],
    'Th': [[230, 230.0331341, 0], [232, 232.0380558, 1]],
    'Pa': [[231, 231.0358842, 1]],
    'U': [[233, 233.0396355, 0], [234, 234.0409523, 0.000054], [235, 235.0439301, 0.007204],
          [236, 236.0455682, 0], [238, 238.0507884, 0.992742]],
    'Np': [[236, 236.04657], [237, 237.0481736]],
    'Pu': [[238, 238.0495601], [239, 239.0521636], [240, 240.0538138], [241, 241.0568517],
           [242, 242.0587428], [244, 244.0642053]],
    'Am': [[241, 241.0568293], [243, 243.0613813]],
    'Cm': [[243, 243.0613893], [244, 244.0627528], [245, 245.0654915], [246, 246.0672238],
           [247, 247.0703541], [248, 248.0723499]],
    'Bk': [[247, 247.0703073], [249, 249.0749877]],
    'Cf': [[249, 249.0748539], [250, 250.0764062], [251, 251.0795886], [252, 252.0816272]],
    'Es': [[252, 252.08298]],
    'Fm': [[257, 257.0951061]],
    'Md': [[258, 258.0984315], [260, 260.10365]],
    'No': [[259, 259.10103]],
    'Lr': [[262, 262.10961]],
    'Rf': [[267, 267.12179]],
    'Db': [[268, 268.12567]],
    'Sg': [[271, 271.13393]],
    'Bh': [[272, 272.13826]],
    'Hs': [[270, 270.13429]],
    'Mt': [[276, 276.15159]],
    'Ds': [[281, 281.16451]],
    'Rg': [[280, 280.16514]],
    'Cn': [[285, 285.17712]],
    'Nh': [[284, 284.17873]],
    'Fl': [[289, 289.19042]],
    'Mc': [[288, 288.19274]],
    'Lv': [[293, 293.20449]],
    'Ts': [[292, 292.20746]],
    'Og': [[294, 294.21392]]}


def get_center_of_mass(coords, numbers=None, symbols=None):
    """
    Calculate and return the 3D position of the center of mass of the current geometry.
    Either ``numbers`` or ``symbols`` must be given.

    Args:
        coords (np.array): Entries are 3-length lists of xyz coordinates for an atom.
        numbers (np.array, list): Entries are atomic numbers corresponding to coords.
        symbols (list): Entries are atom symbols corresponding to coords.

    Returns:
        np.array: The center of mass coordinates.
    """
    if symbols is None and numbers is None:
        raise IndexError('Either symbols or numbers must be given.')
    if numbers is not None:
        symbols = [symbol_by_number[number] for number in numbers]
    center, total_mass = np.zeros(3, np.float64), 0
    for coord, symbol in zip(coords, symbols):
        mass = get_element_mass(symbol)[0]
        center += mass * coord
        total_mass += mass
    center /= total_mass
    return center


def get_moment_of_inertia_tensor(coords, numbers=None, symbols=None):
    """
    Calculate and return the moment of inertia tensor for the current
    geometry in amu*angstrom^2. If the coordinates are not at the center of mass,
    they are temporarily shifted there for the purposes of this calculation.
    Adapted from J.W. Allen: https://github.com/jwallen/ChemPy/blob/master/chempy/geometry.py

    Args:
        coords (np.array): Entries are 3-length lists of xyz coordinates for an atom.
        numbers (np.array, list): Entries are atomic numbers corresponding to coords.
        symbols (list): Entries are atom symbols corresponding to coords.

    Returns:
        np.array: The 3x3 moment of inertia tensor.
    Raises:
        InputError: If neither ``symbols`` nor ``numbers`` are given, or if they have a different length than ``coords``
    """
    if symbols is None and numbers is None:
        raise InputError('Either symbols or numbers must be given.')
    if numbers is not None:
        symbols = [symbol_by_number[number] for number in numbers]
    if len(coords) != len(symbols):
        raise InputError(f'The number of atoms ({len(symbols)}) is not equal to the number of '
                         f'atomic coordinates ({len(list(coords))})')
    tensor = np.zeros((3, 3), np.float64)
    center_of_mass = get_center_of_mass(coords=coords, numbers=numbers, symbols=symbols)
    for symbol, coord in zip(symbols, coords):
        mass = get_element_mass(symbol)[0]
        cm_coord = coord - center_of_mass
        tensor[0, 0] += mass * (cm_coord[1] * cm_coord[1] + cm_coord[2] * cm_coord[2])
        tensor[1, 1] += mass * (cm_coord[0] * cm_coord[0] + cm_coord[2] * cm_coord[2])
        tensor[2, 2] += mass * (cm_coord[0] * cm_coord[0] + cm_coord[1] * cm_coord[1])
        tensor[0, 1] -= mass * cm_coord[0] * cm_coord[1]
        tensor[0, 2] -= mass * cm_coord[0] * cm_coord[2]
        tensor[1, 2] -= mass * cm_coord[1] * cm_coord[2]
    tensor[1, 0] = tensor[0, 1]
    tensor[2, 0] = tensor[0, 2]
    tensor[2, 1] = tensor[1, 2]
    return tensor


def get_principal_moments_of_inertia(coords, numbers=None, symbols=None):
    """
    Calculate and return the principal moments of inertia in amu*angstrom^2 in decending order
    and the corresponding principal axes for the current geometry.
    The moments of inertia are in translated to the center of mass. The principal axes have unit lengths.
    Adapted from J.W. Allen: https://github.com/jwallen/ChemPy/blob/master/chempy/geometry.py

    Args:
        coords (np.array): Entries are 3-length lists of xyz coordinates for an atom.
        numbers (np.array, list): Entries are atomic numbers corresponding to coords.
        symbols (list): Entries are atom symbols corresponding to coords.

    Returns:
        tuple: The principal moments of inertia.
        tuple: The corresponding principal axes.
    """
    tensor0 = get_moment_of_inertia_tensor(coords=coords, numbers=numbers, symbols=symbols)
    # Since tensor0 is real and symmetric, diagonalization is always possible
    principal_moments_of_inertia, axes = np.linalg.eig(tensor0)
    principal_moments_of_inertia, axes = zip(*sorted(zip(np.ndarray.tolist(principal_moments_of_inertia),
                                                         np.ndarray.tolist(axes)), reverse=True))
    return principal_moments_of_inertia, axes


def clean_dir(base_dir_path: str = '',
              files_to_delete: List[str] = None,
              file_extensions_to_delete: List[str] = None,
              files_to_keep: List[str] = None,
              sub_dir_to_keep: List[str] = None,
              ) -> None:
    """
    Clean up a directory. Commonly used for removing unwanted files after unit tests.

    Args:
        base_dir_path (str): absolute path of the directory to clean up.
        files_to_delete (list[str]): full name of the file (includes extension) to delete.
        file_extensions_to_delete: extensions of files to delete.
        files_to_keep: full name of the file (includes extension) to keep, files specified here will NOT be deleted even
                       if its extension is also in file_extensions_to_delete.
        sub_dir_to_keep: name of the subdirectories in the base directory to keep.
    """
    for item in os.listdir(base_dir_path):
        item_path = os.path.join(base_dir_path, item)
        if os.path.isfile(item_path):
            item_extension = os.path.splitext(item_path)[-1]
            if item in files_to_delete or (item_extension in file_extensions_to_delete and item not in files_to_keep):
                os.remove(item_path)
        else:
            # item is sub-directory
            if os.path.split(item_path)[-1] in sub_dir_to_keep:
                continue
            shutil.rmtree(item_path)


def convert_imaginary_freq_to_negative_float(freq: Union[str, float, int]):
    """
    Convert a string representation of an imaginary frequency into a negative float representation, e.g.:
        '635.0i' -> -635.0
        '500.0' -> 500.0

    Args:
        freq (str): The imaginary frequency representation.

    Returns:
        float: A float representation of the frequency value.
    """
    if isinstance(freq, str) and freq.endswith('i'):
        freq = float(freq[:-1]) * -1
    return float(freq)
