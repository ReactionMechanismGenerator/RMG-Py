#!/usr/bin/env python3
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
This module provides the :class:`ErrorCancelingScheme` and related classes for the automatic generation of error
canceling reactions (e.g. isodesmic reactions). This code is heavily based on algorithms and ideas found in the existing
literature, including the following:

Buerger, P., Akroyd, J., Mosbach, S., & Kraft, M. (2018). A systematic method to estimate and validate enthalpies of
formation using error-cancelling balanced reactions. Combustion and Flame (Vol. 187).
https://doi.org/10.1016/j.combustflame.2017.08.013

Dobek, F. J., Ranasinghe, D. S., Throssell, K., & Petersson, G. A. (2013). Evaluation of the heats of formation of
corannulene and C60 by means of inexpensive theoretical procedures. Journal of Physical Chemistry A, 117(22), 4726–4730.
https://doi.org/10.1021/jp404158v
"""

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.quantity import ScalarQuantity


class ErrorCancelingSpecies:
    """Class for target and known (reference) species participating in an error canceling reaction"""

    def __init__(self, molecule, low_level_hf298, model_chemistry, high_level_hf298=None, source=None):
        """

        Args:
            molecule (Molecule): The RMG Molecule object with connectivity information
            low_level_hf298 (ScalarQuantity): evaluated using a lower level of theory (e.g. DFT)
            model_chemistry (str): Level of theory used to calculate the low level thermo
            high_level_hf298 (ScalarQuantity, optional): evaluated using experimental data
                or a high level of theory that is serving as the "reference" for the isodesmic calculation
            source (str): Literature source from which the high level data was taken
        """
        if isinstance(molecule, Molecule):
            self.molecule = molecule
        else:
            raise ValueError(f'ErrorCancelingSpecies molecule attribute must be an rmgpy Molecule object. Instead a '
                             f'{type(molecule)} object was given')

        if isinstance(model_chemistry, str):
            self.model_chemistry = model_chemistry
        else:
            raise ValueError(f'The model chemistry string used to calculate the low level Hf298 must be provided '
                             f'consistency checks. Instead, a {type(model_chemistry)} object was given')

        if not isinstance(low_level_hf298, ScalarQuantity):
            if isinstance(low_level_hf298, tuple):
                low_level_hf298 = ScalarQuantity(*low_level_hf298)
            else:
                raise TypeError(f'Low level Hf298 should be a ScalarQuantity object or its tuple representation, but '
                                f'received {low_level_hf298} instead.')
        self.low_level_hf298 = low_level_hf298

        # If the species is a reference species, then the high level data is already known
        if high_level_hf298 is not None and not isinstance(high_level_hf298, ScalarQuantity):
            if isinstance(high_level_hf298, tuple):
                high_level_hf298 = ScalarQuantity(*high_level_hf298)
            else:
                raise TypeError(f'High level Hf298 should be a ScalarQuantity object or its tuple representation, but '
                                f'received {high_level_hf298} instead.')
        self.high_level_hf298 = high_level_hf298
        self.source = source

    def __repr__(self):
        return f'<ErrorCancelingSpecies {self.molecule.to_smiles()}>'


class ErrorCancelingReaction:
    """Class for representing an error canceling reaction, with the target species being an implicit reactant"""

    def __init__(self, target, species):
        """
        Initialize an error canceling reaction from ErrorCancelingSpecies objects

        The reaction must be written with the target species participating as a reactant with stoichiometric coefficient
        v=-1

        The species dictionary should be given as species/stochiometric coefficient pairs in the following format:
        {ErrorCancelingSpecies_1: v_1, ErrorCancelingSpecies_2: v_2, ... }

        Args:
            target (ErrorCancelingSpecies): high level H_f(298 K) will be estimated for this species
            species (dict): species taking place in the reaction (excluding the target)
        """

        self.target = target
        self.model_chemistry = self.target.model_chemistry

        # Perform a consistency check that all species are using the same model chemistry
        for spcs in species.keys():
            if spcs.model_chemistry != self.model_chemistry:
                raise ValueError(f'Species {spcs} has model chemistry {spcs.model_chemistry}, which does not match the '
                                 f'model chemistry of the reaction of {self.model_chemistry}')

        # Does not include the target, which is handled separately.
        self.species = species

    def __repr__(self):
        reactant_string = f'1*{self.target.molecule.to_smiles()}'
        product_string = ''
        for spcs, coeff in self.species.items():
            if coeff > 0:
                product_string += f' + {int(coeff)}*{spcs.molecule.to_smiles()}'
            else:
                reactant_string += f' + {-1*int(coeff)}*{spcs.molecule.to_smiles()}'

        return f'<ErrorCancelingReaction {reactant_string} <=> {product_string[3:]} >'

    def calculate_target_thermo(self):
        """
        Estimate the high level thermochemistry for the target species using the error canceling scheme

        Returns:
            rmgpy.quantity.ScalarQuantity: Hf298 in 'J/mol' estimated for the target species
        """
        low_level_h_rxn = sum(spec[0].low_level_hf298.value_si*spec[1] for spec in self.species.items()) - \
            self.target.low_level_hf298.value_si

        target_thermo = sum(spec[0].high_level_hf298.value_si*spec[1] for spec in self.species.items()) - \
            low_level_h_rxn
        return ScalarQuantity(target_thermo, 'J/mol')


class SpeciesConstraints:
    """
    A class for defining and enumerating constraints to ReferenceSpecies objects for error canceling reactions
    """

    def __init__(self, target, reference_list, conserve_bonds=True, conserve_ring_size=True):
        """
        Define the constraints that will be enforced, and determine the mapping of indices in the constraint vector to
        the labels for these constraints.

        To reduce the size of the linear programming problem that will try to find error canceling reactions of the
        target and subsets of the reference species, the `reference_species` list is automatically pruned to remove
        species that have additional atom, bond, and/or ring attributes not found in the target molecule.

        Args:
            target (ErrorCancelingSpecies): The target species of the error canceling reaction scheme
            reference_list(list): A list of ErrorCancelingSpecies objects for the reference
                species that can participate in the error canceling reaction scheme
            conserve_bonds (bool, optional): Enforce the number of each bond type be conserved
            conserve_ring_size (bool, optional): Enforce that the number of each ring size be conserved
        """

        self.target = target
        self.all_reference_species = reference_list
        self.reference_species = []
        self.conserve_bonds = conserve_bonds
        self.conserve_ring_size = conserve_ring_size
        self.constraint_map = self._get_constraint_map()

    def _get_constraint_map(self):
        # Enumerate all of the constraints in the target molecule to initialize the constraint mapping
        constraint_map = {label: i for i, label in enumerate(self.target.molecule.get_element_count().keys())}
        if self.conserve_bonds:
            j = len(constraint_map)
            constraint_map.update(
                {label: j + i for i, label in enumerate(self.target.molecule.enumerate_bonds().keys())})
        if self.conserve_ring_size:
            j = len(constraint_map)
            possible_rings_sizes = set(f'{len(sssr)}_ring' for sssr in
                                       self.target.molecule.get_smallest_set_of_smallest_rings())
            constraint_map.update({label: j + i for i, label in enumerate(possible_rings_sizes)})

        return constraint_map

    def _enumerate_constraints(self, species):
        """
        Determine the constraint vector for a molecule given the enforced constraints

        Args:
            species (ErrorCancelingSpecies): Species whose constraints are to be enumerated

        Returns:
            np.ndarray: vector of the number of instances of each constraining feature e.g. number of carbon atoms
        """
        constraint_vector = np.zeros(len(self.constraint_map))
        molecule = species.molecule

        try:
            atoms = molecule.get_element_count()
            for atom_label, count in atoms.items():
                constraint_vector[self.constraint_map[atom_label]] += count

            if self.conserve_bonds:
                bonds = molecule.enumerate_bonds()
                for bond_label, count in bonds.items():
                    constraint_vector[self.constraint_map[bond_label]] += count

            if self.conserve_ring_size:
                rings = molecule.get_smallest_set_of_smallest_rings()
                for ring in rings:
                    constraint_vector[self.constraint_map[f'{len(ring)}_ring']] += 1
        except KeyError:  # This molecule has a feature not found in the target molecule. Return None to exclude this
            return None

        return constraint_vector

    def calculate_constraints(self):
        """
        Calculate the constraint vector for the target and the constraint matrix for all allowable reference species

        Returns:
            tuple(np.ndarray, np.ndarray)
            - target constraint vector (1 x len(constraints))
            - constraint matrix for allowable reference species (len(self.reference_species) x len(constraints))
        """
        target_constraints = self._enumerate_constraints(self.target)
        constraint_matrix = []

        for spcs in self.all_reference_species:
            spcs_constraints = self._enumerate_constraints(spcs)
            if spcs_constraints is not None:  # This species is allowed
                self.reference_species.append(spcs)
                constraint_matrix.append(spcs_constraints)

        return target_constraints, np.array(constraint_matrix, dtype=int)


class ErrorCancelingScheme:
    """
    A Base class for calculating target species thermochemistry using error canceling reactions
    """

    def __init__(self, target, reference_set, conserve_bonds, conserve_ring_size):
        """

        Args:
            target (ErrorCancelingSpecies): Species whose Hf298 will be calculated using error canceling reactions
            reference_set (list): list of reference species (as ErrorCancelingSpecies) that can participate
                in error canceling reactions to help calculate the thermochemistry of the target
            conserve_bonds (bool): Flag to determine if the number and type of each bond must be conserved in each error
                canceling reaction
            conserve_ring_size (bool): Flag to determine if the number of each ring size must be conserved in each error
                canceling reaction
        """

        self.target = target
        self.constraints = SpeciesConstraints(target, reference_set, conserve_bonds=conserve_bonds,
                                              conserve_ring_size=conserve_ring_size)

        self.target_constraint, self.constraint_matrix = self.constraints.calculate_constraints()
        self.reference_species = self.constraints.reference_species


class IsodesmicScheme(ErrorCancelingScheme):
    """
    An error canceling reaction where the number and type of both atoms and bonds are conserved
    """
    def __init__(self, target, reference_set):
        super().__init__(target, reference_set, conserve_bonds=True, conserve_ring_size=False)


class IsodesmicRingScheme(ErrorCancelingScheme):
    """
    A stricter form of the traditional isodesmic reaction scheme where the number of each ring size is also conserved
    """
    def __init__(self, target, reference_set):
        super().__init__(target, reference_set, conserve_bonds=True, conserve_ring_size=True)


if __name__ == '__main__':
    pass
