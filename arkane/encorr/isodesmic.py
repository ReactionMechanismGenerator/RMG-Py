#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2024 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
from copy import deepcopy
from typing import List, Union

import numpy as np
from scipy.optimize import LinearConstraint, milp

from arkane.modelchem import LOT
from rmgpy.molecule import Bond, Molecule
from rmgpy.quantity import ScalarQuantity


class ErrorCancelingSpecies:
    """Class for target and known (reference) species participating in an error canceling reaction"""

    def __init__(
        self,
        molecule,
        low_level_hf298,
        level_of_theory,
        high_level_hf298=None,
        source=None,
    ):
        """

        Args:
            molecule (Molecule): The RMG Molecule object with connectivity information
            low_level_hf298 (ScalarQuantity): evaluated using a lower level of theory (e.g. DFT)
            level_of_theory ((Composite)LevelOfTheory): Level of theory used to calculate the low level thermo
            high_level_hf298 (ScalarQuantity, optional): evaluated using experimental data
                or a high level of theory that is serving as the "reference" for the isodesmic calculation
            source (str): Literature source from which the high level data was taken
        """
        if isinstance(molecule, Molecule):
            self.molecule = molecule
        else:
            raise ValueError(
                f"ErrorCancelingSpecies molecule attribute must be an rmgpy Molecule object. Instead a "
                f"{type(molecule)} object was given"
            )

        if isinstance(level_of_theory, LOT):
            self.level_of_theory = level_of_theory
        else:
            raise ValueError(
                f"The level of theory used to calculate the low level Hf298 must be provided "
                f"for consistency checks. Instead, a {type(level_of_theory)} object was given"
            )

        if not isinstance(low_level_hf298, ScalarQuantity):
            if isinstance(low_level_hf298, tuple):
                low_level_hf298 = ScalarQuantity(*low_level_hf298)
            else:
                raise TypeError(
                    f"Low level Hf298 should be a ScalarQuantity object or its tuple representation, but "
                    f"received {low_level_hf298} instead."
                )
        self.low_level_hf298 = low_level_hf298

        # If the species is a reference species, then the high level data is already known
        if high_level_hf298 is not None and not isinstance(
            high_level_hf298, ScalarQuantity
        ):
            if isinstance(high_level_hf298, tuple):
                high_level_hf298 = ScalarQuantity(*high_level_hf298)
            else:
                raise TypeError(
                    f"High level Hf298 should be a ScalarQuantity object or its tuple representation, but "
                    f"received {high_level_hf298} instead."
                )
        self.high_level_hf298 = high_level_hf298
        self.source = source

    def __repr__(self):
        return f"<ErrorCancelingSpecies {self.molecule.to_smiles()}>"


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
        self.level_of_theory = self.target.level_of_theory

        # Perform a consistency check that all species are using the same level of theory
        for spcs in species.keys():
            if spcs.level_of_theory != self.level_of_theory:
                raise ValueError(
                    f"Species {spcs} has level of theory {spcs.level_of_theory}, which does not match the "
                    f"level of theory of the reaction of {self.level_of_theory}"
                )

        # Does not include the target, which is handled separately.
        self.species = species

    def __repr__(self):
        reactant_string = f"1*{self.target.molecule.to_smiles()}"
        product_string = ""
        for spcs, coeff in self.species.items():
            if coeff > 0:
                product_string += f" + {int(coeff)}*{spcs.molecule.to_smiles()}"
            else:
                reactant_string += f" + {-1*int(coeff)}*{spcs.molecule.to_smiles()}"

        return f"<ErrorCancelingReaction {reactant_string} <=> {product_string[3:]} >"

    def calculate_target_thermo(self):
        """
        Estimate the high level thermochemistry for the target species using the error canceling scheme

        Returns:
            rmgpy.quantity.ScalarQuantity: Hf298 in 'J/mol' estimated for the target species
        """
        low_level_h_rxn = (
            sum(
                spec.low_level_hf298.value_si * coeff
                for spec, coeff in self.species.items()
            )
            - self.target.low_level_hf298.value_si
        )

        target_thermo = (
            sum(
                spec.high_level_hf298.value_si * coeff
                for spec, coeff in self.species.items()
            )
            - low_level_h_rxn
        )
        return ScalarQuantity(target_thermo, "J/mol")


class AtomConstraint:
    def __init__(self, label, connections=None):
        self.label = label
        self.connections = connections if connections is not None else []

    def __eq__(self, other):
        if isinstance(other, AtomConstraint):
            if self.label == other.label:
                return self.match_connections(other)

            return False

        else:
            raise NotImplementedError(
                f"AtomConstraint object has no __eq__ defined for other object of type "
                f"{type(other)}"
            )
        
    def match_connections(self, other):
        if len(self.connections) != len(other.connections):
            return False

        connections = deepcopy(other.connections)
        for c in self.connections:
            for i, c_other in enumerate(connections):
                if c == c_other:
                    break
            else:
                return False
            connections.pop(i)

        return True

    def __repr__(self):
        return f"{self.label}" + "".join([f"({c})" for c in self.connections])


class BondConstraint:
    def __init__(self, atom1, atom2, bond_order):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_order = bond_order

    def __eq__(self, other):
        if isinstance(other, BondConstraint):
            if self.bond_order == other.bond_order:
                if (self.atom1 == other.atom1 and self.atom2 == other.atom2) or (
                    self.atom1 == other.atom2 and self.atom2 == other.atom1
                ):
                    return True
            return False

        if isinstance(other, GenericConstraint):
            return False

        else:
            raise NotImplementedError(
                f"BondConstraint object has no __eq__ defined for other object of type "
                f"{type(other)}"
            )

    def __repr__(self):
        symbols = ["", "-", "=", "#"]
        return f"{self.atom1}{symbols[self.bond_order]}{self.atom2}"


class Connection:
    def __init__(self, atom, bond_order):
        self.atom = atom
        self.bond_order = bond_order

    def __eq__(self, other):
        if isinstance(other, Connection):
            if self.bond_order == other.bond_order:
                if self.atom == other.atom:
                    return True
            return False

        else:
            raise NotImplementedError(
                f"Connection object has no __eq__ defined for other object of type {type(other)}"
            )

    def __repr__(self):
        symbols = ["", "-", "=", "#"]
        return f"{symbols[self.bond_order]}{self.atom}"


class GenericConstraint:
    def __init__(self, constraint_str: str):
        self.constraint_str = constraint_str

    def __eq__(self, other):
        if isinstance(other, BondConstraint):
            return False
        elif isinstance(other, GenericConstraint):
            return self.constraint_str == other.constraint_str
        else:
            raise NotImplementedError(
                f"GenericConstraint object has no __eq__ defined for other object of "
                f"type {type(other)}"
            )
    def __repr__(self):
        return self.constraint_str


def bond_centric_constraints(
    species: ErrorCancelingSpecies, constraint_class: str
) -> List[BondConstraint]:
    constraints = []
    contraint_func = CONSTRAINT_CLASSES[constraint_class]
    molecule = species.molecule

    for bond in molecule.get_all_edges():
        constraints.append(contraint_func(bond))

    return constraints


def _buerger_rc2(bond: Bond) -> BondConstraint:
    atom1 = bond.atom1
    atom2 = bond.atom2

    if atom1.symbol > atom2.symbol:
        atom1, atom2 = atom2, atom1

    atom1 = AtomConstraint(label=atom1.symbol)
    atom2 = AtomConstraint(label=atom2.symbol)

    return BondConstraint(atom1=atom1, atom2=atom2, bond_order=int(bond.order))


def _buerger_rc3(bond: Bond) -> BondConstraint:
    atom1 = bond.atom1
    atom2 = bond.atom2

    if atom1.symbol > atom2.symbol:
        atom1, atom2 = atom2, atom1

    atom1 = AtomConstraint(label=f"{atom1.symbol}{atom1.connectivity1}")
    atom2 = AtomConstraint(label=f"{atom2.symbol}{atom2.connectivity1}")

    return BondConstraint(atom1=atom1, atom2=atom2, bond_order=int(bond.order))


def _buerger_rc4(bond: Bond) -> BondConstraint:
    atom1 = bond.atom1
    atom2 = bond.atom2

    if atom1.symbol > atom2.symbol:
        atom1, atom2 = atom2, atom1

    atoms = []

    for atom in [atom1, atom2]:
        connections = []
        for a, b in atom.bonds.items():
            ac = AtomConstraint(label=f"{a.symbol}{a.connectivity1}")
            bond_order = b.order
            connections.append(Connection(atom=ac, bond_order=bond_order))
        atoms.append(
            AtomConstraint(
                label=f"{atom.symbol}{atom.connectivity1}", connections=connections
            )
        )

    return BondConstraint(atom1=atoms[0], atom2=atoms[1], bond_order=int(bond.order))


class SpeciesConstraints:
    """
    A class for defining and enumerating constraints to ReferenceSpecies objects for error canceling reactions
    """

    def __init__(
        self,
        target,
        reference_list,
        isodesmic_class="rc2",
        conserve_ring_size=True,
        limit_charges=True,
        limit_scope=True,
    ):
        """
        Define the constraints that will be enforced, and determine the mapping of indices in the constraint vector to
        the labels for these constraints.

        Notes:
            To reduce the size of the linear programming problem that will try to find error canceling reactions of the
            target and subsets of the reference species, the `reference_species` list is automatically pruned to remove
            species that have additional atom, bond, and/or ring attributes not found in the target molecule.

            Charge is also explicitly conserved, as there are charged species in the reference database

        Args:
            target (ErrorCancelingSpecies): The target species of the error canceling reaction scheme
            reference_list(list): A list of ErrorCancelingSpecies objects for the reference
                species that can participate in the error canceling reaction scheme
            isodesmic_class (str, optional): Reaction classes as defined by Buerger et al. that determine how specific
                the constraints are.
            conserve_ring_size (bool, optional): Enforce that the number of each ring size be conserved.
            limit_charges (bool, optional): Only allow species in the reaction that are within the range [C, 0] for
                anions or [0, C] for cations where "C" is the charge of the target
            limit_scope (bool, optional): Exclude any molecules from the reference set that have features not contained
                in the target molecule. This will reduce the size of the linear programing problem being solved to yield
                faster, possibly more accurate results
        """

        self.target = target
        self.all_reference_species = reference_list
        self.reference_species = []
        self.isodesmic_class = isodesmic_class
        self.conserve_ring_size = conserve_ring_size
        self.limit_charges = limit_charges
        self.limit_scope = limit_scope

    def _get_constraint_lists(self):
        full_constraints_list = [self._get_all_constraints(self.target)]
        for ref_spcs in self.all_reference_species:
            full_constraints_list.append(self._get_all_constraints(ref_spcs))

        return full_constraints_list

    def _get_ring_constraints(
        self, species: ErrorCancelingSpecies
    ) -> List[GenericConstraint]:
        ring_features = []
        rings = species.molecule.get_symmetrized_smallest_set_of_smallest_rings()
        for ring in rings:
            ring_features.append(GenericConstraint(constraint_str=f"{len(ring)}_ring"))

        return ring_features

    def _get_all_constraints(
        self, species: ErrorCancelingSpecies
    ) -> List[Union[BondConstraint, GenericConstraint]]:
        features = bond_centric_constraints(species, self.isodesmic_class)
        if self.conserve_ring_size:
            features += self._get_ring_constraints(species)

        features.sort(key=lambda x: x.__repr__())
        return features

    def _enumerate_constraints(self, full_constraints_list):
        """
        Return the target constraint counts and the reference constraint counts.
        """

        target_constraints = full_constraints_list[0]
        reference_constraintss = full_constraints_list[1:]

        # Enumerate through the constraints of reference species and keep only those that are present in the target
        enumerated_reference_constraintss = []

        if self.limit_scope:

            for i, ref_spc in enumerate(self.all_reference_species):

                if not all(constraint in target_constraints for constraint in reference_constraintss[i]):
                    continue

                self.reference_species.append(ref_spc)
                enumerated_reference_constraintss.append(reference_constraintss[i])

        else:
            self.reference_species = self.all_reference_species
            enumerated_reference_constraintss = reference_constraintss
        
        # Get a list of the unique constraints sorted by their string representation
        if self.limit_scope:

            # The scope of constraints to consider is the target constraints
            unique_constraints = self._get_unique_constraints(target_constraints)
            unique_constraints.sort(key=lambda x: x.__repr__())

        else:
            all_constraints = target_constraints + [constraint for constraints in enumerated_reference_constraintss for constraint in constraints]
            unique_constraints = self._get_unique_constraints(all_constraints)
            unique_constraints.sort(key=lambda x: x.__repr__())
        
        # Get the counts of each unique constraint in the target and reference constraints
        target_constraint_counts = [target_constraints.count(c) for c in unique_constraints]
        reference_constraint_counts = []

        for i, ref_constraints in enumerate(enumerated_reference_constraintss):
            reference_constraint_counts.append([ref_constraints.count(c) for c in unique_constraints])

        return target_constraint_counts, reference_constraint_counts
    
    def _get_unique_constraints(self, constraints):
        # Constraints are unhashable, so we need to use some workaround to get unique constraints
        constraint_dict = {constraint.__repr__(): constraint for constraint in constraints}
        return list(constraint_dict.values())

    def _enumerate_charge_constraints(self, target_constraints, reference_constraints):
        charge = self.target.molecule.get_net_charge()
        target_constraints.append(charge)

        for i, spcs in enumerate(self.reference_species):
            reference_constraints[i].append(spcs.molecule.get_net_charge())

        if self.limit_charges:
            allowed_reference_species = []
            new_constraints = []

            if charge < 0:
                allowable_charges = list(range(charge, 0))
            elif charge > 0:
                allowable_charges = list(range(1, charge + 1))
            else:
                allowable_charges = [0]
            for i, spcs in enumerate(self.reference_species):
                if reference_constraints[i][-1] in allowable_charges:
                    allowed_reference_species.append(spcs)
                    new_constraints.append(reference_constraints[i])

            reference_constraints = new_constraints
            self.reference_species = allowed_reference_species

        return target_constraints, reference_constraints

    def _enumerate_element_constraints(self, target_constraints, reference_constraints):
        all_elements = set()
        for spc in self.reference_species:
            all_elements.update(spc.molecule.get_element_count().keys())

        # Check that the target and reference species have the same elements to be able to satisfy mass conservation
        if set(self.target.molecule.get_element_count().keys()) != all_elements:
            logging.warning(
                "Target species and reference species do not have the same elements:"
                f"Target: {' '.join(self.target.molecule.get_element_count().keys())}"
                f"Reference: {all_elements}"
            )

        all_elements.update(self.target.molecule.get_element_count().keys())
        all_elements = sorted(list(all_elements))

        element_count = self.target.molecule.get_element_count()
        new_constraints = [element_count.get(element, 0) for element in all_elements]
        target_constraints.extend(new_constraints)

        for i, spc in enumerate(self.reference_species):
            element_count = spc.molecule.get_element_count()
            new_constraints = [
                element_count.get(element, 0) for element in all_elements
            ]
            reference_constraints[i].extend(new_constraints)

        return target_constraints, reference_constraints

    def calculate_constraints(self):
        """
        Calculate the constraint vector for the target and the constraint matrix for all allowable reference species

        Returns:
            tuple(np.ndarray, np.ndarray)
            - target constraint vector (1 x len(constraints))
            - constraint matrix for allowable reference species (len(self.reference_species) x len(constraints))
        """
        full_constraint_list = self._get_constraint_lists()
        target_constraints, reference_constraints = self._enumerate_constraints(
            full_constraint_list
        )
        target_constraints, reference_constraints = self._enumerate_charge_constraints(
            target_constraints, reference_constraints
        )
        target_constraints, reference_constraints = self._enumerate_element_constraints(
            target_constraints, reference_constraints
        )

        target_constraints = np.array(target_constraints, dtype=int)
        constraint_matrix = np.array(reference_constraints, dtype=int)

        return target_constraints, constraint_matrix


def _clean_up_constraints(target_constraints, constraint_matrix):
    # make sure that the constraint matrix is 2d
    if len(constraint_matrix.shape) == 1:
        constraint_matrix = np.array([constraint_matrix], dtype=int)

    # Remove any columns that are all zeros
    zero_indices = np.where(~constraint_matrix.any(axis=0))[0]
    # Check that this wouldn't eliminate a non-zero target entry
    for z in zero_indices:
        if (
            target_constraints[z] != 0
        ):  # This problem is not solvable. Return None to signal this
            return None, None
    indices = [i for i in range(constraint_matrix.shape[1]) if i not in zero_indices]
    constraint_matrix = np.take(constraint_matrix, indices=indices, axis=1)
    target_constraints = np.take(target_constraints, indices=indices)

    return target_constraints, constraint_matrix


class ErrorCancelingScheme:
    """
    A Base class for calculating target species thermochemistry using error canceling reactions
    """

    def __init__(
        self,
        target,
        reference_set,
        isodesmic_class,
        conserve_ring_size,
        limit_charges,
        limit_scope,
    ):
        """

        Args:
            target (ErrorCancelingSpecies): Species whose Hf298 will be calculated using error canceling reactions
            reference_set (list): list of reference species (as ErrorCancelingSpecies) that can participate
                in error canceling reactions to help calculate the thermochemistry of the target
            conserve_ring_size (bool): Flag to determine if the number of each ring size must be conserved in each error
                canceling reaction
        """

        self.target = target
        self.constraints = SpeciesConstraints(
            target,
            reference_set,
            isodesmic_class=isodesmic_class,
            conserve_ring_size=conserve_ring_size,
            limit_charges=limit_charges,
            limit_scope=limit_scope,
        )

        (
            self.target_constraint,
            self.constraint_matrix,
        ) = self.constraints.calculate_constraints()
        self.reference_species = self.constraints.reference_species

    def _find_error_canceling_reaction(self, reference_subset):
        """
        Automatically find a valid error canceling reaction given a subset of the available benchmark species. This
        is done by solving a mixed integer linear programming (MILP) problem similiar to
        Buerger et al. (https://doi.org/10.1016/j.combustflame.2017.08.013)

        Args:
            reference_subset (list): A list of indices from self.reference_species that can participate in the reaction

        Returns:
            tuple(ErrorCancelingReaction, np.ndarray)
            - Reaction with the target species (if a valid reaction is found, else ``None``)
            - Indices (of the subset) for the species that participated in the return reaction
        """

        # Define the constraints based on the provided subset
        c_matrix = np.take(self.constraint_matrix, reference_subset, axis=0)

        # Remove unnecessary constraints
        target_constraint, c_matrix = _clean_up_constraints(
            self.target_constraint, c_matrix
        )
        if target_constraint is None or c_matrix is None:  # The problem is not solvable
            return None, None

        # Setup MILP problem
        c_matrix = np.tile(c_matrix, (2, 1))
        sum_constraints = np.sum(c_matrix, 1, dtype=int)
        targets = -1 * target_constraint
        m = c_matrix.shape[0]
        n = c_matrix.shape[1]
        split = int(m / 2)

        constraints = tuple((LinearConstraint(A=np.concatenate((c_matrix[:split, j], -1 * c_matrix[split:, j])), lb=targets[j], ub=targets[j]) for j in range(n)))

        result = milp(
            sum_constraints,
            integrality=1,
            constraints=constraints,
            options={"time_limit": 10},
        )

        if result.status != 0:
            logging.warning("Optimization could not find a valid solution.")
            return None, None
        
        solution = result.x

        reaction = ErrorCancelingReaction(self.target, dict())
        subset_indices = []
        for index, v in enumerate(solution):
            if v > 0:
                subset_indices.append(index % split)
                if index < split:
                    reaction.species.update(
                        {self.reference_species[reference_subset[index]]: -v}
                    )
                else:
                    reaction.species.update(
                        {self.reference_species[reference_subset[index % split]]: v}
                    )

        return reaction, np.array(subset_indices)

    def multiple_error_canceling_reaction_search(
        self, n_reactions_max=20,
    ):
        """
        Generate multiple error canceling reactions involving the target and a subset of the reference species.

        To do this, a rudimentary search is implemented whereby all possible combinations of the species participating
        in the previously found reaction are excluded from the reference species subset for the next generation process.
        This is implemented using a FIFO queue structure.

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will be returned, after which no
                further searching will occur even if there are possible subsets left in the queue.

        Returns:
            list: A list of the found error canceling reactions
        """
        subset_queue = [np.arange(0, len(self.reference_species))]
        reaction_list = []

        while (len(subset_queue) != 0) and (len(reaction_list) < n_reactions_max):
            subset = subset_queue.pop()
            if len(subset) == 0:
                continue
            reaction, subset_indices = self._find_error_canceling_reaction(
                subset
            )
            if reaction is None:
                continue
            else:
                reaction_list.append(reaction)

                for index in subset_indices:
                    subset_queue.append(np.delete(subset, index))

                # Clean up the queue to remove subsets that would allow for already found reactions
                new_queue = []
                reaction_indices = [subset[i] for i in subset_indices]
                for s in subset_queue:
                    if not all([i in s for i in reaction_indices]):
                        new_queue.append(s)

                subset_queue = new_queue

        return reaction_list

    def calculate_target_enthalpy(self, n_reactions_max=5):
        """
        Perform a multiple error canceling reactions search and calculate hf298 for the target species by taking the
        median hf298 value from among the error canceling reactions found

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will returned, after which no
                further searching will occur even if there are possible subsets left in the queue.

        Returns:
            tuple(ScalarQuantity, list)
            - Standard heat of formation at 298 K calculated for the target species
            - reaction list containing all error canceling reactions found
        """
        reaction_list = self.multiple_error_canceling_reaction_search(n_reactions_max)
        if len(reaction_list) == 0:  # No reactions found
            return None, reaction_list
        h298_list = np.zeros(len(reaction_list))

        for i, rxn in enumerate(reaction_list):
            h298_list[i] = rxn.calculate_target_thermo().value_si

        return ScalarQuantity(np.median(h298_list), "J/mol"), reaction_list


class IsodesmicScheme(ErrorCancelingScheme):
    """
    An error canceling reaction where the number and type of both atoms and bonds are conserved
    """

    def __init__(self, target, reference_set):
        super().__init__(
            target,
            reference_set,
            isodesmic_class="rc2",
            conserve_ring_size=False,
            limit_charges=True,
            limit_scope=True,
        )


class IsodesmicRingScheme(ErrorCancelingScheme):
    """
    A stricter form of the traditional isodesmic reaction scheme where the number of each ring size is also conserved
    """

    def __init__(self, target, reference_set):
        super().__init__(
            target,
            reference_set,
            isodesmic_class="rc2",
            conserve_ring_size=True,
            limit_charges=True,
            limit_scope=True,
        )


CONSTRAINT_CLASSES = {"rc2": _buerger_rc2, "rc3": _buerger_rc3, "rc4": _buerger_rc4}


if __name__ == "__main__":
    pass
