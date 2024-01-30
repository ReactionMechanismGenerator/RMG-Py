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
import signal
from copy import deepcopy
from typing import List, Union

import numpy as np
from lpsolve55 import EQ, LE, lpsolve
from pyutilib.common import ApplicationError

from arkane.modelchem import LOT
from rmgpy.molecule import Bond, Molecule
from rmgpy.quantity import ScalarQuantity

# Optional Imports
try:
    import pyomo.environ as pyo
except ImportError:
    pyo = None


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
                spec[0].low_level_hf298.value_si * spec[1]
                for spec in self.species.items()
            )
            - self.target.low_level_hf298.value_si
        )

        target_thermo = (
            sum(
                spec[0].high_level_hf298.value_si * spec[1]
                for spec in self.species.items()
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
                if len(self.connections) == len(other.connections):
                    connections = deepcopy(other.connections)
                    for c in self.connections:
                        for i, c_other in enumerate(connections):
                            if c == c_other:
                                break
                        else:
                            return False
                        connections.pop(i)

                    return True

            return False

        else:
            raise NotImplementedError(
                f"AtomConstraint object has no __eq__ defined for other object of type "
                f"{type(other)}"
            )

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
        rings = species.molecule.get_smallest_set_of_smallest_rings()
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
        enumerated_constraints = []

        # Initialize list of empty lists. Be careful to avoid making references to a singular empty list
        for _ in range(len(full_constraints_list)):
            enumerated_constraints.append([])

        # Begin enumerating constraints
        while True:
            if (
                len(full_constraints_list[0]) == 0
            ):  # We have exhausted target constraints
                if self.limit_scope:  # No need to enumerate any further
                    break  # Out of the while loop

            # Find a constraint to search for next
            for spcs_constraints in full_constraints_list:
                if len(spcs_constraints) != 0:
                    constraint = spcs_constraints[0]
                    break  # Out of the for loop

            else:
                break  # No more constraints, so break out of the while loop

            # enumerate for each species
            new_constraints_list = []
            for i, spcs_constraints in enumerate(full_constraints_list):
                new_spcs_constraints = [c for c in spcs_constraints if c != constraint]
                matching_constraints = len(spcs_constraints) - len(new_spcs_constraints)
                enumerated_constraints[i].append(matching_constraints)
                new_constraints_list.append(new_spcs_constraints)

            # Finally, update the new list
            full_constraints_list = new_constraints_list

        # Finalize list of reference species and corresponding constraints
        reference_constraints = []
        target_constraints = enumerated_constraints[0]
        if self.limit_scope:
            for i, spcs in enumerate(self.all_reference_species):
                # Add 1 to index to offset for the target
                if (
                    len(full_constraints_list[i + 1]) == 0
                ):  # This species does not have extra constraints
                    self.reference_species.append(spcs)
                    reference_constraints.append(enumerated_constraints[i + 1])

        else:
            self.reference_species = self.all_reference_species
            reference_constraints = enumerated_constraints[1:]

        return target_constraints, reference_constraints

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
            else:
                allowable_charges = list(range(0, charge + 1))
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

    def _find_error_canceling_reaction(self, reference_subset, milp_software=None):
        """
        Automatically find a valid error canceling reaction given a subset of the available benchmark species. This
        is done by solving a mixed integer linear programming (MILP) problem similiar to
        Buerger et al. (https://doi.org/10.1016/j.combustflame.2017.08.013)

        Args:
            reference_subset (list): A list of indices from self.reference_species that can participate in the reaction
            milp_software (list, optional): Solvers to try in order. Defaults to ['lpsolve'] or if pyomo is available
                defaults to ['lpsolve', 'pyomo']. lpsolve is usually faster.

        Returns:
            tuple(ErrorCancelingReaction, np.ndarray)
            - Reaction with the target species (if a valid reaction is found, else ``None``)
            - Indices (of the subset) for the species that participated in the return reaction
        """
        if milp_software is None:
            milp_software = ["lpsolve"]
            if pyo is not None:
                milp_software.append("pyomo")

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

        for solver in milp_software:
            if solver == "pyomo":
                # Check that pyomo is available
                if pyo is None:
                    raise ImportError(
                        "Cannot import optional package pyomo. Either install this dependency with "
                        "`conda install -c conda-forge pyomo glpk` or set milp_software to `lpsolve`"
                    )

                # Diable logging, pyomo outputs too often
                logging.disable()

                # Setup the MILP problem using pyomo
                lp_model = pyo.ConcreteModel()
                lp_model.i = pyo.RangeSet(0, m - 1)
                lp_model.j = pyo.RangeSet(0, n - 1)
                lp_model.r = pyo.RangeSet(
                    0, split - 1
                )  # indices before the split correspond to reactants
                lp_model.p = pyo.RangeSet(
                    split, m - 1
                )  # indices after the split correspond to products
                lp_model.v = pyo.Var(
                    lp_model.i, domain=pyo.NonNegativeIntegers
                )  # The stoich. coef. we are solving for
                lp_model.c = pyo.Param(
                    lp_model.i,
                    lp_model.j,
                    initialize=lambda _, i_ind, j_ind: c_matrix[i_ind, j_ind],
                )
                lp_model.s = pyo.Param(
                    lp_model.i, initialize=lambda _, i_ind: sum_constraints[i_ind]
                )
                lp_model.t = pyo.Param(
                    lp_model.j, initialize=lambda _, j_ind: targets[j_ind]
                )

                lp_model.obj = pyo.Objective(rule=_pyo_obj_expression)
                lp_model.constraints = pyo.Constraint(
                    lp_model.j, rule=_pyo_constraint_rule
                )

                # Solve the MILP problem using the GLPK MILP solver (https://www.gnu.org/software/glpk/)
                opt = pyo.SolverFactory("glpk")
                try:
                    results = opt.solve(lp_model, timelimit=10)
                except ApplicationError:
                    continue

                # Return the solution if a valid reaction is found. Otherwise continue to next solver
                if (
                    results.solver.termination_condition
                    == pyo.TerminationCondition.optimal
                ):
                    # Extract the solution and find the species with non-zero stoichiometric coefficients
                    solution = lp_model.v.extract_values().values()
                    break

                # Re-enable logging
                logging.disable(logging.NOTSET)

            elif solver == "lpsolve":
                # Save the current signal handler
                sig = signal.getsignal(signal.SIGINT)

                # Setup the MILP problem using lpsolve
                lp = lpsolve("make_lp", 0, m)
                lpsolve("set_verbose", lp, 2)  # Reduce the logging from lpsolve
                lpsolve("set_obj_fn", lp, sum_constraints)
                lpsolve("set_minim", lp)

                for j in range(n):
                    lpsolve(
                        "add_constraint",
                        lp,
                        np.concatenate((c_matrix[:split, j], -1 * c_matrix[split:, j])),
                        EQ,
                        targets[j],
                    )

                lpsolve(
                    "set_timeout", lp, 10
                )  # Move on if lpsolve can't find a solution quickly

                # All v_i must be integers
                lpsolve("set_int", lp, [True] * m)

                status = lpsolve("solve", lp)

                # Reset signal handling since lpsolve changed it
                try:
                    signal.signal(signal.SIGINT, sig)
                except TypeError:
                    # This is not being run in the main thread, so we cannot reset signal
                    pass

                # Return the solution if a valid reaction is found. Otherwise continue to next solver
                if status == 0:
                    _, solution = lpsolve("get_solution", lp)[:2]
                    break

            else:
                raise ValueError(
                    f"Unrecognized MILP solver {solver} for isodesmic reaction generation"
                )

        else:
            return None, None

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
        self, n_reactions_max=20, milp_software=None
    ):
        """
        Generate multiple error canceling reactions involving the target and a subset of the reference species.

        To do this, a rudimentary search is implemented whereby all possible combinations of the species participating
        in the previously found reaction are excluded from the reference species subset for the next generation process.
        This is implemented using a FIFO queue structure.

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will be returned, after which no
                further searching will occur even if there are possible subsets left in the queue.
            milp_software (list, optional): Solvers to try in order. Defaults to ['lpsolve'] or if pyomo is available
                defaults to ['lpsolve', 'pyomo']. lpsolve is usually faster.

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
                subset, milp_software=milp_software
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

    def calculate_target_enthalpy(self, n_reactions_max=5, milp_software=None):
        """
        Perform a multiple error canceling reactions search and calculate hf298 for the target species by taking the
        median hf298 value from among the error canceling reactions found

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will returned, after which no
                further searching will occur even if there are possible subsets left in the queue.
            milp_software (list, optional): Solvers to try in order. Defaults to ['lpsolve'] or if pyomo is available
                defaults to ['lpsolve', 'pyomo']. lpsolve is usually faster.

        Returns:
            tuple(ScalarQuantity, list)
            - Standard heat of formation at 298 K calculated for the target species
            - reaction list containing all error canceling reactions found
        """
        reaction_list = self.multiple_error_canceling_reaction_search(
            n_reactions_max, milp_software
        )
        if len(reaction_list) == 0:  # No reactions found
            return None, reaction_list
        h298_list = np.zeros(len(reaction_list))

        for i, rxn in enumerate(reaction_list):
            h298_list[i] = rxn.calculate_target_thermo().value_si

        return ScalarQuantity(np.median(h298_list), "J/mol"), reaction_list


def _pyo_obj_expression(model):
    return pyo.summation(model.v, model.s, index=model.i)


def _pyo_constraint_rule(model, col):
    return (
        sum(model.v[row] * model.c[row, col] for row in model.r)
        - sum(model.v[row] * model.c[row, col] for row in model.p)
        == model.t[col]
    )


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
