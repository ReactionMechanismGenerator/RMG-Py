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

"""
This module contains classes and functions for working with polymers.


Polymer chain-length moments (mu0, mu1, mu2)
--------------------------------------------

This module tracks the evolving polymer chain-length distribution using a small
set of “moments”. Consider a mixture of polymer chains, and let N_n denote the
amount (moles) of chains with degree of polymerization n (i.e., n repeat units).

We define the uncentered moments of the distribution as:

```
μ_k = sum_{n>=1} (n^k * N_n)
```

Interpreting the first three moments:

```
μ0 = sum N_n
    Total moles of polymer chains in the mixture.

μ1 = sum (n * N_n)
    Total moles of repeat units across all chains (“monomer units in chains”).
    Dimensionally it is still moles, but it counts repeat units rather than
    chain molecules.

μ2 = sum (n^2 * N_n)
    The second moment. In this module it is sometimes referred to as the
    “width” moment because it increases when the distribution places more
    weight at larger n and/or becomes more broadly spread. Note that μ2 is
    not itself a standard deviation; it is an intermediate quantity from
    which common breadth metrics are computed.
```

Common derived quantities from (μ0, μ1, μ2):

```
Number-average degree of polymerization:
    M_n = μ1 / μ0

Weight-average degree of polymerization:
    M_w = μ2 / μ1

Polydispersity index (distribution breadth):
    Đ = X_w / X_n = (μ2 * μ0) / (μ1^2)

Variance (and standard deviation) of chain length:
    sigma_n^2 = (μ2 / μ0) - (μ1 / μ0)^2
    sigma_n    = sqrt(sigma_n^2)
```

How “width” may be used:
The raw second moment μ2 is useful for conserving information about the
distribution without storing N_n for all n. In practice, users typically
convert μ2 into one of the breadth measures above. For example, Đ is a
widely used scalar summary of how broad the molecular-weight distribution is,
while sigma_n provides an intuitive “spread in n” measure.

In short:
- μ0 tells you how many chains you have (in moles).
- μ1 tells you how many repeat units are contained in those chains.
- μ2 (“width” moment) captures how strongly the population is weighted
  toward larger chains and supports computing Đ and sigma_n.

"""

import datetime
import json
import logging
import os

import numpy as np
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum, IntEnum
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple, Union

from rmgpy.exceptions import InputError
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.fragment import Fragment
from rmgpy.molecule.group import Group, GroupAtom, GroupBond
from rmgpy.molecule.resonance import generate_resonance_structures
from rmgpy.species import Species
from rmgpy.thermo import Wilhoit, NASA, ThermoData


LABELS_1, LABELS_2 = ('1', '*1'), ('2', '*2')


class PolymerCrosslinkError(Exception):
    """
    Raised when a polymer reaction product is a crosslink / chain-coupling
    structure (>2 intact monomer wings, i.e. two chains joined into one).

    The method-of-moments solver tracks a single chain-length distribution per
    pool and has no representation for chain-chain coupling (which removes a
    chain from the population and merges two distributions). Such reactions are
    therefore *rejected* rather than letting the coupled product fall through to
    a spurious gas-phase species, which would silently break the polymer mass
    balance. Caught in :meth:`CoreEdgeReactionModel.make_new_reaction`, which
    discards the reaction.

    Note: this deliberately does NOT subclass ValueError/RuntimeError so it is
    not swallowed by the ``except (RuntimeError, ValueError)`` guard in
    :func:`rmgpy.data.kinetics.family._handshake_structures`; it must propagate
    up to ``make_new_reaction``.
    """
    pass


class Polymer(Species):
    """
    A class representing a polymer distribution (Polymer Pool).
    input initial_mass is in kg, but attributes here are stored in gr units (Mn, Mw in g/mol, initial_mass in gr).
    This object can accept both adjacency list and SMILES to define the monomer,
    yet it is recommended to use adjacency lists for clarity of radical (connectivity) positions
    (to avoid ambiguity due to resonance) and label the sites with '*1' and '*2' to explicitly define head and tail.

    This class handles the definition of a polymer system, including its microstructure
    (monomer/end-groups) and its macroscopic statistical distribution (moments).

    Adopts the "Distinct Species" architecture:
        1. A 'Baseline' homopolymer represents the bulk statistics.
        2. A 'Feature' polymer represents a specific chemical modification (radical, unsaturation, etc.)
           embedded in the bulk chain.

    The thermodynamics and kinetics are derived using "Proxy Species" (Representative Small Molecules)
    constructed by stitching repeating units together.

    Args:
        label (str): Unique identifier for the polymer pool.
        monomer (str): SMILES or adjacency list of the open shell repeating unit with labeled ends.
                       (It refers to the repeating unit mass, not the reactant mass, so for condensation polymers
                       like PET, the monomer is the repeat unit after loss of small molecules like water.)
                       If an adjacency list is used, then '*1' and '*2' labels will bear the meaning of the polymer's
                       head and tail, respectively. Examples:
                       - PE SMILES: "[CH2][CH2]"  (NOT "C=C")
                       - PS SMILES: "[CH2][CH]c1ccccc1" (NOT "C=Cc1ccccc1")
                       - PS adjacency list:
                            multiplicity 3
                            1 *1 C u1 p0 c0 {2,S} {9,S} {10,S}
                            2 *2 C u1 p0 c0 {1,S} {3,S} {11,S}
                            3    C u0 p0 c0 {2,S} {4,S} {8,D}
                            4    C u0 p0 c0 {3,S} {5,D} {12,S}
                            5    C u0 p0 c0 {4,D} {6,S} {13,S}
                            6    C u0 p0 c0 {5,S} {7,D} {14,S}
                            7    C u0 p0 c0 {6,D} {8,S} {15,S}
                            8    C u0 p0 c0 {3,D} {7,S} {16,S}
                            9    H u0 p0 c0 {1,S}
                            10   H u0 p0 c0 {1,S}
                            11   H u0 p0 c0 {2,S}
                            12   H u0 p0 c0 {4,S}
                            13   H u0 p0 c0 {5,S}
                            14   H u0 p0 c0 {6,S}
                            15   H u0 p0 c0 {7,S}
                            16   H u0 p0 c0 {8,S}
        feature_monomer (Optional[Union[Molecule, str]]): The modified repeating unit graph.
                                                          If None, this is a "baseline" homopolymer.
        end_groups (list): List of 2 SMILES strings representing the chain terminals.
                           These must be open-shell species (radicals) that 'cap' the monomer's head and tail.
                           Format: [Initiator_End, Terminator_End]
                           Examples:
                               1. Polystyrene: ['C[C](C)c1ccccc1', '[H]']
                               2. Polyethylene: ['[CH3]', '[CH3]']
                           If None, defaults to ['[H]', '[H]'] (hydrogen-capped).
        cutoff (int): The hybrid threshold (x_s), the chain length (degree of polymerization) where explicit tracking stops.
                       Chains with DP <= cutoff are tracked explicitly, while longer chains are represented statistically.
                       Must be >= 2.
        Mn (float): Number average molecular weight (g/mol).
        Mw (float): Weight average molecular weight (g/mol).
        initial_mass (float): Initial mass in the reactor (kg).
        moments (list[float]): Moments of the chain-length distribution, [μ0, μ1, μ2],
                               representing the distribution: [Moles of chains, moles of units, and width].

    Attributes:
        label (str): Unique identifier.
        monomer (Molecule): The repeating unit structure.
        feature_monomer (Optional[Molecule]): The *modified* repeating unit graph.
        end_groups (List[Molecule]): List of 2 Molecules representing the chain terminals.
        cutoff (int): The chain length where explicit tracking stops.
        Mn (float): Number average molecular weight (g/mol).
        Mw (float): Weight average molecular weight (g/mol).
        initial_mass_g (float): Initial mass in the reactor (gr).
        monomer_mw_g_mol (float): Molecular weight of the monomer in g/mol.
        moments (np.array): [μ0, μ1, μ2] representing the distribution (in Moles).
        k_unzip (float): Rate constant for unzipping reactions (1/s).
        k_scission (float): Rate constant for random scission reactions (1/s).
    """

    def __init__(self,
                 label: str,
                 monomer: Union[Molecule, str],
                 feature_monomer: Optional[Union[Molecule, str]] = None,
                 end_groups: Optional[List[Union[str, Molecule]]] = None,
                 cutoff: int = 4,
                 Mn: Optional[float] = None,
                 Mw: Optional[float] = None,
                 initial_mass: float = 1.0,
                 moments: Optional[List[float]] = None,
                 **kwargs,
                 ):
        # Pop polymer-specific kwargs before passing the rest to Species.__init__
        # — Species does not accept k_unzip/k_scission and would raise TypeError.
        self.k_unzip = kwargs.pop('k_unzip', 0.0)
        self.k_scission = kwargs.pop('k_scission', 0.0)
        # Discreteness threshold (spec 2026-06-10 §6, D7/D8): chains with
        # literal DP < threshold are candidates for discrete tracking. Default
        # 4 = monomer..trimer explicit. DORMANT under the fixed trimer proxy:
        # the backbone gate is proxy-relative (D2), scission routing keys on
        # is_end_group_reaction (D3), and the conditional DP backstop (D8)
        # activates only when the proxy repeat-count exceeds this threshold
        # (a 3-unit proxy never does). Defined-but-documented beats undefined
        # intent; no behavioral use yet.
        self.discrete_dp_threshold = kwargs.pop('discrete_dp_threshold', 4)

        super(Polymer, self).__init__(label=label, **kwargs)

        self.monomer = self._validate_monomer(monomer, label)
        if feature_monomer:
            self.feature_monomer = self._validate_monomer(feature_monomer, label)
        else:
            self.feature_monomer = None
        self._validate_end_groups(end_groups, label)
        self.cutoff = self._validate_cutoff(cutoff, label)
        self.Mn, self.Mw, self.moments = None, None, None

        self.initial_mass_g = initial_mass * 1000.0  # convert to grams
        self.monomer_mw_g_mol = self.monomer.get_molecular_weight() * 1000.0

        self._baseline_proxy = None
        self._feature_proxy = None
        self._fingerprint = None
        self.thermo = None
        self.is_polymer = True
        self.reactive = True
        self._cached_backbone_group = None

        if moments is not None:
            self.moments = np.array(moments, dtype=np.float64)
            self.Mn, self.Mw = self._calculate_distribution_from_moments()
        else:
            self.Mn = Mn
            self.Mw = Mw
            if self.Mn is not None and self.Mw is not None:
                if self.Mn > self.Mw:
                    raise ValueError(f"Polymer '{self.label}': Physically impossible distribution (Mn > Mw).\n"
                                     f"Got Mn={self.Mn}, Mw={self.Mw}.")
                if self.Mn <= 0 or self.Mw <= 0:
                    raise ValueError(f"Polymer '{self.label}': Molecular weights must be positive.\n"
                                     f"Got Mn={self.Mn}, Mw={self.Mw}.")
                self.moments = self._calculate_moments_from_distribution()
            else:
                raise InputError(f"Polymer '{label}': Must provide either 'moments' OR ('Mn' and 'Mw').\n"
                                 f"Got moments={moments}, Mn={Mn}, Mw={Mw}.")
        reactive_proxy = self.feature_proxy or self.baseline_proxy
        self.molecule = reactive_proxy.molecule

    def __repr__(self):
        return f"<Polymer '{self.label}' Mn={self.Mn:.1f} Mw={self.Mw:.1f} Cutoff={self.cutoff}>"

    @property
    def multiplicity(self):
        """
        Return the multiplicity of the reactive proxy.
        (e.g., if the polymer has a radical end, multiplicity is 2).
        """
        return self.get_proxy_species().multiplicity

    @property
    def molecular_weight(self):
        """
        Returns the Number Average Molecular Weight (Mn) of the distribution.
        We deliberately do NOT return Mn here. Returning Mn would cause
        inconsistency with the 'Per-Site' thermo values.
        """
        return self.get_proxy_species().molecular_weight

    @molecular_weight.setter
    def molecular_weight(self, value):
        """
        Setter required for compatibility with Species.__init__.
        The value is ignored because Polymer MW is dynamically derived from the Proxy.
        """
        pass

    def generate_resonance_structures(self, keep_isomorphic=True, filter_structures=True, save_order=False):
        """
        Override to keep Polymer.molecule in sync with the proxy.

        The base Species.generate_resonance_structures can replace self.molecule
        with a new list, breaking the shared reference with the proxy.  By
        delegating to the proxy and re-binding afterwards, the proxy's internal
        consistency is preserved.
        """
        proxy = self.get_proxy_species()
        if proxy is not None:
            proxy.generate_resonance_structures(
                keep_isomorphic=keep_isomorphic,
                filter_structures=filter_structures,
                save_order=save_order,
            )
            self.molecule = proxy.molecule
        else:
            super().generate_resonance_structures(
                keep_isomorphic=keep_isomorphic,
                filter_structures=filter_structures,
                save_order=save_order,
            )

    def get_symmetry_number(self):
        """Delegates symmetry calculation to the proxy species."""
        proxy = self.get_proxy_species()
        try:
            return proxy.get_symmetry_number()
        except (ValueError, KeyError):
            # Proxy resonance structures may be inconsistent (e.g. stale atom
            # IDs after shared-reference breakage).  Regenerate from the first
            # molecule and retry.
            mol = proxy.molecule[0].copy(deep=True)
            mol.assign_atom_ids()
            proxy.molecule = mol.generate_resonance_structures(
                keep_isomorphic=True, filter_structures=True, save_order=True,
            )
            self.molecule = proxy.molecule
            proxy.symmetry_number = -1  # reset cached value
            return proxy.get_symmetry_number()

    def get_net_charge(self):
        """Delegates charge calculation to the proxy."""
        return self.get_proxy_species().get_net_charge()

    @property
    def fingerprint(self):
        """Fingerprint of this polymer, taken from molecule attribute. Read-only."""
        if self._fingerprint is None:
            if self.monomer:
                feat = f'_Feat-{self.feature_monomer.fingerprint}' if self.feature_monomer else ''
                eg = '_'.join(eg.fingerprint for eg in self.end_groups) if self.end_groups else ''
                self._fingerprint = f'Polymer_{self.monomer.fingerprint}{feat}_EG-{eg}_{self.cutoff}'
        return self._fingerprint

    @property
    def baseline_proxy(self) -> Species:
        """Returns the cached '<Head-Baseline-Tail>' trimer species."""
        if self._baseline_proxy is None:
            self._baseline_proxy = self._stitch_trimer(self.monomer)
        return self._baseline_proxy

    @property
    def feature_proxy(self) -> Species:
        """Returns the cached '<Head-Feature-Tail>' trimer species. (None if homopolymer)"""
        if self.feature_monomer and self._feature_proxy is None:
            self._feature_proxy = self._stitch_trimer(self.feature_monomer)
        return self._feature_proxy

    @property
    def backbone_group(self) -> 'Group':
        """
        Generates and caches a relaxed 'Backbone Pattern' monomer Group for subgraph searching.

        This group is designed to identify the monomer's structural skeleton within a
        larger, potentially reacted polymer chain. It purposefully relaxes electronic
        and bonding constraints to ensure the backbone matches regardless of local
        chemical environment changes (e.g., resonance, charging, or neighbor effects).

        The generation process:
        1.  **Strip Identity:** Removes monomer labels (*1, *2) and radical electrons
            to represent the 'bound' state of the monomer in a chain.
        2.  **Heal Perception:** Updates atom types to prevent RMG from misinterpreting
            the stripped radicals as carbocations or lone pairs.
        3.  **Relax Electronic Constraints:** Sets wildcards for atom types, charges,
            and lone pairs. This ensures a match based purely on element connectivity.
        4.  **Strict Radical Constraint:** Enforces `radical_electrons=[0]`. This is
            critical: it ensures the pattern matches *only* unreacted buffer monomers
            and fails to match monomers that have become active centers (radicals).
        5.  **Relax Bond Orders:** Allows Single, Double, Triple, or Benzene bonds
            between backbone atoms to account for potential resonance delocalization.

        Returns:
            Group: A relaxed Group object instance representing the monomer backbone.
        """
        if getattr(self, '_cached_backbone_group', None) is not None:
            return self._cached_backbone_group

        pat_mol = self.monomer.copy(deep=True)
        pat_mol.clear_labeled_atoms()

        # 1. Strip Radicals
        for atom in pat_mol.atoms:
            atom.radical_electrons = 0

        # 2. Assign base atomtypes (e.g., Cs for backbone, Cd for Kekule ring)
        pat_mol.update_atomtypes()

        # 3. Strip hydrogens to create a pure heavy-atom skeleton matcher
        for atom in pat_mol.atoms[:]:
            if atom.is_hydrogen():
                pat_mol.remove_atom(atom)

        # 4. Convert to Group
        group = pat_mol.to_group()
        group.multiplicity = None

        # 5. Relax Constraints
        for g_atom in group.atoms:

            # Make Kekulé and Aromatic types interchangeable, but strictly preserve aliphatic (Cs) types.
            expanded_types = set(g_atom.atomtype)
            for at in g_atom.atomtype:
                if at.label in ('Cd', 'Cb', 'Cbf', 'Cdd'):
                    expanded_types.update([ATOMTYPES['Cd'], ATOMTYPES['Cb'], ATOMTYPES['Cbf'], ATOMTYPES['Cs']])
            g_atom.atomtype = list(expanded_types)

            g_atom.radical_electrons = [0]  # Strict unreacted
            g_atom.charge = []
            g_atom.lone_pairs = []

        # 6. Relax Bond Orders
        all_orders = [1, 1.5, 2, 3]
        seen_bonds = set()
        for g_atom in group.atoms:
            for neighbor, g_bond in g_atom.bonds.items():
                key = tuple(sorted((id(g_atom), id(neighbor))))
                if key in seen_bonds:
                    continue
                seen_bonds.add(key)
                g_bond.order = all_orders

        self._cached_backbone_group = group
        return group

    def copy(self, deep=True):
        """
        Create a copy of the current species. If the
        kw argument 'deep' is True, then a deep copy will be made of the
        Molecule objects in self.molecule.

        For other complex attributes, a deep copy will always be made.
        """
        other = Polymer.__new__(Polymer)
        other.index = self.index
        other.label = self.label
        other.thermo = deepcopy(self.thermo)
        other.monomer = self.monomer.copy(deep=deep)
        other.feature_monomer = self.feature_monomer.copy(deep=deep) if self.feature_monomer else None
        other.end_groups = list()
        for eg in self.end_groups:
            other.end_groups.append(eg.copy(deep=deep))
        other.molecule = [m.copy(deep=deep) for m in self.molecule]
        other.conformer = deepcopy(self.conformer)
        other.transport_data = deepcopy(self.transport_data)
        other.energy_transfer_model = deepcopy(self.energy_transfer_model)
        other.reactive = self.reactive
        other.props = deepcopy(self.props)
        other.cutoff = self.cutoff
        other.Mn = self.Mn
        other.Mw = self.Mw
        other.moments = self.moments
        other.initial_mass_g = self.initial_mass_g
        other.monomer_mw_g_mol = self.monomer_mw_g_mol
        other._baseline_proxy = self._baseline_proxy.copy(deep=True) if self._baseline_proxy else None
        other._feature_proxy = self._feature_proxy.copy(deep=True) if self._feature_proxy else None
        other._fingerprint = self._fingerprint
        # Attributes set in __init__ that __new__ bypasses — must be carried over,
        # else a copied Polymer loses its identity flag (is_polymer) and, worse,
        # its degradation kinetics (k_scission/k_unzip would silently reset to 0).
        other.k_scission = self.k_scission
        other.k_unzip = self.k_unzip
        other.discrete_dp_threshold = getattr(self, 'discrete_dp_threshold', 4)
        other.is_polymer = True
        other._cached_backbone_group = None
        return other

    @staticmethod
    def _validate_monomer(monomer: Union[Molecule, str],
                          label: str,
                          ) -> Molecule:
        """
        Ensures monomer has labels '*1' and '*2' for connectivity.
        Assigns random labels to radical sites if labels are missing.

        Args:
            monomer (Union[Molecule, str]): The monomer to validate.
            label (str): The polymer label (for error messages).
        """
        if isinstance(monomer, str):
            if len(monomer.splitlines()) > 1:
                mol = Molecule().from_adjacency_list(monomer)
            else:
                mol = Molecule(smiles=monomer)
        elif isinstance(monomer, Molecule):
            mol = monomer
        else:
            raise InputError(f"Polymer '{label}': Monomer must be a SMILES string or Molecule object.\n"
                             f"Got {monomer} of type {type(monomer)}.")
        mol.is_polymer_proxy = True
        has_1 = any(mol.contains_labeled_atom(x) for x in LABELS_1)
        has_2 = any(mol.contains_labeled_atom(x) for x in LABELS_2)
        if not (has_1 and has_2) and mol.get_radical_count() < 2:
            raise InputError(f"Polymer '{label}': Monomer must have 2 reactive sites. Please use labels *1 and *2.")
        if not has_1 or not has_2:
            for atom in mol.atoms:
                if atom.radical_electrons >= 1:
                    if not has_1:
                        atom.label = '*1'
                        has_1 = True
                    elif not has_2:
                        atom.label = '*2'
                        has_2 = True
                if has_1 and has_2:
                    break
        i_1, i_2 = find_labeled_atom(mol, LABELS_1), find_labeled_atom(mol, LABELS_2)
        if i_1 is None or i_2 is None or i_1 == i_2:
            raise InputError(f"Polymer '{label}': failed to define distinct '*1' and '*2' sites on monomer.")
        for idx in [i_1, i_2]:
            atom = mol.atoms[idx]
            if atom.radical_electrons == 0:
                atom.increment_radical()
        mol.update_multiplicity()
        return mol

    def _validate_end_groups(self,
                             end_groups: Optional[List[Union[Molecule, str]]],
                             label: str,
                             ):
        """
        Ensures end groups (head and tail) are valid radicals with proper labels (*1 for head, *2 for tail).
        """
        if not end_groups:
            end_groups = ['[H]', '[H]']
        if not isinstance(end_groups, (list, tuple)) or len(end_groups) != 2:
            raise InputError(f"Polymer '{label}': Must provide exactly 2 end groups.\nGot: {end_groups}")
        validated_groups = list()
        for i, eg in enumerate(end_groups):
            if isinstance(eg, str):
                mol = Molecule(smiles=eg)
            elif isinstance(eg, Molecule):
                mol = eg
            else:
                raise InputError(f"Polymer '{label}': End group #{i + 1} is invalid, got: {eg} of type {type(eg)}.")
            if mol.get_radical_count() != 1:
                raise InputError(f"Polymer '{label}': End group #{i + 1} ('{eg}') is chemically inert. "
                                 "End groups must be mono-radicals.")
            bad = [a.label for a in mol.atoms if a.label and a.label not in (*LABELS_1, *LABELS_2)]
            if bad:
                raise InputError(f"Polymer '{label}': end-group has invalid labels {bad}; only '*1'/'*2' are allowed.")
            validated_groups.append(mol)
        head_mol, tail_mol = validated_groups
        if find_labeled_atom(head_mol, LABELS_1) is None:
            for atom in head_mol.atoms:
                if atom.radical_electrons == 1:
                    atom.label = '*1'
                    break
            if find_labeled_atom(head_mol, LABELS_1) is None:
                raise InputError(f"Polymer '{label}': could not assign '*1' label on head end-group.")
        if find_labeled_atom(tail_mol, LABELS_2) is None:
            for atom in tail_mol.atoms:
                if atom.radical_electrons == 1:
                    atom.label = '*2'
                    break
            if find_labeled_atom(tail_mol, LABELS_2) is None:
                raise InputError(f"Polymer '{label}': could not assign '*2' label on tail end-group.")
        self.end_groups = validated_groups

    def _validate_cutoff(self, cutoff: int, label: str) -> int:
        """Validates that cutoff is an integer >= 2."""
        try:
            cutoff_i = int(cutoff)
        except (ValueError, TypeError):
            raise InputError(f"Polymer '{label}': Cutoff must be an integer. Got cutoff={cutoff}.")
        if cutoff_i < 2:
            raise InputError(f"Polymer '{label}': Cutoff must be at least 2. Got cutoff={cutoff_i}.")
        return int(cutoff_i)

    def _calculate_moments_from_distribution(self) -> np.ndarray:
        """
        Calculates raw moments [μ0, μ1, μ2] in MOLES.

        BASIS NOTE: These moments are tracked on a 'Monomer Unit' basis (Degree of Polymerization),
        not a 'Carbon Atom' basis as done in Vermeire et al., 2025 (https://doi.org/10.1016/j.cej.2025.159455).
        This is consistent with RMG's graph representation.

        Moment Definitions:
        - μ1 (Total Moles of Monomer Units): Derived from initial mass.
             μ1 = Mass(g) / MonomerMW(g/mol)
        - μ0 (Total Moles of Chains): Derived from Number Average Degree of Polymerization (DPn).
             μ0 = μ1 / DPn
        - μ2 (Width Parameter): Derived from Weight Average Degree of Polymerization (DPw).
             μ2 = μ1 * DPw
        """
        DPn = self.Mn / self.monomer_mw_g_mol
        DPw = self.Mw / self.monomer_mw_g_mol
        mu1 = self.initial_mass_g / self.monomer_mw_g_mol
        mu0 = mu1 / DPn
        mu2 = mu1 * DPw
        return np.array([mu0, mu1, mu2])

    def get_closing_moment(self, moments: Optional[Union[List[float], np.ndarray]] = None) -> float:
        """
        Calculates the 3rd moment (μ3) needed to close the moment equations.
        This method isn't used in the solver directly, it is kept here for post-processing.

        Uses Log-Lagrange extrapolation (Vermeire et al., 2025, 2025, https://doi.org/10.1016/j.cej.2025.159455):
            ln(μ3) = 3 ln(μ2) - 3 ln(μ1) + ln(μ0)
            which implies: μ3 = μ0 * (μ2 / μ1)^3

        Lagrange extrapolation allows the data (μ0, μ1, μ2) to define the curvature of the distribution tail.
        Operating in Log-Space (ln μ_k) is critical. Moments span large orders of magnitude.

        Args:
            moments (np.ndarray): Current [μ0, μ1, μ2].
                                  If None, uses self.moments (initial state).

        Returns:
            float: The estimated μ3. Returns 0.0 if input moments are non-positive (e.g. empty reactor or solver noise).
        """
        if moments is None:
            moments = self.moments
        moments = np.asarray(moments, dtype=float)
        if moments.shape[0] < 3:
            raise ValueError(f"get_closing_moment expected array of length >= 3, got shape {moments.shape}")
        mu0, mu1, mu2 = moments[:3]

        if mu0 <= 1e-20 or mu1 <= 1e-20 or mu2 <= 1e-20:
            return 0.0

        if mu1 < mu0:
            # Unrealizable state (μ1 ≥ μ0 always holds for a k≥1 distribution).
            # Kept consistent with the solver's _safe_mu3_from_mu012 guard so
            # post-processing never amplifies an out-of-cone moment vector.
            return 0.0

        # Log-Lagrange Extrapolation: ln(mu3) = 3*ln(mu2) - 3*ln(mu1) + ln(mu0),
        # assuming that the 'curvature' of the distribution in log-space is constant.
        try:
            with np.errstate(divide='raise', invalid='raise', over='raise'):
                ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
                if ln_mu3 > 709.0:  # float64 exp overflow threshold ~709
                    return float('inf')
                if ln_mu3 < -745.0:  # exp underflow to 0
                    return 0.0
                return float(np.exp(ln_mu3))
        except FloatingPointError:
            return 0.0

    def _calculate_distribution_from_moments(self) -> Tuple[Optional[float], Optional[float]]:
        """
        Back-calculates (Mn, Mw) in g/mol from moments.
        """
        if self.moments is None:
            return None, None
        mu0, mu1, mu2 = self.moments
        if mu0 == 0 or mu1 == 0:
            return 0.0, 0.0
        DPn = mu1 / mu0
        DPw = mu2 / mu1
        Mn = DPn * self.monomer_mw_g_mol
        Mw = DPw * self.monomer_mw_g_mol
        return Mn, Mw

    def get_polydispersity(self):
        """Returns PDI (Mw/Mn). Useful for sanity checks and reporting."""
        if self.Mn is None or self.Mn == 0:
            return 0.0
        return self.Mw / self.Mn

    def is_identical(self, other, strict=True):
        """
        Return ``True`` if `other` is a Polymer object and is chemically identical
        to this one (based on their proxy species).

        If ``strict=False``, performs the check ignoring electrons and resonance structures.
        """
        if isinstance(other, Polymer):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.is_identical(molecule2, strict=strict):
                        return True
        return False

    def is_isomorphic(self, other: Union['Polymer', Species, Molecule, Fragment],
                      generate_initial_map=False,
                      save_order=True,
                      strict=True,
                      ) -> bool:
        """
        Return ``True`` if the species is isomorphic to `other`, which can be
        either a :class:`Molecule` object or a :class:`Species` object.

        Args:
            other (Union[Species, Molecule, Fragment]): The other species or molecule to compare to.
            generate_initial_map (bool, optional): If ``True``, make initial map by matching labeled atoms
            save_order (bool, optional):           if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):               If ``False``, perform isomorphism ignoring electrons.

        Returns:
            bool: ``True`` if the species is isomorphic to `other`, ``False`` otherwise.
        """
        if isinstance(other, (Polymer, Species)):
            for mol_1 in self.molecule:
                for mol_2 in other.molecule:
                    if mol_1.copy(deep=True).is_isomorphic(mol_2.copy(deep=True),
                                                           generate_initial_map=generate_initial_map,
                                                           save_order=save_order,
                                                           strict=strict):
                        return True
        elif isinstance(other, Molecule):
            for mol in self.molecule:
                mol_clean = mol.copy(deep=True)
                mol_clean.clear_labeled_atoms()
                if mol_clean.copy(deep=True).is_isomorphic(other.copy(deep=True),
                                                           generate_initial_map=generate_initial_map,
                                                           save_order=save_order,
                                                           strict=strict):
                    return True
            return False
        elif isinstance(other, Fragment):
            return False
        else:
            raise ValueError(f'Unexpected value "{other}" of type {type(other)} for other parameter;'
                             ' should be a Polymer/Species/Molecule object.')
        return False

    def get_proxy_species(self, mode: str = 'auto') -> Optional[Species]:
        """
        Public accessor to get the representative small molecule for this polymer.

        Args:
            mode: 'auto' (returns feature if exists, else baseline),
                  'baseline' (forces baseline trimer),
                  'feature' (forces feature trimer if available, else None).
        """
        if mode == 'baseline':
            return self.baseline_proxy
        elif mode == 'feature':
            return self.feature_proxy
        else:
            return self.feature_proxy or self.baseline_proxy

    def _stitch_trimer(self, center_unit: Molecule) -> Optional[Species]:
        """
        Constructs a capped trimer proxy (3 repeat units) with end-groups:
            [HeadCap *1] -– [*2 Baseline *1] -– [*2 Center *1] –- [*2 Baseline *1] –- [*2 TailCap]
        where Baseline/Center are repeat units and HeadCap/TailCap are terminal end-groups.

        Args:
            center_unit (Molecule): The monomer unit to place in the center.
                                    Either baseline (original) or feature (reacted).

        Returns:
            Optional[Species]: The stitched trimer species.
        """
        baseline = self.monomer.copy(deep=True)
        center = center_unit.copy(deep=True)
        head = self.end_groups[0].copy(deep=True)
        tail = self.end_groups[1].copy(deep=True)

        trimer = stitch_molecules_by_labeled_atoms(head, baseline)
        if trimer is None: return None
        trimer = stitch_molecules_by_labeled_atoms(trimer, center)
        if trimer is None: return None
        trimer = stitch_molecules_by_labeled_atoms(trimer, baseline)
        if trimer is None: return None
        trimer = stitch_molecules_by_labeled_atoms(trimer, tail)
        if trimer is None: return None

        trimer.update()
        trimer.identify_ring_membership()
        spc = Species(molecule=[trimer])
        mol_0 = spc.molecule[0].copy(deep=True)
        mol_0.clear_labeled_atoms()
        mol_0.assign_atom_ids()
        spc.molecule = generate_resonance_structures(mol_0,
                                                     clar_structures=False,
                                                     keep_isomorphic=False,
                                                     filter_structures=True,
                                                     save_order=True)
        spc.is_polymer_proxy = True
        return spc

    def create_reacted_copy(self, reacted_proxy: Molecule) -> Optional['Polymer']:
        """
        Wrapper that ensures any generated polymer fragment is sanitized
        (labels stripped, proxy tagged) before returning to the RMG engine.

        Raises:
            PolymerCrosslinkError: if the product is a crosslink / chain-coupling
                structure (>2 intact wings). Chain-chain coupling is not
                representable in the method-of-moments model, so the caller
                discards the whole reaction instead of leaking the coupled
                product into the gas phase as a spurious small molecule.
        """
        # Guard: reject chain-chain coupling (crosslink) products up front.
        # Without this they fall through _create_reacted_copy_logic to None and
        # get silently registered as gas-phase species, breaking mass balance.
        probe = reacted_proxy.copy(deep=True)
        probe.clear_labeled_atoms()
        probe.update()
        klass, _ = classify_structure(Species(molecule=[probe]), self)
        if klass == PolymerClass.CROSSLINK:
            raise PolymerCrosslinkError(
                f"Reaction product is a crosslink/chain-coupling structure for "
                f"pool '{self.label}'; chain-chain coupling is not representable "
                f"in the method-of-moments model, so the reaction is rejected."
            )

        if klass == PolymerClass.END_MOD:
            # End-group modification (e.g. terminal radical activation, CH3->CH2.)
            # leaves the chain intact: the degree of polymerization, and therefore
            # the moment-tracked chain-length distribution, is unchanged. The
            # method-of-moments model abstracts chain-end activation into k_unzip
            # (dmu1/dt = -k_unzip*mu0, scaling with the chain-end count mu0), so we
            # fold the product back into the parent pool with moments and mass
            # preserved rather than spawning a distinct activated-chain population
            # (which would double-count k_unzip and cannot be represented anyway,
            # since an activated end-cap is di-radical: stitch site + activation).
            #
            # Without this, _create_reacted_copy_logic's raw wing-matching diverges
            # from classify_structure's heavy-view matcher, mis-routes the END_MOD
            # product into a scission branch where it fails the mono-radical
            # end-group assertion, returns None, and the product leaks to the gas
            # phase as a spurious small molecule (a mass-balance leak).
            new_poly = self.copy(deep=True)
        else:
            new_poly = self._create_reacted_copy_logic(reacted_proxy)
        if new_poly is None:
            return None
        # Stamp the classification verdict so the polymer handshake can flag
        # END_MOD reactions for chain-end (mu0) scaling in the solver. Read by
        # is_end_group_reaction(products); a transient generation-time marker.
        new_poly._reacted_class = klass
        # Keep the sanitized reacted fragment so chip product surgery
        # (surge_chip_products, spec 2026-06-10 §4.2) can demote a SCISSION
        # chip back to a discrete Molecule and size chips by MW. Transient
        # generation-time marker like _reacted_class (deliberately not carried
        # by Polymer.copy()).
        new_poly._source_molecule = probe
        proxy_spec = new_poly.get_proxy_species()
        for mol in proxy_spec.molecule:
            mol.clear_labeled_atoms()
            mol.is_polymer_proxy = True
            mol.reactive = True
        return new_poly

    def _create_reacted_copy_logic(self, reacted_proxy: Molecule) -> Optional['Polymer']:
        """
        Creates a new Polymer species from a reacted proxy fragment.
        Handles both modification (intact chain) and scission (broken chain).

        Args:
            reacted_proxy (Molecule): A single product fragment from the reaction.

        Returns:
            Optional['Polymer']: A new Polymer species with updated feature or end-groups.
                                 Returns None if for an invalid/unsupported polymer fragment.
        """
        product = reacted_proxy.copy(deep=True)
        product.clear_labeled_atoms()
        product.update()
        if self.baseline_proxy.is_isomorphic(product):
            return self.copy(deep=True)

        def _count_boundary_edges(m):
            """Counts bonds leaving the matched subgraph."""
            atom_set = get_target_atoms(m)
            cuts = 0
            for a in atom_set:
                for nbr in a.bonds:
                    if nbr not in atom_set:
                        cuts += 1
            return cuts

        head_groups = self._wing_groups("head")
        tail_groups = self._wing_groups("tail")
        head_matches, tail_matches = list(), list()
        for g in head_groups:
            head_matches.extend(product.find_subgraph_isomorphisms(g, save_order=True))
        for g in tail_groups:
            tail_matches.extend(product.find_subgraph_isomorphisms(g, save_order=True))
        head_atoms, tail_atoms = set(), set()
        if head_matches and tail_matches:
            best_pair = None
            best_score = (999, 999, 999)
            for hm in head_matches:
                ha = get_target_atoms(hm)
                h_bnd = _count_boundary_edges(hm)
                for tm in tail_matches:
                    ta = get_target_atoms(tm)
                    t_bnd = _count_boundary_edges(tm)
                    if ha.isdisjoint(ta):
                        score = (abs(h_bnd - 1) + abs(t_bnd - 1), h_bnd + t_bnd, -(len(ha) + len(ta)))
                        if score < best_score:
                            best_score = score
                            best_pair = (ha, ta)

            if best_pair:
                head_atoms, tail_atoms = best_pair
            else:
                def _score_single(m):
                    b = _count_boundary_edges(m)
                    return (abs(b - 1), b, -len(get_target_atoms(m)))
                best_head = min(head_matches, key=_score_single)
                best_tail = min(tail_matches, key=_score_single)
                if _score_single(best_head) <= _score_single(best_tail):
                    head_atoms = get_target_atoms(best_head)
                else:
                    tail_atoms = get_target_atoms(best_tail)

        elif head_matches:
            def _score_single(m):
                b = _count_boundary_edges(m)
                return abs(b - 1), b, -len(get_target_atoms(m))
            best_head = min(head_matches, key=_score_single)
            head_atoms = get_target_atoms(best_head)

        elif tail_matches:
            def _score_single(m):
                b = _count_boundary_edges(m)
                return abs(b - 1), b, -len(get_target_atoms(m))
            best_tail = min(tail_matches, key=_score_single)
            tail_atoms = get_target_atoms(best_tail)

        if head_atoms and tail_atoms:
            if self.baseline_proxy.is_isomorphic(product):
                return self.copy(deep=True)
            atoms_to_remove = head_atoms | tail_atoms
            try:
                new_feature_graph = self._extract_remainder(product, atoms_to_remove)
            except ValueError:
                return None
            self._restore_labels(new_feature_graph,
                                 original_mol=product,
                                 removed_atoms=atoms_to_remove,
                                 head_match_atoms=head_atoms,
                                 tail_match_atoms=tail_atoms)
            try:
                self._assert_feature_unit(new_feature_graph)
            except ValueError:
                return None
            return Polymer(
                label=f"{self.label}_mod",
                monomer=self.monomer,
                feature_monomer=new_feature_graph,
                end_groups=[eg.copy(deep=True) for eg in self.end_groups],
                cutoff=self.cutoff,
                Mn=None if self.moments is not None else self.Mn,
                Mw=None if self.moments is not None else self.Mw,
                moments=self.moments.tolist() if self.moments is not None else None,
                initial_mass=self.initial_mass_g / 1000.0,
            )

        if head_atoms:
            try:
                new_tail = self._extract_remainder(product, head_atoms)
            except ValueError:
                return None
            self._restore_labels(new_tail,
                                 original_mol=product,
                                 removed_atoms=head_atoms,
                                 head_match_atoms=head_atoms,
                                 tail_match_atoms=None)
            try:
                self._assert_end_group(new_tail, want_label='*2')
            except ValueError:
                return None
            # A scission product is a NEW, shorter chain population that starts
            # EMPTY (initial_mass=0 -> zero moments) and is filled by reaction
            # flux during the run; random scission ~halves the chain length, so
            # its Mn/Mw are halved. Seed identically to _scission_head below.
            # (Previously this copied the PARENT's full moments, which both
            # duplicated mass — contradicting initial_mass=0 — and discarded the
            # Mn/Mw halving, since passing `moments` makes __init__ derive Mn/Mw
            # from it and ignore the halved values.)
            new_Mn = self.Mn / 2.0 if self.Mn else None
            new_Mw = self.Mw / 2.0 if self.Mw else None
            return Polymer(label=f"{self.label}_scission_tail",
                           monomer=self.monomer,
                           feature_monomer=None,
                           end_groups=[self.end_groups[0].copy(deep=True), new_tail],
                           cutoff=self.cutoff,
                           Mn=new_Mn,
                           Mw=new_Mw,
                           initial_mass=0.0,
                           moments=None,
                           )

        if tail_atoms:
            try:
                new_head = self._extract_remainder(product, tail_atoms)
            except ValueError:
                return None
            self._restore_labels(new_head,
                                 original_mol=product,
                                 removed_atoms=tail_atoms,
                                 head_match_atoms=None,
                                 tail_match_atoms=tail_atoms)
            try:
                self._assert_end_group(new_head, want_label='*1')
            except ValueError:
                return None
            new_Mn = self.Mn / 2.0 if self.Mn else None
            new_Mw = self.Mw / 2.0 if self.Mw else None
            return Polymer(label=f"{self.label}_scission_head",
                           monomer=self.monomer,
                           feature_monomer=None,
                           end_groups=[new_head, self.end_groups[1].copy(deep=True)],
                           cutoff=self.cutoff,
                           Mn=new_Mn,
                           Mw=new_Mw,
                           initial_mass=0.0,
                           moments=None)

        return None

    def _stitch_wing(self, side: str) -> Molecule:
        """
        Helper to build [Head-Base] or [Base-Tail] for subtraction.

        Args:
            side (str): 'head' or 'tail' wing.

        Returns:
            Molecule: The stitched wing molecule.
        """
        baseline = self.monomer.copy(deep=True)
        if side == 'head':
            end = self.end_groups[0].copy(deep=True)
            return stitch_molecules_by_labeled_atoms(end, baseline)
        else:
            end = self.end_groups[1].copy(deep=True)
            return stitch_molecules_by_labeled_atoms(baseline, end)

    def _wing_groups(self, side: str) -> list['Group']:
        """
        Return wing patterns (an End group stitched to one Monomer) as Groups,
        using resonance structures to get Clar & Kekulé (S/D/B) variants.
        Prune duplicate Group patterns by adjacency-list string.
        Then relax radical constraints so an open-shell wing can match a closed-shell proxy.

        Args:
            side (str): 'head' or 'tail' wing.

        Returns:
            list[Group]: A list of Group objects representing the wing patterns.
        """
        wing = self._stitch_wing(side=side)
        for a in wing.atoms:
            a.label = ''
        wing.update_multiplicity()
        spc = Species(molecule=[wing])
        mol_0 = spc.molecule[0].copy(deep=True)
        molecules = generate_resonance_structures(mol_0,
                                                  clar_structures=False,
                                                  keep_isomorphic=False,
                                                  filter_structures=True,
                                                  save_order=True)
        uniq = dict()
        for m in molecules:
            g = m.to_group()
            g.multiplicity = []
            for ga in g.atoms:
                expanded_types = set(ga.atomtype)
                for at in ga.atomtype:
                    if at.label in ('Cd', 'Cb', 'Cbf', 'Cdd'):
                        expanded_types.update([ATOMTYPES['Cd'], ATOMTYPES['Cb'], ATOMTYPES['Cbf']])
                ga.atomtype = list(expanded_types)
                ga.radical_electrons = []
                ga.charge = []
                ga.lone_pairs = []
                for gb in ga.bonds.values():
                    gb.order = [1, 1.5, 2, 3]
            g.update()
            key = g.to_adjacency_list()
            if key not in uniq:
                uniq[key] = g
        return list(uniq.values())

    @staticmethod
    def _extract_remainder(complex_mol: Molecule, atoms_to_remove) -> Molecule:
        """
        Creates a new Molecule containing only the atoms NOT in 'atoms_to_remove'.
        Strips labels on the copied atoms so downstream label restoration is deterministic.

        Args:
            complex_mol (Molecule): The original complex molecule.
            atoms_to_remove (set): Set of atoms to exclude.

        Returns:
            Molecule: The extracted remainder molecule.
        """
        atoms_to_remove = set(atoms_to_remove)
        if any(a not in complex_mol.atoms for a in atoms_to_remove):
            raise ValueError("atoms_to_remove contains atoms not in complex_mol (wrong object identity / wrong copy).")
        remainder = Molecule()
        old_to_new_map = dict()
        for atom in complex_mol.atoms:
            if atom not in atoms_to_remove:
                new_atom = atom.copy()
                new_atom.label = ''
                remainder.add_atom(new_atom)
                old_to_new_map[atom] = new_atom
        if not remainder.atoms:
            raise ValueError("Polymer extraction failed: No atoms remained after wing removal.")
        added_bonds = set()
        for old_atom in complex_mol.atoms:
            if old_atom not in old_to_new_map:
                continue
            for bonded_neighbor, bond in old_atom.bonds.items():
                if bonded_neighbor not in old_to_new_map:
                    continue
                bond_key = frozenset((old_atom, bonded_neighbor))
                if bond_key in added_bonds:
                    continue
                remainder.add_bond(Bond(old_to_new_map[old_atom],
                                        old_to_new_map[bonded_neighbor],
                                        order=bond.order))
                added_bonds.add(bond_key)
        remainder.update_multiplicity()
        return remainder

    @staticmethod
    def _restore_labels(new_mol: Molecule,
                        original_mol: Molecule,
                        removed_atoms: set,
                        head_match_atoms=None,
                        tail_match_atoms=None,
                        ):
        """
        Identifies cut sites and restores *1/*2 labels and radical character.
        Updates multiplicity after changes.

        Args:
            new_mol (Molecule): The extracted remainder molecule to modify.
            original_mol (Molecule): The original complex molecule.
            removed_atoms (set): Set of atoms that were removed from original_mol.
            head_match_atoms (set, optional): Set of atoms in original_mol that matched the head wing.
            tail_match_atoms (set, optional): Set of atoms in original_mol that matched the tail wing.
        """
        orig_to_new_map = dict()
        new_atom_iter = iter(new_mol.atoms)
        kept_count = 0
        for atom in original_mol.atoms:
            if atom not in removed_atoms:
                kept_count += 1
                try:
                    orig_to_new_map[atom] = next(new_atom_iter)
                except StopIteration:
                    raise ValueError("Mapping failure: new_mol has fewer atoms than expected.")
        if kept_count != len(new_mol.atoms):
            raise ValueError(f"Mapping failure: new_mol has {len(new_mol.atoms)} atoms, expected {kept_count}.")
        if len(orig_to_new_map) != len(new_mol.atoms):
            raise ValueError("Mapping failure: could not map all kept atoms into new_mol.")
        for atom in original_mol.atoms:
            if atom in removed_atoms:
                continue
            for neighbor in atom.bonds:
                if neighbor not in removed_atoms:
                    continue
                target_atom = orig_to_new_map[atom]
                if head_match_atoms and neighbor in head_match_atoms:
                    if target_atom.label and target_atom.label != '*2':
                        raise ValueError(f"Label conflict: Atom already {target_atom.label}, wants *2")
                    target_atom.label = '*2'
                    _ensure_open_site(target_atom)
                elif tail_match_atoms and neighbor in tail_match_atoms:
                    if target_atom.label and target_atom.label != '*1':
                        raise ValueError(f"Label conflict: Atom already {target_atom.label}, wants *1")
                    target_atom.label = '*1'
                    _ensure_open_site(target_atom)
        new_mol.update_multiplicity()

    @staticmethod
    def _assert_end_group(mol: Molecule, want_label: str):
        """
        Validate that a scission fragment can serve as an end-group.
        Contract we enforce:
        - Must have exactly one labeled atom, and it must be `want_label` ('*1' or '*2').
        - Must be a mono-radical overall (get_radical_count() == 1).
        - That single radical must be located on the labeled atom (strongly recommended for stitching).
        - No stray labels ('1','2', etc.) allowed.

        Args:
            mol (Molecule): The molecule to validate.
            want_label (str): The expected label for the end-group ('*1' or '*2').
        """
        if want_label not in ('*1', '*2'):
            raise ValueError(f"want_label must be '*1' or '*2', got {want_label!r}")
        labels = _labels_present(mol)
        if labels.count(want_label) != 1:
            raise ValueError(f"End-group must contain exactly one {want_label} label. Got labels={labels}")
        allowed = {want_label}
        extras = [lab for lab in labels if lab not in allowed]
        if extras:
            raise ValueError(f"End-group has invalid extra labels {extras}; expected only {want_label}.")
        rad_count = mol.get_radical_count()
        if rad_count != 1:
            raise ValueError(f"End-group must be mono-radical (radical_count==1). Got {rad_count}.")
        labeled_atom = next(a for a in mol.atoms if a.label == want_label)
        if labeled_atom.radical_electrons < 1:
            raise ValueError(f"End-group labeled atom {want_label} must carry a radical electron for stitching.")

    @staticmethod
    def _assert_feature_unit(mol: Molecule):
        """
        Validate that an extracted 'feature monomer' is stitchable as a repeat unit.
        Contract we enforce:
        - Exactly one '*1' and exactly one '*2' label.
        - No other labels (including '1','2') survive.
        - Must have at least two radical electrons overall (typically exactly 2).
          We enforce exactly 2 by default because stitching consumes one at each end.
        - The labeled atoms must each have >=1 radical electron (open sites).

        Args:
            mol (Molecule): The molecule to validate.
        """
        labels = _labels_present(mol)
        n1 = _count_label(mol, '*1')
        n2 = _count_label(mol, '*2')
        if n1 != 1 or n2 != 1:
            raise ValueError(f"Feature unit must have exactly one '*1' and one '*2'. "
                             f"Got counts: '*1'={n1}, '*2'={n2}. Labels={labels}")
        extras = [lab for lab in labels if lab not in {'*1', '*2'}]
        if extras:
            raise ValueError(f"Feature unit has invalid extra labels {extras}; only '*1'/'*2' are allowed.")
        rad_count = mol.get_radical_count()
        if rad_count != 2:
            raise ValueError(f"Feature unit must have exactly 2 radical electrons total (one at each end). "
                             f"Got radical_count={rad_count}.")
        atom_1 = next(a for a in mol.atoms if a.label == '*1')
        atom_2 = next(a for a in mol.atoms if a.label == '*2')
        if atom_1.radical_electrons < 1 or atom_2.radical_electrons < 1:
            raise ValueError("Feature unit labeled atoms '*1' and '*2' must each carry a radical electron.")

    def get_thermo_data(self, solvent_name='', mode='auto'):
        """
        Returns the thermodynamic data of the polymer's proxy species.
        Delegates generation to the proxy if it doesn't exist.

        CRITICAL NOTE: This returns thermo on a "Per Reaction Site" basis (Proxy),
        not a "Per Chain" basis. This is required for RMG to calculate
        chemically meaningful reaction enthalpies (dH_rxn) and rates.

        For bulk heat capacity calculations, the solver should treat the concentration
        of this species as 'concentration of monomer units'.

        Args:
            solvent_name (str): Solvent for liquid phase corrections.
            mode (str): 'auto' (default), 'baseline', or 'feature'. Controls which proxy is used.
        """
        from rmgpy.thermo.thermoengine import submit
        proxy = self.get_proxy_species()
        if proxy is None:
            raise RuntimeError(f"Polymer '{self.label}': Could not determine proxy species for mode='{mode}'.")
        if not proxy.thermo:
            submit(proxy, solvent_name)
        if proxy.thermo is None:
            raise RuntimeError(f"Polymer '{self.label}': Thermo generation failed for proxy '{proxy.label}'.")
        if hasattr(proxy.thermo, 'result') and not isinstance(proxy.thermo, (NASA, Wilhoit, ThermoData)):
            proxy.thermo = proxy.thermo.result()
        self.thermo = proxy.thermo
        if self.thermo.comment and not self.thermo.comment.endswith(', Polymer'):
            self.thermo.comment += ', Polymer'
        elif not self.thermo.comment:
            self.thermo.comment = 'Polymer'
        return self.thermo

    def get_enthalpy(self, T):
        """Return enthalpy of the proxy (per-site basis) in J/mol at the specified temperature `T` in K."""
        return self.get_thermo_data().get_enthalpy(T)

    def get_entropy(self, T):
        """Return entropy of the proxy (per-site basis) in J/mol*K at the specified temperature `T` in K."""
        return self.get_thermo_data().get_entropy(T)

    def get_free_energy(self, T):
        """Return Gibbs free energy of the proxy (per-site basis) in J/mol at the specified temperature `T` in K."""
        return self.get_thermo_data().get_free_energy(T)

    def get_heat_capacity(self, T):
        """
        Return heat capacity of the proxy (per-site basis) in J/mol*K at the specified temperature `T` in K.

        To get the bulk heat capacity of a chain with degree of polymerization DP,
        you would calculate: DP * get_heat_capacity(T).
        """
        return self.get_thermo_data().get_heat_capacity(T)

    def get_bulk_heat_capacity(self, T, DP: float) -> float:
        """
        Helper to calculate the total heat capacity for a chain of length DP.
        Cp_bulk(T) ≈ DP * Cp_proxy(T)

        Args:
            T (float): Temperature in Kelvin.
            DP (float): Degree of Polymerization (number of monomer units).
        """
        return DP * self.get_heat_capacity(T)

    def calculate_cp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        Delegates Cp0 calculation to the proxy molecule.
        """
        proxy = self.get_proxy_species()
        if not proxy or not proxy.molecule:
            return 0.0
        return proxy.molecule[0].calculate_cp0()

    def calculate_cpinf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        Delegates CpInf calculation to the proxy molecule.
        """
        proxy = self.get_proxy_species()
        if not proxy or not proxy.molecule:
            return 0.0
        return proxy.molecule[0].calculate_cpinf()

    def generate_transport_data(self):
        """
        Generates transport data for the proxy species.
        (Future improvement: We can scale the diffusivity D based on polymer Mn).
        """
        proxy = self.get_proxy_species()
        if not proxy:
            raise RuntimeError(f"Polymer '{self.label}': No proxy available for transport generation.")
        if not proxy.transport_data:
            proxy.generate_transport_data()
        self.transport_data = proxy.transport_data
        return self.transport_data

    def generate_statmech(self):
        """
        Generates statistical mechanics data (frequencies, modes) for the proxy species.
        Used for Master Equation calculations.
        """
        proxy = self.get_proxy_species()
        if not proxy:
            raise RuntimeError(f"Polymer '{self.label}': No proxy available for statmech generation.")
        if not proxy.has_statmech():
            proxy.generate_statmech()
        self.conformer = proxy.conformer
        return self.conformer


class PolymerClass(str, Enum):
    """
    Classification of reaction products relative to the original polymer proxy structure.
    Notation: X = head, Y = tail, O = monomer, . = feature, [ ] = wing element, W - scission.
    """
    # --- Non-Polymer or Error States ---
    GAS = 'GAS'              # No intact wings found (Small molecule byproduct/gas)
    UNKNOWN = 'UNKNOWN'      # Classifier failed to confidently parse the topology

    # --- Intact Backbone States (2 Wings Found) ---
    BASELINE = 'BASELINE'    # Exactly matches the unreacted starting proxy [X-O]-O-[O-Y]
    FEATURE = 'FEATURE'      # Center monomer modified [X-O]-O.-[O-Y]
    END_MOD = 'END_MOD'      # Terminal head or tail modified [X.-O]-O-[O-Y] or [X-O]-O-[O-Y.]
    DISCARD = 'DISCARD'      # Buffer monomer modified (ignore to prevent double-counting) [X-O.]-O-[O-Y] or [X-O]-O-[O.-Y]

    # --- Chain Breaking/Linking States ---
    SCISSION = 'SCISSION'    # Only 1 wing found; chain broke (e.g., [X-O]-W)
    CROSSLINK = 'CROSSLINK'  # >2 wings found; bi-molecular polymer recombination
    CHIP = 'CHIP'            # Fold-back parent copy left by DISCRETE_CHIP product
                             # surgery (surge_chip_products, spec 2026-06-10 §4.2).
                             # Never produced by classify_structure itself.


def is_end_group_reaction(products) -> bool:
    """
    True iff this reaction's products mark it as an end-group (terminal)
    modification — i.e. some product :class:`Polymer` was classified ``END_MOD``
    by :meth:`Polymer.create_reacted_copy` (which stamps ``_reacted_class``).

    End-group reactions occur at chain ends, so the polymer hybrid solver scales
    them by chain-end density (mu0) rather than monomer-unit density (mu1).
    Non-Polymer products and polymers without a stamped verdict are ignored,
    leaving the default mu1 scaling.
    """
    return any(getattr(p, '_reacted_class', None) == PolymerClass.END_MOD
               for p in products if isinstance(p, Polymer))


class PolymerFluxArchetype(IntEnum):
    """
    Per-reaction pool moment-flux archetype, stamped at generation time on
    ``Reaction.polymer_flux_archetype`` and dispatched by the polymer hybrid
    solver's residual. The solver reads the STORED reaction flags; nothing
    downstream may recompute them from product stamps (chip surgery re-stamps
    products). See
    docs/superpowers/specs/2026-06-09-proxy-moment-flux-apportionment-design.md
    and docs/superpowers/specs/2026-06-10-discreteness-gate-discrete-chip-design.md
    (DISCRETE_CHIP).
    """
    NONE = 0               # no proxy involvement
    SAME_POOL = 1          # product folds back into the reactant pool (net-zero)
    MIGRATION = 2          # whole chain migrates to a different pool
    SCISSION_FRAGMENT = 3  # chain cut; fragment to daughter, complement stays
    UNRESOLVED = 4         # ambiguous/unstamped: solver applies legacy mu1 flux
    DISCRETE_CHIP = 5      # end-anchored cut ejects a stamped a-unit discrete
                           # chip; complement folds back to the SAME pool
                           # (spec 2026-06-10). Mirror: solver FLUX_DISCRETE_CHIP.


_flux_archetype_warned = set()


def _warn_unresolved_archetype(reason: str, detail: tuple) -> None:
    """Log each distinct UNRESOLVED-archetype cause once (not per reaction)."""
    key = (reason, detail)
    if key not in _flux_archetype_warned:
        _flux_archetype_warned.add(key)
        logging.warning(
            "Polymer flux archetype UNRESOLVED (%s): %s -- the solver will "
            "apply legacy mu1-only moment flux for this reaction shape.",
            reason, detail)


_chip_tripwire_warned = set()


def _warn_probable_end_cut(detail) -> None:
    """
    Diagnostics-only census (spec 2026-06-10 §4.4): warn once per distinct
    piece when a mu1-scaled scission's represented piece is end-confined.
    Never affects routing. The accumulated count on real decks measures how
    much chemistry waits on the end-anchor detector follow-up item.
    """
    if detail not in _chip_tripwire_warned:
        _chip_tripwire_warned.add(detail)
        logging.warning(
            "Polymer scission piece %s is end-confined (wing + <=1 repeat "
            "unit) but the reaction is mu1-scaled: probable mis-scaled "
            "end-anchored cut; routed SCISSION_FRAGMENT pending the "
            "end-anchor detector item.", detail)


def classify_reaction_flux_archetype(reactants, products) -> PolymerFluxArchetype:
    """
    Classify a reaction's pool moment-flux archetype from its (handshaked)
    reactant and product lists. Product Polymers carry ``_reacted_class``
    stamped by :meth:`Polymer.create_reacted_copy` (or re-stamped by chip
    product surgery, spec 2026-06-10 §4.2); pool identity is the Polymer
    label (the same key the solver's ``initialize_model`` uses to map
    species to pools).
    """
    reactant_pools = {r.label for r in reactants if isinstance(r, Polymer)}
    product_polymers = [p for p in products if isinstance(p, Polymer)]
    if not reactant_pools and not product_polymers:
        return PolymerFluxArchetype.NONE

    if any(getattr(p, '_reacted_class', None) == PolymerClass.CHIP
           for p in product_polymers):
        # Chip product surgery (surge_chip_products) already rewrote this
        # product list to [discrete chip, CHIP-stamped fold-back]. This check
        # MUST precede the SCISSION branch: after the (b)-surgery there is no
        # END_MOD member left, so is_end_group_reaction(products) would
        # recompute False and misroute. The solver reads the STORED Reaction
        # flag; nothing downstream may recompute it from product stamps.
        return PolymerFluxArchetype.DISCRETE_CHIP

    if any(getattr(p, '_reacted_class', None) == PolymerClass.SCISSION
           for p in product_polymers):
        if is_end_group_reaction(products):
            # Unsurged end-initiated scission: chip product surgery
            # (surge_chip_products, spec 2026-06-10 §4.2) was either not
            # attempted or infeasible for this shape. Surged shapes never
            # reach here (the CHIP branch above short-circuits). UNRESOLVED,
            # never SCISSION_FRAGMENT: uniform-cut statistics near a chain
            # end are wrong AND the chip mass would go unaccounted.
            _warn_unresolved_archetype(
                "unsurged end-initiated scission",
                tuple(sorted(p.label for p in product_polymers)))
            return PolymerFluxArchetype.UNRESOLVED
        # Tripwire diagnostic (spec 2026-06-10 §4.4): structure is used for
        # DIAGNOSTICS only, never routing. "wing + at most 1 repeat unit"
        # (not "wing only") so cap+1-unit pieces -- plausibly the most common
        # single-step end cuts -- are counted. Heavy-atom bound:
        # max(cap heavy) + 2 * monomer heavy (wing = cap + 1 unit, plus 1).
        parent = next((r for r in reactants if isinstance(r, Polymer)), None)
        piece = next(
            (p for p in product_polymers
             if getattr(p, '_reacted_class', None) == PolymerClass.SCISSION),
            None)
        src_mol = getattr(piece, '_source_molecule', None) if piece is not None else None
        if parent is not None and src_mol is not None:
            piece_heavy = sum(1 for a in src_mol.atoms if not a.is_hydrogen())
            mon_heavy = sum(1 for a in parent.monomer.atoms if not a.is_hydrogen())
            cap_heavy = max(
                sum(1 for a in eg.atoms if not a.is_hydrogen())
                for eg in parent.end_groups)
            if piece_heavy <= cap_heavy + 2 * mon_heavy:
                _warn_probable_end_cut(piece.label)
        return PolymerFluxArchetype.SCISSION_FRAGMENT

    if not product_polymers:
        # Polymer reactant but no polymer product (e.g. full conversion to
        # gas). No flux rule exists for this shape; the solver-level phase
        # check skips such core reactions anyway, so flag it loudly.
        _warn_unresolved_archetype(
            "polymer reactant with no polymer product",
            tuple(sorted(reactant_pools)))
        return PolymerFluxArchetype.UNRESOLVED

    cross_pool = [p for p in product_polymers if p.label not in reactant_pools]
    if not cross_pool:
        return PolymerFluxArchetype.SAME_POOL
    if (len(product_polymers) == 1 and len(cross_pool) == 1
            and len(reactant_pools) == 1):
        return PolymerFluxArchetype.MIGRATION
    _warn_unresolved_archetype(
        "ambiguous cross-pool shape",
        tuple(sorted(p.label for p in product_polymers)))
    return PolymerFluxArchetype.UNRESOLVED


def surge_chip_products(products, parent: 'Polymer') -> Optional[int]:
    """
    Chip product surgery (spec 2026-06-10 §4.2): rewrite a flag-true
    (end-group-scaled) scission-shaped product list IN PLACE into the
    canonical DISCRETE_CHIP end state [discrete chip Molecule, CHIP-stamped
    fold-back Polymer], and return the chip's repeat-unit count
    ``a = round(chip_MW / parent.monomer_mw_g_mol)``. ``a == 0`` is legal
    (bare cap ejection) and distinct from the infeasible return ``None``.

    Sub-shapes (chip identification = the smaller piece; the cut POSITION is
    already known from the reaction's stored is_end_group_reaction flag):

    (b) END_MOD fold-back present -- the only flag-true shape live today:
        the SCISSION-stamped Polymer is the chip. Demote it back to its
        sanitized source Molecule (undoing the handshake conversion) and
        re-stamp the existing END_MOD fold-back as CHIP. Applying (a) here
        would replace the chip with a second fold-back -- losing the chip
        and double-folding the parent.
    (a) no END_MOD product (dormant until the end-anchor detector item):
        the SCISSION-stamped Polymer is the MACRO daughter; the chip is the
        single non-Polymer co-product (already discrete -- the handshake left
        it a Molecule). Replace the daughter with ``parent.copy(deep=True)``
        stamped ``PolymerClass.CHIP``.

    Returns ``None`` WITHOUT modifying ``products`` when the shape is not a
    feasible chip shape (no/ambiguous scission piece, chip unrepresentable or
    not the smaller piece, missing source molecule). The caller
    (stamp_polymer_flux_archetype) then stamps UNRESOLVED -- never
    SCISSION_FRAGMENT.
    """
    product_polymers = [p for p in products if isinstance(p, Polymer)]
    scissions = [p for p in product_polymers
                 if getattr(p, '_reacted_class', None) == PolymerClass.SCISSION]
    end_mods = [p for p in product_polymers
                if getattr(p, '_reacted_class', None) == PolymerClass.END_MOD]
    if len(scissions) != 1:
        return None  # not a chip shape, or ambiguous piece identification
    daughter = scissions[0]
    proxy_mw = parent.baseline_proxy.molecule[0].get_molecular_weight()

    if end_mods:
        # --- sub-shape (b): SCISSION piece = chip, END_MOD = fold-back ---
        if len(end_mods) != 1:
            return None
        chip_src = getattr(daughter, '_source_molecule', None)
        if chip_src is None:
            return None  # cannot demote: handshake source unavailable
        if chip_src.get_molecular_weight() >= proxy_mw:
            return None  # "chip" is not the smaller piece
        chip_mol = chip_src.copy(deep=True)
        chip_mol.clear_labeled_atoms()
        for i, p in enumerate(products):
            if p is daughter:
                products[i] = chip_mol
                break
        end_mods[0]._reacted_class = PolymerClass.CHIP
        chip_mw_g = chip_mol.get_molecular_weight() * 1000.0
        return int(round(chip_mw_g / parent.monomer_mw_g_mol))

    # --- sub-shape (a): SCISSION piece = macro daughter, chip discrete ---
    chips = [p for p in products
             if not isinstance(p, Polymer) and isinstance(p, Molecule)]
    if len(chips) != 1:
        return None  # chip absent or ambiguous
    chip_mol = chips[0]
    daughter_src = getattr(daughter, '_source_molecule', None)
    ref_mw = (daughter_src.get_molecular_weight()
              if daughter_src is not None else proxy_mw)
    if chip_mol.get_molecular_weight() >= ref_mw:
        return None  # the discrete co-product is not the smaller piece
    fold = parent.copy(deep=True)
    fold._reacted_class = PolymerClass.CHIP
    for i, p in enumerate(products):
        if p is daughter:
            products[i] = fold
            break
    chip_mw_g = chip_mol.get_molecular_weight() * 1000.0
    return int(round(chip_mw_g / parent.monomer_mw_g_mol))


def stamp_polymer_flux_archetype(forward, reactants, polymer_reactants) -> None:
    """
    Stamp ``forward.polymer_flux_archetype`` (and ``polymer_chip_units``)
    AFTER the handshake and AFTER ``forward.is_end_group_reaction`` is stored.
    Called from both make_new_reaction stamping branches (rmgpy/rmg/model.py).

    Flag-true shapes run chip product surgery FIRST so the classifier sees
    the surged product list (its CHIP branch precedes the SCISSION
    recompute). An infeasible flag-true scission shape stamps UNRESOLVED +
    warn-once -- never SCISSION_FRAGMENT (spec 2026-06-10 §4.2).
    """
    chip_a = None
    if forward.is_end_group_reaction and len(polymer_reactants) == 1:
        chip_a = surge_chip_products(forward.products, polymer_reactants[0])
    if chip_a is not None:
        forward.polymer_chip_units = chip_a
    elif forward.is_end_group_reaction and any(
            getattr(p, '_reacted_class', None) == PolymerClass.SCISSION
            for p in forward.products if isinstance(p, Polymer)):
        _warn_unresolved_archetype(
            "infeasible chip surgery",
            tuple(sorted(getattr(p, 'label', '?') for p in forward.products
                         if isinstance(p, Polymer))))
        forward.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
        return
    forward.polymer_flux_archetype = int(
        classify_reaction_flux_archetype(reactants, forward.products))


MatchMapping = Mapping[Any, Any]


@dataclass(frozen=True)
class MatchSummary:
    raw: int
    disjoint: int
    best_matches: List[MatchMapping]


def stitch_molecules_by_labeled_atoms(mol_1: Optional[Molecule],
                                      mol_2: Optional[Molecule],
                                      left_labels: Optional[Tuple[str, ...]] = None,
                                      right_labels: Optional[Tuple[str, ...]] = None,
                                      ) -> Optional[Molecule]:
    """
    Stitches two molecules together at their labeled '*1', '*2' atoms.

    Args:
        mol_1 (Optional[Molecule]): The first molecule (with '*1' label).
        mol_2 (Optional[Molecule]): The second molecule (with '*2' label).
        left_labels (Optional[Tuple[str, ...]]): Labels to search for in the first molecule.
        right_labels (Optional[Tuple[str, ...]]): Labels to search for in the second molecule.

    Returns:
        Optional[Molecule]: The stitched molecule.
    """
    left_labels = left_labels or LABELS_1
    right_labels = right_labels or LABELS_2

    if not set(left_labels).isdisjoint(right_labels):
        raise ValueError("Stitch error: left_labels and right_labels overlap; ambiguous stitching.")

    if mol_1 is None or mol_2 is None:
        return None
    m1 = mol_1.copy(deep=True)
    m2 = mol_2.copy(deep=True)

    if sum(1 for a in m1.atoms if a.label in left_labels) != 1:
        raise ValueError("Stitch error: mol_1 must have exactly one left label.")
    if sum(1 for a in m2.atoms if a.label in right_labels) != 1:
        raise ValueError("Stitch error: mol_2 must have exactly one right label.")

    merged = m1.merge(m2)
    idx_1 = find_labeled_atom(merged, left_labels)
    idx_2 = find_labeled_atom(merged, right_labels)

    if idx_1 is None or idx_2 is None:
        raise ValueError("Stitch error: Could not locate labels after merge.")
    if idx_1 == idx_2:
        raise ValueError("Stitch error: Stitch sites resolved to the same atom.")

    atom_1 = merged.atoms[idx_1]
    atom_2 = merged.atoms[idx_2]

    if atom_1.radical_electrons < 1:
        raise ValueError(f"Stitch site 1 must have at least one radical electron (got {atom_1.radical_electrons}).")
    if atom_2.radical_electrons < 1:
        raise ValueError(f"Stitch site 2 must have at least one radical electron (got {atom_2.radical_electrons}).")

    bond = Bond(atom_1, atom_2, order=1)
    merged.add_bond(bond)
    atom_1.decrement_radical()
    atom_2.decrement_radical()
    atom_1.label = ''
    atom_2.label = ''
    merged.update_multiplicity()
    merged.update_atomtypes()
    return merged


def get_target_atoms(match: Mapping[Any, Any]) -> set:
    """
    Extracts ONLY the target molecule atoms from an RMG match mapping.
    Ensures we don't accidentally include GroupAtoms from the search pattern.
    """
    if not match:
        return set()
    atoms = set()
    for k, v in match.items():
        for item in (k, v):
            if hasattr(item, 'element') and not isinstance(item, type):
                atoms.add(item)
    return atoms


def find_labeled_atom(mol: Molecule, labels: Optional[tuple[str, ...]] = None) -> Optional[int]:
    """
    Finds the first atom in the molecule with any of the specified labels.

    Args:
        mol (Molecule): The molecule to search.
        labels (tuple[str, ...]): The labels to look for. If None, defaults to ('1', '*1', '2', '*2').

    Returns:
        Optional[int]: The index of the labeled atom, or None if not found.
    """
    labels = labels or (*LABELS_1, *LABELS_2)
    return next((i for i, a in enumerate(mol.atoms) if a.label in labels), None)


def _labels_present(mol: Molecule) -> list[str]:
    """Return all non-empty atom labels in the molecule (preserves duplicates)."""
    return [a.label for a in mol.atoms if a.label]


def _count_label(mol: Molecule, label: str) -> int:
    """Count the number of atoms with a specific label in the molecule."""
    return sum(1 for a in mol.atoms if a.label == label)


def _ensure_open_site(atom: 'Atom') -> None:
    """
    Ensure `atom` has at least one radical electron so it can serve as a stitch site.

    This is used after we identify a cut bond between the kept remainder and a removed wing:
    we label the kept-side atom (*1 or *2) and must ensure it is open-shell.

    Args:
        atom (Atom): Atom in the remainder molecule to make radical if needed.
    """
    if atom.radical_electrons >= 1:
        return
    atom.increment_radical()


def classify_structure(species: 'Species',
                       original_polymer,
                       *,
                       monomer_group: Optional['Group'] = None,
                       ) -> Tuple[PolymerClass, Dict[str, Any]]:
    """
    Classifies a reaction product structurally ("Topological Partitioning") by counting intact monomer subgraphs.
    """
    base_details = {"raw_matches": 0, "disjoint_matches": 0}

    if not species.molecule:
        return PolymerClass.UNKNOWN, {**base_details, "reason": "no_molecule"}
    if monomer_group is None:
        monomer_group = original_polymer.backbone_group
    if monomer_group is None:
        return PolymerClass.GAS, {**base_details, "reason": "no_backbone_group"}
    if len(species.molecule[0].atoms) < len(monomer_group.atoms):
        return PolymerClass.GAS, {**base_details, "reason": "too_few_atoms_for_monomer"}

    head_wing = original_polymer._wing_groups('head')
    tail_wing = original_polymer._wing_groups('tail')
    wing_count, wing_match_details = _analyze_wing_matches(species.molecule[0], head_wing, tail_wing, monomer_group)

    # BRANCH A: INTACT BACKBONE (2 or more wings) ---
    if wing_count >= 2:
        base_details.update(wing_match_details)

        # 1. Crosslinking Check (>2 wings means bi-molecular polymer combination)
        if wing_count > 2:
            return PolymerClass.CROSSLINK, {**base_details, "reason": "more_than_two_wings_found"}

        # 1.5 Discreteness gate: backbone impostor rejection (spec 2026-06-10
        # docs/superpowers/specs/2026-06-10-discreteness-gate-discrete-chip-design.md
        # §3.1). A genuine 2-wing candidate is an image of the baseline proxy
        # REPRESENTATION; a small molecule that merely contains both wing
        # subgraphs (e.g. bibenzyl against a PS proxy) sits far below it in
        # heavy atoms. One-sided on purpose: legitimate images can be larger
        # (FEATURE side groups) or modestly smaller (H loss, side-group
        # elimination). Heavy atoms, not MW: H-insensitive, and a lost side
        # group is a known heavy-atom delta.
        # No upper ceiling (spec §10-V3, verified 2026-06-10): polymer+polymer
        # coupling cannot deliver a >=2x-proxy candidate here -- a coupling
        # product carries >2 wings and returns CROSSLINK above, and
        # create_reacted_copy raises PolymerCrosslinkError for it upstream
        # (polymer.py create_reacted_copy crosslink guard) and make_new_reaction
        # discards the reaction, so coupled shapes never reach this gate.
        proxy_heavy = sum(
            1 for a in original_polymer.baseline_proxy.molecule[0].atoms
            if not a.is_hydrogen())
        cand_heavy = sum(
            1 for a in species.molecule[0].atoms if not a.is_hydrogen())
        if cand_heavy < proxy_heavy - round(0.35 * proxy_heavy):
            return PolymerClass.GAS, {**base_details, "reason": "backbone_impostor"}

        # 2. Baseline Check (Is it the exact unreacted proxy?)
        # Evaluates: [X-O]-O-[O-Y]
        if original_polymer.baseline_proxy.is_isomorphic(species):
            return PolymerClass.BASELINE, {**base_details, "reason": "unreacted_proxy"}

        # 3. End-Group Modification Check
        # Evaluates: [X.-O]-O-[O-Y] or [X-O]-O-[O-Y.]
        if _is_end_group_modified(wing_match_details, original_polymer):
            return PolymerClass.END_MOD, {**base_details, "reason": "terminal_end_modified"}

        # 4. Buffer Monomer Modification Check
        # Evaluates: [X-O.]-O-[O-Y] or [X-O]-O-[O.-Y]
        if _is_buffer_monomer_modified(wing_match_details, original_polymer):
            return PolymerClass.DISCARD, {**base_details, "reason": "buffer_monomer_modified"}

        # 5. Feature Check
        # Evaluates: [X-O]-O.-[O-Y]
        # If it has 2 wings, is not the baseline, not an end mod, and not a buffer mod,
        # it MUST be a valid central feature modification. We verify just to be safe.
        if _is_center_feature_modified(species.molecule[0], wing_match_details):
            return PolymerClass.FEATURE, {**base_details, "reason": "center_monomer_modified"}

        # 6. Fallback for anomalous intact backbones (e.g. graph matching weirdness)
        return PolymerClass.UNKNOWN, {**base_details, "reason": "unclassified_intact_backbone"}

    # BRANCH B: BROKEN BACKBONE OR GAS (0 or 1 wing) ---
    elif wing_count == 1:
        base_details.update(wing_match_details)
        # Evaluates: [X-O]-W or W-[O-Y]
        return PolymerClass.SCISSION, {**base_details, "reason": "single_terminal_wing"}

    # BRANCH C: ANOMALOUS CASE ---
    else:
        return PolymerClass.GAS, {**base_details, "reason": "no_intact_wings"}


def _analyze_wing_matches(product_mol: Molecule,
                          head_wings: List['Group'],
                          tail_wings: List['Group'],
                          monomer_group: 'Group',
                          ) -> Tuple[int, Dict[str, Any]]:
    """
    Performs a single, pooled subgraph search for all wing patterns.
    Uses Maximum Set Packing with a 25% monomer overlap threshold to find the optimal wings.

    Args:
        product_mol (Molecule): The molecule to analyze.
        head_wings (List[Group]): List of head wing patterns to search for.
        tail_wings (List[Group]): List of tail wing patterns to search for.
        monomer_group (Group): The monomer group used to set the overlap threshold.

    Returns:
        Tuple[int, Dict[str, Any]]: The number of valid wings found and detailed match
    """
    mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())
    overlap_threshold = int(0.25 * mon_heavy_count)

    # Normalize product molecule to Kekulé (localized) Heavy View
    full_mol = product_mol.copy(deep=True)
    if full_mol.is_aromatic():
        full_mol.kekulize()
    full_mol.clear_labeled_atoms()
    full_mol.update_multiplicity()

    heavy_mol, copied_heavy_to_full = get_heavy_view_with_maps(full_mol)
    heavy_to_full = {}
    for heavy_atom, copied_full_atom in copied_heavy_to_full.items():
        idx = full_mol.atoms.index(copied_full_atom)
        heavy_to_full[heavy_atom] = product_mol.atoms[idx]

    raw_matches = []
    for side, wings in [('head', head_wings), ('tail', tail_wings)]:
        for g in wings:
            g_heavy_group = Group()
            g_mapping = {}

            for ga in g.atoms:
                symbol = get_element_symbol(ga)
                if symbol == 'H' or getattr(ga, 'is_hydrogen', lambda: False)():
                    continue

                new_ga = GroupAtom(atomtype=[ATOMTYPES[symbol]],
                                   radical_electrons=[], charge=[], lone_pairs=[])
                g_heavy_group.add_atom(new_ga)
                g_mapping[ga] = new_ga

            for ga1 in g.atoms:
                if ga1 not in g_mapping: continue
                for ga2, g_bond in ga1.edges.items():
                    if ga2 not in g_mapping: continue
                    if id(ga1) < id(ga2):
                        new_bond = GroupBond(g_mapping[ga1], g_mapping[ga2], order=[1, 1.5, 2, 3, 4])
                        g_heavy_group.add_bond(new_bond)

            matches = heavy_mol.find_subgraph_isomorphisms(g_heavy_group, save_order=True)
            for m in matches:
                atom_set = set(m.keys()) if isinstance(m, dict) else set(m)
                cut_edges = len(get_heavy_cut_edges(atom_set))

                if cut_edges > 1:
                    continue

                raw_matches.append({
                    'side': side,
                    'atoms': atom_set,
                    'cut_edges': cut_edges,
                    'size': len(atom_set)
                })

    raw_matches = sorted(raw_matches, key=lambda x: (x['cut_edges'], -x['size']))

    def find_max_disjoint_set(candidates: List[Dict]) -> List[Dict]:
        if not candidates:
            return []

        first = candidates[0]
        compatible_with_first = [
            c for c in candidates[1:]
            if len(first['atoms'].intersection(c['atoms'])) <= overlap_threshold
        ]
        universe_take = [first] + find_max_disjoint_set(compatible_with_first)
        universe_leave = find_max_disjoint_set(candidates[1:])

        return universe_take if len(universe_take) >= len(universe_leave) else universe_leave

    best_wings = find_max_disjoint_set(raw_matches)
    wing_count = len(best_wings)

    head_matches = [w for w in best_wings if w['side'] == 'head']
    tail_matches = [w for w in best_wings if w['side'] == 'tail']

    if len(head_matches) >= 2 and not tail_matches:
        tail_matches = [head_matches.pop()]
    elif len(tail_matches) >= 2 and not head_matches:
        head_matches = [tail_matches.pop()]

    match_details = {
        "num_disjoint_wings": wing_count,
        "head_match": head_matches[0] if head_matches else None,
        "tail_match": tail_matches[0] if tail_matches else None,
        "all_optimal_wings": best_wings,
        "raw_search_hits": len(raw_matches),
        "heavy_to_full_map": heavy_to_full
    }

    return wing_count, match_details


def get_heavy_view_with_maps(full_mol: Molecule) -> Tuple[Molecule, Dict[Atom, Atom]]:
    """
    Creates a heavy-atom-only view of a molecule while preserving a map back to the original.

    Args:
        full_mol (Molecule): The original, full Molecule (with hydrogens).

    Returns:
        Tuple[Molecule, Dict[Atom, Atom]]:
            - The hydrogen-stripped Molecule.
            - A dictionary mapping the new heavy atoms back to the exact Atom objects in full_mol.
    """
    heavy_mol = full_mol.copy(deep=True)
    heavy_to_full = {}
    for heavy_atom, original_atom in zip(heavy_mol.atoms, full_mol.atoms):
        if not heavy_atom.is_hydrogen():
            heavy_to_full[heavy_atom] = original_atom
    # Iterate over a slice [:] so we don't mutate the list while iterating over it
    for atom in heavy_mol.atoms[:]:
        if atom.is_hydrogen():
            heavy_mol.remove_atom(atom)
    return heavy_mol, heavy_to_full


def get_heavy_cut_edges(atom_set: Set[Atom]) -> List[Tuple[Atom, Atom]]:
    """
    Identifies the bonds that connect a subset of atoms (the wing) to the rest of the molecule (the Cut Set).

    Args:
        atom_set (Set[Atom]): A set of Atom objects representing the matched subgraph.

    Returns:
        List[Tuple[Atom, Atom]]: A list of tuples, where each tuple is (atom_inside_wing, atom_outside_wing).
    """
    cut_edges = []
    for atom in atom_set:
        for neighbor in atom.bonds.keys():
            if neighbor not in atom_set:
                cut_edges.append((atom, neighbor))
    return cut_edges


def _is_end_group_modified(wing_match_details: Dict[str, Any],
                           original_polymer) -> bool:
    """
    Determines if the structural modification in the product is located strictly
    within the terminal end-caps (Head or Tail).
    """
    heavy_to_full = wing_match_details['heavy_to_full_map']
    monomer_group = original_polymer.backbone_group
    mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())
    for side in ['head_match', 'tail_match']:
        match = wing_match_details.get(side)
        if not match:
            continue
        heavy_wing_atoms = match['atoms']
        end_group_heavy, buffer_heavy = _slice_wing(heavy_wing_atoms, mon_heavy_count)
        end_group_full = {heavy_to_full[ha] for ha in end_group_heavy}
        if any(atom.radical_electrons > 0 for atom in end_group_full):
            return True
        if any(atom.label for atom in end_group_full if atom.label.startswith('*')):
            return True
    return False


def _slice_wing(heavy_wing_atoms: Set['Atom'], mon_heavy_count: int) -> Tuple[Set['Atom'], Set['Atom']]:
    """
    Uses Breadth-First Search (BFS) starting from the cut edge to slice a wing
    into its Buffer Monomer and true End-Group.

    Returns:
        Tuple[Set[Atom], Set[Atom]]: (end_group_atoms, buffer_atoms)
    """
    cut_edges = get_heavy_cut_edges(heavy_wing_atoms)
    if not cut_edges:
        return heavy_wing_atoms, set()
    wing_root_atom = cut_edges[0][0]
    buffer_atoms = set()
    queue = [wing_root_atom]
    visited = {wing_root_atom}
    while queue and len(buffer_atoms) < mon_heavy_count:
        curr = queue.pop(0)
        buffer_atoms.add(curr)
        for neighbor in curr.bonds.keys():
            if neighbor in heavy_wing_atoms and neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    end_group_atoms = heavy_wing_atoms - buffer_atoms
    return end_group_atoms, buffer_atoms


def _is_buffer_monomer_modified(wing_match_details: Dict[str, Any],
                                original_polymer) -> bool:
    """
    Determines if the structural modification in the product is located strictly
    within the buffer monomer section of the wings.

    Args:
        wing_match_details (Dict): The dictionary returned by _analyze_wing_matches.
        original_polymer (Polymer): The original polymer object.

    Returns:
        bool: True if a modification is found in the buffer zone, False otherwise.
    """
    heavy_to_full = wing_match_details.get('heavy_to_full_map', {})
    if not heavy_to_full:
        return False
    monomer_group = original_polymer.backbone_group
    if not monomer_group:
        return False
    mon_heavy_count = sum(1 for ga in monomer_group.atoms if not ga.is_hydrogen())
    for side in ['head_match', 'tail_match']:
        match = wing_match_details.get(side)
        if not match:
            continue
        heavy_wing_atoms = match['atoms']
        _, buffer_heavy = _slice_wing(heavy_wing_atoms, mon_heavy_count)
        buffer_full = {heavy_to_full[ha] for ha in buffer_heavy}

        for atom in buffer_full:
            if atom.radical_electrons > 0:
                return True
            if atom.label and str(atom.label).startswith('*'):
                return True
    return False


def _is_center_feature_modified(product_mol: 'Molecule',
                                wing_match_details: Dict[str, Any],
                                ) -> bool:
    """
    Determines if the structural modification in the product is located strictly
    within the central repeating units (the feature) of the polymer.

    Args:
        product_mol (Molecule): The generated reaction product molecule.
        wing_match_details (Dict): The dictionary returned by _analyze_wing_matches.

    Returns:
        bool: True if a modification is found in the central feature, False otherwise.
    """
    heavy_to_full = wing_match_details.get('heavy_to_full_map', {})
    if not heavy_to_full:
        return False
    wing_full_atoms = set()
    for side in ['head_match', 'tail_match']:
        match = wing_match_details.get(side)
        if match and 'atoms' in match:
            for heavy_atom in match['atoms']:
                if heavy_atom in heavy_to_full:
                    full_atom = heavy_to_full[heavy_atom]
                    wing_full_atoms.add(full_atom)
                    for neighbor in full_atom.bonds.keys():
                        if neighbor.is_hydrogen():
                            wing_full_atoms.add(neighbor)
    center_full_atoms = set(product_mol.atoms) - wing_full_atoms
    for atom in center_full_atoms:
        if atom.radical_electrons > 0:
            return True
        if atom.label and str(atom.label).startswith('*'):
            return True
    return False


def process_polymer_candidates(candidates: List[Species],
                               _reaction_model,
                               original_polymer,
                               ) -> List[Species]:
    """
    Handshake function to convert generic Species into Polymer objects.
    """
    processed_list: List['Species'] = []
    stats = {k: 0 for k in PolymerClass}
    monomer_group = original_polymer.backbone_group
    for cand in candidates:
        classification, details = classify_structure(cand, original_polymer, monomer_group=monomer_group)

        if not isinstance(classification, PolymerClass):
            raise TypeError(f"Expected PolymerClass enum, got {type(classification)}")

        stats[classification] += 1
        is_proxy = bool(classification != PolymerClass.GAS)

        if not hasattr(cand, "props"):
            cand.props = {}
        cand.props["is_polymer_proxy"] = is_proxy
        cand.is_polymer_proxy = is_proxy

        if getattr(cand, "molecule", None):
            for m in cand.molecule:
                if not hasattr(m, "props"):
                    m.props = {}
                m.props["is_polymer_proxy"] = is_proxy
                m.is_polymer_proxy = is_proxy

        if classification == PolymerClass.DISCARD:
            continue
        processed_list.append(cand)
    return processed_list


def get_element_symbol(atom: Union[Atom, GroupAtom]) -> str:
    """
    Extracts the element symbol from an Atom or GroupAtom.
    For GroupAtoms, searches the atomtype.generic list for the
    shortest string that doesn't contain 'R' or 'Val'.
    """
    if hasattr(atom, 'element') and atom.element is not None:
        return atom.element.symbol
    if hasattr(atom, 'atomtype') and atom.atomtype:
        at = atom.atomtype[0]
        candidates = [g.label for g in at.generic if 'R' not in g.label and 'Val' not in g.label]
        if candidates:
            return min(candidates, key=len)
        return at.label
    raise ValueError(f"Could not extract element from {type(atom)}: {atom}")


# ---------------------------------------------------------------------------
# Dynamic multi-pool spawning — see docs/multi_pool_design.md
# ---------------------------------------------------------------------------

@dataclass
class SpawnIntent:
    """A queued request to spawn a new polymer pool.

    Mirrors the ``polymer_pools.json`` sidecar schema entries (design doc §6).
    Created during product classification when a structurally novel chain
    population is detected; drained between RMG iterations to grow the pool
    registry (the solver is rebuilt on the next iteration — no in-place resize).
    """

    parent_pool: 'Polymer'
    monomer: Union[Group, Molecule]
    end_groups: List[str]
    triggering_product: Optional['Species'] = None
    triggering_dp: int = 0
    triggering_moles: float = 0.0
    triggering_reaction_index: Optional[int] = None
    mass_flux_at_spawn: float = 0.0


class MassFluxAccumulator:
    """Trailing-window accumulator for mass produced into a candidate motif.

    Used by the spawn gate (design doc §4.4): a motif must accumulate at
    least ``threshold`` fraction of total polymer-derived mass over the
    last ``window`` RMG iterations before it is allowed to spawn its own
    pool. Single transient peaks therefore cannot trigger spawning.
    """

    def __init__(self, window: int = 3):
        if window < 1:
            raise ValueError(f"window must be >= 1; got {window}")
        self.window = window
        # {motif_key: list of (iteration, mass)}
        self._records: Dict[str, List[Tuple[int, float]]] = defaultdict(list)

    def record(self, motif_key: str, mass: float, iteration: int) -> None:
        """Record ``mass`` produced into ``motif_key`` at iteration ``iteration``.

        Entries older than ``window`` iterations relative to the recording
        iteration are evicted on each call.
        """
        cutoff = iteration - self.window + 1
        self._records[motif_key] = [
            (i, m) for (i, m) in self._records[motif_key] if i >= cutoff
        ]
        self._records[motif_key].append((iteration, float(mass)))

    def flux(self, motif_key: str) -> float:
        """Sum of masses currently in the rolling window for ``motif_key``."""
        if motif_key not in self._records:
            return 0.0
        return sum(m for (_, m) in self._records[motif_key])


def _bfs_grow_heavy_subset(
    start: 'Atom',
    size: int,
) -> Optional[Set['Atom']]:
    """BFS-grow a connected subset of ``size`` heavy atoms starting at ``start``.

    Returns the set, or ``None`` if a connected subset of that size is not
    reachable from ``start`` via heavy-atom-only edges.
    """
    if start.element.symbol == 'H':
        return None
    visited: Set[Atom] = {start}
    if size == 1:
        return visited
    queue: List[Atom] = [start]
    while queue and len(visited) < size:
        a = queue.pop(0)
        for nbr in a.edges.keys():
            if nbr.element.symbol != 'H' and nbr not in visited:
                visited.add(nbr)
                if len(visited) == size:
                    break
                queue.append(nbr)
    return visited if len(visited) == size else None


def _atom_subset_to_group(sub_atoms: Set['Atom']) -> Optional[Group]:
    """Build a :class:`Group` pattern from an arbitrary subset of Atoms.

    Mirrors :meth:`Molecule.to_group` but iterates a subset; bonds are added
    only for edges whose both endpoints are in the subset.
    """
    if not sub_atoms:
        return None
    atom_to_group: Dict['Atom', GroupAtom] = {}
    for atom in sub_atoms:
        ga = GroupAtom(
            atomtype=[atom.atomtype],
            radical_electrons=[atom.radical_electrons],
            charge=[atom.charge],
            lone_pairs=[atom.lone_pairs],
            label=getattr(atom, 'label', '') or '',
        )
        atom_to_group[atom] = ga
    group = Group(atoms=list(atom_to_group.values()))
    seen: Set[Tuple[int, int]] = set()
    for atom in sub_atoms:
        for bonded, bond in atom.edges.items():
            if bonded in atom_to_group:
                key = tuple(sorted((id(atom), id(bonded))))
                if key in seen:
                    continue
                seen.add(key)
                group.add_bond(
                    GroupBond(atom_to_group[atom], atom_to_group[bonded],
                              order=[bond.order])
                )
    group.update()
    return group


def count_disjoint_subgraph_isomorphisms(
    mol: Molecule,
    group: Group,
) -> int:
    """Count how many mutually disjoint occurrences of ``group`` appear in ``mol``.

    Greedy Maximum Set Packing over :meth:`Molecule.find_subgraph_isomorphisms`.
    Mirrors the pattern used inside :func:`_analyze_wing_matches`.
    """
    if group is None:
        return 0
    try:
        mappings = mol.find_subgraph_isomorphisms(group, save_order=True)
    except (NotImplementedError, AttributeError, ValueError):
        return 0
    if not mappings:
        return 0
    occupied: Set[int] = set()
    count = 0
    for mapping in mappings:
        atoms_used = {id(a) for a in mapping.keys()}
        if not atoms_used & occupied:
            count += 1
            occupied |= atoms_used
    return count


def discover_repeat_motif(
    mol: Molecule,
    *,
    min_motif_size: int = 2,
) -> Optional[Group]:
    """Auto-detect a repeat motif within ``mol`` (design doc §4.2).

    Returns a :class:`Group` pattern that occurs at least twice as a
    disjoint subgraph in ``mol``, or ``None`` if no such motif exists.

    Algorithm: enumerate connected heavy-atom subsets of varying sizes,
    test each for ≥2 disjoint isomorphisms, prefer the motif that
    maximises the disjoint occurrence count (tie-break: smaller motif —
    more "basic" repeat unit). Phase-1 implementation; pathologically
    large products may be slow. See design doc §10.
    """
    heavy = [a for a in mol.atoms if a.element.symbol != 'H']
    if len(heavy) < 2 * min_motif_size:
        return None
    max_size = len(heavy) // 2

    seen_signatures: Set[frozenset] = set()
    best: Optional[Tuple[int, int, Group]] = None  # (n_occ, size, group)

    for start in heavy:
        for size in range(min_motif_size, max_size + 1):
            sub_atoms = _bfs_grow_heavy_subset(start, size)
            if sub_atoms is None or len(sub_atoms) != size:
                continue
            sig = frozenset(id(a) for a in sub_atoms)
            if sig in seen_signatures:
                continue
            seen_signatures.add(sig)

            group = _atom_subset_to_group(sub_atoms)
            if group is None:
                continue
            n_occ = count_disjoint_subgraph_isomorphisms(mol, group)
            if n_occ < 2:
                continue
            # Selection: max n_occ, tie-break smaller size (more "basic" motif).
            score = (n_occ, -size)
            if best is None or score > (best[0], -best[1]):
                best = (n_occ, size, group)
    return best[2] if best else None


def _estimate_relative_flux(
    candidate: 'Species',
    pool_registry: List['Polymer'],
    reaction_model: Any,
) -> float:
    """Estimate the fraction of polymer-derived mass flowing into ``candidate``.

    Phase-1 simplification: when ``reaction_model`` is ``None`` (unit-test
    path) return 0.5, so the spawn gate is exercised by the threshold knob
    without needing a full reaction-rate integrator. The real implementation,
    used during an RMG run, will integrate reaction rates and species
    molecular weights over the trailing window via :class:`MassFluxAccumulator`
    (design doc §4.4).
    """
    if reaction_model is None:
        return 0.5
    # TODO(multi-pool §4.4): real flux calculation against reaction_model
    return 0.5


def process_polymer_candidates_multipool(
    candidates: List['Species'],
    reaction_model: Any,
    pool_registry: List['Polymer'],
    *,
    max_pools: int = 5,
    mass_flux_threshold: float = 0.01,
    iteration: int = 0,
    flux_accumulator: Optional[MassFluxAccumulator] = None,
) -> Tuple[List['Species'], List[SpawnIntent]]:
    """Multi-pool aware product classification + spawn-intent generation.

    Extends :func:`process_polymer_candidates` (single-pool) per design doc §4.1
    by:

    * Classifying each candidate against EVERY pool in ``pool_registry``.
    * Running :func:`discover_repeat_motif` when no existing pool classifies.
    * Similarity-merging the discovered motif against existing pool patterns.
    * Gating spawns on a mass-flux threshold and a ``max_pools`` cap.

    Returns
    -------
    processed : list of Species
        Candidates that survived classification (i.e. were not dropped as
        :attr:`PolymerClass.DISCARD`). All survivors are tagged with
        ``is_polymer_proxy``.
    spawn_intents : list of :class:`SpawnIntent`
        Queued spawn requests to drain between RMG iterations (the daughter
        pools register `_muN` dummy species and the solver is rebuilt on the
        next iteration — no in-place resize; design doc §4.5/§7).
    """
    processed: List['Species'] = []
    spawn_intents: List[SpawnIntent] = []

    for cand in candidates:
        # Phase A: classify against every existing pool, take the first non-trivial hit.
        matched_pool: Optional['Polymer'] = None
        matched_class: Optional[PolymerClass] = None
        saw_unknown = False  # intact backbone (>=2 wings) but no clean pool match
        for pool in pool_registry:
            try:
                klass, _ = classify_structure(cand, pool)
            except Exception:
                continue
            if klass not in (PolymerClass.GAS, PolymerClass.UNKNOWN):
                matched_pool = pool
                matched_class = klass
                break
            if klass == PolymerClass.UNKNOWN:
                saw_unknown = True

        if matched_pool is not None:
            _tag_polymer_proxy(cand, is_proxy=(matched_class != PolymerClass.GAS))
            if matched_class != PolymerClass.DISCARD:
                processed.append(cand)
            continue

        # No clean pool match. An UNKNOWN classification means an intact backbone
        # (>=2 wings) that simply didn't match a pool's exact structure — it is
        # still structurally a polymer, so keep it in the polymer phase instead of
        # risking a leak to gas via Phase B. This unifies the semantics with the
        # single-pool process_polymer_candidates (classification != GAS -> proxy).
        if saw_unknown:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue

        # Phase B: novel-monomer discovery (only for candidates that were GAS
        # against every pool, i.e. not chain-like w.r.t. any existing pool).
        mol = cand.molecule[0] if getattr(cand, "molecule", None) else None
        if mol is None:
            continue
        motif = discover_repeat_motif(mol)
        if motif is None:
            _tag_polymer_proxy(cand, is_proxy=False)
            continue

        # Phase C: similarity-merge against existing pools.
        merged_pool = similarity_merge(motif, pool_registry)
        if merged_pool is not None:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue

        # Phase D: gates — relative flux and max_pools cap.
        relative_flux = _estimate_relative_flux(cand, pool_registry, reaction_model)
        if relative_flux < mass_flux_threshold:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue
        if len(pool_registry) >= max_pools:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue

        # Phase E: queue the spawn intent.
        triggering_dp = count_disjoint_subgraph_isomorphisms(mol, motif)
        parent_for_intent = pool_registry[0] if pool_registry else None
        if parent_for_intent is None:
            continue
        spawn_intents.append(
            SpawnIntent(
                parent_pool=parent_for_intent,
                monomer=motif,
                end_groups=list(parent_for_intent.end_groups),
                triggering_product=cand,
                triggering_dp=triggering_dp,
                triggering_moles=float(getattr(cand, "amount", 1.0)),
                mass_flux_at_spawn=relative_flux,
            )
        )
        _tag_polymer_proxy(cand, is_proxy=True)
        processed.append(cand)

    return processed, spawn_intents


POLYMER_POOLS_SIDECAR_SCHEMA_VERSION = "2.0"
POLYMER_POOLS_SIDECAR_FILENAME = "polymer_pools.json"


def _artifact_species_label(spc) -> str:
    """chem.yaml species name for ``spc`` — must match rmgpy.cantera.get_label:
    ``label(index)`` when index > 0, bare label otherwise (µ-dummies have
    index = -1 and appear bare in chem.yaml)."""
    index = getattr(spc, "index", -1)
    label = getattr(spc, "label", "")
    return f"{label}({index})" if index > 0 else label


def _species_base_label(spc) -> str:
    """Strip the RMG ``(N)`` index suffix — the solver's pool-membership rule
    (rmgpy/solver/polymer.pyx:478-485 uses label.partition('(')[0])."""
    return getattr(spc, "label", "").partition('(')[0]


def _get_rmg_commit():
    """Best-effort git SHA of the emitting RMG-Py checkout (envelope field)."""
    try:
        import subprocess
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None


def _serialize_pool_for_sidecar(pool: 'Polymer',
                                core_species: Optional[List['Species']] = None,
                                monomer_routing: Optional[str] = None) -> Dict[str, Any]:
    """Convert a :class:`Polymer` instance to a JSON-serialisable dict.

    Schema 2.0: 1.0 fields (docs/multi_pool_design.md §6) preserved verbatim;
    additions per docs/polymer_moments_format.md §2.
    """
    monomer_smiles = ""
    monomer_adj_list = ""
    try:
        if getattr(pool, "monomer", None) is not None:
            monomer_smiles = pool.monomer.to_smiles() if hasattr(pool.monomer, "to_smiles") else ""
            monomer_adj_list = (
                pool.monomer.to_adjacency_list() if hasattr(pool.monomer, "to_adjacency_list") else ""
            )
    except Exception:
        pass

    feature_smiles: List[str] = []
    feature_attr = getattr(pool, "feature_monomers", None) or (
        [pool.feature_monomer] if getattr(pool, "feature_monomer", None) else []
    )
    for fm in feature_attr:
        try:
            if hasattr(fm, "to_smiles"):
                feature_smiles.append(fm.to_smiles())
        except Exception:
            continue

    spawn_metadata = getattr(pool, "spawn_metadata", None) or {"source": "input"}
    mu_indices = getattr(pool, "mu_indices", None)
    if mu_indices is not None and not isinstance(mu_indices, dict):
        try:
            mu0_idx, mu1_idx, mu2_idx = mu_indices
            mu_indices = {"mu0_idx": mu0_idx, "mu1_idx": mu1_idx, "mu2_idx": mu2_idx}
        except Exception:
            mu_indices = None

    d = {
        "label": getattr(pool, "label", ""),
        "monomer_smiles": monomer_smiles,
        "monomer_adj_list": monomer_adj_list,
        "feature_monomers_smiles": feature_smiles,
        "end_groups": [
            eg.to_smiles() if hasattr(eg, "to_smiles") else str(eg)
            for eg in (getattr(pool, "end_groups", []) or [])
        ],
        "cutoff": getattr(pool, "cutoff", None),
        "parent_pool": getattr(pool, "parent_pool_label", None),
        "spawn_iteration": getattr(pool, "spawn_iteration", 0),
        "spawn_event_metadata": spawn_metadata,
        "mu_indices": mu_indices,
    }

    # --- schema 2.0 additions (field names pinned by TA's loader,
    #     ~/Code/TA/ta/mechanism.py) ---
    moments = getattr(pool, "moments", None)
    d["moments"] = [float(m) for m in moments] if moments is not None else None
    d["monomer_mw_g_mol"] = (float(pool.monomer_mw_g_mol)
                             if getattr(pool, "monomer_mw_g_mol", None) is not None else None)
    d["mn_g_mol"] = float(pool.Mn) if getattr(pool, "Mn", None) is not None else None
    d["mw_g_mol"] = float(pool.Mw) if getattr(pool, "Mw", None) is not None else None
    d["initial_mass_g"] = (float(pool.initial_mass_g)
                           if getattr(pool, "initial_mass_g", None) is not None else None)
    d["channels"] = {
        "scission": {"A": float(getattr(pool, "k_scission", 0.0)), "n": 0.0, "Ea": 0.0,
                     "units": {"A": "s^-1", "Ea": "J/mol"}},
        "unzip": {"A": float(getattr(pool, "k_unzip", 0.0)), "n": 0.0, "Ea": 0.0,
                  "units": {"A": "s^-1", "Ea": "J/mol"}},
    }
    phase_species: List[str] = []
    if core_species:
        member_bases = {pool.label, f"{pool.label}_mu0", f"{pool.label}_mu1",
                        f"{pool.label}_mu2"}
        for spc in core_species:
            if _species_base_label(spc) in member_bases:
                phase_species.append(_artifact_species_label(spc))
    if monomer_routing and monomer_routing not in phase_species:
        phase_species.append(monomer_routing)
    d["phase_species"] = phase_species
    d["monomer_routing"] = monomer_routing
    d["mu3_closure"] = "log_lagrange/1"
    return d


def write_polymer_pools_sidecar(
    pool_registry: List['Polymer'],
    output_dir: str,
    iteration: int = 0,
    filename: str = POLYMER_POOLS_SIDECAR_FILENAME,
    core_species=None,
    core_reactions=None,
    configured_pool_labels=None,
    condensed_species=None,
    monomer_routing_by_pool=None,
    cantera_index_map=None,
    rmg_commit=None,
) -> str:
    """Emit ``polymer_pools.json`` alongside ``chem.yaml`` (design doc §6).

    The TA-side mechanism loader (``~/Code/TA``) consumes this file in
    lock-step with the cantera YAML to recover pool semantics that cannot
    be reverse-inferred from ``<label>_muN`` pseudo-species names.

    Parameters
    ----------
    pool_registry : list of Polymer
        All live pools at the time of writing.
    output_dir : str
        Directory where ``chem.yaml`` is being written.
    iteration : int
        RMG iteration number. Recorded in the sidecar header.
    filename : str
        Override the default basename. Defaults to ``polymer_pools.json``.
    core_species : list of Species, optional
        Live core species list — used to populate ``phase_species`` on each
        pool block and ``condensed_species`` in the conventions envelope.
        ``core_species``/``monomer_routing_by_pool`` are populated by the
        live-run hook in ``save_everything``; legacy callers omit them, which
        yields ``phase_species: []`` and ``monomer_routing: null`` per pool.
    core_reactions : list of Reaction, optional
        Live core reactions list — compiled into the ``reactions[]`` block.
    configured_pool_labels : list of str, optional
        Solver-configured pool labels (may be a subset of pool_registry).
        Defaults to the label of every pool in pool_registry.
    condensed_species : list of Species, optional
        Condensed-phase core species (gas_species_mask == False). Used to
        populate ``conventions.condensed_species``.
    monomer_routing_by_pool : dict, optional
        ``{pool_label: routing_label_string}`` for monomer-in-poly routing.
        The caller must pass CONDENSED-phase labels (appended unchecked).
    cantera_index_map : dict, optional
        ``{id(rxn): [cantera indices]}`` from
        ``generate_cantera_data(..., return_reaction_index_map=True)``.
    rmg_commit : str, optional
        Override the auto-detected git SHA.

    Returns
    -------
    str
        Absolute path of the written file.
    """
    payload = build_polymer_moments_artifact(
        pool_registry,
        core_species=core_species,
        core_reactions=core_reactions,
        configured_pool_labels=configured_pool_labels,
        condensed_species=condensed_species,
        monomer_routing_by_pool=monomer_routing_by_pool,
        cantera_index_map=cantera_index_map,
        iteration=iteration,
        rmg_commit=rmg_commit,
    )
    path = os.path.join(output_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return path


# Archetype int -> versioned term-type name (docs/polymer_moments_format.md §3).
ARCHETYPE_TERM_NAMES = {
    int(PolymerFluxArchetype.SAME_POOL): "same_pool/1",
    int(PolymerFluxArchetype.MIGRATION): "migration/1",
    int(PolymerFluxArchetype.SCISSION_FRAGMENT): "scission_fragment/1",
    int(PolymerFluxArchetype.UNRESOLVED): "legacy_mu1/1",
    int(PolymerFluxArchetype.DISCRETE_CHIP): "discrete_chip/1",
}

_ARRHENIUS_A_UNITS = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}


def _resolve_reaction_pools(rxn, pool_set):
    """Mirror the solver's src/dst pool resolution (polymer.pyx:535-556):
    src = first reactant slot in a configured pool; dst = first cross-pool
    product, falling back to the same-pool fold-back product."""
    src = None
    for s in rxn.reactants:
        b = _species_base_label(s)
        if b in pool_set:
            src = b
            break
    dst = None
    for s in rxn.products:
        b = _species_base_label(s)
        if b not in pool_set:
            continue
        if b != src:
            dst = b
            break
        if dst is None:
            dst = b
    return src, dst


def compile_polymer_reaction_entries(core_reactions, core_species,
                                     configured_pool_labels,
                                     cantera_index_map=None):
    """Compile stamped proxy-touching core reactions into schema-2.0
    ``reactions[]`` entries (docs/polymer_moments_format.md §3).

    Mirrors the solver's init-time pool resolution and demotions
    (rmgpy/solver/polymer.pyx:527-588): the artifact describes what the
    oracle DOES, including legacy/unresolved fallbacks (design spec Q3).

    Parameters
    ----------
    cantera_index_map : dict, optional
        ``{id(rxn): [entry indices in chem.yaml's reactions list]}`` from
        ``rmgpy.cantera.generate_cantera_data(..., return_reaction_index_map=True)``.
        Reactions absent from the map are emitted with ``cantera: null``
        (the unbalanced-proxy filter dropped them) and MUST carry kinetics.
    """
    from rmgpy.cantera import get_reaction_equation
    from rmgpy.kinetics import Arrhenius as _Arrhenius

    cantera_index_map = cantera_index_map or {}
    pool_set = set(configured_pool_labels)
    entries = []
    dropped_counters: Dict[tuple, int] = {}

    NONE_ = int(PolymerFluxArchetype.NONE)
    MIG = int(PolymerFluxArchetype.MIGRATION)
    SCI = int(PolymerFluxArchetype.SCISSION_FRAGMENT)
    UNR = int(PolymerFluxArchetype.UNRESOLVED)
    CHIP = int(PolymerFluxArchetype.DISCRETE_CHIP)

    for rxn in core_reactions:
        arch = int(getattr(rxn, "polymer_flux_archetype", 0))
        src, dst = _resolve_reaction_pools(rxn, pool_set)
        if arch == NONE_ and src is None and dst is None:
            continue  # ordinary chemistry — Cantera handles it untouched

        # Mirror solver demotions (polymer.pyx:557-578).
        unresolved = False
        if arch == NONE_:
            arch, unresolved = UNR, True
        elif arch == UNR:
            unresolved = True
        elif arch in (MIG, SCI) and (src is None or dst is None):
            arch, unresolved = UNR, True
        elif arch == CHIP and src is None:
            arch, unresolved = UNR, True

        equation = get_reaction_equation(rxn, core_species)
        indices = cantera_index_map.get(id(rxn))
        if indices:
            if len(indices) > 1:
                logging.warning(
                    "Polymer artifact: reaction %s maps to %d Cantera entries; "
                    "emitting the first index. The consumer must zero ALL "
                    "duplicate entries (see format doc §4 step 0).",
                    equation, len(indices))
            cantera = {"index": int(indices[0]), "equation": equation}
            entry_id = f"r{int(indices[0])}"
        else:
            cantera = None
            family = str(getattr(rxn, "family", None)
                         or getattr(rxn, "library", None) or "rxn")
            key = (family, equation)
            occurrence = dropped_counters.get(key, 0)
            dropped_counters[key] = occurrence + 1
            entry_id = f"d{family}:{equation}:{occurrence}"

        kin = getattr(rxn, "kinetics", None)
        # Fix 2: use exact type check so SurfaceArrhenius (and other
        # subclasses) fall through to kinetics=None instead of emitting
        # volumetrically wrong units via _ARRHENIUS_A_UNITS.
        if type(kin) is _Arrhenius:
            # Fix 1: fold T0 into A to match the T0=1 convention used by
            # both Cantera's ArrheniusRate and the artifact consumer.
            # Mirrors Arrhenius.to_cantera_kinetics (arrhenius.pyx:259-262):
            #   A_folded = A.value_si / T0.value_si ** n.value_si
            kinetics = {
                "A": float(kin.A.value_si / (kin.T0.value_si ** kin.n.value_si)),
                "n": float(kin.n.value_si),
                "Ea": float(kin.Ea.value_si),
                "units": {"A": _ARRHENIUS_A_UNITS.get(len(rxn.reactants), "SI"),
                          "Ea": "J/mol"},
                "reversible": bool(rxn.reversible),
            }
        else:
            kinetics = None
            # Fix 3: warn for ALL non-Arrhenius entries (retained or dropped).
            # Retained entries still carry kinetics=null; consumers needing
            # reversibility must treat them as reversible because chem.yaml
            # always writes <=> even for irreversible reactions.
            logging.warning(
                "Polymer artifact: reaction %s has non-Arrhenius kinetics "
                "(%s); entry carries kinetics=null (no A/n/Ea and no "
                "reversible flag — consumers needing reversibility must "
                "treat it as reversible, matching chem.yaml).",
                equation, type(kin).__name__)

        entry = {
            "id": entry_id,
            "cantera": cantera,
            "kinetics": kinetics,
            "reactants": [_artifact_species_label(s) for s in rxn.reactants],
            "products": [_artifact_species_label(s) for s in rxn.products],
            "proxy_reactants": [_artifact_species_label(s) for s in rxn.reactants
                                if _species_base_label(s) in pool_set],
            "proxy_products": [_artifact_species_label(s) for s in rxn.products
                               if _species_base_label(s) in pool_set],
            "scaling": "mu0" if getattr(rxn, "is_end_group_reaction", False) else "mu1",
            "src_pool": src,
            "dst_pool": dst,
            "archetype": ARCHETYPE_TERM_NAMES[arch],
            "unresolved": unresolved,
        }
        if arch == CHIP:
            entry["params"] = {"a": int(getattr(rxn, "polymer_chip_units", 0))}
        entries.append(entry)
    return entries


def build_polymer_moments_artifact(pool_registry,
                                   core_species=None,
                                   core_reactions=None,
                                   configured_pool_labels=None,
                                   condensed_species=None,
                                   monomer_routing_by_pool=None,
                                   cantera_index_map=None,
                                   iteration=0,
                                   rmg_commit=None):
    """Assemble the full schema-2.0 polymer moments artifact payload.

    Normative contract: docs/polymer_moments_format.md. The payload mirrors
    the HybridPolymerSystem oracle, including its init-time demotions —
    ``configured_pool_labels`` must be the SOLVER-configured pools (which can
    be a subset of ``pool_registry``: spawned daughters are registry pools
    without solver configs and run as ordinary species).

    ``core_species``/``monomer_routing`` are populated by
    ``build_polymer_moments_artifact`` in the live run (legacy callers omit
    them → ``phase_species: []``, ``monomer_routing: null``). The caller must
    pass a CONDENSED-phase routing label (it is appended to phase_species
    unchecked by ``_serialize_pool_for_sidecar``).
    """
    if configured_pool_labels is None:
        configured_pool_labels = [getattr(p, "label", "") for p in pool_registry]
    monomer_routing_by_pool = monomer_routing_by_pool or {}

    pools = [
        _serialize_pool_for_sidecar(
            p,
            core_species=core_species,
            monomer_routing=monomer_routing_by_pool.get(getattr(p, "label", "")),
        )
        for p in pool_registry
    ]
    reactions = compile_polymer_reaction_entries(
        core_reactions or [], core_species or [],
        configured_pool_labels, cantera_index_map)

    conventions = {
        "format_doc": "docs/polymer_moments_format.md (polymer_moments_format/2.0)",
        "moment_basis": "extensive mol, DP basis (mu1 = moles of repeat units)",
        "volumes": {
            "V_poly": "constant, consumer-supplied [m^3]",
            "V_gas": "ideal gas, dynamic: V_gas = n_gas*R*T/P (1.0 m^3 floor when n_gas <= 0)",
        },
        "configured_pools": list(configured_pool_labels),
        "condensed_species": sorted(_artifact_species_label(s)
                                    for s in (condensed_species or [])),
        "site_scaling": ("site = max(0, mu_scaling)/V_poly read from the first proxy "
                         "reactant's pool; multiplies ONCE; scales rf AND rr"),
        "chip_site_throttle": ("site = min(max(0,mu0), max(0,mu1)/a)/V_poly when "
                               "archetype=discrete_chip/1 and scaling=mu0 and a>0"),
        "kb_recipe": ("kb = kf/Keq; Keq(T) = (P0/(R*T))^dn_gas * exp(-dG0/(R*T)), "
                      "P0 = 1e5 Pa, dG0 from chem.yaml NASA thermo"),
        "mu3_closure": "log_lagrange/1",
        "invariants": {
            "discrete_subset": ("sum_pools(mu1) + sum_chip_species(a_i * n_i) is "
                                "invariant over the discrete-reaction subset only"),
            "with_unzip": ("add + n(monomer_routing) per pool with an active unzip "
                           "channel (unzip moves units from mu1 into that species)"),
        },
    }

    return {
        "schema_version": POLYMER_POOLS_SIDECAR_SCHEMA_VERSION,
        "generated_at": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "rmg_commit": rmg_commit if rmg_commit is not None else _get_rmg_commit(),
        "rmg_iteration": int(iteration),
        "conventions": conventions,
        "pools": pools,
        "reactions": reactions,
    }


def schulz_flory_mu2(mu0: float, mu1: float) -> float:
    """Predict μ₂ from (μ₀, μ₁) under a Schulz-Flory chain-length distribution.

    For P(n) = p^(n-1)·(1-p), the analytic relations are
    μ₀=N, μ₁=N/(1-p), μ₂=N·(1+p)/(1-p)². Eliminating ``p`` gives the
    closed form ``μ₂ = 2·μ₁²/μ₀ − μ₁``.

    Used in inter-pool transfer-reaction moment effects (design doc §5)
    when the second-moment source term cannot be read off directly.
    """
    if mu0 == 0:
        return 0.0
    return 2.0 * mu1 * mu1 / mu0 - mu1


def drain_spawn_intents(
    intents: List[SpawnIntent],
    iteration: int,
    existing_pools: Optional[List['Polymer']] = None,
) -> List['Polymer']:
    """Materialise queued :class:`SpawnIntent`s into new :class:`Polymer` pools.

    Iteration-boundary hook (design doc §4.5). Each returned Polymer carries
    ``parent_pool_label``, ``spawn_iteration``, and ``end_groups_str`` so the
    sidecar writer can serialise lineage without re-deriving it.

    ``existing_pools`` is consulted to namespace daughter labels — the
    n-th daughter of a given parent across all calls gets ``<parent>_d{n}``,
    preventing collisions when the registry grows across iterations.

    NOTE: there is no in-place state-vector resize or ``CVodeReInit``. The
    daughter pools' ``_muN`` dummy species are registered with the reaction
    model; the next RMG iteration rebuilds the solver and
    ``HybridPolymerSystem.initialize_model`` resolves their moment indices by
    label. (The ``mu_indices`` computed here are only used for the sidecar JSON.)
    """
    taken: Set[str] = {p.label for p in (existing_pools or [])}
    # μ_index allocator: next free slot is one past the max index any
    # existing pool already claims. Pools without explicit mu_indices are
    # assumed to occupy [0, 1, 2] (root-pool convention).
    max_idx = -1
    for p in (existing_pools or []):
        mi = getattr(p, "mu_indices", None)
        if mi is None:
            max_idx = max(max_idx, 2)
        else:
            try:
                max_idx = max(max_idx, max(int(x) for x in (
                    mi.values() if isinstance(mi, dict) else mi
                )))
            except (TypeError, ValueError):
                continue
    next_idx = max_idx + 1
    new_pools: List['Polymer'] = []
    for intent in intents:
        parent = intent.parent_pool
        n = 1
        while f"{parent.label}_d{n}" in taken:
            n += 1
        new_label = f"{parent.label}_d{n}"
        taken.add(new_label)
        # Daughter monomer: kept as the parent's monomer (a labeled-radical
        # Molecule passes Polymer.__init__ validation). The true detected
        # motif lives on intent.monomer (a Group) and is referenced via
        # spawn_metadata for downstream chemistry; the placeholder here is
        # just a valid Molecule the constructor will accept.
        new_pool = Polymer(
            label=new_label,
            monomer=parent.monomer,
            end_groups=intent.end_groups,
            cutoff=parent.cutoff,
            Mn=parent.Mn,
            Mw=parent.Mw,
            initial_mass=0.001,
        )
        # Override fingerprint so _register_polymer's dedup sees the daughter
        # as distinct from the parent (which shares monomer + end_groups +
        # cutoff and would otherwise hash to the same fingerprint).
        new_pool._fingerprint = f"{parent.fingerprint}_daughter-{new_label}"
        # Event-spawn initialisation (design doc §5): μ_k = N · DP^k.
        N = float(intent.triggering_moles)
        DP = float(intent.triggering_dp)
        new_pool.moments = np.array([N, N * DP, N * DP * DP], dtype=np.float64)
        new_pool.parent_pool_label = parent.label
        new_pool.spawn_iteration = iteration
        new_pool.end_groups_str = list(intent.end_groups)
        new_pool.mu_indices = (next_idx, next_idx + 1, next_idx + 2)
        next_idx += 3
        new_pool.spawn_metadata = {
            "triggering_dp": int(intent.triggering_dp),
            "triggering_moles": N,
            "mass_flux_at_spawn": float(intent.mass_flux_at_spawn),
        }
        tp = intent.triggering_product
        if tp is not None:
            try:
                if getattr(tp, "molecule", None) and tp.molecule:
                    new_pool.spawn_metadata["triggering_product_smiles"] = (
                        tp.molecule[0].to_smiles()
                    )
            except Exception:
                pass
        new_pools.append(new_pool)
    return new_pools


def apply_spawn_intents(
    reaction_model: Any,
    intents: List[SpawnIntent],
    iteration: int,
    existing_pools: Optional[List['Polymer']] = None,
) -> List['Polymer']:
    """Iteration-boundary entry point.

    Drains queued :class:`SpawnIntent`s into new :class:`Polymer` pools and
    registers each with ``reaction_model``. Returns the materialised pools
    so the caller can extend its pool registry.

    This is the single hook the RMG main loop will call between iterations
    (design doc §4.5). The next :meth:`HybridPolymerSystem.initialize_model`
    sees the expanded core species list and grows the state vector
    automatically.
    """
    new_pools = drain_spawn_intents(intents, iteration, existing_pools=existing_pools)
    register_spawned_pools(reaction_model, new_pools)
    return new_pools


def register_spawned_pools(
    reaction_model: Any,
    new_pools: List['Polymer'],
) -> None:
    """Register each daughter ``Polymer`` with the reaction model.

    Iteration-boundary glue (design doc §4.5). ``reaction_model.make_new_species``
    (which dispatches to ``_register_polymer`` for ``Polymer`` arguments)
    handles label disambiguation and the automatic creation of the three
    ``_mu0`` / ``_mu1`` / ``_mu2`` moment-dummy core species. The next call
    to :meth:`HybridPolymerSystem.initialize_model` then picks them up via
    the standard polymer-pool resolution path — no Cython hot-reinit needed.
    """
    for poly in new_pools:
        reaction_model.make_new_species(poly)


def _tag_polymer_proxy(cand: 'Species', *, is_proxy: bool) -> None:
    """Stamp an ``is_polymer_proxy`` flag on a Species and its Molecules.

    Mirrors the tagging block inside :func:`process_polymer_candidates`.
    """
    if not hasattr(cand, "props"):
        cand.props = {}
    cand.props["is_polymer_proxy"] = is_proxy
    cand.is_polymer_proxy = is_proxy
    if getattr(cand, "molecule", None):
        for m in cand.molecule:
            if not hasattr(m, "props"):
                m.props = {}
            m.props["is_polymer_proxy"] = is_proxy
            m.is_polymer_proxy = is_proxy


def similarity_merge(
    candidate: Union[Group, Molecule],
    pool_registry: List['Polymer'],
) -> Optional['Polymer']:
    """Return an existing pool whose monomer pattern matches ``candidate``.

    Used as the first phase of the spawn-trigger pipeline (design doc §4.3):
    before treating a detected motif as novel, check whether it is already
    represented by a live pool. If so, the caller should extend that pool's
    feature_monomer set rather than spawning a new pool.

    Returns ``None`` if no pool matches.
    """
    if not pool_registry:
        return None

    for pool in pool_registry:
        try:
            pool_pattern = pool.backbone_group
        except Exception:
            continue
        if pool_pattern is None:
            continue

        if isinstance(candidate, Group):
            try:
                if candidate.is_isomorphic(pool_pattern):
                    return pool
            except (NotImplementedError, AttributeError, ValueError):
                pass
        elif isinstance(candidate, Molecule):
            # A whole molecule "merges" if it contains the pool's backbone as
            # a subgraph. Best-effort: gives the right answer for the synthetic
            # tests; the production caller passes Group from discover_repeat_motif.
            try:
                if candidate.is_subgraph_isomorphic(pool_pattern):
                    return pool
            except (NotImplementedError, AttributeError, ValueError):
                pass
    return None
