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

import numpy as np
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Tuple, Union

from rmgpy.exceptions import InputError
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.fragment import Fragment
from rmgpy.molecule.resonance import generate_resonance_structures
from rmgpy.species import Species
from rmgpy.thermo import Wilhoit, NASA, ThermoData

if TYPE_CHECKING:
    from rmgpy.molecule.group import Group


LABELS_1, LABELS_2 = ('1', '*1'), ('2', '*2')


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
        super(Polymer, self).__init__(label=label, **kwargs)

        self.monomer = self._validate_monomer(monomer, label)
        if feature_monomer:
            self.feature_monomer = self._validate_monomer(feature_monomer, label)
        else:
            self.feature_monomer = None
        self._validate_end_groups(end_groups, label)
        self.cutoff = self._validate_cutoff(cutoff, label)
        self.Mn, self.Mw, self.moments = None, None, None

        self.k_unzip = kwargs.get('k_unzip', 0.0)
        self.k_scission = kwargs.get('k_scission', 0.0)

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

    def get_symmetry_number(self):
        """Delegates symmetry calculation to the proxy species."""
        return self.get_proxy_species().get_symmetry_number()

    def get_net_charge(self):
        """Delegates charge calculation to the proxy."""
        return self.get_proxy_species().get_net_charge()

    @property
    def fingerprint(self):
        """Fingerprint of this polymer, taken from molecule attribute. Read-only."""
        if self._fingerprint is None:
            if self.monomer:
                feat = f'_Feat-{self.feature_monomer.fingerprint}' if self.feature_monomer else ''
                self._fingerprint = f'Polymer_{self.monomer.fingerprint}{feat}_{self.cutoff}'
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
        """
        new_poly = self._create_reacted_copy_logic(reacted_proxy)
        if new_poly is None:
            return None
        proxy_spec = new_poly.get_proxy_species()
        for mol in proxy_spec.molecule:
            for atom in mol.atoms:
                atom.label = ''
            mol.update()
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
            new_Mn = self.Mn / 2.0 if self.Mn else None
            new_Mw = self.Mw / 2.0 if self.Mw else None
            new_moments = self.moments.copy() if self.moments is not None else None
            return Polymer(label=f"{self.label}_scission_tail",
                           monomer=self.monomer,
                           feature_monomer=None,
                           end_groups=[self.end_groups[0].copy(deep=True), new_tail],
                           cutoff=self.cutoff,
                           Mn=new_Mn,
                           Mw=new_Mw,
                           initial_mass=0.0,
                           moments=new_moments,
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

    def _wing_pattern(self, side: str) -> Molecule:
        """
        Generates a purely structural "wing" template, [Head-Base] or [Base-Tail],
        by stripped labels for isomorphism checks.

        Args:
            side (str): 'head' or 'tail' wing.

        Returns:
            Molecule: The wing pattern for subgraph matching.
        """
        wing = self._stitch_wing(side=side)
        for atom in wing.atoms:
            atom.label = ''
            atom.charge = 0
            while atom.radical_electrons > 0:
                atom.decrement_radical()
        wing.update_multiplicity()
        return wing

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

    def _wing_groups(self, side: str):
        """
        Return wing patterns as Groups, using resonance structures to get Clar & Kekulé (S/D/B) variants.
        Prune duplicate Group patterns by adjacency-list string.
        Then relax radical constraints so an open-shell wing can match a closed-shell proxy.
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

    def _pick_disjoint_pair(self,
                            head_matches: list[dict],
                            tail_matches: list[dict],
                            target: Molecule,
                            ) -> tuple[set, set]:
        """
        Pick a head-wing match and a tail-wing match that are disjoint in the *target*.
        Selection heuristic:
          1) boundary penalty: prefer each wing having ~1 boundary edge to remainder
          2) remainder connectedness (prefer connected)
          3) maximize total covered atoms |H|+|T|
          4) maximize min(|H|,|T|)
          5) minimize remainder size (prefer the smallest plausible feature)

        Args:
            head_matches (list): List of head wing subgraph matches.
            tail_matches (list): List of tail wing subgraph matches.
            target (Molecule): The molecule being analyzed.

        Returns:
            Tuple[set, set]: The atom sets for the head and tail wings.
        """
        if not head_matches or not tail_matches:
            raise ValueError("Need both head and tail matches to pick a disjoint pair.")

        def boundary_edges(atom_set: set) -> int:
            """Count unique undirected edges from atom_set to outside."""
            edges = set()
            for a in atom_set:
                for b in a.bonds:
                    if b not in atom_set:
                        edges.add(frozenset((a, b)))
            return len(edges)

        def is_connected(atom_set: set) -> bool:
            """Connectivity in the induced subgraph on atom_set using Molecule bonds."""
            if not atom_set:
                return False
            start = next(iter(atom_set))
            stack = [start]
            seen = {start}
            while stack:
                a = stack.pop()
                for b in a.bonds:
                    if b in atom_set and b not in seen:
                        seen.add(b)
                        stack.append(b)
            return len(seen) == len(atom_set)

        best_pair, best_key = None, None

        for hm in head_matches:
            H = get_target_atoms(hm)
            if any(a not in target.atoms for a in H):
                raise ValueError("Unexpected mapping direction in head match: keys are not target atoms.")

            for tm in tail_matches:
                T = get_target_atoms(tm)
                if any(a not in target.atoms for a in T):
                    raise ValueError("Unexpected mapping direction in tail match: keys are not target atoms.")

                if H & T:
                    continue

                remainder = set(target.atoms) - (H | T)
                if not remainder:
                    continue

                bh, bt = boundary_edges(H), boundary_edges(T)
                boundary_penalty = abs(bh - 1) + abs(bt - 1)
                boundary_score = -boundary_penalty
                rem_connected = is_connected(remainder)
                total = len(H) + len(T)
                min_side = min(len(H), len(T))
                rem_size = len(remainder)
                key = (boundary_score, rem_connected, total, min_side, -rem_size)
                if best_key is None or key > best_key:
                    best_key = key
                    best_pair = (H, T)

        if best_pair is None:
            raise ValueError("Could not find disjoint head/tail wing matches.")
        return best_pair

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

    def _restore_labels(self,
                        new_mol: Molecule,
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

    def _assert_feature_unit(self, mol: Molecule):
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
    Classification of reaction products relative to the original polymer structure.
    """
    GAS = 'GAS'           # Not a polymer species
    FEATURE = 'FEATURE'   # Polymer with modified center monomer
    DISCARD = 'DISCARD'   # Polymer with modified buffer monomer (ignore)
    SCISSION = 'SCISSION' # Scission fragment (requires moment update)
    END_MOD = 'END_MOD'   # Polymer with modified end-group (modified head or tail, backbone intact)


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
    Safely extracts the target molecule atoms from an RMG match mapping.
    Uses '.element' because both Atom and GroupAtom possess a '.bonds' attribute.
    """
    if not match:
        return set()
    vals = set(match.values())
    if vals and isinstance(next(iter(vals)), Atom):
        return vals
    keys = set(match.keys())
    if keys and isinstance(next(iter(keys)), Atom):
        return keys
    return {x for x in vals | keys if isinstance(x, Atom)}


def find_max_disjoint_matches(matches: List[Dict[Any, Any]]) -> List[Dict[Any, Any]]:
    """
    Return the largest subset of pairwise-disjoint isomorphism matches.

    A "match" is a dict whose values are target atoms (or atom-like objects). Two matches
    overlap if they share any value. This function finds the maximum-cardinality subset
    with no overlaps (avoids greedy undercounting for symmetric/overlapping matches).
    """
    if not matches:
        return []
    matches = sorted(matches, key=lambda m: len(m), reverse=True)

    def solve(candidates: List[Dict[Any, Any]]) -> List[Dict[Any, Any]]:
        if not candidates:
            return []
        first = candidates[0]
        first_atoms = get_target_atoms(first)
        remaining_compatible = [m for m in candidates[1:] if first_atoms.isdisjoint(get_target_atoms(m))]
        with_first = [first] + solve(remaining_compatible)
        without_first = solve(candidates[1:])
        return with_first if len(with_first) >= len(without_first) else without_first

    return solve(matches)


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


def count_backbone_matches(product_mol: Molecule, monomer_group: 'Group') -> MatchSummary:
    """
    Abstracts the subgraph isomorphism and disjoint filtering logic.
    """
    matches = product_mol.find_subgraph_isomorphisms(monomer_group)
    raw_count = len(matches)

    best_disjoint_set = find_max_disjoint_matches(matches)
    return MatchSummary(
        raw=raw_count,
        disjoint=len(best_disjoint_set),
        best_matches=best_disjoint_set
    )


def classify_structure(species: 'Species',
                       original_polymer,
                       *,
                       monomer_group: Optional['Group'] = None,
                       use_proxy_as_gate: bool = False,
                       ) -> Tuple[PolymerClass, Dict[str, Any]]:
    """
    Classifies a reaction product structurally by counting intact monomer subgraphs.
    """
    base_details = {"raw_matches": 0, "disjoint_matches": 0}

    # 1. Gates and Fail-Fasts
    if use_proxy_as_gate:
        spc_flag = True
        if hasattr(species, "props") and "is_polymer_proxy" in species.props:
            spc_flag = species.props["is_polymer_proxy"]
        else:
            spc_flag = getattr(species, "is_polymer_proxy", True)
        mol_flag = True
        if getattr(species, "molecule", None) and species.molecule:
            mol = species.molecule[0]
            if hasattr(mol, "props") and "is_polymer_proxy" in mol.props:
                mol_flag = mol.props["is_polymer_proxy"]
            else:
                mol_flag = getattr(mol, "is_polymer_proxy", True)

        if spc_flag == False or mol_flag == False:
            return PolymerClass.GAS, {**base_details, "reason": "proxy_flag_false"}

    if not getattr(species, "molecule", None) or not species.molecule[0]:
        return PolymerClass.GAS, {**base_details, "reason": "no_molecule"}
    if monomer_group is None:
        monomer_group = getattr(original_polymer, "backbone_group", None)
    if monomer_group is None:
        return PolymerClass.GAS, {**base_details, "reason": "no_backbone_group"}

    # 2. Minimal Safe Normalization
    product_mol = species.molecule[0].copy(deep=True)
    product_mol.clear_labeled_atoms()
    product_mol.update_multiplicity()

    # 3. Search and Count
    summary = count_backbone_matches(product_mol, monomer_group)
    details: Dict[str, Any] = {
        "raw_matches": summary.raw,
        "disjoint_matches": summary.disjoint
    }

    # 4. Type Guard
    if summary.disjoint > 0:
        sample_atoms = get_target_atoms(summary.best_matches[0])
        sample = next(iter(sample_atoms))
        if not hasattr(sample, "bonds"):
            raise TypeError(f"Expected Atom objects in match mapping, got {type(sample).__name__}")
        details["match_value_type"] = type(sample).__name__

    # 5. Classification
    if summary.disjoint >= 3:
        if summary.disjoint > 3:
            details["note"] = "more_than_3_matches"
        return PolymerClass.END_MOD, details

    if summary.disjoint == 2:
        m1_atoms = get_target_atoms(summary.best_matches[0])
        m2_atoms = frozenset(get_target_atoms(summary.best_matches[1]))
        is_connected = any(
            nbr in m2_atoms
            for atom in m1_atoms
            for nbr in atom.bonds.keys()
        )
        details["connected"] = is_connected
        return (PolymerClass.DISCARD if is_connected else PolymerClass.FEATURE), details

    if summary.disjoint == 1:
        return PolymerClass.SCISSION, details

    return PolymerClass.GAS, {**details, "reason": "no_backbone_matches"}


def process_polymer_candidates(candidates: List[Species],
                               _reaction_model,
                               original_polymer,
                               ) -> List[Species]:
    """
    Handshake function to convert generic Species into Polymer objects.
    """
    processed_list: List['Species'] = []
    stats = {k: 0 for k in PolymerClass}

    monomer_group = getattr(original_polymer, "backbone_group", None)

    for cand in candidates:
        classification, details = classify_structure(
            cand,
            original_polymer,
            monomer_group=monomer_group,
            use_proxy_as_gate=False,
        )

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

def _find_disjoint_pair(head_matches, tail_matches):
    """
    Returns (head_atoms, tail_atoms) for the first disjoint head/tail match pair.
    Uses get_target_atoms to ensure we are comparing target Molecule Atom objects.
    """
    for hm in head_matches:
        ha = get_target_atoms(hm)
        for tm in tail_matches:
            ta = get_target_atoms(tm)
            if ha.isdisjoint(ta):
                return ha, ta
    return set(), set()
