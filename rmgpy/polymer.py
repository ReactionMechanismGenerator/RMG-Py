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
import math
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
# Temporary label used by _stitch_end_radical_oligomer to keep an open chain
# terminus out of the stitcher's *1/*2 label resolution (restored before the
# strict end-radical assertion runs).
_END_RADICAL_PARKED_LABEL = '*er'


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
        radical_qssa_unzip (Optional[dict]): Radical QSSA unzip channel config
            (initiation/depropagation/termination Arrhenius triplets {A, n, Ea},
            optional transfer triplet, efficiency, monomer_yield, basis).
            Passive storage on the Polymer object: validation lives in the
            polymer() deck helper, PolymerPool.to_config and the solver
            (validate_radical_qssa_unzip in rmgpy/solver/polymer.pyx documents
            the full contract). Mutually exclusive with k_unzip > 0. M1: the
            stored config is inert (no RHS reads it until the M2 rate law).
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
                 k_unzip: float = 0.0,
                 k_scission: float = 0.0,
                 radical_qssa_unzip: Optional[dict] = None,
                 k_homolysis: Optional[dict] = None,
                 k_depropagation: Optional[dict] = None,
                 end_radical_site: Optional[str] = None,
                 side_group_homolysis: Optional[list] = None,
                 side_loss_channel: Optional[str] = None,
                 **kwargs,
                 ):
        # k_unzip/k_scission/radical_qssa_unzip/k_homolysis/end_radical_site
        # are named parameters, so they never
        # appear in kwargs; assign them directly. discrete_dp_threshold (below) is
        # popped from kwargs before Species.__init__, which does not accept it and
        # would raise TypeError.
        self.k_unzip = k_unzip
        self.k_scission = k_scission
        self.radical_qssa_unzip = radical_qssa_unzip
        # Radical-homolysis initiation kernel (Stage 1, adjudicated round 66):
        # Arrhenius triplet {A, n, Ea} (SI: A [s^-1], Ea [J/mol]) or None.
        # Passive storage here (like radical_qssa_unzip); validated by
        # rmgpy.solver.polymer.validate_k_homolysis at the deck / config /
        # solver layers. NOT a bare scalar by design: the solver evaluates
        # k(T) at its runtime temperature.
        self.k_homolysis = k_homolysis
        # End-radical depropagation kernel (adjudicated round 74 SS2):
        # Arrhenius triplet {A, n, Ea} (SI: A [s^-1], Ea [J/mol]) or None.
        # Passive storage here (like k_homolysis); validated by
        # rmgpy.solver.polymer.validate_k_depropagation at the deck /
        # config / solver layers. Declared on the PARENT pool's k_homolysis
        # context; the producer (generate_end_radical_daughters) copies it
        # onto both spawned end-radical daughter pools, where the kernel
        # actually runs.
        self.k_depropagation = k_depropagation
        # end_radical_site marks this Polymer as one of the two homolysis
        # end-radical daughter pools: 'primary' = the open-*1 chain terminus,
        # 'secondary' = the open-*2 terminus of the backbone C-C cut. When
        # set, the reactive proxy is the end-radical oligomer (mono-radical
        # on the uncapped terminal heavy atom) and the fingerprint carries an
        # _EndRad segment so the pools stay distinct from each other, from
        # the parent, and from mid-chain _mod feature pools.
        if end_radical_site not in (None, 'primary', 'secondary'):
            raise InputError(
                f"Polymer '{label}': end_radical_site must be None, "
                f"'primary' or 'secondary', got {end_radical_site!r}.")
        self.end_radical_site = end_radical_site
        # Side-group homolysis initiation kernel (FR1-K1, adjudicated round
        # 70): LIST of channel dicts (each exactly {label, A, n, Ea,
        # site_selector, sites_per_unit, gas_product}; site_selector is the
        # round-72 REQUIRED structural site selector) or None. Passive
        # storage here (like k_homolysis); validated by
        # rmgpy.solver.polymer.validate_side_group_homolysis at the deck /
        # config / solver layers. The producer
        # (generate_side_loss_daughters) spawns ONE X-loss feature pool per
        # channel at model setup.
        self.side_group_homolysis = side_group_homolysis
        # side_loss_channel marks this Polymer as the X-loss feature pool of
        # one (parent, channel) pair: the raw channel label (e.g.
        # 'aliphatic_C-Br'). When set, the fingerprint carries a _SideLoss
        # segment keyed on the sanitized channel label + the feature unit's
        # radical environment, so channels NEVER lump (round-70 ruling) and
        # the pool stays distinct from the parent, from _EndRad daughters
        # and from S2 H-loss _Feat/_mod pools.
        if side_loss_channel is not None and (
                not isinstance(side_loss_channel, str)
                or not side_loss_channel.strip()):
            raise InputError(
                f"Polymer '{label}': side_loss_channel must be None or a "
                f"non-empty channel-label string, got {side_loss_channel!r}.")
        self.side_loss_channel = side_loss_channel
        # Exact per-chain mass defect [g/mol] of an X-loss feature pool
        # (FR1-K1 mass contract): every feature chain lost exactly ONE X
        # (v1, no multi-loss cascade), so the pool's exact condensed mass is
        # mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol. Pinned to M_X
        # of the spawning channel's gas_product by
        # generate_side_loss_daughters; 0.0 on ordinary pools.
        self.chain_mass_defect_g_mol = 0.0
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
        self._end_radical_proxy = None
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
        if self.end_radical_site is not None:
            # End-radical daughter pool: the reactive proxy IS the
            # end-radical oligomer (built + strictly asserted in
            # end_radical_proxy), never the closed-shell baseline trimer --
            # otherwise RMG would dedup the daughter against the parent's
            # species graph.
            reactive_proxy = self.end_radical_proxy
        else:
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
        """Fingerprint of this polymer, taken from molecule attribute. Read-only.

        The ``_Feat-`` segment is the feature monomer's element-count
        fingerprint; for RADICAL feature units (radical-feature producer
        path: more than the two stitch radicals) a ``-rad...`` radical-site
        descriptor is appended, because the element-count string alone
        cannot distinguish H-loss units of the same formula (e.g. the ~3 PP
        C3H5 units) and ``_register_polymer`` dedups pools BY fingerprint --
        distinct abstraction environments must stay distinct pools while
        positional twins (identical unit graph) still collapse to one.
        Classic 2-radical features are byte-identical to before.
        """
        if self._fingerprint is None:
            if self.monomer:
                feat = ''
                side_ch = getattr(self, 'side_loss_channel', None)
                if self.feature_monomer and side_ch:
                    # Side-group X-loss feature pools (FR1-K1): a DISTINCT
                    # fingerprint class keyed on the sanitized CHANNEL label
                    # (round-70 ruling: channels never lump -- two channels
                    # of the same element with identical unit graphs must
                    # stay distinct pools) plus the unit's radical-site
                    # environment. The _SideLoss literal never collides with
                    # _Feat (S2 H-loss) or _EndRad (k_homolysis) segments.
                    from rmgpy.solver.polymer import \
                        sanitize_side_group_channel_label
                    feat = (f'_SideLoss-'
                            f'{sanitize_side_group_channel_label(side_ch)}'
                            f'-rad{_radical_site_descriptor(self.feature_monomer)}')
                elif self.feature_monomer:
                    feat = f'_Feat-{self.feature_monomer.fingerprint}'
                    if self.feature_monomer.get_radical_count() > 2:
                        feat += f'-rad{_radical_site_descriptor(self.feature_monomer)}'
                if getattr(self, 'end_radical_site', None):
                    # End-radical daughter pools (k_homolysis conduit): the
                    # _EndRad segment keys the pool identity on WHICH chain
                    # terminus carries the radical (site name) plus the
                    # radical atom's environment descriptor -- distinct from
                    # the parent (no segment), from each other (different
                    # site + environment) and from mid-chain _mod pools
                    # (_Feat segment).
                    er_mol = self.end_radical_proxy.molecule[0]
                    feat += (f'_EndRad-{self.end_radical_site}'
                             f'-rad{_radical_site_descriptor(er_mol)}')
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
    def end_radical_proxy(self) -> Optional[Species]:
        """The cached end-radical oligomer proxy (k_homolysis daughter pools
        only; None unless end_radical_site is set). Built with ONE uncapped
        stitch terminus so the surviving stitch radical IS the chain-end
        radical; strictly asserted by _assert_end_radical_proxy."""
        site = getattr(self, 'end_radical_site', None)
        if site and getattr(self, '_end_radical_proxy', None) is None:
            self._end_radical_proxy = self._stitch_end_radical_oligomer(site)
        return getattr(self, '_end_radical_proxy', None)

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
        other._end_radical_proxy = (self._end_radical_proxy.copy(deep=True)
                                    if getattr(self, '_end_radical_proxy', None) else None)
        other._fingerprint = self._fingerprint
        # Attributes set in __init__ that __new__ bypasses — must be carried over,
        # else a copied Polymer loses its identity flag (is_polymer) and, worse,
        # its degradation kinetics (k_scission/k_unzip would silently reset to 0).
        other.k_scission = self.k_scission
        other.k_unzip = self.k_unzip
        # deepcopy, not shallow-assign (review round 21, finding 3): shallow
        # assignment aliases the nested channel dict across copies, so mutating
        # one Polymer's channel would silently rewrite every copy's.
        other.radical_qssa_unzip = deepcopy(getattr(self, 'radical_qssa_unzip', None))
        # k_homolysis kernel triplet + end-radical pool identity (Stage 1,
        # round 66): losing the triplet on copy would silently disable the
        # initiation kernel; losing end_radical_site would collapse a
        # daughter pool's proxy/fingerprint back onto the parent's.
        other.k_homolysis = deepcopy(getattr(self, 'k_homolysis', None))
        # k_depropagation kernel triplet (r74 SS2): losing it on copy would
        # silently disable the radical-end consumption channel (the run-6
        # no-outlet wall would reopen on the copied pool).
        other.k_depropagation = deepcopy(getattr(self, 'k_depropagation',
                                                 None))
        other.end_radical_site = getattr(self, 'end_radical_site', None)
        # Side-group homolysis kernel channels + X-loss feature-pool
        # identity + exact mass-defect contract (FR1-K1): losing the
        # channel list on copy would silently disable the kernel; losing
        # side_loss_channel would collapse a feature pool's fingerprint
        # back onto an S2 _Feat pool; losing chain_mass_defect_g_mol would
        # re-open the round-70 P1 mass-minting trap.
        other.side_group_homolysis = deepcopy(
            getattr(self, 'side_group_homolysis', None))
        other.side_loss_channel = getattr(self, 'side_loss_channel', None)
        other.chain_mass_defect_g_mol = getattr(
            self, 'chain_mass_defect_g_mol', 0.0)
        # Released-monomer routing target (input.py:432). Shared BY REFERENCE
        # deliberately, NOT deep-copied: routing resolution downstream
        # (derive_daughter_pool_configs' spc_map, PolymerPool.to_config) is
        # object-keyed, so identity with the core Species is load-bearing --
        # a deep copy would silently resolve to no core index at all.
        other.monomer_product_species = getattr(self, 'monomer_product_species', None)
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

        FR1-K2 mass audit note (round-72 P2): the Mass/MonomerMW inversion
        is deliberately NOT defect-aware. It seeds t=0 moments from a
        deck-declared initial mass, and only INPUT pools declare one --
        X-loss feature pools (chain_mass_defect_g_mol > 0) are
        producer-spawned born-at-zero (initial_mass=0.0 pinned by
        generate_side_loss_daughters), so this path never runs with a
        nonzero defect. Inverting the defect-aware formula would need mu0,
        which mass alone cannot determine.
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
        elif getattr(self, 'end_radical_site', None):
            # End-radical daughter pool: its reactive identity IS the
            # end-radical oligomer (mono-radical chain terminus).
            return self.end_radical_proxy
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

        # The repeat units stitched between the end-caps. This list IS the
        # proxy size: PROXY_STITCH_REPEAT_UNITS (and the impostor threshold
        # derived from it, r85 P2(b)) must move with it.
        repeat_units = [baseline, center, baseline]
        if len(repeat_units) != PROXY_STITCH_REPEAT_UNITS:
            raise ValueError(
                "Polymer._stitch_trimer stitches "
                f"{len(repeat_units)} repeat units but "
                f"PROXY_STITCH_REPEAT_UNITS = {PROXY_STITCH_REPEAT_UNITS}; "
                "update the constant together with the construction "
                "(the impostor threshold is derived from it).")

        trimer = head
        for unit in repeat_units:
            trimer = stitch_molecules_by_labeled_atoms(trimer, unit)
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

    def _stitch_end_radical_oligomer(self, site: str) -> Optional[Species]:
        """
        Constructs the end-radical oligomer proxy for a k_homolysis daughter
        pool (Stage 1, adjudicated round 66): a 3-unit chain with ONE
        uncapped stitch terminus, so the surviving stitch radical IS the
        chain-end radical (no H bookkeeping -- homolysis of the inter-unit
        *1-*2 backbone bond leaves exactly these two ends):

            primary   [HeadCap *1]-[*2 u *1]-[*2 u *1]-[*2 u *1(RADICAL)]
            secondary [(RADICAL)*2 u *1]-[*2 u *1]-[*2 u *1]-[*2 TailCap]

        The surviving stitch label is the TERMINALITY WITNESS: the strict
        assertion (_assert_end_radical_proxy -- exactly one radical, ON the
        labeled terminal heavy atom) runs BEFORE labels are cleared. Finishing
        recipe (update / resonance / is_polymer_proxy) matches _stitch_trimer
        so the proxy passes the same downstream validity expectations as the
        S2 feature proxies.
        """
        if site not in ('primary', 'secondary'):
            raise ValueError(
                f"Polymer '{self.label}': unknown end-radical site "
                f"{site!r} (expected 'primary' or 'secondary').")
        unit = self.monomer
        if site == 'primary':
            chain = stitch_molecules_by_labeled_atoms(
                self.end_groups[0].copy(deep=True), unit.copy(deep=True))
            if chain is None:
                return None
            chain = stitch_molecules_by_labeled_atoms(chain, unit.copy(deep=True))
            if chain is None:
                return None
            chain = stitch_molecules_by_labeled_atoms(chain, unit.copy(deep=True))
        else:
            # The open *2 terminus must survive three stitches, but
            # stitch_molecules_by_labeled_atoms resolves labels on the MERGED
            # graph (first match wins), so an unconsumed *2 on the left
            # operand would mis-bond the first unit onto itself. Park the
            # open-site label out of the stitcher's label sets for the
            # duration and restore it afterwards (the atom keeps its radical
            # electron throughout -- it is simply never bonded).
            u1 = unit.copy(deep=True)
            parked = next((a for a in u1.atoms if a.label in LABELS_2), None)
            if parked is None:
                return None
            parked.label = _END_RADICAL_PARKED_LABEL
            chain = stitch_molecules_by_labeled_atoms(u1, unit.copy(deep=True))
            if chain is None:
                return None
            chain = stitch_molecules_by_labeled_atoms(chain, unit.copy(deep=True))
            if chain is None:
                return None
            chain = stitch_molecules_by_labeled_atoms(
                chain, self.end_groups[1].copy(deep=True))
            if chain is None:
                return None
            restored = next((a for a in chain.atoms
                             if a.label == _END_RADICAL_PARKED_LABEL), None)
            if restored is None:
                return None
            restored.label = LABELS_2[-1]
        if chain is None:
            return None
        chain.update()
        chain.identify_ring_membership()
        # STRICT assertion path (round 66: do NOT relax _assert_feature_unit
        # / _assert_end_group) -- before the labels are cleared.
        _assert_end_radical_proxy(chain, site)
        spc = Species(molecule=[chain])
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

    def generate_end_radical_daughters(self) -> Tuple['Polymer', 'Polymer']:
        """
        Producer for the radical-homolysis initiation conduit (Stage 1,
        adjudicated round 66): build BOTH end-radical daughter pools of this
        parent, for registration at MODEL SETUP (not lazily).

        Returns ``(primary, secondary)``:
        - ``{label}_rad_primary_end``   -- open-*1 end radical,
        - ``{label}_rad_secondary_end`` -- open-*2 end radical
        of the backbone C-C cut. The suffix names are POSITIONAL (round 67
        ruling a): they identify WHICH stitch terminus carries the radical,
        not primary/secondary radical character -- that chemical reading
        holds only for PP under its head-to-tail repeat convention, is
        orientation-dependent, and is false for other polymers.

        Each daughter: born-at-zero moments (spawned_empty pattern, mirrors
        _born_at_zero_mod_daughter), parent's monomer (monomer_mw_g_mol
        pinned by __init__), spawned provenance markers, and an end-radical
        reactive proxy with its own strict assertion path.
        """
        daughters = []
        for site, suffix in (('primary', '_rad_primary_end'),
                             ('secondary', '_rad_secondary_end')):
            if self.Mn and self.Mw:
                dist_kwargs = dict(Mn=self.Mn, Mw=self.Mw, moments=None)
            else:
                dist_kwargs = dict(moments=[0.0, 0.0, 0.0])
            daughter = Polymer(
                label=f"{self.label}{suffix}",
                monomer=self.monomer,
                end_groups=[eg.copy(deep=True) for eg in self.end_groups],
                cutoff=self.cutoff,
                initial_mass=0.0,
                end_radical_site=site,
                **dist_kwargs,
            )
            daughter.parent_pool_label = self.label
            # Literal pinned as module constant HOMOLYSIS_SPAWN_SOURCE
            # (defined below; used at call time) -- both closure guards
            # (producer + runner _HOMOLYSIS_SPAWN_SOURCE) pin it exactly.
            daughter.spawn_metadata = {"source": HOMOLYSIS_SPAWN_SOURCE}
            # End-radical depropagation kernel (r74 SS2): the deck declares
            # k_depropagation on the PARENT's k_homolysis context (daughters
            # are never deck-declared); apply it to the spawned radical-end
            # pools where the kernel actually runs. The triplet is
            # deep-copied (post-hoc mutation of the parent's dict must not
            # reach the daughters), while the released-monomer routing
            # Species is held BY REFERENCE -- identity with the core
            # Species is load-bearing for the object-keyed spc_map
            # resolution in derive_daughter_pool_configs (same contract as
            # the QSSA _inherit_unzip_channel path).
            kdep = getattr(self, 'k_depropagation', None)
            if kdep is not None:
                daughter.k_depropagation = deepcopy(kdep)
                daughter.monomer_product_species = getattr(
                    self, 'monomer_product_species', None)
            daughters.append(daughter)
        return tuple(daughters)

    def _side_loss_feature_unit(self, channel_label: str,
                                gas_product: str,
                                site_selector: str) -> Molecule:
        """
        Build the X-loss feature repeat unit for one side_group_homolysis
        channel (FR1-K1): the parent monomer with ONE X atom removed and a
        radical electron left on its ex-neighbor (the mid-chain radical
        defect of homolytic C-X cleavage).

        DETERMINISTIC transformation (documented; the channel LABEL keys
        pool identity, the unit is the representative structure): the
        removed atom is the FIRST atom, in the unit's atom order, matched
        by the channel's REQUIRED structural site selector
        (rmgpy.solver.polymer.side_group_site_atom_indices; round-72 P1:
        the kinetics label alone must never pick the site -- on a
        mixed-site monomer a channel labeled 'aliphatic_C-Br' would
        otherwise remove the first X in atom order, an ARYL one, minting
        the wrong defect structure under that channel's rate). The
        removable-X predicate is unchanged: terminal (exactly one bond),
        single-bonded, closed-shell and not a stitch-labeled atom. No
        selector-matched atom is a hard error naming the channel (the
        validator enforces the same law with the full three-part
        structural check).

        Validated by its own STRICT assertion path
        (_assert_side_loss_unit: element census = monomer minus one X,
        exactly 3 radical electrons, radical on the ex-neighbor) plus the
        existing feature-unit machinery (_assert_feature_unit with the
        3-radical budget, mirroring the S2 H-loss units).
        """
        from rmgpy.solver.polymer import side_group_site_atom_indices
        gas = Molecule().from_smiles(gas_product)
        if len(gas.atoms) != 1 or gas.get_radical_count() != 1:
            raise ValueError(
                f"Polymer '{self.label}': side_group_homolysis channel "
                f"'{channel_label}' gas_product={gas_product!r} must be a "
                f"monoatomic mono-radical (e.g. '[Br]') in v1.")
        sym = gas.atoms[0].symbol
        unit = self.monomer.copy(deep=True)
        # Selector-directed selection, resolved ON THE COPY the removal
        # operates on (self-consistent indices whatever the copy's atom
        # ordering).
        idxs = side_group_site_atom_indices(unit, sym, site_selector)
        if not idxs:
            raise ValueError(
                f"Polymer '{self.label}': side_group_homolysis channel "
                f"'{channel_label}' (site_selector='{site_selector}') "
                f"needs a removable side-group {sym} atom (terminal, "
                f"single-bonded, closed-shell, unlabeled) on a matching "
                f"carbon environment in the repeat unit, but the monomer "
                f"has none. The X-loss feature pool cannot be built.")
        x_atom = unit.atoms[idxs[0]]
        neighbor = next(iter(x_atom.bonds.keys()))
        unit.remove_atom(x_atom)
        neighbor.increment_radical()
        unit.update()
        _assert_side_loss_unit(unit, sym, self.monomer, channel_label,
                               radical_atom=neighbor)
        self._assert_feature_unit(unit, allow_h_loss_radical=True)
        return unit

    def generate_side_loss_daughters(self) -> Tuple['Polymer', ...]:
        """
        Producer for the side-group homolysis conduit (FR1-K1, adjudicated
        round 70): build ONE X-loss feature pool per configured channel of
        this parent, for registration at MODEL SETUP (not lazily) --
        mirrors generate_end_radical_daughters.

        Each daughter '{label}_sidegrp_{sanitized channel label}':
        born-at-zero moments (spawned_empty pattern), the parent's monomer
        (monomer_mw_g_mol pinned by __init__ -- the chain transfers INTACT,
        no chain cut), an X-loss feature_monomer built by the documented
        deterministic transformation (_side_loss_feature_unit) with its own
        strict assertion path, a distinct _SideLoss fingerprint class keyed
        on the channel label, spawn provenance {'source':
        'side_group_homolysis', 'channel': <label>}, and the EXACT
        mass-defect contract chain_mass_defect_g_mol = M_X of the channel's
        gas_product (round-70 #1 P1 trap: the feature pool keeps the
        parent's monomer_mw, so its condensed mass is mu1*MW - mu0*M_X).

        v1 LIMITATION (adjudicated): daughters carry NO side_group_homolysis
        of their own -- the kernel acts on the PARENT pool only, so feature
        pools saturate as terminal X-loss sinks (no multi-loss cascade).
        """
        from rmgpy.solver.polymer import (side_group_daughter_pool_label,
                                          validate_side_group_homolysis)
        channels = getattr(self, 'side_group_homolysis', None)
        if not channels:
            return tuple()
        # Round-72 P1: structural validation against the parsed monomer --
        # selector required + matching >= 1 X atom, sites_per_unit checked
        # against the match count, no two channels on one atom set.
        channels = validate_side_group_homolysis(self.label, channels,
                                                 monomer=self.monomer)
        daughters = []
        for ch in channels:
            unit = self._side_loss_feature_unit(ch["label"],
                                                ch["gas_product"],
                                                ch["site_selector"])
            if self.Mn and self.Mw:
                dist_kwargs = dict(Mn=self.Mn, Mw=self.Mw, moments=None)
            else:
                dist_kwargs = dict(moments=[0.0, 0.0, 0.0])
            daughter = Polymer(
                label=side_group_daughter_pool_label(self.label,
                                                     ch["label"]),
                monomer=self.monomer,
                feature_monomer=unit,
                end_groups=[eg.copy(deep=True) for eg in self.end_groups],
                cutoff=self.cutoff,
                initial_mass=0.0,
                side_loss_channel=ch["label"],
                **dist_kwargs,
            )
            gas = Molecule().from_smiles(ch["gas_product"])
            daughter.chain_mass_defect_g_mol = (
                gas.get_molecular_weight() * 1000.0)
            daughter.parent_pool_label = self.label
            # Literal pinned as module constant
            # SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE (defined below, near
            # HOMOLYSIS_SPAWN_SOURCE) -- the schema-2.7 consumer (K2) will
            # pin it exactly, like the sibling kernel's provenance.
            daughter.spawn_metadata = {
                "source": SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE,
                "channel": ch["label"],
            }
            daughters.append(daughter)
        return tuple(daughters)

    def _capped_chain_species(self, dp: int) -> Optional[Species]:
        """
        The complete capped chain at the given degree of polymerization for
        this pool: end_groups[0]--(monomer)x``dp``--end_groups[1], as a
        Species with resonance structures (same recipe as _stitch_trimer) so
        kekule/aromatic representations compare equal. Single stitching
        recipe shared by the DP-1 fold-back gate (dp=1) and the explicit-DP
        oligomer generator (dp=cutoff). Returns None when the stitch is not
        constructible.
        """
        if dp < 1:
            raise ValueError(
                f"Polymer '{self.label}': capped chain DP must be >= 1, got {dp}.")
        mol = self.end_groups[0].copy(deep=True)
        for _ in range(int(dp)):
            mol = stitch_molecules_by_labeled_atoms(mol, self.monomer.copy(deep=True))
            if mol is None:
                return None
        mol = stitch_molecules_by_labeled_atoms(mol, self.end_groups[1].copy(deep=True))
        if mol is None:
            return None
        mol.update()
        mol.clear_labeled_atoms()
        mol.assign_atom_ids()
        spc = Species(molecule=[mol])
        spc.molecule = generate_resonance_structures(mol,
                                                     clar_structures=False,
                                                     keep_isomorphic=False,
                                                     filter_structures=True,
                                                     save_order=True)
        return spc

    def _dp1_capped_species(self) -> Optional[Species]:
        """
        The complete capped chain at DP = 1 for this pool:
        end_groups[0]--monomer--end_groups[1] (e.g. propane for a PP pool
        with H/H end caps, ethylbenzene for a PS pool with H/H end caps).
        Built by the shared _capped_chain_species recipe.
        Used by the DP-1 fold-back gate in _create_reacted_copy_logic.
        Returns None when the stitch is not constructible. Cached lazily per
        instance (copies rebuild it on first use).
        """
        cached = getattr(self, '_dp1_capped_cache', None)
        if cached is None:
            spc = self._capped_chain_species(1)
            # False sentinel: not constructible, do not retry.
            cached = spc if spc is not None else False
            self._dp1_capped_cache = cached
        return cached or None

    def generate_explicit_dp_species(self) -> Species:
        """
        Build the explicit boundary oligomer at DP == self.cutoff (xs) for
        the explicit-DP handshake (deck flag ``explicit_dp=True`` on the
        polymer() input block): the complete capped chain
        end_groups[0]--(monomer)x``cutoff``--end_groups[1], labeled
        ``{label}_dp{cutoff}`` (same naming family as the ``_mu{k}`` moment
        dummies).

        Hard errors (never silent -- a missing species would leave the
        handshake structurally inert):

        * ``cutoff < 2``: DP-1 capped chains are forced GAS by the DP-1
          fold-back gate (_create_reacted_copy_logic refuses products
          isomorphic to the DP-1 chain), so an explicit DP=1 species would
          collide with that gate. The deck validator already refuses
          cutoff < 2; this is the generator's own defensive gate against a
          future override mapping DP 1.
        * stitch failure: the capped chain is not constructible from
          monomer + end_groups.
        """
        dp = int(self.cutoff)
        if dp < 2:
            raise ValueError(
                f"Polymer '{self.label}': explicit_dp=True requires cutoff >= 2, "
                f"got cutoff={dp}. DP-1 capped chains are forced GAS by the DP-1 "
                f"fold-back gate (products at DP <= 1 stay gas; "
                f"_create_reacted_copy_logic), so an explicit DP=1 species would "
                f"collide with that gate. Raise the pool's cutoff to >= 2 or "
                f"disable explicit_dp.")
        spc = self._capped_chain_species(dp)
        if spc is None:
            raise ValueError(
                f"Polymer '{self.label}': explicit_dp=True could not construct "
                f"the capped DP={dp} oligomer from monomer + end_groups (stitch "
                f"failed). Check the monomer's [*:1]/[*:2] labels and the "
                f"end_groups, or disable explicit_dp.")
        spc.label = f"{self.label}_dp{dp}"
        return spc

    def create_reacted_copy(self, reacted_proxy: Molecule,
                            h_loss_feature: bool = False) -> Optional['Polymer']:
        """
        Wrapper that ensures any generated polymer fragment is sanitized
        (labels stripped, proxy tagged) before returning to the RMG engine.

        Args:
            reacted_proxy: A single product fragment from the reaction.
            h_loss_feature: HANDSHAKE CONTEXT flag (radical-feature producer
                path, stage S1a of the feature-pool conduit arc). The CALLER
                asserts, from reaction-level knowledge, that the reaction
                removed exactly ONE H atom from the polymer reactant's proxy
                (the H-abstraction shape). Only then may the route-first
                producer path materialize a ``{label}_mod`` RADICAL FEATURE
                daughter whose feature_monomer carries the mid-chain radical
                (see :meth:`_create_h_loss_feature_copy`). The flag is never
                inferred from the product's structure alone; the live
                handshake wiring (stage S2) computes the per-product verdict
                where reactants AND products are both visible
                (:func:`compute_h_loss_feature_verdicts`, threaded through
                ``_handshake_structures`` from ``make_new_reaction``).

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
            # method-of-moments model abstracts chain-end activation into the
            # pool's sanctioned depropagation representation -- either the lumped
            # k_unzip rate (dmu1/dt = -k_unzip*mu0, scaling with the chain-end
            # count mu0) or the radical_qssa_unzip channel (initiation/
            # depropagation/termination QSSA; the two are mutually exclusive per
            # pool, enforced at config validation) -- so we fold the product back
            # into the parent pool with moments and mass preserved rather than
            # spawning a distinct activated-chain population (which would
            # double-count whichever unzip representation the pool carries, and
            # cannot be represented anyway, since an activated end-cap is
            # di-radical: stitch site + activation).
            #
            # Without this, _create_reacted_copy_logic's raw wing-matching diverges
            # from classify_structure's heavy-view matcher, mis-routes the END_MOD
            # product into a scission branch where it fails the mono-radical
            # end-group assertion, returns None, and the product leaks to the gas
            # phase as a spurious small molecule (a mass-balance leak).
            new_poly = self.copy(deep=True)
        else:
            new_poly = self._create_reacted_copy_logic(
                reacted_proxy, h_loss_feature=h_loss_feature)
        if new_poly is None:
            return None
        # Stamp the classification verdict so the polymer handshake can flag
        # END_MOD reactions for chain-end (mu0) scaling in the solver. Read by
        # is_end_group_reaction(products); a transient generation-time marker.
        new_poly._reacted_class = klass
        # Daughter-pool channel inheritance (radical_qssa_unzip milestone 5),
        # applied at the SINGLE exit point and GATED on same repeat chemistry
        # (round-25 P2-1): scission tail/head and END_MOD/isomorphic copies
        # share the parent's monomer chemistry and inherit the QSSA unzip
        # channel (deep-copied) + monomer routing -- without this a PS
        # scission cascade freezes. A _mod product whose feature_monomer
        # CHANGED does NOT share that chemistry: it stays channel-free and
        # the withheld inheritance is disclosed once per pool.
        _inherit_unzip_channel(new_poly, self)
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

    def _create_reacted_copy_logic(self, reacted_proxy: Molecule,
                                   h_loss_feature: bool = False,
                                   ) -> Optional['Polymer']:
        """
        Creates a new Polymer species from a reacted proxy fragment.
        Handles both modification (intact chain) and scission (broken chain)
        via the wing-matching logic. ROUTE-FIRST ordering (stage S2,
        adjudicated after the S1a fallback ordering proved a design P1 in the
        live PP run): when the caller threaded the H-loss handshake context
        (``h_loss_feature``), the radical-feature producer path
        (:meth:`_create_h_loss_feature_copy`, stage S1a) runs BEFORE the wing
        logic, so the scission branches never see a DP-preserving H-loss
        daughter (live defect: the center-tertiary PP daughter spawned
        ``polypropylene_scission_tail`` -- a same-length chain misbooked as a
        half-length population with a malformed +1-cation C15H31 proxy). The
        wing logic independently enforces the scission invariant (a scission
        daughter must be strictly shorter than the parent proxy), so a lying
        or absent flag cannot re-open the scission-spawn hole.

        Args:
            reacted_proxy (Molecule): A single product fragment from the reaction.
            h_loss_feature (bool): Handshake context -- the reaction removed
                exactly ONE H from the unit (see :meth:`create_reacted_copy`).

        Returns:
            Optional['Polymer']: A new Polymer species with updated feature or end-groups.
                                 Returns None if for an invalid/unsupported polymer fragment.
        """
        if h_loss_feature:
            new_poly = self._create_h_loss_feature_copy(reacted_proxy)
            if new_poly is not None:
                return new_poly
        return self._create_reacted_copy_wing_logic(reacted_proxy)

    def _create_reacted_copy_wing_logic(self, reacted_proxy: Molecule) -> Optional['Polymer']:
        """
        Wing-matching core of :meth:`_create_reacted_copy_logic`: classifies
        the fragment by locating intact head/tail wings and builds the
        modification (intact chain) or scission (broken chain) daughter.
        """
        product = reacted_proxy.copy(deep=True)
        product.clear_labeled_atoms()
        product.update()
        if self.baseline_proxy.is_isomorphic(product):
            return self.copy(deep=True)

        # DP-1 fold-back gate (ratified 2026-07-04, PP run-2 gate/conduit
        # diagnosis Phenomenon 2): a product that is a COMPLETE capped chain
        # at DP = 1 -- end_groups[0]-monomer-end_groups[1] exactly, e.g.
        # propane == H-(C3H6)-H under end_groups=['[H]','[H]'] -- is the
        # monomer-hydride/solvent-class VOLATILE co-product of abstraction/
        # disproportionation against the proxy, NOT a chain population. The
        # wing matcher below would otherwise accept its single wing, build a
        # "{label}_scission_tail" whose reconstructed proxy is isomorphic to
        # the PARENT proxy, and species dedup would fold the volatile back
        # into the parent pool -- leaving the reaction's melt reference state
        # unpaired (measured U = 11.63 decades vs 0.33 with the volatile kept
        # gas; the run-2 tripwire refusal). Ratified rule: products at
        # DP <= 1 STAY GAS; there is no automatic DP-1 fold-back (a deck
        # wanting DP-1 condensed bookkeeping must configure it explicitly).
        # DP-0 (end_groups[0]-end_groups[1], e.g. H2) needs no gate: it
        # cannot contain a wing (end + full monomer) and already falls
        # through to None structurally.
        dp1 = self._dp1_capped_species()
        if dp1 is not None and dp1.is_isomorphic(product):
            return None

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
            return self._born_at_zero_mod_daughter(new_feature_graph,
                                                   source="feature_mod")

        if head_atoms or tail_atoms:
            # Scission invariant (stage S2, feature-pool conduit arc): a
            # scission head/tail daughter is a PIECE of the cut chain, so it
            # must be strictly shorter (fewer heavy atoms) than the parent
            # proxy. A same-heavy-skeleton H-loss radical daughter reaches
            # here with the full proxy heavy count (its internal radical
            # coincides with the *2 cut atom, so _ensure_open_site adds
            # nothing) and used to spawn a DP-preserving chain misbooked as a
            # brand-new half-length population (Mn/2, Mw/2) with a malformed
            # +1-cation proxy -- the live PP-run species-25 defect. Refuse
            # regardless of the h_loss_feature flag: routing is the producer
            # path's job (route-first in _create_reacted_copy_logic), and a
            # non-routed same-length daughter must fall to the refuse stamp,
            # never to a scission spawn.
            proxy_heavy = sum(
                1 for a in self.baseline_proxy.molecule[0].atoms
                if not a.is_hydrogen())
            product_heavy = sum(
                1 for a in product.atoms if not a.is_hydrogen())
            if product_heavy >= proxy_heavy:
                return None

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

    def _h_loss_feature_units(self) -> List[Molecule]:
        """
        Enumerate the candidate H-loss feature units of this pool's monomer:
        for every H-bearing heavy atom of the (labeled) monomer, the unit
        obtained by removing ONE H there and adding one radical electron --
        i.e. the repeat unit as it looks after a mid-chain H-abstraction,
        carrying the two stitch radicals on '*1'/'*2' plus exactly one extra
        internal radical. Deterministic monomer-atom order (symmetric-
        equivalent positions yield isomorphic units and the first match wins
        downstream). For PP ('[CH2][CH](C)') this is the documented-v1 set of
        ~3 H-environments: backbone CH2, backbone CH, pendant CH3.
        """
        units = []
        n = len(self.monomer.atoms)
        for i in range(n):
            unit = self.monomer.copy(deep=True)
            if len(unit.atoms) != n:  # defensive: copy must preserve order
                return units
            target = unit.atoms[i]
            if target.is_hydrogen():
                continue
            h_neighbor = next((a for a in target.bonds if a.is_hydrogen()), None)
            if h_neighbor is None:
                continue
            unit.remove_atom(h_neighbor)
            target.increment_radical()
            unit.update(sort_atoms=False)
            units.append(unit)
        return units

    def _h_loss_positional_species(self, unit: Molecule, position: int) -> Optional[Species]:
        """
        The complete capped 3-unit chain with ``unit`` substituted at
        ``position`` (0 = head-side flank, 1 = center, 2 = tail-side flank)
        and the baseline monomer everywhere else -- the positional rendering
        of a single mid-chain H-abstraction on the trimer proxy. Same
        stitching + resonance recipe as :meth:`_capped_chain_species`.
        Returns None when the stitch is not constructible.
        """
        mol = self.end_groups[0].copy(deep=True)
        for i in range(3):
            piece = unit if i == position else self.monomer
            mol = stitch_molecules_by_labeled_atoms(mol, piece.copy(deep=True))
            if mol is None:
                return None
        mol = stitch_molecules_by_labeled_atoms(mol, self.end_groups[1].copy(deep=True))
        if mol is None:
            return None
        mol.update()
        mol.clear_labeled_atoms()
        mol.assign_atom_ids()
        spc = Species(molecule=[mol])
        spc.molecule = generate_resonance_structures(mol,
                                                     clar_structures=False,
                                                     keep_isomorphic=False,
                                                     filter_structures=True,
                                                     save_order=True)
        return spc

    @staticmethod
    def _h_loss_unit_rendered_radical_envs(unit: Molecule) -> List[Tuple[int, int]]:
        """
        (n_H, n_heavy) environments of the unit's EXTRA radical site(s) as
        rendered mid-chain: each labeled stitch atom consumes one radical and
        gains one heavy (stitch-bond) neighbor when the unit sits inside a
        chain. Used only as the ambiguity tie-break in
        :meth:`_create_h_loss_feature_copy`.
        """
        envs = []
        for a in unit.atoms:
            if a.radical_electrons <= 0:
                continue
            rendered_rad = a.radical_electrons - (1 if a.label else 0)
            if rendered_rad <= 0:
                continue
            n_h = sum(1 for nb in a.bonds if nb.is_hydrogen())
            n_heavy = (sum(1 for nb in a.bonds if not nb.is_hydrogen())
                       + (1 if a.label else 0))
            envs.append((n_h, n_heavy))
        return envs

    def _create_h_loss_feature_copy(self, product: Molecule) -> Optional['Polymer']:
        """
        Radical-feature producer path (stage S1a, feature-pool conduit arc,
        adversarially ratified). ONLY reached when the caller threaded the
        H-loss handshake context (``h_loss_feature=True`` -- the reaction
        removed exactly ONE H from the polymer reactant's proxy) AND the
        wing logic refused the product; never triggered from the product's
        structure alone.

        Structural cross-check (so a lying context cannot smuggle garbage):
        the product must be isomorphic to one of the pool's single-H-loss
        POSITIONAL VARIANTS -- the capped trimer with exactly one H-loss
        feature unit at the center or either flank. Di-radical (2-H-loss)
        products, crosslinks, wrong skeletons and any other shape match no
        variant and stay refused exactly as today.

        Positional twins collapse by construction: a center hit (classified
        FEATURE) and its flank twins (classified DISCARD, the positional
        renderings of the SAME abstraction) all map to the SAME feature unit
        graph, hence the same ``{label}_mod`` daughter fingerprint and ONE
        pool after :meth:`CoreEdgeReactionModel._register_polymer` dedup.
        Distinct H-environments (PP documented v1: backbone CH2, backbone
        CH, pendant CH3 -- ~3 pools) map to distinct units and distinct
        pools. A chain-END twin can be positionally ambiguous (e.g. the PP
        C1 methyl daughter renders identically from the capped backbone-CH2
        unit and the pendant-methyl unit); the tie-break prefers the unit
        whose MID-CHAIN radical environment (n_H, n_heavy) matches the
        product's radical atom, falling back to deterministic monomer-atom
        order.
        """
        matched_units = []
        for unit in self._h_loss_feature_units():
            for position in (1, 0, 2):  # center first: canonical representative
                spc = self._h_loss_positional_species(unit, position)
                if spc is not None and spc.is_isomorphic(product):
                    matched_units.append(unit)
                    break
        if not matched_units:
            return None
        unit = matched_units[0]
        if len(matched_units) > 1:
            radical_atoms = [a for a in product.atoms if a.radical_electrons > 0]
            if len(radical_atoms) == 1:
                pa = radical_atoms[0]
                product_env = (sum(1 for nb in pa.bonds if nb.is_hydrogen()),
                               sum(1 for nb in pa.bonds if not nb.is_hydrogen()))
                exact = [u for u in matched_units
                         if product_env in self._h_loss_unit_rendered_radical_envs(u)]
                if exact:
                    unit = exact[0]
        feature = unit.copy(deep=True)
        try:
            self._assert_feature_unit(feature, allow_h_loss_radical=True)
        except ValueError:
            return None
        return self._born_at_zero_mod_daughter(feature,
                                               source="radical_feature_h_loss")

    def _born_at_zero_mod_daughter(self, feature_monomer: Molecule,
                                   source: str) -> 'Polymer':
        """
        Construct the ``{label}_mod`` feature-pool daughter BORN AT ZERO
        (stage S1b, the ratified P1 mass-duplicator fix): a just-spawned
        daughter genuinely contains nothing, so it seeds ``initial_mass=0``
        -> moments [0, 0, 0], exactly unifying with the scission-tail
        convention and ``drain_spawn_intents``' honest-empty daughters --
        instead of duplicating the PARENT's moments (the parent keeps its
        mass; a verbatim copy re-declared it). H abstraction / feature
        modification does not cut the chain, so the parent's Mn/Mw ride
        along unchanged as lineage/DP metadata (with zero mass they derive
        zero interim moments). Spawn provenance markers
        (``parent_pool_label`` + ``spawn_metadata``) make the sidecar's
        secondary spawned-pool signal classify the pool
        ``moments_provenance: spawned_empty`` even in legacy default-label
        serializer calls. ``monomer_mw_g_mol`` is set by ``Polymer.__init__``
        from the (shared) monomer -- pinned by test: the TA consumer sizes
        sample mass with it.
        """
        daughter = Polymer(
            label=f"{self.label}_mod",
            monomer=self.monomer,
            feature_monomer=feature_monomer,
            end_groups=[eg.copy(deep=True) for eg in self.end_groups],
            cutoff=self.cutoff,
            Mn=self.Mn,
            Mw=self.Mw,
            moments=None,
            initial_mass=0.0,
        )
        daughter.parent_pool_label = self.label
        daughter.spawn_metadata = {"source": source}
        return daughter

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
    def _assert_feature_unit(mol: Molecule, allow_h_loss_radical: bool = False):
        """
        Validate that an extracted 'feature monomer' is stitchable as a repeat unit.
        Contract we enforce:
        - Exactly one '*1' and exactly one '*2' label.
        - No other labels (including '1','2') survive.
        - Radical budget (deliberately NOT a blanket >=2 relaxation, which
          would accept garbage -- crosslinks, multi-radicals, bad resonance
          forms):
          * default: exactly 2 radical electrons total (one per stitch end;
            stitching consumes one at each end) -- unchanged.
          * ``allow_h_loss_radical=True`` (radical-feature producer path
            ONLY; the caller must hold explicit handshake context that the
            reaction removed exactly ONE H from this unit -- never inferred
            from the unit's structure alone): exactly 3 radical electrons
            total -- the two stitch radicals plus exactly ONE extra internal
            radical from the abstracted H. A 2-radical unit (no H actually
            lost) or a 4+-radical unit (di-radical / 2-extra-radical
            garbage) refuses under the flag too.
        - The labeled atoms must each have >=1 radical electron (open sites).

        Args:
            mol (Molecule): The molecule to validate.
            allow_h_loss_radical (bool): See radical budget above.
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
        expected_rad_count = 3 if allow_h_loss_radical else 2
        if rad_count != expected_rad_count:
            if allow_h_loss_radical:
                raise ValueError(
                    f"H-loss radical feature unit must have exactly 3 radical electrons "
                    f"total (two stitch radicals + exactly ONE extra internal radical "
                    f"from the abstracted H). Got radical_count={rad_count}.")
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


def is_qssa_eliminating_radical(mol):
    """True if a radical-bearing molecule is QSSA-ELIMINATING (fast beta-scission,
    not resonance-stabilized) -> conduit-representable; False if ACCUMULATING
    (resonance-stabilized, e.g. allylic) -> QSSA-invalid. Probe F: allylic k_beta
    is 3-4 orders slower; resonance count is the in-codebase signal. Runs once per
    FEATURE reaction at generation, never in the RHS."""
    from rmgpy.molecule.resonance import generate_resonance_structures
    if not mol.is_radical():
        return True  # not a radical: not the accumulating case
    return len(generate_resonance_structures(mol.copy(deep=True))) <= 1


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
    VOLATILE_EJECTION = 6  # depolymerization: polymer A -> discrete volatile(s)
                           # + cross-pool polymer B; mass leaves the chain as
                           # the ejected volatile (a = fractional source-monomer
                           # equivalents ejected). Mass-losing, unlike MIGRATION.


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


_refused_census_warned = set()

_pool_cap_warned = set()


def _warn_pool_cap_exhausted(cand_label: str, n_pools: int, max_pools: int) -> None:
    """Warn once per (registry size, cap) that the max_pools cap declined a
    gate-cleared spawn -- and NAME radical-feature pool exhaustion: the
    registry the cap counts includes radical-feature '_mod' daughter pools
    (one per distinct H-abstraction environment of a pool's monomer, e.g.
    ~3 for polypropylene), so feature pools alone can exhaust the cap and
    silently suppress genuine novel-motif spawns."""
    key = (n_pools, max_pools)
    if key in _pool_cap_warned:
        return
    _pool_cap_warned.add(key)
    logging.warning(
        "Polymer pool cap reached (max_pools=%d, registry holds %d pools): "
        "declining to spawn a new pool for candidate %s despite a cleared "
        "mass-flux gate. NOTE: radical-feature '_mod' daughter pools (one "
        "per distinct H-abstraction environment, e.g. ~3 for polypropylene) "
        "count toward this cap -- if feature-pool exhaustion is suppressing "
        "spawns, raise max_pools.",
        max_pools, n_pools, cand_label)


def _warn_once_refused(entry: dict) -> None:
    """Log each distinct refused-FEATURE-radical reaction once (correct-but-loud)."""
    key = (entry["reaction"], entry["reason"])
    if key not in _refused_census_warned:
        _refused_census_warned.add(key)
        logging.warning(
            "FEATURE-RADICAL REFUSED CENSUS: %s -- %s radical refused (%s); "
            "no flux applied (stamp-but-keep). Deferred to item 20's conduit.",
            entry["reaction"], entry["radical_class"], entry["reason"])


def _reaction_census_label(rxn) -> str:
    """Stable label for a refused reaction's census entry."""
    return str(rxn)


_double_count_warned = set()


def _warn_once_double_count(entry: dict) -> None:
    """Log each distinct double-count pool once (correct-but-loud). Severity is
    census-loud only (no refuse) in item 18; the refuse-vs-warn cliff is
    calibrated by item 19's magnitude probe."""
    key = (entry["pool"],)
    if key not in _double_count_warned:
        _double_count_warned.add(key)
        logging.warning(
            "DOUBLE-COUNT TRIPWIRE: pool '%s' has BOTH explicit beta-scission/chip "
            "reactions AND nonzero phenomenological k_scission=%.3e / k_unzip=%.3e -- "
            "chain degradation may be double-counted. (Severity calibrated by item 19; "
            "census-only in item 18.)",
            entry["pool"], entry["k_scission"], entry["k_unzip"])


_qssa_double_count_warned = set()


def warn_once_qssa_double_count(entry: dict) -> None:
    """Log each distinct QSSA double-count census entry once
    (correct-but-loud, mirroring ``_warn_once_double_count``). QSSA
    initiation is random backbone homolysis, so it can overlap two other
    backbone-cutting channels on the same pool, distinguished by
    ``entry["overlap"]``:

    - ``"generated_scission_ve"``: surviving generated-chemistry SCISSION /
      VOLATILE_EJECTION reactions sourced from the pool.
    - ``"k_scission"``: the pool's own phenomenological k_scission -- the
      most direct initiation double-count (both cut random backbone bonds).

    Census-only, NEVER refuse: the overlapping channel may have been
    parameterized for DIFFERENT/disjoint physics and only the user can know
    (hard exclusion would be wrong); only the k_unzip co-presence is a hard
    error (M1 mutual exclusion). Warn-once keyed per (pool, overlap kind)."""
    overlap = entry.get("overlap", "generated_scission_ve")
    key = (entry["pool"], overlap)
    if key in _qssa_double_count_warned:
        return
    _qssa_double_count_warned.add(key)
    if overlap == "k_scission":
        logging.warning(
            "QSSA DOUBLE-COUNT TRIPWIRE: pool '%s' has BOTH a radical_qssa_unzip "
            "channel (QSSA initiation = random backbone homolysis) AND a nonzero "
            "pool-level k_scission=%.3e [s^-1] -- both cut random backbone bonds, "
            "the most direct initiation double-count. Ensure the two are "
            "parameterized to cover DISJOINT physics. Census-only (warn, never "
            "refuse).",
            entry["pool"], entry["k_scission"])
    else:
        logging.warning(
            "QSSA DOUBLE-COUNT TRIPWIRE: pool '%s' has BOTH a radical_qssa_unzip "
            "channel (backbone-homolysis initiation) AND surviving explicit "
            "scission/volatile-ejection reactions sourced from it -- chain "
            "degradation may be double-counted. Census-only (warn, never refuse: "
            "the generated chemistry may cover different bonds).",
            entry["pool"])


# Back-compat alias: the helper was born underscore-private and the solver
# (and possibly external callers) imported it by that name. The public name
# above is canonical; keep the private binding pointing at the same function
# (the warn-once state set `_qssa_double_count_warned` keeps its name).
_warn_once_qssa_double_count = warn_once_qssa_double_count


_homolysis_supersession_warned = set()


def warn_once_homolysis_supersession(entry: dict) -> None:
    """Log each k_homolysis supersession census entry once (Stage 1,
    adjudicated round 66): refused explicit homolysis/association rows
    co-existing with a live k_homolysis kernel are FINE -- they carry zero
    flux (``reaction_refused``) -- but the pairing must be stated explicitly
    in the run log: the refused conduit-deferred rows are SUPERSEDED by the
    pool's kernel, which now carries the initiation flux they were refused
    for. Census-only, NEVER refuse. Warn-once keyed per pool."""
    key = entry["pool"]
    if key in _homolysis_supersession_warned:
        return
    _homolysis_supersession_warned.add(key)
    rows = entry.get("superseded_rows", [])
    logging.warning(
        "HOMOLYSIS KERNEL CENSUS: pool '%s' has a live k_homolysis initiation "
        "kernel; %d refused conduit-deferred row(s) are superseded by the "
        "kernel (they remain refused and carry zero flux): %s",
        entry["pool"], len(rows), "; ".join(str(r) for r in rows))


_side_group_supersession_warned = set()


def warn_once_side_group_supersession(entry: dict) -> None:
    """Log each side_group_homolysis supersession census entry once
    (FR1-K1; same contract as warn_once_homolysis_supersession): refused
    explicit gas-radical<->condensed rows co-existing with a live
    side-group kernel REMAIN refused and carry zero flux, and the pairing
    is stated explicitly in the run log. Census-only, NEVER refuse.
    Warn-once keyed per pool."""
    key = entry["pool"]
    if key in _side_group_supersession_warned:
        return
    _side_group_supersession_warned.add(key)
    rows = entry.get("superseded_rows", [])
    logging.warning(
        "SIDE-GROUP HOMOLYSIS KERNEL CENSUS: pool '%s' has a live "
        "side_group_homolysis initiation kernel; %d refused "
        "conduit-deferred row(s) are superseded by the kernel (they remain "
        "refused and carry zero flux): %s",
        entry["pool"], len(rows), "; ".join(str(r) for r in rows))


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


def _nonpolymer_product_mw_g_mol(product) -> Optional[float]:
    """Molecular weight (g/mol) of a NON-Polymer reaction product, robust to the
    object type the LIVE pipeline actually passes at stamp time.

    A product may be a ``Molecule`` (has ``get_molecular_weight`` -> kg/mol) OR a
    bare pre-thermo ``Species`` (NO ``get_molecular_weight`` method, and often an
    empty label, but a populated ``.molecule`` structure). ``Polymer`` (a Species
    subclass) is a pool chain, not a discrete volatile -> returns None. Returns
    None when no structure is available (cannot weigh -> caller treats it as "no
    volatile", i.e. a mass-conserving relabel).

    Empirically load-bearing: at ``stamp_polymer_flux_archetype`` time the
    depolymerization volatiles arrive as bare ``Species``, which lack
    ``get_molecular_weight`` -- a ``hasattr(p, 'get_molecular_weight')`` gate
    silently misses them and mis-stamps VOLATILE_EJECTION as MIGRATION
    (mass-conserving), flattening the TGA. Reach through to ``.molecule[0]``.
    """
    if isinstance(product, Polymer):
        return None
    getter = getattr(product, 'get_molecular_weight', None)
    if callable(getter):
        return getter() * 1000.0
    mols = getattr(product, 'molecule', None)
    if mols:
        return mols[0].get_molecular_weight() * 1000.0
    return None


# Signed-VOLATILE_EJECTION contract (Codex round-13). Routing gate excludes only
# exact net-zero fold-backs; the DP-vs-mass THRESHOLD is deferred, so we never
# hard-gate on a large threshold -- instead a diagnostic warn-once census fires
# when |a| is atom-transfer-scale (see _warn_atom_transfer_ve).
_VE_NET_MASS_EPS_G = 1e-6      # g/mol; |net| below this is a mass-conserving relabel
_VE_ATOM_TRANSFER_UNITS = 0.5  # source-monomer-equivalents; below => census-only warn

# Number of REPEAT UNITS spanned by the stitched pool proxy
# (Polymer._stitch_trimer: head-cap + [baseline, center, baseline] + tail-cap
# -- the construction asserts against this constant, so a future proxy-size
# change must move it and everything derived from it moves too).
PROXY_STITCH_REPEAT_UNITS = 3

# r82 impostor-row refusal (FR1 run-2): a discrete (non-Polymer) participant at
# or above this many source-monomer-equivalents on EVERY axis (mass AND
# heavy-atom count, see _discrete_is_polymer_sized -- r85 P2(a): BOTH axes
# must be COMPUTABLE and at/above threshold to refuse; an uncomputable axis
# makes the case undecidable and it is warned, never refused) is
# "polymer-sized" -- proxy-scale, not volatile-scale. Proxy-derivation is NOT
# traceable at the stamp/restamp seam (r82 probe 3), so the conservative
# monomer-scale threshold form is adjudicated instead, DERIVED from the
# actual proxy construction (r85 P2(b)): the stitched proxy spans
# PROXY_STITCH_REPEAT_UNITS repeat units, and the run-2 impostors (proxy
# minus Br2/HBr) sit at 2.76-3.00 monomer-equivalents on both axes, while
# the largest adjudicated-LIVE volatile shape (DP-2 dimer volatiles, e.g.
# hexene off polypropylene -- pinned negative control) sits at exactly 2.0.
# proxy_repeat_units - 0.5 = 2.5 splits the two conservatively: only clearly
# proxy-scale discretes are refused; anything smaller stays live for the r71
# unstamped-proxy hard-fail to catch LOUDLY rather than being silently
# refused.
_IMPOSTOR_DISCRETE_MONOMER_UNITS = PROXY_STITCH_REPEAT_UNITS - 0.5


def _net_nonpolymer_mass_g_mol(reactants, products) -> float:
    """SIGNED net non-polymer mass (g/mol) that leaves the chain in a reaction:

        (Σ non-polymer PRODUCT MW) − (Σ non-polymer REACTANT MW)

    weighed per participant by :func:`_nonpolymer_product_mw_g_mol` (``None`` ->
    skipped, i.e. Polymer pool chains and structureless species do not count).

    - net > 0 : chain LOSES mass (ejection / depropagation).
    - net < 0 : chain GAINS mass (monomer / radical addition, reverse ejection).
    - |net| ~ 0 : mass-conserving relabel / fold-back (NOT VOLATILE_EJECTION).

    Netting the REACTANT non-polymers is REQUIRED so a bimolecular H_Abstraction
    ``chain + R• -> chain + RH`` debits only the ~1 H actually shed (net = MW(RH)
    − MW(R•) ≈ MW(H)), not the full RH mass -- otherwise every radical abstraction
    would fabricate a whole co-reactant's worth of chain mass loss.
    """
    products_g = 0.0
    for p in products:
        mw = _nonpolymer_product_mw_g_mol(p)
        if mw is not None:
            products_g += mw
    reactants_g = 0.0
    for r in reactants:
        mw = _nonpolymer_product_mw_g_mol(r)
        if mw is not None:
            reactants_g += mw
    return products_g - reactants_g


def _warn_atom_transfer_ve(signature, a) -> None:
    """Diagnostic-only warn-once census (signed-VE spec): a VOLATILE_EJECTION whose
    ``|a|`` is below one atom-transfer unit (``_VE_ATOM_TRANSFER_UNITS`` = 0.5
    source-monomer-equivalents) is an atom-scale mass transfer (e.g. the single H
    shed in an H_Abstraction), not a monomer-scale unzip. μ1 is then treated as a
    MASS BUCKET, not a degree-of-polymerization count, for this reaction. Routing
    is UNCHANGED (the DP-vs-mass threshold is deferred); this census only surfaces
    the a-distribution on real decks so the threshold can be set empirically.
    Reuses the ``_flux_archetype_warned`` warn-once set, keyed to avoid spam."""
    key = ("atom_transfer_ve", signature)
    if key not in _flux_archetype_warned:
        _flux_archetype_warned.add(key)
        logging.warning(
            "Polymer VOLATILE_EJECTION atom-transfer census: |a|=%.4g < %.2f "
            "source-monomer-equivalents for %s -- mu1 treated as a mass bucket "
            "(not DP) for this reaction; routing unchanged (threshold deferred).",
            abs(a), _VE_ATOM_TRANSFER_UNITS, signature)


def classify_reaction_flux_archetype(reactants, products,
                                     trust_reacted_class=True) -> PolymerFluxArchetype:
    """
    Classify a reaction's pool moment-flux archetype from its (handshaked)
    reactant and product lists. Product Polymers carry ``_reacted_class``
    stamped by :meth:`Polymer.create_reacted_copy` (or re-stamped by chip
    product surgery, spec 2026-06-10 §4.2); pool identity is the Polymer
    label (the same key the solver's ``initialize_model`` uses to map
    species to pools).

    ``trust_reacted_class=False`` (r92 flip-restamp): classify STAMP-BLIND,
    from pool labels + net non-polymer mass only. Used when re-classifying a
    kinetics-FLIPPED reaction (:func:`restamp_flipped_polymer_archetype`),
    where the participant lists hold REGISTRY objects: a registry daughter
    Polymer persistently carries the ``_reacted_class`` stamped by the
    reaction that CREATED it (``_register_polymer`` registers the reacted
    copy itself), which describes a DIFFERENT reaction -- trusting it would
    spuriously route the flipped row SCISSION_FRAGMENT/DISCRETE_CHIP. The
    stamp-gated branches are simply unreachable then; shapes that would
    need them classify UNRESOLVED (and the caller refuses).
    """
    reactant_pools = {r.label for r in reactants if isinstance(r, Polymer)}
    product_polymers = [p for p in products if isinstance(p, Polymer)]
    if not reactant_pools and not product_polymers:
        return PolymerFluxArchetype.NONE

    if trust_reacted_class and any(
            getattr(p, '_reacted_class', None) == PolymerClass.CHIP
            for p in product_polymers):
        # Chip product surgery (surge_chip_products) already rewrote this
        # product list to [discrete chip, CHIP-stamped fold-back]. This check
        # MUST precede the SCISSION branch: after the (b)-surgery there is no
        # END_MOD member left, so is_end_group_reaction(products) would
        # recompute False and misroute. The solver reads the STORED Reaction
        # flag; nothing downstream may recompute it from product stamps.
        return PolymerFluxArchetype.DISCRETE_CHIP

    if trust_reacted_class and any(
            getattr(p, '_reacted_class', None) == PolymerClass.SCISSION
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
        # All polymer products fold back into a reactant pool. A genuine
        # net-zero fold-back is SAME_POOL -- BUT if a discrete mass-carrying
        # non-polymer product ALSO leaves (radical unzip / depropagation:
        # `pool A -> monomer + pool A`, the chain sheds a monomer and the
        # shorter chain stays in the SAME pool), the chain LOSES that mass.
        # That is VOLATILE_EJECTION with src == dst (same pool), NOT SAME_POOL:
        # SAME_POOL would apply zero moment loss while the gas monomer is still
        # produced, re-creating the flat-TGA / mass-fabrication bug under a new
        # archetype (Codex round-12). The solver's VE dispatch handles src==dst
        # (the full-bundle-out + a-shifted-bundle-in on one pool nets to the
        # chip drain mu0:0, mu1:-a, mu2:-(2a*E[n]-a^2)).
        # Signed-net gate (Codex round-13): route VE only when the NET non-polymer
        # mass moved is non-zero (abs(net) > eps). This nets the reactant side too,
        # so a same-pool bimolecular step (chain + R• -> chain + RH) that only
        # shifts an atom -- and a monomer ADDITION (chain + monomer -> chain, net
        # < 0, chain grows) -- both route VE, while an exact fold-back nets 0 and
        # stays SAME_POOL.
        if (len(reactant_pools) == 1 and len(product_polymers) == 1
                and abs(_net_nonpolymer_mass_g_mol(reactants, products))
                > _VE_NET_MASS_EPS_G):
            return PolymerFluxArchetype.VOLATILE_EJECTION
        return PolymerFluxArchetype.SAME_POOL
    if (len(product_polymers) == 1 and len(cross_pool) == 1
            and len(reactant_pools) == 1):
        # Same shape as MIGRATION (1 reactant pool, 1 cross-pool polymer
        # product), BUT if a discrete non-polymer product ALSO leaves (a
        # gas-phase volatile carrying mass), the chain LOSES that mass -- this
        # is depolymerization, not a mass-conserving whole-chain relabel. Route
        # VOLATILE_EJECTION so the ejected volatile's mass is debited from the
        # chain. A non-polymer product that carries structure (Molecule OR bare
        # pre-thermo Species) is a mass-carrying volatile; _nonpolymer_product_mw
        # is robust to both (a plain hasattr(get_molecular_weight) gate misses
        # the Species the live pipeline actually passes -- see that helper).
        # Signed-net gate (Codex round-13): use abs(net) > eps (netting reactant
        # non-polymers too), so a bimolecular cross-pool step debits only the net
        # mass moved, not a whole co-reactant; an exact whole-chain relabel nets 0
        # and stays MIGRATION.
        if (abs(_net_nonpolymer_mass_g_mol(reactants, products))
                > _VE_NET_MASS_EPS_G):
            return PolymerFluxArchetype.VOLATILE_EJECTION
        return PolymerFluxArchetype.MIGRATION
    _warn_unresolved_archetype(
        "ambiguous cross-pool shape",
        tuple(sorted(p.label for p in product_polymers)))
    return PolymerFluxArchetype.UNRESOLVED


def compute_volatile_ejection_units(reactants, products, source_polymer: 'Polymer') -> float:
    """SIGNED source-monomer-equivalents transferred in a VOLATILE_EJECTION
    reaction: ``a = _net_nonpolymer_mass_g_mol(reactants, products) /
    source_polymer.monomer_mw_g_mol``. NOT rounded (alpha-methylstyrene off a
    styrene pool = 118.18/104.15 ~= 1.135).

    Signed (Codex round-13):
    - a > 0 : chain LOSES mass (net non-polymer PRODUCT mass, e.g. depropagation).
    - a < 0 : chain GAINS mass (net non-polymer REACTANT mass, e.g. monomer /
      radical addition).
    Netting the reactant side is required so a bimolecular H_Abstraction debits
    only the ~1 H shed (a ~= MW(H)/monomer), NOT the full co-reactant.

    Raises on a missing/zero monomer_mw_g_mol rather than silently producing 0
    (which would fabricate mass conservation)."""
    monomer_mw = getattr(source_polymer, 'monomer_mw_g_mol', 0.0) or 0.0
    if monomer_mw <= 0.0:
        raise ValueError(
            f"VOLATILE_EJECTION: source polymer {getattr(source_polymer, 'label', '?')!r} "
            f"has no positive monomer_mw_g_mol ({monomer_mw!r}); cannot compute "
            "eject_units without fabricating mass conservation.")
    return _net_nonpolymer_mass_g_mol(reactants, products) / monomer_mw


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
        # Demoted chip is a discrete gas-phase Molecule, not a Polymer
        # fold-back -- scrub the stale handshake proxy stamp AND stamp the
        # durable gas-volatile veto so it is not solver-visible as a melt
        # participant (see set_polymer_gas_veto).
        clear_polymer_proxy(chip_mol)
        set_polymer_gas_veto(chip_mol)
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
    # The discrete co-product chip stays a Molecule (only the daughter is
    # replaced by a fold-back Polymer below) -- scrub its stale proxy stamp
    # AND stamp the durable gas-volatile veto (see set_polymer_gas_veto).
    clear_polymer_proxy(chip_mol)
    set_polymer_gas_veto(chip_mol)
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
    # VOLATILE_EJECTION: stamp the fractional source-monomer-equivalents that
    # leave as discrete volatile(s) (a = sum non-polymer product MW /
    # source monomer_mw_g_mol). The solver's mass-losing dispatch (later phase)
    # and the sidecar emit read this STORED value.
    if forward.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION):
        source = polymer_reactants[0] if polymer_reactants else next(
            (r for r in reactants if isinstance(r, Polymer)), None)
        a = float(
            compute_volatile_ejection_units(reactants, forward.products, source))
        forward.polymer_eject_units = a
        # Diagnostic-only census (signed-VE spec): an atom-transfer-scale |a|
        # (< 0.5 monomer-equivalents, e.g. an H shed in H_Abstraction) means mu1
        # is being moved as a mass bucket, not a DP count. Warn once; routing is
        # unchanged (DP-vs-mass threshold deferred).
        if abs(a) < _VE_ATOM_TRANSFER_UNITS:
            _warn_atom_transfer_ve(_reaction_census_label(forward), a)
    # Refuse detection (item 18, spec section 4): a reaction stamped UNRESOLVED
    # because a polymer reactant produced no polymer product -- where a
    # non-polymer (gas) product structurally classifies FEATURE against a
    # polymer reactant's pool -- is a mid-chain radical the handshake dropped to
    # gas (mass-fabricating in core). FLAG it refused; keep the reaction and its
    # products untouched (stamp-but-keep). NEVER raise here: raising would
    # trigger the discard path. The infeasible-chip-surgery UNRESOLVED above
    # returns early and never reaches this check (correct: that is an
    # end-scission shape, not a FEATURE-radical-lost-to-gas).
    if forward.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED):
        lost = _chain_radical_lost_to_gas(forward, polymer_reactants)
        if lost is not None:
            forward.polymer_refused = True
            forward.polymer_refused_accumulating = (
                not is_qssa_eliminating_radical(lost))


# r92 flip-restamp census accumulators (module-level, mirror the warn-once
# census sets above). Refusals log a cumulative WARNING every time (they are
# rare and each one is a zeroed row); clean restamps log an INFO line once per
# (family, before -> after) signature and DEBUG otherwise, always carrying the
# cumulative clean-restamp count.
_flip_restamp_census = {
    "restamped": 0,
    "refused": 0,
    "refused_archetypes": set(),
    "refused_families": set(),
    "refused_rows": [],
    "restamp_signatures": set(),
}


def _archetype_name(value) -> str:
    try:
        return PolymerFluxArchetype(int(value)).name
    except ValueError:
        return str(value)


def restamp_flipped_polymer_archetype(forward) -> None:
    """
    r92 flip-restamp-or-refuse (PP run-10 killer, artifact rows r8/r30-32).

    ``apply_kinetics_to_reaction`` flips a reaction's reactants/products in
    place when the kinetics are estimated in reverse. The polymer flux
    archetype stamped by ``make_new_reaction`` encodes the GENERATION
    direction (parent/daughter roles, eject sign), so it is stale after the
    flip. The pre-r92 behavior -- blind demotion to UNRESOLVED -- dispatched
    live legacy mu1-only flux with RESOLVED pools: the r71-banned
    unclassified-pool-flux class through a generation-time door.

    Instead: CLEAR the stale direction-bound fields and RE-RUN the flux
    classification on the FLIPPED direction.

    * Classification is STAMP-BLIND (``trust_reacted_class=False``): the
      participant lists hold registry objects, and a registry daughter
      Polymer persistently carries the ``_reacted_class`` stamped by the
      reaction that CREATED it (``_register_polymer`` registers the reacted
      copy itself) -- species-level residue describing a different reaction.
      Only pool labels + net non-polymer mass are orientation-honest here,
      and they ARE direction-aware: a flipped MIGRATION re-derives its
      swapped src/dst, a flipped VOLATILE_EJECTION recomputes its eject
      units with the flipped SIGN (the chain now gains what it lost).
    * Classification UNRESOLVED, or required metadata uncomputable
      (VOLATILE_EJECTION ``eject_units`` with no positive source monomer
      MW): REFUSE conduit-deferred -- ``polymer_refused=True``,
      ``polymer_refused_accumulating=False`` -- zero whole-row flux via the
      solver's ``reaction_refused``, never a live legacy row.
    * An upstream refusal is preserved (OR semantics; qssa-invalid wins),
      matching :func:`merge_polymer_adjudication_stamps`.
    * ``is_end_group_reaction`` is deliberately kept: it flags WHERE on the
      chain the event happens (chain end -> mu0 scaling), which a direction
      flip does not move; its END_MOD stamp basis is unreconstructable here
      for the same registry-residue reason above.

    Pure-gas rows (no Polymer participant) are left untouched.
    """
    participants = (forward.reactants or []) + (forward.products or [])
    if not any(isinstance(s, Polymer) for s in participants):
        return

    before = int(getattr(forward, 'polymer_flux_archetype', 0))
    # CLEAR the stale direction-bound archetype fields (written for the
    # pre-flip direction) before re-classifying.
    forward.polymer_flux_archetype = int(PolymerFluxArchetype.NONE)
    forward.polymer_chip_units = 0
    forward.polymer_eject_units = 0.0

    new_arch = classify_reaction_flux_archetype(
        forward.reactants or [], forward.products or [],
        trust_reacted_class=False)

    refusal_detail = None
    if new_arch == PolymerFluxArchetype.VOLATILE_EJECTION:
        source = next((r for r in (forward.reactants or [])
                       if isinstance(r, Polymer)), None)
        try:
            forward.polymer_eject_units = float(compute_volatile_ejection_units(
                forward.reactants, forward.products, source))
        except (ValueError, AttributeError) as e:
            refusal_detail = "eject_units uncomputable: %s" % e
            new_arch = PolymerFluxArchetype.UNRESOLVED

    fam = str(getattr(forward, 'family', None)
              or getattr(forward, 'library', None) or '?')
    c = _flip_restamp_census
    if new_arch == PolymerFluxArchetype.UNRESOLVED:
        forward.polymer_flux_archetype = int(PolymerFluxArchetype.UNRESOLVED)
        forward.polymer_refused = True
        # conduit-deferred; an upstream qssa-invalid refusal is preserved.
        forward.polymer_refused_accumulating = bool(
            getattr(forward, 'polymer_refused_accumulating', False))
        c["refused"] += 1
        c["refused_archetypes"].add(_archetype_name(before))
        c["refused_families"].add(fam)
        if len(c["refused_rows"]) < 5:
            c["refused_rows"].append(_reaction_census_label(forward))
        logging.warning(
            "FLIPPED-POLYMER RESTAMP REFUSAL: %d kinetics-flipped polymer "
            "row(s) could not be restamped on the flipped direction and are "
            "refused conduit-deferred (zero whole-row flux; never live "
            "legacy mu1) [latest: %s]; archetypes-before=%s families=%s "
            "first_rows=%s",
            c["refused"],
            refusal_detail or "flipped-direction classification UNRESOLVED",
            ",".join(sorted(c["refused_archetypes"])),
            ",".join(sorted(c["refused_families"])),
            "; ".join(c["refused_rows"]))
        return

    forward.polymer_flux_archetype = int(new_arch)
    c["restamped"] += 1
    sig = (fam, _archetype_name(before), new_arch.name)
    level = logging.DEBUG if sig in c["restamp_signatures"] else logging.INFO
    c["restamp_signatures"].add(sig)
    logging.log(
        level,
        "FLIPPED-POLYMER RESTAMP: %d kinetics-flipped polymer row(s) "
        "restamped cleanly so far (latest %s -> %s, family %s: %s)",
        c["restamped"], _archetype_name(before), new_arch.name, fam,
        _reaction_census_label(forward))


def stamp_gas_association_refusal(forward, pool_registry=None) -> None:
    """PP v1 campaign refusal (adjudicated adversarial round 63, grounded in
    the run-5 rerun): an R_Recombination-style row that bridges PURE gas-phase
    radicals into a condensed polymer proxy -- e.g.
    ``[CH2]C(C)C + C[CH]CCC <=> polypropylene`` -- carries an unpaired
    thermo reference state (U ~ 10-12 decades; species thermo is uniformly
    gas-referenced, and this shape has NO same-mass chain counterparty to
    cancel against), exactly the class the solver's reference-state tripwire
    exists to stop. The correct LONG-TERM representation is a
    pool-moment-credit conduit (out of scope); until it lands the shape is
    predeclared out of the PP v1 campaign scope and REFUSED via the shipped
    refused-row contract (schema 2.4): stamp-but-keep, whole flux zeroed in
    the solver / reference runner / TA, sidecar reason "conduit-deferred"
    (``polymer_refused_accumulating`` stays False -- the closed-vocabulary
    reason for rows deferred pending the conduit).

    SHAPE-specific, deliberately NOT family-specific:

    * gas radical(s) -> condensed proxy (association orientation, the run-5
      generated direction) and condensed proxy -> gas radicals (homolysis,
      the reverse generated orientation) are refused;
    * gas+gas->gas R_Recombination termination has no condensed side and is
      untouched (whitelisted chemistry);
    * H-abstraction / routing shapes carry the proxy on the mixed side --
      neither side is "all gas radicals" -- untouched;
    * beta-scission / volatile-producing shapes carry a closed-shell product
      on the all-gas side -- untouched.

    Phase info at this layer: "condensed polymer proxy" means an
    :class:`Polymer` participant (pool proxies and spawned daughters ARE
    ``Polymer`` Species here -- the association product resolves onto the
    registered pool object via ``species_dict`` isomorphism in
    ``make_new_species``). "Gas" is any non-``Polymer`` participant; the
    sticky ``is_polymer_proxy`` tag is deliberately NOT consulted
    (family.py blanket-tags small radicals on proxy-touching reactions, and
    keying on it would let the exact run-5 radicals slip through).

    Called from ``make_new_reaction`` (rmgpy/rmg/model.py) AFTER the resolved
    reactants/products are assigned onto ``forward`` -- the association
    orientation has no polymer REACTANT, so the ``polymer_reactants``
    handshake/stamping block never sees it; this is the one seam where both
    resolved sides are visible for every branch.

    Unconditional engine behavior, matching the existing item-18 refusal
    stamp (no deck-level gate exists for refusals today); the long-term
    conduit will replace it. An earlier refusal stamp (item-18 detector) is
    never overwritten -- its census reason wins.

    SECOND refused shape (adjudicated adversarial round 74, grounded in the
    PP run-5 TA diagnostic): the SAME-PROXY degenerate row, e.g.
    ``H(1) + rad_end <=> rad_end`` -- the resolved pool proxy Species is
    IDENTICAL on both sides (the H-capped adduct folds back onto the
    rad_end pool in the handshake and resolves to the same registered
    object) while gas species appear/disappear across the row. Because the
    polymer participant cancels, Keq collapses to the gas species' free
    energy alone and the thermo-reversed direction becomes a unimolecular
    gas source: run-5's TA replay volatilized the entire 10 mg PP sample as
    pure H2 at 889 C -- element-impossible (carbon leaving as hydrogen).
    The row is a moment identity wearing a reaction's clothes; its Keq is
    meaningless BY CONSTRUCTION. See :func:`_stamp_same_proxy_refusal` for
    the three-conjunct predicate (incl. the DP-scale carve-out that keeps
    genuine same-pool depropagation/unzip/chip rows live).

    THIRD refused shape (adjudicated adversarial round 82, grounded in FR1
    run-2): IMPOSTOR row, e.g. ``BrBr + <C36-scale unsaturated discrete>
    <=> FR1`` (XY_Addition_MultipleBond, 15 rows in run-2). Exactly one
    side carries a Polymer participant; the other side carries NO Polymer
    but at least one POLYMER-SIZED discrete molecule (the pool proxy minus
    a small closed-shell partner, living in the model as an ordinary gas
    Species). The gas partner is closed-shell, so the r63 all-gas-radicals
    conjunct never fires -- pre-r82 these rows arrived UNSTAMPED at the
    solver rebuild and died at the r71 unstamped-proxy hard-fail. Same
    reference-state bridge, same contract: stamp-but-keep,
    conduit-deferred, both written orientations. "Polymer-sized" is decided
    by the conservative monomer-scale threshold
    ``_IMPOSTOR_DISCRETE_MONOMER_UNITS`` (2.5 source-monomer-equivalents;
    proxy-derivation is not traceable at this seam), which keeps DP-2 dimer
    volatiles (e.g. hexene off polypropylene, 2.0 monomer-equivalents --
    pinned negative control) live.

    FOURTH and FIFTH refused shapes (adjudicated adversarial round 87,
    grounded in FR1 run-3, /home/alon/Projects/polymer/FR1/rmg/run3): rows
    coupling an UNREPRESENTED chain-scale proxy-DERIVED discrete species to
    polymer pool state (Br-ring-addition adducts of the FR1 trimer proxy:
    de-aromatized defect states, neither pool-representable nor genuinely
    gas -- run-3's (110)-(119) cohort, C36H32Br19O3 = proxy + 2 Br).

    * Shape A (:func:`_stamp_chain_scale_defect_refusal`): UNEQUAL-count
      same-pool fold-back, ``adduct + pool <=> pool + pool`` -- run-3
      serialized 20 such rows LIVE as ``same_pool/1`` (the [pool] vs
      [pool, pool] participant-count step escapes both r74's
      identical-pairing conjunct and the classifier's discrete-VE reroute,
      which requires exactly one polymer product), so SAME_POOL's zero
      moment change fabricated/annihilated the whole adduct mass
      (element-unbalanced by 2 Br).
    * Shape B (:func:`_stamp_reference_state_split_refusal`): the
      discrete<=>discrete isomerization between such adducts that SPLITS
      the reference-state classification (run-3 ``(117) <=> (111)``, one
      side melt-classified, the other veto-suppressed, U ~ 13 decades) --
      the row has no Polymer participant at all, so the monomer scale is
      taken from ``pool_registry`` (a list of registered :class:`Polymer`
      pools, or a zero-arg callable returning one, resolved lazily; when
      absent the stamp stays conservative and the solver tripwire remains
      the loud backstop).

    Same contract as r63/r74/r82: stamp-but-keep, conduit-deferred, both
    written orientations; the third shape (BrBr + pool <=> adduct) was
    already refused qssa-invalid by the item-18 detector and the top
    early-return keeps that census reason.
    """
    if getattr(forward, "polymer_refused", False):
        return  # already refused upstream; keep that census reason
    reactants = getattr(forward, "reactants", None) or []
    products = getattr(forward, "products", None) or []
    polymer_r = [s for s in reactants if isinstance(s, Polymer)]
    polymer_p = [s for s in products if isinstance(s, Polymer)]
    r_condensed = bool(polymer_r)
    p_condensed = bool(polymer_p)
    if r_condensed and p_condensed:
        # condensed on BOTH sides: ordinarily routing chemistry, EXCEPT the
        # r74 same-proxy degenerate shape (run-5 H fountain) and the r87
        # unequal-count same-pool fold-back coupling a chain-scale
        # proxy-derived discrete adduct (FR1 run-3; checked second so an
        # r74 stamp -- same census reason -- is never re-derived).
        _stamp_same_proxy_refusal(forward, reactants, products,
                                  polymer_r, polymer_p)
        if not getattr(forward, "polymer_refused", False):
            _stamp_chain_scale_defect_refusal(forward, reactants, products,
                                              polymer_r, polymer_p)
        return
    if not r_condensed and not p_condensed:
        # no condensed participant at all: ordinary gas chemistry, EXCEPT
        # the r87 reference-state-SPLIT isomerization between chain-scale
        # proxy-derived discretes (FR1 run-3 ``(117) <=> (111)``) -- the
        # row carries no Polymer participant, so the monomer scale comes
        # from ``pool_registry``.
        _stamp_reference_state_split_refusal(forward, reactants, products,
                                             pool_registry)
        return

    def _all_gas_radicals(side):
        if not side:
            return False
        for s in side:
            if isinstance(s, Polymer):
                return False
            mol_list = getattr(s, "molecule", None)
            mol = mol_list[0] if mol_list else (
                s if isinstance(s, Molecule) else None)
            if mol is None or mol.get_radical_count() < 1:
                return False
        return True

    if ((p_condensed and _all_gas_radicals(reactants))
            or (r_condensed and _all_gas_radicals(products))):
        forward.polymer_refused = True
        forward.polymer_refused_accumulating = False  # -> "conduit-deferred"
        return

    # THIRD refused shape (adjudicated adversarial round 82, grounded in FR1
    # run-2, /home/alon/Projects/polymer/FR1/rmg/run2/RMG.out ~line 1774):
    # IMPOSTOR row. Exactly one side carries the Polymer participant
    # (guaranteed here -- the both-condensed and no-condensed branches
    # returned above); the other side carries NO Polymer but at least one
    # POLYMER-SIZED discrete molecule (>= _IMPOSTOR_DISCRETE_MONOMER_UNITS
    # source-monomer-equivalents) -- a proxy-scale structure living in the
    # model as an ordinary gas Species. Run-2: XY_Addition_MultipleBond
    # generated 15 rows ``Br/BrBr + <C36-scale unsaturated proxy-minus-XY>
    # <=> FR1``; the gas partner is CLOSED-SHELL, so the r63
    # all-gas-radicals conjunct never fires, and the rows arrived UNSTAMPED
    # at the solver rebuild (r71 hard-fail killed the run). The shape is
    # the same reference-state bridge as r63 -- a discrete gas-referenced
    # participant with no same-mass counterparty crossing into the
    # condensed pool -- and is refused identically: stamp-but-keep,
    # conduit-deferred, orientation-independent by construction.
    poly_side = polymer_p if p_condensed else polymer_r
    discrete_side = reactants if p_condensed else products
    for s in discrete_side:
        if any(_discrete_is_polymer_sized(s, poly)
               for poly in poly_side):
            forward.polymer_refused = True
            forward.polymer_refused_accumulating = False  # "conduit-deferred"
            return


def _heavy_atom_count(mol) -> int:
    """Number of non-hydrogen atoms in ``mol`` (0 when structureless)."""
    if mol is None:
        return 0
    return mol.get_num_atoms() - mol.get_num_atoms('H')


def _warn_impostor_axis_undecidable(species, poly, missing) -> None:
    """Warn-once census (r85 P2(a)): an impostor-sized decision arrived with
    an UNCOMPUTABLE axis and the computable evidence (if any) did not
    already decide False -- the case is undecidable and the predicate
    refuses to refuse. Announced instead of silently degenerating to a
    single-axis refusal. Reuses the ``_flux_archetype_warned`` warn-once
    set, keyed to avoid spam."""
    pool_label = getattr(poly, 'label', None)
    spc_label = getattr(species, 'label', None) or repr(species)
    key = ("impostor_axis_undecidable", pool_label, spc_label, missing)
    if key not in _flux_archetype_warned:
        _flux_archetype_warned.add(key)
        logging.warning(
            "IMPOSTOR AXIS UNDECIDABLE (r85 P2a census): discrete '%s' vs "
            "pool '%s' -- the %s axis is uncomputable, so the r82 "
            "polymer-sized conjunct cannot be established on BOTH axes; "
            "NOT refusing (never refuse blind; the r71 unstamped-proxy "
            "hard-fail stays the loud backstop).",
            spc_label, pool_label, missing)


def _discrete_is_polymer_sized(species, poly) -> bool:
    """True when the non-``Polymer`` participant ``species`` is POLYMER-SIZED
    relative to pool ``poly``'s monomer: at/above
    ``_IMPOSTOR_DISCRETE_MONOMER_UNITS`` monomer-equivalents on BOTH axes,
    and BOTH axes computable (r82 impostor conjunct, hardened r85 P2(a)).

    Two axes, conservative-by-conjunction:

    * mass: ``MW(species) / poly.monomer_mw_g_mol``;
    * structure: heavy-atom count ratio vs the monomer structure. This axis
      is what keeps a HEAVY small gas honest -- Br2 is 3.8 propene-monomer-
      equivalents of MASS but 0.67 of heavy atoms (a diatomic is never
      polymer-sized).

    r85 P2(a): an axis whose denominator is unavailable (no monomer_mw / no
    monomer structure on a defensive copy) does NOT defer to the other --
    a single-axis refusal is exactly the Br2-with-stripped-monomer false
    positive. Both axes must be computable AND at/above threshold to
    refuse; otherwise the answer is False, and the genuinely undecidable
    case (the computable evidence did not already decide False) is
    announced through the warn-once census."""
    if isinstance(species, Polymer):
        return False
    mol_list = getattr(species, 'molecule', None)
    mol = mol_list[0] if mol_list else (
        species if isinstance(species, Molecule) else None)
    if mol is None:
        return False
    # Mass axis.
    monomer_mw = float(getattr(poly, 'monomer_mw_g_mol', 0.0) or 0.0)
    mass_computable = monomer_mw > 0.0
    mass_polymer_sized = (
        mass_computable
        and mol.get_molecular_weight() * 1000.0
        >= _IMPOSTOR_DISCRETE_MONOMER_UNITS * monomer_mw)
    # Structure (heavy-atom) axis.
    monomer_heavy = _heavy_atom_count(getattr(poly, 'monomer', None))
    structure_computable = monomer_heavy > 0
    structure_polymer_sized = (
        structure_computable
        and _heavy_atom_count(mol)
        >= _IMPOSTOR_DISCRETE_MONOMER_UNITS * monomer_heavy)
    if mass_computable and structure_computable:
        return mass_polymer_sized and structure_polymer_sized
    # At least one axis uncomputable: never refuse. Announce the case as
    # undecidable UNLESS a computable axis already decided False (then the
    # answer is a decided False, not an undecidable one).
    if not any([mass_computable and not mass_polymer_sized,
                structure_computable and not structure_polymer_sized]):
        missing = ("mass and structure" if not mass_computable
                   and not structure_computable
                   else ("mass" if not mass_computable else "structure"))
        _warn_impostor_axis_undecidable(species, poly, missing)
    return False


def _discrete_is_chain_scale_proxy_derived(species, pools) -> bool:
    """r87 classifier (FR1 run-3): True when the non-``Polymer`` participant
    ``species`` is a CHAIN-SCALE PROXY-DERIVED discrete -- it carries
    proxy-derivation evidence AND is polymer-sized against at least one pool
    in ``pools``.

    Evidence is the sticky ``is_polymer_proxy`` melt tag OR the durable gas
    veto (:func:`has_polymer_gas_veto`) -- both are stamped only by
    proxy-touching machinery, and the run-3 adduct cohort splits across
    them ((111)/(116) melt-tagged, (117) gas-vetoed), so r87 requires the
    OR: veto status must NOT change the refusal outcome.

    Size is :func:`_discrete_is_polymer_sized` -- the r82/r85 dual-axis
    (mass AND heavy-atom) polymer-sized threshold at
    ``_IMPOSTOR_DISCRETE_MONOMER_UNITS`` (2.5) monomer-equivalents, NOT the
    solver tripwire's per-repeat-unit MW window (monomer + 10 g/mol).
    Deliberate, r87 over-refusal pin: declared-monomer-routing volatiles
    clear the tripwire window (AMS 118.2 vs the styrene window 114.2,
    hexene 84.2 vs the propene window 52.1) and are routinely
    proxy-contaminated and/or gas-vetoed, so a window-based conjunct would
    refuse them -- 'the next over-refusal bug'. The polymer-sized threshold
    keeps all of them (propene 1.0, hexene 2.0, AMS 1.13
    monomer-equivalents) below the bar while the run-3 adducts (3.13
    monomer-equivalents, 2030.8 g/mol) clear it on both axes."""
    if isinstance(species, Polymer):
        return False
    if not (getattr(species, 'is_polymer_proxy', False)
            or has_polymer_gas_veto(species)):
        return False
    return any(_discrete_is_polymer_sized(species, poly) for poly in pools)


def _stamp_chain_scale_defect_refusal(forward, reactants, products,
                                      polymer_r, polymer_p) -> None:
    """r87 shape A (FR1 run-3): refuse (stamp-but-keep, census reason
    "conduit-deferred") the UNEQUAL-count same-pool fold-back coupling a
    chain-scale proxy-derived discrete adduct to pool state, e.g. run-3's
    20 live ``same_pool/1`` Disproportionation-Y rows
    ``adduct + pool <=> pool + pool`` (element-unbalanced by 2 Br). All
    conjuncts required -- dropping any one keeps the row LIVE (pinned):

    1. every Polymer participant folds back (reactant and product pool
       label SETS equal -- label-changing rows are MIGRATION's routing
       business);
    2. Polymer participant COUNTS differ between the sides (the [pool] vs
       [pool, pool] step that escaped r74's identical-pairing conjunct and
       that the classifier's SAME_POOL/VE fork cannot represent: SAME_POOL
       applies zero moment change, and the discrete-VE reroute requires
       exactly one polymer product -- equal-count fold-backs are the
       legitimate VE/deprop/chip shapes and stay live);
    3. some discrete participant on either side is chain-scale
       proxy-derived against a pool of THIS row
       (:func:`_discrete_is_chain_scale_proxy_derived`).

    Orientation-independent by construction (both written orientations
    present the same label sets, count asymmetry and participants)."""
    if ({getattr(p, 'label', None) for p in polymer_r}
            != {getattr(p, 'label', None) for p in polymer_p}):
        return  # cross-pool routing chemistry
    if len(polymer_r) == len(polymer_p):
        return  # equal-count fold-back: representable VE/deprop/chip shapes
    pools = polymer_r + polymer_p
    for s in list(reactants) + list(products):
        if _discrete_is_chain_scale_proxy_derived(s, pools):
            forward.polymer_refused = True
            forward.polymer_refused_accumulating = False  # "conduit-deferred"
            return


def _stamp_reference_state_split_refusal(forward, reactants, products,
                                         pool_registry) -> None:
    """r87 shape B (FR1 run-3): refuse (stamp-but-keep, "conduit-deferred")
    the discrete<=>discrete row -- NO Polymer participant on either side --
    that (a) involves a chain-scale proxy-derived discrete
    (:func:`_discrete_is_chain_scale_proxy_derived`) and (b) SPLITS the
    thermo reference-state classification: the per-side multisets of
    melt-classified participant MWs differ, so the solver tripwire's U
    cannot pair off and the build dies on the cliff-sign ValueError (run-3:
    ``(117) <=> (111)``/``(116)``, U = 13.10 decades, RMG.out:5140-5144).

    The melt classification mirrors the tripwire's tag branch EXACTLY
    (polymer.pyx ``_reference_state_tripwire``, r89 dual-axis amendment of
    spec 5.1 C3): ``is_polymer_proxy`` AND NOT gas-vetoed AND polymer-sized
    against at least one registered pool per
    :func:`_discrete_is_polymer_sized` (mass AND heavy-atom axes both
    computable and at/above ``_IMPOSTOR_DISCRETE_MONOMER_UNITS``
    monomer-equivalents). ONE size predicate shared with the classifier and
    (constant-shared) with the solver gate, so the stamp and the sensor
    cannot drift apart -- the pre-r89 window form (MW >= largest pool
    monomer + slack) melt-classified DP-2 gas volatiles (PP run-9
    1,5-hexadiene, 82.15 g/mol > propene window 52.1 but 2.0
    monomer-equivalents) and disagreed with conjunct (a)'s dual-axis
    membership. A NO-split row (both sides melt-classified isomers) pairs
    off exactly in U and stays LIVE (pinned negative control).

    ``pool_registry`` is a list of registered :class:`Polymer` pools or a
    zero-arg callable returning one (resolved lazily, AFTER the cheap
    evidence pre-gate, so gas-only decks never pay the collection); with no
    registry the stamp answers conservatively (live) and the tripwire
    remains the loud backstop. Both chokepoints supply real ``Polymer``
    objects (make_new_reaction collects core/edge/new species; the r71
    solver rebuild restamp collects core/edge species), keeping the r85
    dual-axis size gate computable -- the solver's ``PolymerPoolConfig``
    sidecars carry no monomer structure and are deliberately NOT used."""
    participants = list(reactants) + list(products)
    if not any(getattr(s, 'is_polymer_proxy', False)
               or has_polymer_gas_veto(s) for s in participants):
        return  # cheap evidence pre-gate: no proxy-derivation anywhere
    if callable(pool_registry):
        pool_registry = pool_registry()
    pools = [p for p in (pool_registry or []) if isinstance(p, Polymer)]
    if not pools:
        return  # no pool context: conservative; tripwire stays the backstop
    if not any(_discrete_is_chain_scale_proxy_derived(s, pools)
               for s in participants):
        return

    def _melt_mws(side):
        # r89 dual-axis melt mirror: membership = proxy-tagged AND unvetoed
        # AND polymer-sized (both axes) against >= 1 pool -- the SAME
        # predicate conjunct (a) used above, and the same gate the solver
        # tripwire applies, so refused/live sets cannot drift. A
        # structureless participant is undecidable to
        # _discrete_is_polymer_sized (mol is None -> False) and stays
        # conservative-gas here exactly as at the solver seam.
        out = []
        for s in side:
            if not getattr(s, 'is_polymer_proxy', False):
                continue
            if has_polymer_gas_veto(s):
                continue
            mol_list = getattr(s, 'molecule', None)
            mol = mol_list[0] if mol_list else (
                s if isinstance(s, Molecule) else None)
            if mol is None:
                continue
            if not any(_discrete_is_polymer_sized(s, poly) for poly in pools):
                continue
            mw = mol.get_molecular_weight() * 1000.0
            # Round away summation-order ulps so isomers pair exactly.
            out.append(round(mw, 6))
        return sorted(out)

    if _melt_mws(reactants) == _melt_mws(products):
        return  # classification pairs off; the tripwire's U cancels
    forward.polymer_refused = True
    forward.polymer_refused_accumulating = False  # -> "conduit-deferred"


def _polymer_participants_identical(polymer_reactants, polymer_products) -> bool:
    """True when the Polymer participants on the two sides of a row pair off
    IDENTICALLY: same resolved pool proxy Species by object identity (the
    post-``make_new_species`` / rebuild state -- ``species_dict`` and the
    solver's core/edge lists hand every side the same registered object),
    with a registered-pool label fallback for copies (pool labels are unique
    by ``_register_polymer``'s first-writer-wins dedup). The fallback
    additionally requires an equal ``_reacted_class``: a handshake
    reacted-copy daughter (FEATURE / SCISSION / CHIP / ...) carries that
    stamp as the RECORD of a polymer state transition, so an unresolved
    same-label daughter never pairs off against its unreacted parent (the
    pinned H-abstraction routing negative control). Identical pairing means
    NO polymer participant undergoes any pool/proxy/state transition
    anywhere in the row -- r74 conjuncts (1) and (3a). Any polymer that
    moves to a distinct pool (real H-capping's target end-state, feature
    '_mod' daughters, scission daughters) breaks the pairing and the row is
    NOT the degenerate shape."""
    if not polymer_reactants or len(polymer_reactants) != len(polymer_products):
        return False
    unmatched = list(polymer_products)
    for r in polymer_reactants:
        for p in unmatched:
            if r is p or (getattr(r, 'label', '')
                          and r.label == getattr(p, 'label', None)
                          and getattr(r, '_reacted_class', None)
                          == getattr(p, '_reacted_class', None)):
                unmatched.remove(p)
                break
        else:
            return False
    return True


def _stamp_same_proxy_refusal(forward, reactants, products,
                              polymer_r, polymer_p) -> None:
    """r74 same-proxy refusal (PP run-5 "H fountain"; forensics
    ~/Projects/polymer/PP/rmg/run5 + ta-diag-run5). Refuse (stamp-but-keep,
    census reason "conduit-deferred", identical to the r63 gas-association
    refusal contract) a condensed-both-sides row when ALL of:

    1. the resolved polymer pool/proxy participants are IDENTICAL on the
       reactant and product side (:func:`_polymer_participants_identical`),
       AND
    2. the non-polymer (gas) sides genuinely DIFFER -- gas species appear or
       disappear across the row (a pure identity row nets nothing and is
       not the laundered shape), AND
    3. there is no distinct polymer state transition the row could be
       booking. Decided at DP scale (the reformulation adjudication asked
       to be verified at the chokepoint): a net gas-side mass change BELOW
       ``_VE_ATOM_TRANSFER_UNITS`` (0.5) source-monomer-equivalents cannot
       be a repeat-unit-count change -- it claims a chemical end-state
       transition (H-capping / H-loss) that same-pool resolution cannot
       represent, so mu1 would be moved as a laundered mass bucket with a
       Keq of the gas species alone (run-5's |a| = MW(H)/monomer = 0.0240,
       the H fountain). A monomer-scale change (|a| >= 0.5) IS a genuine
       within-pool DP transition -- depropagation / unzip / discrete-chip,
       the solver VE/chip dispatch's design domain -- and stays LIVE
       (pinned bitwise by r71's
       test_negative_control_live_ve_row_with_real_polymer_keeps_flux).

    Real H-capping to a DISTINCT target pool state fails conjunct 1 and is
    never refused. Rides every r71 chokepoint for free: this function is
    called from :func:`stamp_gas_association_refusal`, which runs at
    make_new_reaction classification, merges across canonical dedup
    (``merge_polymer_adjudication_stamps``), is idempotently re-run over
    chain(core, edge) at the solver rebuild restamp, and feeds the
    item-18/r71 flux gates (zero core RHS, zero edge enlargement inputs)
    and the unstamped-proxy hard-fail exemption."""
    if not _polymer_participants_identical(polymer_r, polymer_p):
        return  # a polymer moves pool/proxy/state -- real conduit chemistry
    gas_r = defaultdict(int)
    gas_p = defaultdict(int)
    for s in reactants:
        if not isinstance(s, Polymer):
            gas_r[id(s)] += 1
    for s in products:
        if not isinstance(s, Polymer):
            gas_p[id(s)] += 1
    if gas_r == gas_p:
        return  # pure identity row: no net gas change, nothing laundered
    monomer_mw = max((float(getattr(p, 'monomer_mw_g_mol', 0.0) or 0.0)
                      for p in polymer_r), default=0.0)
    net_mass = abs(_net_nonpolymer_mass_g_mol(reactants, products))
    if monomer_mw > 0.0 and net_mass >= _VE_ATOM_TRANSFER_UNITS * monomer_mw:
        # Monomer-scale DP transition (depropagation/unzip/chip):
        # representable within the pool by the VE/chip dispatch -- live.
        return
    forward.polymer_refused = True
    forward.polymer_refused_accumulating = False  # -> "conduit-deferred"


def merge_polymer_adjudication_stamps(source, target) -> None:
    """r71 FIX 1 (PP run-5 stall): ``check_for_existing_reaction`` can discard
    a freshly stamped candidate and return a pre-existing UNSTAMPED canonical
    equivalent -- silently dropping the polymer adjudication from the model
    (the refused rows then run live legacy mu1-only flux against the pool,
    the run-5 DASPK collapse). Called from ``make_new_reaction`` when the
    canonical-dedup path returns an existing reaction: OR/merge the discarded
    ``source`` candidate's adjudication-bearing stamps onto the canonical
    ``target``.

    Merge rules (adjudicated):

    * ``polymer_refused`` ORs; ``polymer_refused_accumulating`` ORs across
      refused stamps, so the census reason "qssa-invalid" (accumulating=True)
      WINS over "conduit-deferred" in both merge directions and an existing
      qssa-invalid canonical is never demoted.
    * ``polymer_flux_archetype`` / ``polymer_chip_units`` /
      ``polymer_eject_units`` / ``is_end_group_reaction`` FILL the canonical's
      unstamped (NONE / zero / False) slots only -- a live canonical stamp is
      never overwritten (it was made with its own full handshake context).

    No-op for ordinary gas chemistry (every getattr defaults falsy).
    """
    if source is target:
        return
    if getattr(source, "polymer_refused", False):
        target.polymer_refused = True
        target.polymer_refused_accumulating = bool(
            getattr(target, "polymer_refused_accumulating", False)
            or getattr(source, "polymer_refused_accumulating", False))
    if (int(getattr(target, "polymer_flux_archetype", 0) or 0) == 0
            and int(getattr(source, "polymer_flux_archetype", 0) or 0) != 0):
        target.polymer_flux_archetype = int(source.polymer_flux_archetype)
    if (int(getattr(target, "polymer_chip_units", 0) or 0) == 0
            and int(getattr(source, "polymer_chip_units", 0) or 0) != 0):
        target.polymer_chip_units = int(source.polymer_chip_units)
    if (float(getattr(target, "polymer_eject_units", 0.0) or 0.0) == 0.0
            and float(getattr(source, "polymer_eject_units", 0.0) or 0.0) != 0.0):
        target.polymer_eject_units = float(source.polymer_eject_units)
    if getattr(source, "is_end_group_reaction", False):
        target.is_end_group_reaction = True


def _chain_radical_lost_to_gas(forward, polymer_reactants):
    """Return the gas-product :class:`Molecule` that a polymer reactant dropped to
    gas (``create_reacted_copy`` returned ``None``) and that carries CHAIN-SCALE
    mass -- i.e. a backbone radical the solver's UNRESOLVED single-monomer-debit
    fabricates mass for (spec section 4), else ``None``. The conservation bug
    fires for ANY chain-scale gas product on the UNRESOLVED "polymer reactant, no
    polymer product" leg, not only ``classify_structure == FEATURE`` -- the
    FEATURE/DISCARD split is a positional artifact of the 3-unit proxy (center vs
    cap-adjacent monomer), not chemistry, so both leak the same MW-211 C15
    backbone radical and fabricate the same mass. Mechanism-keyed, not
    label-keyed: dropping the FEATURE check avoids re-introducing the
    under-inclusive trap one label at a time.

    Size gate: a species is chain-scale iff its MW >= one monomer-unit plus the
    chain-window slack (the SAME ``REFERENCE_STATE_MW_SLACK_G_MOL`` the solver's
    ``chain_window_kg`` uses). Genuinely-small fragments (H2, CH4, small radicals;
    MW within one monomer + slack) fall below threshold and are NOT refused -- the
    single-monomer-debit accounts for those leaks correctly. Polymer products are
    skipped, which correctly excludes the conserving pool->pool case (e.g.
    ``epdm_scission_tail``).

    r89 audit note: this threshold DELIBERATELY stays window-based (one
    monomer + slack), NOT the dual-axis 2.5-monomer-equivalent melt gate.
    It is a per-event MASS-CONSERVATION debit audit, not a phase
    classification: the UNRESOLVED leg debits exactly ONE monomer, so ANY
    gas product materially heavier than one repeat unit (a DP-2 dimer
    volatile included) fabricates the difference and must be caught here.
    Raising the bar to 2.5 units would let ~2-unit products leak un-audited
    mass -- the run-9 DP-2 misclassification was a melt-CLASSIFICATION
    error (a different consumer), not a debit-audit error.

    Imported function-locally: there is a documented solver->polymer import cycle
    (``polymer.pyx`` avoids importing ``polymer.py`` at module level), so the
    constant must be pulled in at call time (generation time), not module import.

    The handshake leaves the leaked radical a plain ``Molecule``; the returned
    value is the ``Molecule`` (the downstream ``is_qssa_eliminating_radical``
    probe needs a ``Molecule``)."""
    from rmgpy.solver.polymer import REFERENCE_STATE_MW_SLACK_G_MOL
    for poly in polymer_reactants:               # Polymer instances
        threshold = poly.monomer_mw_g_mol + REFERENCE_STATE_MW_SLACK_G_MOL
        for prod in forward.products:
            if isinstance(prod, Polymer):
                continue
            mol = prod.molecule[0] if isinstance(prod, Species) else prod
            if mol.get_molecular_weight() * 1000.0 >= threshold:
                return mol
    return None


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


def _assert_end_radical_proxy(mol: Molecule, site: str) -> None:
    """
    STRICT validity assertion for an end-radical oligomer proxy (k_homolysis
    daughter pools, Stage 1 / round 66). Runs on the stitched chain BEFORE
    labels are cleared, so the surviving stitch label is the terminality
    witness. Deliberately its own assertion path: _assert_feature_unit
    (mid-chain, 2-3 radicals + both stitch labels) and _assert_end_group
    (stitch fragment) describe different shapes and are NOT relaxed.

    Requirements:
    - exactly ONE surviving open-site label, from the site's label set
      (*1 for 'primary', *2 for 'secondary'), and no other labels;
    - the molecule is MONO-radical (exactly one radical electron total);
    - the radical sits ON the labeled terminal atom, which is a heavy atom.
    """
    want = LABELS_1 if site == 'primary' else LABELS_2
    labeled = [a for a in mol.atoms if a.label]
    on_site = [a for a in labeled if a.label in want]
    if len(on_site) != 1 or len(labeled) != 1:
        raise ValueError(
            f"End-radical proxy ({site}) must carry exactly one surviving "
            f"open-site label from {want} and no other labels; got "
            f"{[a.label for a in labeled] or 'none'}.")
    rad_count = mol.get_radical_count()
    if rad_count != 1:
        raise ValueError(
            f"End-radical proxy ({site}) must be mono-radical (exactly one "
            f"radical electron total). Got {rad_count}.")
    atom = on_site[0]
    if atom.is_hydrogen():
        raise ValueError(
            f"End-radical proxy ({site}): the open-site atom must be a "
            f"terminal HEAVY atom, got a hydrogen.")
    if atom.radical_electrons != 1:
        raise ValueError(
            f"End-radical proxy ({site}): the radical must sit ON the "
            f"terminal (stitch-labeled) atom; the labeled atom carries "
            f"{atom.radical_electrons} radical electrons.")


def _assert_side_loss_unit(unit: Molecule, element_symbol: str,
                           monomer: Molecule, channel_label: str,
                           radical_atom: Optional['Atom'] = None) -> None:
    """
    STRICT validity assertion for a side-group X-loss feature unit (FR1-K1;
    the analog of _assert_end_radical_proxy for the side-group kernel --
    deliberately its own assertion path: _assert_feature_unit and
    _assert_end_group describe different shapes and are NOT relaxed; the
    caller additionally runs _assert_feature_unit with the 3-radical
    budget).

    Requirements:
    - element census: the unit carries EXACTLY one fewer atom of
      element_symbol than the monomer, and identical counts for every
      other element (one X removed, nothing else touched);
    - the unit is 3-radical total (two stitch radicals + exactly ONE
      mid-chain defect radical from the homolyzed C-X bond);
    - when the builder passes the ex-neighbor handle (radical_atom), that
      atom carries at least one radical electron (the defect sits WHERE
      the X was lost).
    """
    def _census(mol):
        counts = {}
        for a in mol.atoms:
            counts[a.symbol] = counts.get(a.symbol, 0) + 1
        return counts

    unit_counts = _census(unit)
    mono_counts = _census(monomer)
    expected = dict(mono_counts)
    expected[element_symbol] = expected.get(element_symbol, 0) - 1
    if expected.get(element_symbol, 0) <= 0:
        expected.pop(element_symbol, None)
    if unit_counts != expected:
        raise ValueError(
            f"Side-group X-loss unit (channel '{channel_label}') must be "
            f"the monomer minus exactly ONE {element_symbol} atom; got "
            f"element census {unit_counts} vs expected {expected}.")
    rad_count = unit.get_radical_count()
    if rad_count != 3:
        raise ValueError(
            f"Side-group X-loss unit (channel '{channel_label}') must "
            f"carry exactly 3 radical electrons total (two stitch radicals "
            f"+ ONE mid-chain defect radical from the homolyzed C-"
            f"{element_symbol} bond). Got {rad_count}.")
    if radical_atom is not None and radical_atom.radical_electrons < 1:
        raise ValueError(
            f"Side-group X-loss unit (channel '{channel_label}'): the "
            f"defect radical must sit on the atom that lost the "
            f"{element_symbol} substituent; that atom carries "
            f"{radical_atom.radical_electrons} radical electrons.")


def _radical_site_descriptor(mol: Molecule) -> str:
    """
    Deterministic radical-placement descriptor for a feature-unit fingerprint
    (radical-feature producer path): a sorted, joined encoding of every
    radical-bearing atom as ``{label|element}{radicals}h{n_H}c{n_heavy}``.
    Distinguishes same-formula H-loss units by WHERE the extra radical sits
    (e.g. the ~3 PP C3H5 units: backbone CH2 / backbone CH / pendant CH3)
    while staying identical for isomorphic units, so positional twins still
    dedup to one pool. v1 coarseness is documented: exotic units whose
    radical sites share (element, radical count, H count, heavy-neighbor
    count) multisets would collide -- acceptable for the documented ~3-
    environment scope; refine when a real deck needs it.
    """
    sites = []
    for a in mol.atoms:
        if a.radical_electrons > 0:
            n_h = sum(1 for nb in a.bonds if nb.is_hydrogen())
            n_heavy = sum(1 for nb in a.bonds if not nb.is_hydrogen())
            key = a.label if a.label else a.symbol
            sites.append(f"{key}{a.radical_electrons}h{n_h}c{n_heavy}")
    return '-'.join(sorted(sites))


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

    def gate_statistic(self, motif_key: str) -> float:
        """Window sum divided by the FIXED window length (zero-filled
        semantics, spec 2026-06-10 §4.4 step 4): a single-snapshot spike must
        be ``window``x the bar to clear the gate; a channel persisting at
        fraction f for ``window`` iterations reads f."""
        return self.flux(motif_key) / float(self.window)

    def window_occupancy(self, motif_key: str) -> int:
        """Number of records currently in the rolling window."""
        return len(self._records.get(motif_key, []))


@dataclass
class MotifLedgerEntry:
    """Per-motif spawn-gate ledger entry (design doc §4.4; spec 2026-06-10
    §4.3, AMENDED).

    Lives in ``reaction_model.polymer_motif_ledger`` — in-memory ONLY (an RMG
    restart resets windows and deferred motifs re-earn their bar:
    correct-but-loud, same philosophy as unstamped-reaction demotion).
    Lookup is by Group isomorphism (:func:`_ledger_lookup`), never a canonical
    string key. ``representatives`` are ``(species_label, parent_pool_label)``
    pairs — the parent pool recorded at absorption (for Phase-D deferred
    candidates: the pool that would parent the spawn intent, currently
    ``pool_registry[0]``; if multi-parent attribution ever lands, this
    follows it). E[n]/monomer_MW for a representative are read LIVE from that
    pool's stats in the snapshot — live freshness kept, the wrong-mapping
    hole closed. ``accumulator_key`` is the opaque per-entry id used with
    :class:`MassFluxAccumulator`.
    """

    motif: Group
    accumulator_key: str
    representatives: List[Tuple[str, str]] = field(default_factory=list)
    last_recorded_iteration: int = -1
    spawned: bool = False


def _ledger_lookup(
    ledger: List[MotifLedgerEntry],
    motif: Group,
) -> Optional[MotifLedgerEntry]:
    """Find the entry whose motif Group is isomorphic to ``motif``.

    Same matching idiom as :func:`similarity_merge` — sidesteps Group
    canonicalization. The ledger is ``max_pools``-scale, so the O(n)
    isomorphism scan is fine.
    """
    for entry in ledger:
        try:
            if motif.is_isomorphic(entry.motif):
                return entry
        except (NotImplementedError, AttributeError, ValueError):
            continue
    return None


def _snapshot_event_mass(
    snapshot: Tuple[Dict[str, float],
                    Dict[str, Tuple[float, float, float]], float],
    species_label: str,
    parent_pool_label: str,
) -> float:
    """g_i for one representative per spec §3 (amended twice):
    ``gross[label] * max(0, E[n]_parent * monomer_MW_parent -
    chain_mass_defect_parent)``, with the (E[n], MW, defect) triple read
    from the snapshot's ``pool_stats`` for the parent pool RECORDED AT
    ABSORPTION. The defect term is the FR1-K2 mass-consumer audit
    (round-72 P2): an X-loss feature pool's chains each lost exactly one
    X, so the EXACT per-chain mass is the normative
    ``condensed_mass_g(mu0=1, mu1=E[n]) = E[n]*MW - defect`` -- the raw
    E[n]*MW would overstate every event by M_X. Ordinary pools carry
    defect 0 and reduce to the legacy product bit-identically. Labels
    absent from ``gross`` (absorbed this iteration, not yet simulated)
    and parent pools absent from ``pool_stats`` contribute 0 — stated,
    not incidental; the max(0, .) clamp and both absences err toward
    deferral."""
    gross, pool_stats, _ = snapshot
    e_n, mw, defect = pool_stats.get(parent_pool_label, (0.0, 0.0, 0.0))
    return gross.get(species_label, 0.0) * max(0.0, e_n * mw - defect)


def _spawn_gate_fraction(
    entry: MotifLedgerEntry,
    ledger: List[MotifLedgerEntry],
    snapshot: Optional[Tuple[Dict[str, float], Dict[str, Tuple[float, float]], float]],
) -> float:
    """fraction(motif) per spec §3 (AMENDED), from a stashed engine snapshot.

    ``snapshot`` is the 3-tuple from
    ``HybridPolymerSystem.spawn_gate_flux_snapshot()``:
    ``(gross, pool_stats, proxy_event_mass_total)``.

        numerator   = sum of g_i over THIS entry's (label, parent_pool) pairs
        denominator = proxy_event_mass_total
                      + sum of g_i over DEDUPED representatives across the
                        WHOLE ledger (a species in multiple motif entries
                        counts ONCE here, keyed by species label)

    Numerators are subsets of denominator terms, so each motif's fraction is
    in [0,1]; the SUM across motifs may exceed 1 (the stated multi-motif
    double-counting decision — competing for different pool slots). No
    snapshot (iteration 0, or no polymer reaction system) or an empty
    denominator -> 0.0: honest degradation, the gate defers; no production
    code path fakes a number.
    """
    if not snapshot:
        return 0.0
    try:
        _, _, proxy_event_mass_total = snapshot
    except (TypeError, ValueError):
        return 0.0
    numerator = sum(
        _snapshot_event_mass(snapshot, lbl, pool_lbl)
        for (lbl, pool_lbl) in entry.representatives
    )
    # First-seen-pool-wins is a DELIBERATE simplification, not an accident:
    # per spec §3, a species is absorbed into exactly ONE pool and carries
    # that pool's E[n] calibration everywhere it appears, so all of a
    # label's ledger records share one parent_pool_label today and
    # setdefault is exact. If multi-parent attribution ever lands, this
    # dedup must pick per-pool terms consistently with the numerator or the
    # fraction's [0,1] pin breaks (a numerator g_i(P2) could exceed its
    # denominator counterpart g_i(P1)).
    deduped: Dict[str, str] = {}
    for e in ledger:
        for (lbl, pool_lbl) in e.representatives:
            deduped.setdefault(lbl, pool_lbl)
    representative_total = sum(
        _snapshot_event_mass(snapshot, lbl, pool_lbl)
        for lbl, pool_lbl in deduped.items()
    )
    denominator = float(proxy_event_mass_total) + representative_total
    if denominator <= 0.0:
        return 0.0
    return numerator / denominator


def _evaluate_spawn_gate(
    cand: 'Species',
    motif: Group,
    reaction_model: Any,
    iteration: int,
    flux_accumulator: Optional[MassFluxAccumulator],
    mass_flux_threshold: float,
) -> Tuple[bool, float, Optional[MotifLedgerEntry]]:
    """Mass-flux spawn gate (design doc §4.4; spec 2026-06-10 §4.4).

    Arrival-driven: runs only when a NEW candidate carrying ``motif`` reaches
    Phase D (deferred candidates are never re-presented). Records at most ONE
    snapshot-attributed gross mass-flux fraction per motif per RMG iteration
    — same-iteration burst arrivals re-check the existing window only,
    otherwise "windowed over N iterations" silently becomes "windowed over N
    arrivals". The gate statistic is the window sum divided by the FIXED
    window length (zero-filled): a single-snapshot spike must be window x the
    bar to clear it.

    Returns ``(spawn, statistic, entry)``. ``entry`` is ``None`` when there
    is no reaction model / accumulator to hold gate state (bare unit-test
    path) — then the candidate defers with statistic 0.0.
    """
    cand_label = getattr(cand, "label", "") or repr(cand)
    if reaction_model is None or flux_accumulator is None:
        logging.info(
            "Polymer spawn gate: no reaction-model gate state available; "
            "deferring spawn for %s (statistic 0.0 < bar %.4g).",
            cand_label, mass_flux_threshold,
        )
        return False, 0.0, None

    ledger = getattr(reaction_model, "polymer_motif_ledger", None)
    if ledger is None:
        ledger = []
        reaction_model.polymer_motif_ledger = ledger

    entry = _ledger_lookup(ledger, motif)
    if entry is None:
        entry = MotifLedgerEntry(motif=motif, accumulator_key=f"motif-{len(ledger)}")
        ledger.append(entry)

    if entry.spawned:
        # Should be unreachable: Phase A classifies arrivals against the
        # spawned pool first. Assert-log; never re-run the gate.
        logging.warning(
            "Polymer spawn gate: arrival %s hit already-spawned ledger entry "
            "%s — Phase A should have classified it against the spawned pool.",
            cand_label, entry.accumulator_key,
        )
        return False, 0.0, entry

    snapshot = getattr(reaction_model, "polymer_flux_snapshot", None)
    if entry.last_recorded_iteration < iteration:
        fraction = _spawn_gate_fraction(entry, ledger, snapshot)
        flux_accumulator.record(entry.accumulator_key, fraction, iteration)
        entry.last_recorded_iteration = iteration

    statistic = flux_accumulator.gate_statistic(entry.accumulator_key)
    spawn = statistic >= mass_flux_threshold
    if not spawn:
        logging.info(
            "Polymer spawn gate: deferring spawn for %s — statistic %.4g < "
            "bar %.4g (window %d/%d records%s).",
            cand_label, statistic, mass_flux_threshold,
            flux_accumulator.window_occupancy(entry.accumulator_key),
            flux_accumulator.window,
            "" if snapshot is not None else "; no flux snapshot stashed",
        )
    return spawn, statistic, entry


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
    # The shared accumulator lives on the reaction model (alongside the
    # motif ledger and the stashed snapshot, spec 2026-06-10 §4.3); an
    # explicitly-passed accumulator (test injection) wins.
    if flux_accumulator is None and reaction_model is not None:
        flux_accumulator = getattr(reaction_model, "polymer_flux_accumulator", None)

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

        # Phase D: gates — mass-flux spawn gate (design doc §4.4, spec
        # 2026-06-10) and max_pools cap.
        spawn_ok, gate_statistic, ledger_entry = _evaluate_spawn_gate(
            cand=cand,
            motif=motif,
            reaction_model=reaction_model,
            iteration=iteration,
            flux_accumulator=flux_accumulator,
            mass_flux_threshold=mass_flux_threshold,
        )
        # Either way, the arriving candidate becomes a representative of its
        # motif (spec §4.4 step 6): it is absorbed as a proxy variant and is
        # the handle a future snapshot attributes flux to. Recorded as a
        # (species_label, parent_pool_label) pair — the parent pool that
        # would parent its spawn intent (currently pool_registry[0]; if
        # multi-parent attribution ever lands, this follows it), per the
        # spec-§3 attribution rules.
        cand_label = getattr(cand, "label", "") or ""
        parent_pool_label = pool_registry[0].label if pool_registry else ""
        if (ledger_entry is not None and cand_label
                and all(lbl != cand_label for (lbl, _) in ledger_entry.representatives)):
            ledger_entry.representatives.append((cand_label, parent_pool_label))
        if not spawn_ok:
            _tag_polymer_proxy(cand, is_proxy=True)
            processed.append(cand)
            continue
        if len(pool_registry) >= max_pools:
            _warn_pool_cap_exhausted(cand_label, len(pool_registry), max_pools)
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
                mass_flux_at_spawn=gate_statistic,
            )
        )
        if ledger_entry is not None:
            ledger_entry.spawned = True
        _tag_polymer_proxy(cand, is_proxy=True)
        processed.append(cand)

    return processed, spawn_intents


POLYMER_POOLS_SIDECAR_SCHEMA_VERSION = "2.0"
# Schema 2.1 = 2.0 + the channels.radical_qssa_unzip vocabulary. Channel-
# vocabulary growth is a MINOR shape bump per the versioning policy
# (docs/polymer_moments_format.md), so the emitter stamps 2.1 exactly when
# at least one serialized pool carries the QSSA block — and keeps stamping
# 2.0 otherwise, so legacy artifacts stay byte-identical (pinned by test).
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_QSSA = "2.1"
# Rate-recipe revision marker, emitted as conventions["recipe_revision"].
# Bump rule: assign a NEW date whenever the RATE RECIPE changes — site
# scaling, the chip exhaustion throttle, the kb/Keq recipe, or channel /
# flux-archetype algebra — independent of
# POLYMER_POOLS_SIDECAR_SCHEMA_VERSION, which governs artifact SHAPE only.
# Downstream consumers (TA) hard-fail on unknown values, so bumping this is
# a consumer-coordination event, not a cosmetic edit.
# 2026-07-03-monomer-gas (incident 2026-07-03, design B-prime): the
# monomer_routing target is a GAS-phase species -- it is no longer listed in
# pools[].phase_species or conventions.condensed_species, and the unzip/QSSA
# release deposits into the gas amount basis. Consumers that treated the
# routed monomer as condensed (V_rxn selection, headspace/V_gas inventory,
# phase reporting) MUST be updated; the revision bump forces that
# coordination (docs/polymer_moments_format.md revision note).
POLYMER_RATE_RECIPE_REVISION = "2026-07-03-monomer-gas"
# The radical_qssa_unzip channel is new channel/flux algebra, so a QSSA
# artifact carries a NEW recipe revision. Conditional for the same reason as
# the 2.1 schema stamp: no QSSA anywhere -> legacy revision, byte-identical.
POLYMER_RATE_RECIPE_REVISION_QSSA = "2026-07-03-qssa-monomer-gas"
# Schema 2.2 = 2.1 + the weak-link allyl/U-state vocabulary INSIDE the
# radical_qssa_unzip block (initiation_allyl / termination_recombination /
# termination_disproportionation / unsaturated_tail_ends_initial). Same
# minor-bump policy: the emitter stamps 2.2 exactly when at least one
# serialized pool carries the weak-link vocabulary (the artifact-level stamp
# is governed by the STRONGEST vocabulary present, so a mixed artifact is
# 2.2); no weak-link anywhere -> the 2.1/2.0 stamps apply byte-identically.
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_WEAKLINK = "2.2"
# The weak-link U-state channel is new rate algebra (allyl initiation with
# nu=1, split termination, dU/dt law), so a weak-link artifact carries a NEW
# recipe revision — conditional for the same reason as the schema stamps.
POLYMER_RATE_RECIPE_REVISION_WEAKLINK = "2026-07-03-weaklink-u-monomer-gas"
# Schema 2.3 = 2.2 + the POOL-LEVEL explicit_dp block (stage B of the
# explicit-DP arc): the capped DP=cutoff handshake oligomer generated by the
# deck flag explicit_dp=True (rmgpy/rmg/input.py polymer() step 4c) is
# serialized as pool state/topology — deliberately NOT inside channels
# (it is not kinetics; the boundary-flux law it participates in is carried
# by the block's normative recipe). Same minor-bump policy: the emitter
# stamps 2.3 exactly when at least one serialized pool carries the block;
# no explicit-DP anywhere -> the 2.2/2.1/2.0 stamps apply byte-identically
# (golden-pinned by test). STRICT-MINOR acceptance means a 2.2 consumer
# hard-rejects a 2.3 artifact instead of silently dropping the handshake
# target (whose flux drain would then be unaccounted).
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_EXPLICIT_DP = "2.3"
# The explicit-DP boundary flux (gamma-conditional handshake at DP == xs,
# polymer.pyx:3475-3524) is rate algebra that was previously inert in every
# emitted artifact (run-path pools carried no explicit map), so an
# explicit-DP artifact carries a NEW recipe revision — one token per channel
# vocabulary, exactly the 3-way composition rule the 2026-07-03-monomer-gas
# family established (base / qssa / weaklink-u each get a re-stamped token;
# the artifact-level token is picked by the strongest channel vocabulary
# present). All three implement the monomer-gas contract (suffix retained).
POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP = "2026-07-04-explicit-dp-monomer-gas"
POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP_QSSA = (
    "2026-07-04-explicit-dp-qssa-monomer-gas")
POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP_WEAKLINK = (
    "2026-07-04-explicit-dp-weaklink-u-monomer-gas")
# The pool block's own revision token (block-local, like the QSSA recipe
# sub-block): names the explicit-DP recipe itself, independent of which
# channel family the artifact-level token composes with.
EXPLICIT_DP_BLOCK_RECIPE_REVISION = "2026-07-04-explicit-dp"
# Schema 2.4 = 2.3 + the per-row refused marker on reactions[] entries
# (refused-row sidecar contract): a reaction stamped polymer_refused=True
# (item 18 stamp-but-keep) is zeroed WHOLESALE by the generating solver
# (polymer.pyx reaction_refused), so its row carries "refused": true +
# "refused_reason" (the census reason available at stamp time:
# conduit-deferred / qssa-invalid). The row STAYS listed — consumers need it
# to zero the mapped Cantera reaction — but MUST skip its moment flux AND
# its manual species dispatch (docs/polymer_moments_format.md §12). Same
# presence-based minor-bump policy as 2.1/2.2/2.3: the emitter stamps 2.4
# exactly when at least one row carries the marker; no refused row anywhere
# -> the 2.3/2.2/2.1/2.0 stamps apply byte-identically (non-refused rows
# never gain the keys: absent, not false). NO recipe_revision change:
# refused is SHAPE vocabulary with consumption semantics, not new rate
# algebra, and STRICT-MINOR acceptance already stops older consumers at the
# envelope. Erratum: pre-2.4 sidecars generated while refused stamps
# existed carry such rows UNMARKED and consumers over-integrate them.
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_REFUSED = "2.4"
# Schema 2.5 = 2.4 + the spawned-pool closure on the conventions block
# (S4 serializer closure): conventions.spawned_pools lists every registry
# pool that is NOT solver-configured (runtime-spawned scission daughters and
# S2 feature pools <parent>_mod), disjoint from configured_pools by
# construction, and conventions.condensed_species is CLOSED over those
# pools' phase_species (proxy + mu-dummies, already declared condensed
# row-side) so the TA consumer classifies them CONDENSED instead of
# defaulting them GAS (the item-16 mass-balance hazard). Same presence-based
# minor-bump policy as 2.1-2.4: the emitter stamps 2.5 exactly when at
# least one spawned pool is present in the registry; no spawned pool
# anywhere -> the key is ABSENT (never an empty list) and the 2.4/2.3/...
# stamps apply byte-identically (golden-pinned). NO recipe_revision change:
# spawned_pools is SHAPE vocabulary with classification semantics, not new
# rate algebra (spawned pools still carry no solver-config channels
# integration), and STRICT-MINOR acceptance already stops older consumers
# at the envelope — the exact 2.4 refused-row precedent.
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_SPAWNED = "2.5"
# Schema 2.6 = 2.5 + the POOL-LEVEL homolysis_initiation block (Stage 2 of
# the radical-homolysis initiation arc, adjudicated rounds 66/67): a pool
# with the k_homolysis kernel configured serializes its Arrhenius triplet
# (structured, explicit SI units), the two end-radical daughter pool labels
# under the open-*1/open-*2 POSITIONAL field names (round-67 ruling (a):
# never a primary/secondary radical-character claim), the solver kernel name
# (ruling (c): machine-checkable supersession -- refused explicit homolysis
# rows stay refused/zero-flux, consumers must never infer the kernel's flux
# from them), a block-local recipe_revision, and the machine-pinned moment
# law in the STABLE product forms actually implemented (round-67 P2). Same
# presence-based minor-bump policy as 2.1-2.5: the emitter stamps 2.6
# exactly when at least one serialized pool carries the block; no homolysis
# pool anywhere -> the older stamps apply byte-identically (pinned by test).
# 2.6 sits BESIDE 2.5, never subsumes it (round-67 §Stage 2 Scope): the
# spawned-pool closure + condensed closure keep their exact 2.5 semantics,
# and both daughter pools must appear in pools[], be classified in the
# configured/spawned closure, and be condensed per the condensed closure
# (the producer refuses to emit any other shape; the consumer hard-rejects
# it). NO artifact-level recipe_revision change: the block carries its own
# block-local revision token (the explicit-dp 2.3 precedent), and
# STRICT-MINOR acceptance already stops pre-2.6 consumers at the envelope.
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_HOMOLYSIS = "2.6"
# Schema 2.7 = 2.6 + the POOL-LEVEL side_group_homolysis block (FR1-K2,
# adjudicated rounds 72/73 -- the sibling of the 2.6 homolysis_initiation
# surface): a pool with the side-group C-X homolysis kernel configured
# serializes its FULL channel list (label, structured SI kinetics, the
# round-72 REQUIRED site_selector from the CLOSED vocabulary
# {'aryl','benzylic','aliphatic'}, sites_per_unit PLUS the selector's
# resolved site_atom_indices in monomer_adj_list atom order -- so the
# loader validates selector/feature closure from serialized data alone,
# never re-deriving from a monomer graph it does not have (round-73: K2
# must not inherit the solver-backstop structural gap), the gas_product
# SMILES with its routed artifact species label + pinned M_X, and the
# ratified X-loss feature pool label), the solver kernel name, a
# block-local recipe_revision, the machine-pinned moment law AND the
# NORMATIVE mass formula
#     condensed_mass_g = mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol
# (round-70 #1 P1 trap: the feature pool keeps the parent's monomer_mw, so
# the one-lost-X-per-chain defect must be explicit). Each X-loss FEATURE
# pool entry carries chain_mass_defect_g_mol (a POOL property, = M_X of
# the spawning channel's gas_product) -- vocabulary schema 2.6 cannot
# express. Same presence-based minor-bump policy as 2.1-2.6: the emitter
# stamps 2.7 exactly when at least one serialized pool carries the block
# or the defect field; a kernel-free artifact (including a 2.6-only
# homolysis artifact) keeps its older stamp byte-identically (pinned).
# 2.7 sits BESIDE 2.5/2.6, never subsumes them. NO artifact-level
# recipe_revision change (the 2.3/2.6 block-local precedent).
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_SIDE_GROUP = "2.7"
# Schema 2.8 = 2.7 + the POOL-LEVEL end_radical_depropagation block (r74
# SS2 kernel, r78 serialization rulings -- the real serialization that
# replaced the 'schema 2.8 pending' producer hard-refusal): each spawned
# end-radical DAUGHTER pool carrying the k_depropagation kernel serializes
# its full Arrhenius triplet (structured kinetics with pinned SI units),
# the released gas monomer identity + MW (cross-pinned to the pool
# surface's monomer_routing / monomer_mw_g_mol -- each unzip event moves
# exactly ONE repeat unit from the condensed basis into ONE mole of
# gas_species), the machine-readable smooth-gate width (bitwise ==
# KDEP_GATE_WIDTH, the 1e-12 TA-twin contract), the solver kernel name, a
# block-local recipe_revision, and the machine-pinned AS-IMPLEMENTED
# moment law (gated smooth-pos mu2 under-drain disclosed; UN-gated dmu0
# half-bin N1 gamma closure with its terminal smoothstep floor). The
# PARENT pool's k_depropagation is deck DECLARATION surface only (the
# kernel never integrates there: validate_configuration excludes
# k_depropagation + k_homolysis on one pool and parent configs never
# carry the triplet), so the parent entry emits NO block. Same
# presence-based minor-bump policy as 2.1-2.7: emitter stamps 2.8 exactly
# when at least one serialized pool carries the block; a kernel-free
# artifact keeps its older stamp byte-identically (pinned). 2.8 sits
# BESIDE 2.5/2.6/2.7, never subsumes them. NO artifact-level
# recipe_revision change (the 2.3/2.6/2.7 block-local precedent).
POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_DEPROPAGATION = "2.8"
POLYMER_POOLS_SIDECAR_FILENAME = "polymer_pools.json"


def collect_polymer_pool_registry(*species_lists) -> List['Polymer']:
    """Build the sidecar/spawn-pass pool registry from species lists
    (typically ``core.species``, ``edge.species``, ``new_species_list``):
    every :class:`Polymer`, identity-deduped, first occurrence wins,
    order preserved.

    Identity dedup is load-bearing: a freshly-promoted daughter Polymer sits
    in BOTH ``core.species`` and ``new_species_list`` until the next enlarge
    clears the latter, so a plain concatenation serializes the same pool
    object twice (observed live: pools [PS, tail, tail_2, tail, tail_2]).
    Dedup is by object identity only — two DISTINCT Polymer objects are both
    kept even if they compare equal."""
    seen = set()
    registry = []
    for species_list in species_lists:
        for spc in species_list:
            if isinstance(spc, Polymer) and id(spc) not in seen:
                seen.add(id(spc))
                registry.append(spc)
    return registry


def _artifact_species_label(spc) -> str:
    """chem.yaml species name for ``spc`` — must match rmgpy.cantera.get_label:
    ``label(index)`` when index > 0, bare label otherwise (µ-dummies have
    index = -1 and appear bare in chem.yaml)."""
    index = getattr(spc, "index", -1)
    label = getattr(spc, "label", "")
    return f"{label}({index})" if index > 0 else label


def _species_base_label(spc) -> str:
    """Strip the RMG ``(N)`` index suffix — the solver's pool-membership rule
    (``_apply_pool_phase_overrides`` in rmgpy/solver/polymer.pyx; ONE
    convention: :func:`strip_rmg_index_suffix`)."""
    return strip_rmg_index_suffix(getattr(spc, "label", ""))


def derive_condensed_species(core_species, pools_cfg, mask=None):
    """Resolve the condensed-phase subset of ``core_species`` for the artifact.

    ``conventions.condensed_species`` (docs/polymer_moments_format.md §4 step 2,
    §8) MUST mirror the live solver's final-core phase membership (the oracle).
    The solver only marks a species condensed for pools it has CONFIGURED
    (``HybridPolymerSystem.polymer_pools``); a daughter pool spawned mid-run but
    never solver-configured (e.g. ``epdm_scission_tail``) runs its proxies as
    ordinary GAS species. So ``pools_cfg`` must be the CONFIGURED pools — the
    caller (the save_everything hook) passes the live engine's
    ``polymer_pools``, falling back to the full registry only when no engine is
    resolvable (direct/test invocation).

    Two complementary sources are UNIONED:

    1. The solver's ``gas_species_mask`` (True=gas, False=condensed), honored
       VERBATIM when it matches ``len(core_species)`` — the authoritative oracle
       phase verdict (rmgpy/solver/polymer.pyx:478-516). When the mask is absent
       or length-mismatched (e.g. a constructor-era mask sized to a smaller
       core, or the blueprint surfacing ``None``), only source 2 is used.
    2. Label/index membership derived from ``pools_cfg`` the same way the solver
       builds its mask: a core species is condensed iff its base label (RMG
       ``(N)`` suffix stripped) matches a pool's proxy label or
       ``{label}_mu0/_mu1/_mu2``, or its index is a pool's explicit-oligomer /
       moment index. The routed-monomer index is NOT condensed (recipe
       revision 2026-07-03-monomer-gas: the release target is a gas
       volatile). Keyed on the CONFIGURED pools, this
       reproduces exactly what the solver mask would have marked — so when the
       mask IS present they agree, and when it is missing this is a faithful
       reconstruction (NOT an over-broad registry union).

    Parameters
    ----------
    core_species : list of Species
        The FINAL core species list.
    pools_cfg : iterable, optional
        The solver-CONFIGURED pool objects, exposing ``label`` (and optionally
        ``mu_indices``, ``explicit_dp_to_species_index``, ``monomer_poly_index``).
        Accepts both solver ``PolymerPoolConfig`` and ``Polymer`` registry
        entries. Do NOT pass spawned-but-unconfigured daughter pools here: the
        oracle does not phase-resolve them as condensed.
    mask : sequence of bool, optional
        The solver's ``gas_species_mask``. Honored verbatim when
        ``len(mask) == len(core_species)``.

    Returns
    -------
    list of Species
        The condensed-phase core species, in core order.
    """
    n = len(core_species)
    condensed_idx = set()

    # Source 1: authoritative final-core mask (when length-matched).
    if mask is not None and len(mask) == n:
        condensed_idx.update(i for i in range(n) if not mask[i])

    # Source 2: pool membership — mirror polymer.pyx:478-516. Catches pools
    # spawned mid-run that the solver mask never marked condensed.
    for pool in (pools_cfg or []):
        label = getattr(pool, "label", None)
        if label:
            member_bases = {label, f"{label}_mu0", f"{label}_mu1", f"{label}_mu2"}
            for i in range(n):
                if _species_base_label(core_species[i]) in member_bases:
                    condensed_idx.add(i)
        explicit_map = getattr(pool, "explicit_dp_to_species_index", None) or {}
        for idx in explicit_map.values():
            if isinstance(idx, int) and 0 <= idx < n:
                condensed_idx.add(idx)
        # mu_indices may be a (mu0,mu1,mu2) tuple of core indices (solver config
        # or Polymer.mu_indices) — a dict form carries no core indices, skip it.
        mu_indices = getattr(pool, "mu_indices", None)
        if isinstance(mu_indices, (tuple, list)):
            for idx in mu_indices:
                if isinstance(idx, int) and 0 <= idx < n:
                    condensed_idx.add(idx)
        # monomer_poly_index (the unzip/QSSA release target) is deliberately
        # NOT added: since recipe revision 2026-07-03-monomer-gas (incident
        # 2026-07-03, design B-prime) the routed monomer is a GAS-phase
        # species -- the solver validates it gas and the release deposits
        # into the gas amount basis. Mirroring the oracle means leaving it
        # out here too.

    return [core_species[i] for i in range(n) if i in condensed_idx]


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


# Pinned sidecar units convention for the radical_qssa_unzip channel
# (single source for the emitter; the artifact loader in
# rmgpy/tools/polymer_moments_runner.py pins the same strings independently
# on its side of the boundary). initiation/depropagation are unimolecular
# [s^-1]; termination is the ONLY bimolecular block [m^3 mol^-1 s^-1];
# transfer is PSEUDO-first-order [s^-1] -- ktr multiplies the active-end
# concentration R directly in the M2 rate law, so a literature bimolecular
# k_tr must be premultiplied by the substrate concentration [mol/m^3]
# BEFORE it enters the config. Ea always [J/mol]. A consumer seeing a
# different string must ERROR, never convert.
RADICAL_QSSA_SIDECAR_A_UNITS = {
    "initiation": "s^-1",
    "depropagation": "s^-1",
    "termination": "m^3/(mol*s)",
    "transfer": "s^-1",
    # Weak-link vocabulary (schema 2.2): allylic initiation is unimolecular
    # in the active unsaturated end; the split termination triplets carry
    # the bimolecular units of the summed block they replace.
    "initiation_allyl": "s^-1",
    "termination_recombination": "m^3/(mol*s)",
    "termination_disproportionation": "m^3/(mol*s)",
}

# Pinned units note for the serialized U-state initial condition (schema
# 2.2). U is STATE, not a rate constant: same amount basis as mu0 [mol]; the
# consumer divides by V_poly for the concentration form. The loader pins the
# same string byte-for-byte (reject, never convert).
RADICAL_QSSA_SIDECAR_U0_UNITS = (
    "mol — tail-distribution state; consumer divides by V_poly")

# Normative machine-readable QSSA rate recipe, emitted as the channel block's
# ``recipe`` sub-block (schema 2.1). Unlike the human-readable ``provenance``
# strings, these are CONTRACT: the artifact loader
# (rmgpy/tools/polymer_moments_runner.py _QSSA_PINNED_RECIPE) pins the same
# values independently and errors on any mismatch or absence — reject, never
# adapt, exactly like the units pin above. Every entry is verified against
# the implemented RHS in rmgpy/solver/polymer.pyx (M2 rate law: B_qssa /
# R_ss / r_qssa in the pool loop; moment signature dmu1/dmu2 drains;
# SMALL_EPS = 1e-30 at polymer.pyx:71) — editing the solver algebra without
# re-verifying and re-dating these strings breaks the contract.
RADICAL_QSSA_SIDECAR_RECIPE = {
    "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                   "(mol/m^3 condensed)"),
    "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                         "sqrt(efficiency * ki * B / kt)"),
    "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                           "(sqrt(ktr^2 + 8*kt*(2*efficiency*ki*B)) "
                           "- ktr) / (4*kt)"),
    "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                         "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
    "small_eps": 1e-30,
    "volume_note": ("kt is bimolecular: rates depend on condensed "
                    "volume V_poly; consumers MUST evaluate on "
                    "concentration moments mu_k = n_k / V_poly and "
                    "convert emitted rate back with *V_poly"),
}

# Normative machine-readable recipe for the WEAK-LINK allyl/U-state variant
# of the QSSA channel (schema 2.2). Same contract as
# RADICAL_QSSA_SIDECAR_RECIPE: the loader (polymer_moments_runner
# _QSSA_PINNED_RECIPE_WEAKLINK) pins the same values independently and
# errors on mismatch or absence — reject, never adapt. Every entry is
# transcribed from the implemented weak-link RHS in
# rmgpy/solver/polymer.pyx (weak-link branch of the pool loop; U slot
# layout in _flatten_radical_qssa_state; U0 census in
# set_initial_conditions) — editing means re-verifying and re-dating.
RADICAL_QSSA_SIDECAR_RECIPE_WEAKLINK = {
    "bond_basis": ("B = max(mu1 - mu0, 0) on concentration moments "
                   "(mol/m^3 condensed)"),
    "channel_gate": ("channel gates on B > 0: at B = 0 it is INERT even at "
                     "U > 0 (no monomer drain, no U production, no U sink; "
                     "the DP->1 self-termination floor survives)"),
    "u_state": ("U is a per-pool amount state [mol], tail-distribution "
                "basis, MASSLESS (it never enters condensed mass or any "
                "Mn/Mw/PDI consumer); daughter pools spawn with U0 = 0 "
                "(constants inherit, state resets)"),
    "u_active": ("u_active = min(max(U, 0)/V_poly, B) (active-site clamp, "
                 "mol/m^3)"),
    "radical_generation": ("G_R = 2*efficiency*ki*B + "
                           "1*efficiency*ki_allyl*u_active (nu = 1: one "
                           "unzipping radical per weak-link fission; the "
                           "allylic co-fragment does not unzip)"),
    "kt_split": ("kt_total = kt_rec + kt_disp replaces the legacy summed kt "
                 "everywhere (halved radical-disappearance convention "
                 "unchanged)"),
    "rate_no_transfer": ("r_mono = monomer_yield * kdp * "
                         "sqrt(G_R / (2*kt_total))"),
    "rate_with_transfer": ("r_mono = monomer_yield * kdp * "
                           "(sqrt(ktr^2 + 8*kt_total*G_R) - ktr) / "
                           "(4*kt_total)"),
    "du_dt": ("dU/dt = kt_disp*R_ss^2*max(0, 1 - U/(2*mu0))*V_poly - "
              "efficiency*ki_allyl*u_active*V_poly; production is the "
              "disproportionation EVENT rate with NO efficiency factor "
              "(R_ss already carries the escape efficiency) under a linear "
              "capacity throttle exactly zero at the TAIL chain-end "
              "capacity 2*mu0 [mol], where mu0 is the INSTANTANEOUS "
              "tail-distribution mu0(t) amount read from the current "
              "state at each RHS evaluation (NOT the frozen initial mu0 "
              "of the u0_census bound); at mu0(t) <= 0 the throttle is "
              "EXACTLY 0 by explicit branch (no ends, no capacity, no U "
              "production; no division is performed and no small_eps "
              "floor is applied); the sink is efficiency-SYMMETRIC "
              "with the G_R allyl term (a caged recombination restores "
              "the allylic bond)"),
    "u0_census": ("U0 = unsaturated_tail_ends_initial [mol] is rejected at "
                  "solver init when U0 > 2*mu0_tail evaluated on the "
                  "INITIAL (t = 0) tail mu0 (TAIL-only chain-end capacity, "
                  "amount basis); the loader passes the value through"),
    "u_transport": ("U is NEVER advected or transferred by any inter-pool "
                    "or ejection flux archetype (migration, "
                    "scission_fragment, discrete_chip and "
                    "volatile_ejection move pool moments and species "
                    "amounts only); the ONLY writers of U are this "
                    "recipe's du_dt law and the t = 0 initial condition, "
                    "and daughter pools spawn with U0 = 0 (constants "
                    "inherit, state resets)"),
    "moment_signature": ("dmu0 = 0; dmu1 -= r_mono; dmu2 -= r_mono * "
                         "max(2*mu1/max(mu0, small_eps) - 1, 0)"),
    "small_eps": 1e-30,
    "volume_note": ("kt is bimolecular: rates depend on condensed "
                    "volume V_poly; consumers MUST evaluate on "
                    "concentration moments mu_k = n_k / V_poly and "
                    "convert emitted rate back with *V_poly"),
}

# Normative machine-readable explicit-DP handshake recipe, emitted as the
# pool-level ``explicit_dp`` block's ``recipe`` sub-block (schema 2.3). Same
# contract as RADICAL_QSSA_SIDECAR_RECIPE: consumers pin these strings
# independently and reject on mismatch or absence — never adapt. Every entry
# is transcribed from the implemented oracle:
#   * boundary_flux — the hybrid handshake in rmgpy/solver/polymer.pyx
#     (:3475-3524 at the stage-B commit): gamma-conditional boundary
#     population with the triangular fallback and the mu0/mu1/mu2 flux
#     clamps; gamma moment-matching in _gamma_params_from_mu012 /
#     _gamma_prob_conditional_hybrid (polymer.pyx:491-525).
#   * k_chain — the per-chain handshake frequency arm selection
#     (polymer.pyx:3484-3488; k_unzip and QSSA are mutually exclusive M1).
#   * tail_split — the t=0 accounting. NORMATIVE RULE: declared
#     pools[].moments are TOTAL-INCLUSIVE (explicit chains counted in). The
#     solver seeds the explicit species' initial_moles as species amounts
#     (set_initial_conditions step 2, clamped >= 0), seeds the declared
#     moments (step 3), and then step 6 — the Moment Consistency Check
#     (polymer.pyx:2109-2132) — subtracts each mapped DP's contribution,
#     _explicit_moment_contributions (polymer.pyx:471-488): (N_dp, dp*N_dp,
#     dp^2*N_dp) on concentration moments, from mu0/mu1/mu2, clamped >= 0
#     with a clamp WARNING. The integrated tail moments are therefore
#     total - explicit (behaviorally pinned: declared [1,5,30] with 0.25 mol
#     at DP=3 seeds mu = [0.75, 4.25, 27.75]). The generation-side V_poly
#     mass split (PolymerPhase.calculate_volume, rmgpy/rmg/polymer_input.py)
#     applies the same subtraction and HARD-ERRORS when explicit mu1 exceeds
#     the declared mu1 beyond -1e-12.
#   * transport — one-way by construction (no reverse flux term exists).
# Editing the solver algebra without re-verifying and re-dating these
# strings (and EXPLICIT_DP_BLOCK_RECIPE_REVISION) breaks the contract.
EXPLICIT_DP_SIDECAR_RECIPE = {
    "tail_split": ("declared pools[].moments are TOTAL-INCLUSIVE (explicit "
                   "chains counted in); the solver seeds the explicit "
                   "species' initial_moles as species amounts (clamped "
                   ">= 0), then subtracts each mapped DP's contribution "
                   "(N_dp, dp*N_dp, dp^2*N_dp on concentration moments) "
                   "from the seeded mu0/mu1/mu2, clamped >= 0 with a clamp "
                   "warning (set_initial_conditions step 6) -- the "
                   "integrated tail moments are total - explicit; the "
                   "generation-side V_poly mass split applies the same "
                   "subtraction and hard-errors when explicit mu1 exceeds "
                   "declared mu1 beyond -1e-12"),
    "boundary_flux": ("gated on mu0 > 1e-9 mol/m^3 AND mu1/mu0 > xs + 1e-9; "
                      "p_cond = P(DP = xs+1 | DP > xs) from the gamma "
                      "distribution moment-matched to (mu0, mu1, mu2) "
                      "(k = 1/(PDI - 1), theta = mean/k; half-integer bins: "
                      "[F((xs+1.5)/theta) - F((xs+0.5)/theta)] / "
                      "[1 - F((xs+0.5)/theta)] with F the regularized lower "
                      "incomplete gamma); triangular fallback on tail_mean "
                      "in (xs+1, xs+2) peaking 1.0 at xs+1.5 when the gamma "
                      "is unrealizable (any moment <= 1e-30, PDI <= 1+1e-6, "
                      "or non-finite params); p_cond clamped to [0, 1]; "
                      "N_boundary = min(mu0*p_cond, mu0, mu1/xs, mu2/xs^2) "
                      "[mol/m^3]; F_flux = k_chain * N_boundary; "
                      "dn(species[xs])/dt += F_flux*V_poly; dmu0 -= F_flux; "
                      "dmu1 -= xs*F_flux; dmu2 -= xs^2*F_flux"),
    "k_chain": ("k_unzip when k_unzip > 0, else r_qssa/max(mu0, 1e-30) when "
                "the radical_qssa_unzip channel is active (mutually "
                "exclusive upstream; at most one arm fires)"),
    "transport": ("one-way: statistical moment tail -> explicit real "
                  "condensed species at DP == handshake_target_dp; no "
                  "reverse flux"),
}


def _serialize_explicit_dp_block(pool: 'Polymer',
                                 core_species: Optional[List['Species']] = None,
                                 initial_explicit_moles: Optional[Dict[int, float]] = None
                                 ) -> Optional[Dict[str, Any]]:
    """Serialize the pool's explicit-DP handshake target (schema 2.3), or
    return None when the pool carries none (legacy pools emit nothing —
    flag-OFF artifacts stay byte-identical, golden-pinned).

    The block is POOL-LEVEL state/topology, not a channel: the species map
    and initial loading describe what the handshake deposits INTO; the rate
    algebra is carried by the normative ``recipe`` sub-block (pinned by
    consumers, reject-never-adapt).

    Hard errors (the producer never emits a silently-inert or
    universe-inconsistent shape):

    * ``explicit_dp=True`` with no attached species — would recreate the
      structurally inert handshake (same posture as compile_polymer_phase,
      rmgpy/rmg/polymer_input.py).
    * a core universe is supplied and the species is not in it BY IDENTITY —
      artifact labels must come from the same universe as every other
      species reference (same posture as _derive_qssa_monomer_routing).
    """
    dp_spc = getattr(pool, "explicit_dp_species", None)
    label = getattr(pool, "label", "")
    if dp_spc is None:
        if getattr(pool, "explicit_dp", False):
            raise ValueError(
                f"Pool '{label}': explicit_dp=True but no auto-generated "
                f"DP={getattr(pool, 'cutoff', '?')} oligomer species is "
                f"attached (explicit_dp_species is None) -- serializing the "
                f"flag without its species would emit a structurally inert "
                f"handshake. The species is created by the polymer() input "
                f"block when explicit_dp=True; refusing to serialize.")
        return None
    if core_species and not any(spc is dp_spc for spc in core_species):
        raise ValueError(
            f"Pool '{label}': explicit-DP species "
            f"('{getattr(dp_spc, 'label', '')}') is not in the core species "
            f"universe the artifact labels come from (identity check). "
            f"Refusing to serialize an explicit_dp block whose species the "
            f"artifact cannot name consistently.")
    dp = int(getattr(pool, "cutoff", 0))
    moles = float((initial_explicit_moles or {}).get(dp, 0.0))
    return {
        "enabled": True,
        "species": {str(dp): _artifact_species_label(dp_spc)},
        "initial_moles": {str(dp): moles},
        "handshake_target_dp": dp,
        "recipe_revision": EXPLICIT_DP_BLOCK_RECIPE_REVISION,
        "recipe": dict(EXPLICIT_DP_SIDECAR_RECIPE),
    }


# Block-local revision token of the pool-level ``homolysis_initiation``
# block (schema 2.6). Same contract as EXPLICIT_DP_BLOCK_RECIPE_REVISION:
# names the kernel's rate recipe itself, independent of the artifact-level
# token; consumers pin it by exact match and reject on mismatch or absence.
HOMOLYSIS_BLOCK_RECIPE_REVISION = "2026-07-05-radical-homolysis"
# The solver kernel this block supersedes refused explicit rows FOR (round
# 67 ruling (c)): versioned term-type style, same convention as
# ARCHETYPE_TERM_NAMES / mu3_closure. Consumers pin it by exact match --
# an unknown kernel is flux they cannot reproduce.
HOMOLYSIS_KERNEL_NAME = "radical_homolysis_initiation/1"
# Stage-1 spawn provenance pin stamped on every homolysis end-radical
# daughter (generate_end_radical_daughters). Consumers pin
# spawn_event_metadata.source by exact match (runner
# _HOMOLYSIS_SPAWN_SOURCE) and reject any other provenance; the producer
# closure guard mirrors that rejection (r68).
HOMOLYSIS_SPAWN_SOURCE = "k_homolysis_end_radical"
# FR1-K1 spawn provenance pin stamped on every side-group X-loss feature
# daughter (generate_side_loss_daughters), together with the spawning
# channel label: spawn_metadata = {"source": SIDE_GROUP_HOMOLYSIS_SPAWN_
# SOURCE, "channel": <label>}. The schema-2.7 consumer (K2, pending) will
# pin it exactly, mirroring the sibling kernel's HOMOLYSIS_SPAWN_SOURCE
# contract.
SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE = "side_group_homolysis"
# Normative moment law of the radical-homolysis initiation kernel, pinned
# in the STABLE product forms actually implemented (round-67 P2;
# rmgpy/solver/polymer.pyx, the khom_* RHS section). Same contract as
# EXPLICIT_DP_SIDECAR_RECIPE: consumers pin these strings independently and
# reject on mismatch -- editing the solver algebra without re-verifying and
# re-dating these strings (and HOMOLYSIS_BLOCK_RECIPE_REVISION) breaks the
# contract.
HOMOLYSIS_SIDECAR_RECIPE = {
    "event_rate": ("R = k(T)*max(mu1 - mu0, 0) [mol/(m^3 s)]; "
                   "k(T) = A*T^n*exp(-Ea/(R_gas*T)) evaluated at the "
                   "RUNTIME temperature (round 66: never a precomputed "
                   "scalar); a chain of length n has n-1 breakable bonds"),
    "parent_debit": ("dmu0 -= R; dmu1 -= R*B1 computed as k*(mu2 - mu1); "
                     "dmu2 -= R*B2 computed as k*(mu3 - mu2); mu3 from the "
                     "log_lagrange/1 closure"),
    "daughter_credit": ("EACH of the two end-radical daughter pools: "
                        "dmu0 += R; dmu1 += k*(mu2 - mu1)/2; "
                        "dmu2 += k*(2*mu3 - 3*mu2 + mu1)/6 -- STABLE "
                        "direct forms (round-67 P2): never B1/B2 ratios "
                        "re-multiplied by R (catastrophic cancellation "
                        "near DP -> 1 exhaustion)"),
    "totals": ("dmu0_total = +R = k(T)*max(mu1 - mu0, 0); dmu1_total = 0 "
               "(mass conserved, machine precision); dmu2_total = "
               "k*(mu1 - mu3)/3 = -k*(mu3 - mu1)/3 (the legacy "
               "Ziff-McGrady random-scission second-moment source)"),
    "out_of_domain": ("zero flux (warn once per pool per rebuild) when "
                      "mu2 < mu1, mu3 is nonfinite, or "
                      "2*mu3 - 3*mu2 + mu1 < 0 (B1, B2 >= 0 does NOT imply "
                      "a nonnegative daughter mu2 credit; round-67 P1)"),
    "reversibility": ("one-way, NO gas release: homolysis releases no "
                      "volatiles, and recombination arrives via the "
                      "discovered chemistry conduit, never this kernel"),
}


def _serialize_homolysis_initiation_block(pool: 'Polymer'
                                          ) -> Optional[Dict[str, Any]]:
    """Serialize the pool's k_homolysis kernel as the pool-level
    ``homolysis_initiation`` block (schema 2.6), or return None when the
    kernel is absent (kernel-free pools emit nothing -- presence-gated
    artifacts stay byte-identical, pinned by test).

    The triplet is re-normalized through ``validate_k_homolysis`` (the
    shared single source of truth, same posture as
    ``_serialize_radical_qssa_channel``) so a directly-constructed Polymer
    carrying a garbage dict raises here instead of poisoning the artifact.
    The daughter fields are named for the open-*1/open-*2 POSITIONAL
    contract (round-67 ruling (a)); their values follow the ratified label
    convention (K_HOMOLYSIS_DAUGHTER_SUFFIXES), which consumers pin."""
    khom = getattr(pool, "k_homolysis", None)
    if khom is None:
        return None
    # Lazy import (the QSSA serializer's cycle-avoidance idiom).
    from rmgpy.solver.polymer import (K_HOMOLYSIS_DAUGHTER_SUFFIXES,
                                      validate_k_homolysis)
    label = getattr(pool, "label", "")
    trip = validate_k_homolysis(label, khom)
    return {
        "enabled": True,
        "kinetics": {
            "A": float(trip["A"]), "n": float(trip["n"]),
            "Ea": float(trip["Ea"]),
            "units": {"A": "s^-1", "Ea": "J/mol"},
        },
        "open_site_1_radical_pool":
            f"{label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]}",
        "open_site_2_radical_pool":
            f"{label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]}",
        "kernel": HOMOLYSIS_KERNEL_NAME,
        "recipe_revision": HOMOLYSIS_BLOCK_RECIPE_REVISION,
        "recipe": dict(HOMOLYSIS_SIDECAR_RECIPE),
    }


def _assert_homolysis_serialization_closure(pools: List[Dict[str, Any]],
                                            carriers: List[Dict[str, Any]],
                                            condensed_labels,
                                            configured_labels,
                                            spawned_labels) -> None:
    """Producer-side closure guard for the schema-2.6 block, the exact
    mirror of the consumer's _check_homolysis_initiation (r68 adjudication:
    build_polymer_moments_artifact must never emit an artifact the
    reference loader rejects): the kernel-carrying pool itself must be
    solver-configured (the loader skips unconfigured pools -- a block
    there is a silently dropped kernel), and each of its two end-radical
    daughter pools must be present in pools[], solver-configured (in
    ``configured_labels``) and NOT spawned-classified (``spawned_labels``
    is the configured set's complement per schema 2.5), condensed (a
    non-empty ``phase_species`` fully inside the artifact's condensed
    closure), and carry the Stage-1 spawn provenance pin
    (spawn_event_metadata.source == HOMOLYSIS_SPAWN_SOURCE). Daughters are
    eagerly registered at model setup (rmgpy/rmg/input.py step 4d) and
    solver-configured (polymer.pyx _flatten_homolysis_state hard-errors
    otherwise), so a miss here means a broken caller, never a valid
    artifact."""
    by_label = {p.get("label"): p for p in pools}
    configured_labels = set(configured_labels)
    spawned_labels = set(spawned_labels)
    for carrier in carriers:
        lab = carrier.get("label", "")
        if lab not in configured_labels:
            raise ValueError(
                f"Pool '{lab}': carries a homolysis_initiation block but "
                f"is not among the solver-configured pools "
                f"(configured_pool_labels). The consumer only builds "
                f"configured pools, so the kernel would be silently "
                f"dropped on load; refusing to serialize (consumers "
                f"hard-reject it).")
        block = carrier["homolysis_initiation"]
        for field in ("open_site_1_radical_pool",
                      "open_site_2_radical_pool"):
            d_label = block[field]
            d_entry = by_label.get(d_label)
            if d_entry is None:
                raise ValueError(
                    f"Pool '{lab}': k_homolysis is enabled but its "
                    f"end-radical daughter pool '{d_label}' ({field}) is "
                    f"missing from the serialized pool registry. The "
                    f"producer spawns both daughters at model setup "
                    f"(rmgpy/rmg/input.py step 4d); refusing to serialize "
                    f"a homolysis_initiation block whose daughter the "
                    f"artifact cannot name (consumers hard-reject it).")
            if d_label not in configured_labels:
                raise ValueError(
                    f"Pool '{lab}': end-radical daughter pool '{d_label}' "
                    f"is not among the solver-configured pools "
                    f"(configured_pool_labels). Homolysis daughters are "
                    f"eagerly solver-configured by design (rmgpy/rmg/"
                    f"input.py step 4d; polymer.pyx "
                    f"_flatten_homolysis_state hard-errors otherwise) and "
                    f"the consumer only builds configured pools -- the "
                    f"kernel's credits would have no solver home; "
                    f"refusing to serialize (consumers hard-reject it).")
            if d_label in spawned_labels:
                raise ValueError(
                    f"Pool '{lab}': end-radical daughter pool '{d_label}' "
                    f"is classified in the spawned-pool closure "
                    f"(conventions.spawned_pools, the configured set's "
                    f"complement per schema 2.5). A spawned-classified "
                    f"homolysis daughter contradicts the eager-configured "
                    f"daughter design; refusing to serialize (consumers "
                    f"hard-reject it).")
            members = d_entry.get("phase_species") or []
            missing = sorted(m for m in members if m not in condensed_labels)
            if not members or missing:
                raise ValueError(
                    f"Pool '{lab}': end-radical daughter pool '{d_label}' "
                    f"is not condensed per the artifact's condensed closure "
                    f"(phase_species={members}, not condensed={missing}). "
                    f"The kernel's fragment credits land in the daughter's "
                    f"moment slots, so an un-condensed daughter is the "
                    f"item-16 mass-balance hazard shape -- pass the live "
                    f"core_species/condensed_species (the engine mask marks "
                    f"configured daughter pools condensed); refusing to "
                    f"serialize.")
            meta = d_entry.get("spawn_event_metadata")
            if not isinstance(meta, dict) or \
                    meta.get("source") != HOMOLYSIS_SPAWN_SOURCE:
                raise ValueError(
                    f"Pool '{lab}': end-radical daughter pool '{d_label}' "
                    f"lacks the Stage-1 spawn provenance pin (expected "
                    f"spawn_event_metadata.source == "
                    f"{HOMOLYSIS_SPAWN_SOURCE!r}, got {meta!r}). Homolysis "
                    f"daughters are producer-spawned pools, never user "
                    f"pools, and the consumer hard-rejects any other "
                    f"provenance; refusing to serialize.")


# Block-local revision token of the pool-level ``side_group_homolysis``
# block (schema 2.7, FR1-K2). Same contract as
# HOMOLYSIS_BLOCK_RECIPE_REVISION: names the kernel's rate recipe itself,
# independent of the artifact-level token; consumers pin it by exact match
# and reject on mismatch or absence.
SIDE_GROUP_BLOCK_RECIPE_REVISION = "2026-07-06-side-group-homolysis"
# The solver kernel the block serializes (and whose refused explicit
# gas-radical<->condensed rows it supersedes, the sibling of round-67
# ruling (c)): versioned term-type style. Consumers pin it by exact match
# -- an unknown kernel is flux they cannot reproduce.
SIDE_GROUP_KERNEL_NAME = "side_group_homolysis/1"
# Pinned A units: the kernel triplet is SI per SITE (the RHS multiplies by
# sites_per_unit*mu1), unlike the sibling kernel's plain s^-1.
SIDE_GROUP_SIDECAR_A_UNITS = "s^-1 per site"
# Normative moment law + EXACT mass contract of the side-group homolysis
# kernel, pinned in the forms actually implemented
# (rmgpy/solver/polymer.pyx, the sgh_* RHS section + PolymerPoolConfig.
# condensed_mass_g). Same contract as HOMOLYSIS_SIDECAR_RECIPE: consumers
# pin these strings independently and reject on mismatch -- editing the
# solver algebra or the mass accessor without re-verifying and re-dating
# these strings (and SIDE_GROUP_BLOCK_RECIPE_REVISION) breaks the
# contract. The "mass" entry is the NORMATIVE mass formula (round-72 K2
# scope): every consumer of a feature pool's condensed mass MUST use it.
SIDE_GROUP_SIDECAR_RECIPE = {
    "event_rate": ("per channel: R = k(T)*s*mu1 [mol/(m^3 s)]; "
                   "k(T) = A*T^n*exp(-Ea/(R_gas*T)) per site evaluated at "
                   "the RUNTIME temperature; s = sites_per_unit -- every "
                   "repeat unit carries s X sites, so a chain of length n "
                   "reacts at k*s*n and the reacting chain is picked with "
                   "probability ~ n (site-weighted)"),
    "parent_debit": ("dmu_j -= k*s*mu_{j+1} for j = 0, 1, 2; mu3 from the "
                     "log_lagrange/1 closure -- STABLE direct product "
                     "forms, never picked-chain ratios re-multiplied "
                     "by R"),
    "feature_credit": ("the channel's X-loss feature pool "
                       "('{parent}_sidegrp_{sanitized channel label}') is "
                       "credited EXACTLY the parent debit -- the chain "
                       "transfers INTACT (no chain cut), two-pool "
                       "mu0/mu1/mu2 totals conserved bitwise, arriving-"
                       "flux mean length = the parent's length-biased "
                       "mean mu2/mu1"),
    "gas_credit": ("dn(gas_product) += R via the small-source -> "
                   "dn_dt*V_poly path (bit-identical to the feature mu0 "
                   "credit); channels sharing one gas species are "
                   "additive"),
    "out_of_domain": ("zero flux (warn once per (pool, channel) per "
                      "rebuild) when mu2 < mu1, mu3 is nonfinite, or "
                      "mu3 < mu2; mu1 <= 0 / mu0 <= 0 is the silent "
                      "in-domain zero (no sites, no events); a degenerate "
                      "k(T) raises"),
    "mass": ("condensed_mass_g = mu1*monomer_mw_g_mol - "
             "mu0*chain_mass_defect_g_mol -- NORMATIVE feature-pool mass "
             "formula: the feature pool keeps the parent's monomer_mw "
             "(intact-chain transfer) and every feature chain lost "
             "exactly ONE X (v1), so chain_mass_defect_g_mol = M_X of the "
             "spawning channel's gas_product and d(condensed mass)/dt + "
             "d(gas X mass)/dt = 0 exactly on the kernel's "
             "contributions"),
    "saturation": ("v1 LIMITATION: the kernel acts on the PARENT pool "
                   "only; feature pools carry no side_group_homolysis of "
                   "their own and saturate as terminal X-loss sinks (no "
                   "multi-loss cascade)"),
}


def _serialize_side_group_homolysis_block(
        pool: 'Polymer',
        core_species: Optional[List['Species']] = None
        ) -> Optional[Dict[str, Any]]:
    """Serialize the pool's side_group_homolysis kernel as the pool-level
    ``side_group_homolysis`` block (schema 2.7, FR1-K2), or return None
    when the kernel is absent (kernel-free pools emit nothing --
    presence-gated artifacts stay byte-identical, pinned by test).

    The channel list is re-normalized through
    ``validate_side_group_homolysis`` WITH the pool's parsed monomer (the
    emitter holds the structure, so the full round-72 structural selector
    law runs: selector matches >= 1 removable X atom, sites_per_unit
    equals the match count, no two channels on one atom set) -- a
    directly-constructed Polymer carrying a garbage channel list raises
    here instead of poisoning the artifact.

    Serialized per channel, beyond the deck vocabulary:

    * ``site_atom_indices`` -- the selector's resolved match indices in
      ``pool.monomer.atoms`` order, the SAME atom order
      ``monomer_adj_list`` is written in. This is what lets the loader
      validate selector/feature closure from serialized data alone
      (round-73: K2 must not inherit the solver-backstop structural gap):
      match count == sites_per_unit and channel-pairwise atom-set
      disjointness are both checkable without re-deriving from a monomer
      graph the consumer does not have.
    * ``gas_species`` -- the artifact species label of the channel's
      registered gas-product Species (rmgpy/rmg/input.py step 4e), the
      routing the consumer wires the kernel's +R gas credit to (the
      monomer_routing idiom). Resolved from
      ``pool.side_group_gas_species`` by IDENTITY against the core
      universe when one is supplied; an unresolvable gas target is a
      DEFINED-MALFORMED shape the producer refuses to emit (the ejected X
      would silently vanish -- un-conserved mass).
    * ``gas_mw_g_mol`` -- M_X of the gas_product, the value the feature
      pool's ``chain_mass_defect_g_mol`` must pin exactly (checked by the
      producer closure guard, the consumer, and the solver's
      _flatten_side_group_state).
    * ``feature_pool`` -- the ratified X-loss feature pool label
      ('{parent}_sidegrp_{sanitized channel label}').
    """
    channels = getattr(pool, "side_group_homolysis", None)
    if not channels:
        return None
    # Lazy import (the QSSA serializer's cycle-avoidance idiom).
    from rmgpy.solver.polymer import (_side_group_gas_mw_g_mol,
                                      side_group_daughter_pool_label,
                                      side_group_site_atom_indices,
                                      validate_side_group_homolysis)
    label = getattr(pool, "label", "")
    monomer = getattr(pool, "monomer", None)
    if monomer is None or not hasattr(monomer, "atoms"):
        raise ValueError(
            f"Pool '{label}': carries the side_group_homolysis kernel but "
            f"no parsed monomer structure -- the schema-2.7 block "
            f"serializes the selector's resolved site_atom_indices "
            f"(round-73 loader-side structural closure), which cannot be "
            f"derived without the repeat unit. Refusing to serialize.")
    channels = validate_side_group_homolysis(label, channels,
                                             monomer=monomer)
    gas_list = getattr(pool, "side_group_gas_species", None) or []
    if len(gas_list) != len(channels):
        raise ValueError(
            f"Pool '{label}': side_group_homolysis has {len(channels)} "
            f"channel(s) but {len(gas_list)} registered gas-product "
            f"species (rmgpy/rmg/input.py step 4e registers one live gas "
            f"Species per channel). Refusing to serialize a block whose "
            f"gas X-radical routing the artifact cannot name (the ejected "
            f"X would silently vanish -- consumers hard-reject it).")
    out_channels = []
    for ci, ch in enumerate(channels):
        g_spc = gas_list[ci]
        if core_species and not any(spc is g_spc for spc in core_species):
            raise ValueError(
                f"Pool '{label}': side_group_homolysis channel "
                f"'{ch['label']}' gas_product species "
                f"('{getattr(g_spc, 'label', '')}') is not in the core "
                f"species universe the artifact labels come from "
                f"(identity check). Refusing to serialize a gas routing "
                f"the consumer cannot resolve.")
        gas_mol = Molecule().from_smiles(ch["gas_product"])
        idxs = side_group_site_atom_indices(
            monomer, gas_mol.atoms[0].symbol, ch["site_selector"])
        out_channels.append({
            "label": ch["label"],
            "kinetics": {
                "A": float(ch["A"]), "n": float(ch["n"]),
                "Ea": float(ch["Ea"]),
                "units": {"A": SIDE_GROUP_SIDECAR_A_UNITS,
                          "Ea": "J/mol"},
            },
            "site_selector": ch["site_selector"],
            "sites_per_unit": float(ch["sites_per_unit"]),
            "site_atom_indices": [int(i) for i in idxs],
            "gas_product": ch["gas_product"],
            "gas_species": _artifact_species_label(g_spc),
            "gas_mw_g_mol": float(
                _side_group_gas_mw_g_mol(ch["gas_product"])),
            "feature_pool": side_group_daughter_pool_label(label,
                                                           ch["label"]),
        })
    return {
        "enabled": True,
        "channels": out_channels,
        "kernel": SIDE_GROUP_KERNEL_NAME,
        "recipe_revision": SIDE_GROUP_BLOCK_RECIPE_REVISION,
        "recipe": dict(SIDE_GROUP_SIDECAR_RECIPE),
    }


def _adj_list_atom_count(adj_list_text: Any) -> int:
    """Number of numbered atom lines in an adjacency-list TEXT (one line
    per atom, the first token being the 1-based atom number) -- counted
    without building a molecule graph, the same way the consumer
    bounds-anchors the serialized site_atom_indices (r75 P1-2)."""
    if not isinstance(adj_list_text, str):
        return 0
    return sum(1 for line in adj_list_text.splitlines()
               if line.split() and line.split()[0].isdigit())


def _assert_side_group_adj_list_alignment(label: str, monomer: Any,
                                          adj_list_text: str) -> None:
    """Producer-side order-stability assertion (r75 P1-2): the schema-2.7
    site_atom_indices are computed against ``monomer.atoms`` order
    (side_group_site_atom_indices), while the loader bounds-anchors them
    in ``monomer_adj_list`` atom order -- ONE index space only if the
    adj-list serializer walked the same atoms sequence. Checked from the
    emitted text (element symbol per numbered atom line, in order);
    round-trip order preservation (to_adjacency_list walks
    molecule.atoms; from_adjacency_list appends atoms in read order) is
    pinned by test_adjacency_list_round_trip_preserves_atom_order."""
    lines = [line.split() for line in (adj_list_text or "").splitlines()]
    # Atom lines are '<number> [<*label>] <element> ...' -- labeled atoms
    # (the pool monomer's stitch atoms, e.g. '*1') carry the label as the
    # second token and the element third.
    walked = [tok[2] if tok[1].startswith("*") and len(tok) > 2 else tok[1]
              for tok in lines
              if tok and tok[0].isdigit() and len(tok) > 1]
    if not walked:
        raise ValueError(
            f"Pool '{label}': carries the side_group_homolysis kernel "
            f"but its serialized monomer_adj_list is empty/unparseable "
            f"-- the consumer bounds-anchors site_atom_indices in its "
            f"atom order and hard-rejects the block without it (r75 "
            f"P1-2). Refusing to serialize.")
    expected = [a.symbol for a in monomer.atoms]
    if walked != expected:
        raise ValueError(
            f"Pool '{label}': the serialized monomer_adj_list walks the "
            f"atom order {walked} but site_atom_indices were computed "
            f"against monomer.atoms order {expected} -- the two must be "
            f"the SAME index space or the serialized indices point at "
            f"the wrong atoms on load (r75 P1-2 order stability). "
            f"Refusing to serialize.")


def _assert_side_group_serialization_closure(pools: List[Dict[str, Any]],
                                             carriers: List[Dict[str, Any]],
                                             condensed_labels,
                                             configured_labels,
                                             spawned_labels) -> None:
    """Producer-side closure guard for the schema-2.7 block, the exact
    mirror of the consumer's _check_side_group_homolysis (the r68
    producer/consumer mirror-property, inherited from the 2.6 sibling:
    build_polymer_moments_artifact must never emit an artifact the
    reference loader rejects): the kernel-carrying pool itself must be
    solver-configured, and each channel's X-loss feature pool must be
    present in pools[], solver-configured and NOT spawned-classified,
    condensed (non-empty ``phase_species`` fully inside the condensed
    closure), provenance-pinned (spawn_event_metadata == {'source':
    SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE, 'channel': <label>}), carrying the
    parent's monomer_mw_g_mol EXACTLY (intact-chain transfer) and a
    chain_mass_defect_g_mol pinning the channel's gas M_X. The reverse
    direction is guarded too: a pool entry carrying the
    SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE provenance that NO carrier channel
    claims is an orphan the consumer hard-rejects, so the producer
    refuses it (a copy-carried defect WITHOUT the provenance is legal --
    see the guard body). Feature daughters are eagerly registered at
    model setup (rmgpy/rmg/input.py steps 4e/4f) and solver-configured
    (polymer.pyx _flatten_side_group_state hard-errors otherwise), so a
    miss here means a broken caller, never a valid artifact."""
    by_label = {p.get("label"): p for p in pools}
    configured_labels = set(configured_labels)
    spawned_labels = set(spawned_labels)
    claimed: Dict[str, tuple] = {}
    for carrier in carriers:
        lab = carrier.get("label", "")
        if lab not in configured_labels:
            raise ValueError(
                f"Pool '{lab}': carries a side_group_homolysis block but "
                f"is not among the solver-configured pools "
                f"(configured_pool_labels). The consumer only builds "
                f"configured pools, so the kernel would be silently "
                f"dropped on load; refusing to serialize (consumers "
                f"hard-reject it).")
        if "chain_mass_defect_g_mol" in carrier:
            raise ValueError(
                f"Pool '{lab}': serializes BOTH a side_group_homolysis "
                f"block and a chain_mass_defect_g_mol field "
                f"({carrier.get('chain_mass_defect_g_mol')!r}). v1 "
                f"saturation: the parent pool owns the kernel and X-loss "
                f"feature pools own the defect, NEVER both -- a defected "
                f"carrier claims its chains already lost an X while the "
                f"live kernel debits them for losing another (no "
                f"multi-loss cascade in v1; r75 P1); refusing to "
                f"serialize (consumers hard-reject it).")
        # r75 P1-2: the consumer bounds-anchors site_atom_indices in the
        # carrier's monomer_adj_list atom order and hard-rejects a block
        # without it -- the producer refuses the same shapes.
        n_atoms = _adj_list_atom_count(carrier.get("monomer_adj_list"))
        if n_atoms <= 0:
            raise ValueError(
                f"Pool '{lab}': carries a side_group_homolysis block but "
                f"its serialized monomer_adj_list names no atoms -- the "
                f"consumer bounds-anchors site_atom_indices in its atom "
                f"order and hard-rejects the block without it (r75 "
                f"P1-2); refusing to serialize.")
        block = carrier["side_group_homolysis"]
        carrier_mw = carrier.get("monomer_mw_g_mol")
        seen_sites: Dict[str, list] = {}
        for ch in block["channels"]:
            oob = sorted(int(i) for i in ch["site_atom_indices"]
                         if int(i) >= n_atoms)
            if oob:
                raise ValueError(
                    f"Pool '{lab}': side_group_homolysis channel "
                    f"'{ch['label']}' site_atom_indices {oob} are out of "
                    f"range for the serialized monomer_adj_list, which "
                    f"names {n_atoms} atom(s) (the indices are 0-based "
                    f"positions in its atom order; r75 P1-2); refusing "
                    f"to serialize (consumers hard-reject it).")
            # r75 P1-1 (mirror of the consumer's double-carry guard):
            # EMPTY pairwise intersection between channels' atom sets per
            # gas element -- overlapping (subset/superset) sets
            # double-carry the shared site exactly like identical sets.
            # Unreachable through the builder (the structural selectors
            # partition the atom classes), guarded anyway (the r68
            # broken-caller posture).
            gas_sym = Molecule().from_smiles(
                ch["gas_product"]).atoms[0].symbol
            new_set = {int(i) for i in ch["site_atom_indices"]}
            for prev_label, prev_set in seen_sites.get(gas_sym, []):
                shared = sorted(prev_set & new_set)
                if shared:
                    raise ValueError(
                        f"Pool '{lab}': side_group_homolysis channels "
                        f"'{prev_label}' and '{ch['label']}' resolve to "
                        f"the SAME {gas_sym} atom set or overlap on atom "
                        f"indices {shared} -- two rate channels claiming "
                        f"one structural site double-carry the loss "
                        f"(rounds 72/75 P1: disjointness is EMPTY "
                        f"pairwise intersection, not merely "
                        f"non-identical sets); refusing to serialize "
                        f"(consumers hard-reject it).")
            seen_sites.setdefault(gas_sym, []).append(
                (ch["label"], new_set))
            d_label = ch["feature_pool"]
            claimed[d_label] = (lab, ch["label"])
            d_entry = by_label.get(d_label)
            if d_entry is None:
                raise ValueError(
                    f"Pool '{lab}': side_group_homolysis channel "
                    f"'{ch['label']}' is enabled but its X-loss feature "
                    f"pool '{d_label}' is missing from the serialized "
                    f"pool registry. The producer spawns it at model "
                    f"setup (rmgpy/rmg/input.py step 4f); refusing to "
                    f"serialize a block whose feature pool the artifact "
                    f"cannot name (consumers hard-reject it).")
            if d_label not in configured_labels:
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' "
                    f"(channel '{ch['label']}') is not among the "
                    f"solver-configured pools (configured_pool_labels). "
                    f"Feature daughters are eagerly solver-configured by "
                    f"design (rmgpy/rmg/input.py step 4f; polymer.pyx "
                    f"_flatten_side_group_state hard-errors otherwise) "
                    f"and the consumer only builds configured pools -- "
                    f"the kernel's transferred chains would have no "
                    f"solver home; refusing to serialize (consumers "
                    f"hard-reject it).")
            if d_label in spawned_labels:
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' "
                    f"(channel '{ch['label']}') is classified in the "
                    f"spawned-pool closure (conventions.spawned_pools, "
                    f"the configured set's complement per schema 2.5). A "
                    f"spawned-classified feature daughter contradicts the "
                    f"eager-configured design; refusing to serialize "
                    f"(consumers hard-reject it).")
            members = d_entry.get("phase_species") or []
            missing = sorted(m for m in members
                             if m not in condensed_labels)
            if not members or missing:
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' "
                    f"(channel '{ch['label']}') is not condensed per the "
                    f"artifact's condensed closure "
                    f"(phase_species={members}, not condensed={missing}). "
                    f"The kernel's intact-chain credits land in the "
                    f"feature pool's moment slots, so an un-condensed "
                    f"feature pool is the item-16 mass-balance hazard "
                    f"shape; refusing to serialize.")
            meta = d_entry.get("spawn_event_metadata")
            if (not isinstance(meta, dict)
                    or meta.get("source") !=
                    SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE
                    or meta.get("channel") != ch["label"]):
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' "
                    f"lacks the FR1-K1 spawn provenance pin (expected "
                    f"spawn_event_metadata == {{'source': "
                    f"{SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE!r}, 'channel': "
                    f"{ch['label']!r}}}, got {meta!r}). Feature daughters "
                    f"are producer-spawned pools, never user pools, and "
                    f"the consumer hard-rejects any other provenance; "
                    f"refusing to serialize.")
            d_mw = d_entry.get("monomer_mw_g_mol")
            if (not isinstance(carrier_mw, float) or carrier_mw <= 0.0
                    or not isinstance(d_mw, float)
                    or abs(d_mw - carrier_mw) >
                    1.0e-9 * max(abs(carrier_mw), 1.0)):
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' "
                    f"monomer_mw_g_mol={d_mw!r} does not pin the "
                    f"parent's ({carrier_mw!r}). The chain transfers "
                    f"INTACT (the X loss is carried by "
                    f"chain_mass_defect_g_mol, never by the monomer "
                    f"mass); a diverging repeat-unit mass "
                    f"fabricates/destroys condensed mass -- refusing to "
                    f"serialize (consumers hard-reject it).")
            defect = d_entry.get("chain_mass_defect_g_mol")
            mw_x = ch["gas_mw_g_mol"]
            if (not isinstance(defect, float) or not (defect > 0.0)
                    or abs(defect - mw_x) > 1.0e-6 * mw_x):
                raise ValueError(
                    f"Pool '{lab}': X-loss feature pool '{d_label}' has "
                    f"chain_mass_defect_g_mol={defect!r} but channel "
                    f"'{ch['label']}' ejects "
                    f"gas_product={ch['gas_product']!r} with "
                    f"M_X={mw_x:g} g/mol. The NORMATIVE mass formula "
                    f"(condensed_mass_g = mu1*monomer_mw_g_mol - "
                    f"mu0*chain_mass_defect_g_mol) requires the defect "
                    f"to pin M_X exactly -- anything else mints or "
                    f"destroys condensed mass while gas X appears; "
                    f"refusing to serialize (consumers hard-reject it).")
    # Reverse closure (orphan guard), keyed on the spawn PROVENANCE only:
    # a pool stamped as a producer-spawned X-loss feature daughter that NO
    # serialized channel claims is corrupted (its label is a deterministic
    # function of (parent, channel), so unclaimed provenance means the
    # carrier block was lost) -- the consumer hard-rejects it, so the
    # producer refuses to emit it. A defect-carrying pool WITHOUT the
    # provenance is deliberately legal: Polymer.copy() carries
    # chain_mass_defect_g_mol to downstream daughters (e.g. an S2 _mod
    # child of a feature pool -- its chains each lost one X too), and the
    # NORMATIVE mass formula stays exact there with no live channel.
    for p in pools:
        lbl = p.get("label", "")
        meta = p.get("spawn_event_metadata")
        if (isinstance(meta, dict) and meta.get("source") ==
                SIDE_GROUP_HOMOLYSIS_SPAWN_SOURCE and lbl not in claimed):
            raise ValueError(
                f"Pool '{lbl}': carries the side_group_homolysis spawn "
                f"provenance (a producer-spawned X-loss feature daughter) "
                f"but NO serialized side_group_homolysis channel claims "
                f"it as its feature pool. The feature-pool label is a "
                f"deterministic function of (parent, channel), so an "
                f"unclaimed provenance pin means the carrier's block was "
                f"lost -- the kernel's flux into this pool would be "
                f"unreproducible; refusing to serialize (consumers "
                f"hard-reject it).")


# Block-local revision token of the pool-level ``end_radical_depropagation``
# block (schema 2.8, r74 SS2). Same contract as
# HOMOLYSIS_BLOCK_RECIPE_REVISION / SIDE_GROUP_BLOCK_RECIPE_REVISION: names
# the kernel's rate recipe itself, independent of the artifact-level token;
# consumers pin it by exact match and reject on mismatch or absence.
DEPROPAGATION_BLOCK_RECIPE_REVISION = "2026-07-06-end-radical-depropagation"
# The solver kernel the block serializes: versioned term-type style.
# Consumers pin it by exact match -- an unknown kernel is flux they cannot
# reproduce.
DEPROPAGATION_KERNEL_NAME = "end_radical_depropagation/1"
# Normative AS-IMPLEMENTED moment law of the end-radical depropagation
# kernel, pinned in the forms the RHS actually integrates
# (rmgpy/solver/polymer.pyx, the kdep_* RHS section + _smooth_pos +
# _deprop_dp1_fraction). r78 ruling: the block pins the ARTIFACT's actual
# integrated behavior, not the idealized r74 law -- everything a TA twin
# needs to replicate the RHS bitwise (1e-12 cross-code parity target) from
# the artifact + this block alone. Same contract as
# HOMOLYSIS_SIDECAR_RECIPE: consumers pin these strings independently and
# reject on mismatch -- editing the solver algebra without re-verifying and
# re-dating these strings (and DEPROPAGATION_BLOCK_RECIPE_REVISION) breaks
# the contract.
DEPROPAGATION_SIDECAR_RECIPE = {
    "event_rate": ("R = k(T)*mu0*g [mol/(m^3 s)] per radical-end pool "
                   "(ONE active radical end per chain); "
                   "k(T) = A*T^n*exp(-Ea/(R_gas*T)) with "
                   "R_gas = 8.314 J/(mol K), evaluated at the RUNTIME "
                   "temperature (round 66: never a precomputed scalar); "
                   "the kernel contributes only while mu0 > 0"),
    "gate": ("g = 1 - sp(1 - mu1/mu0) with sp(x) = x^3/(x^2 + W^2) for "
             "x > 0 else exactly 0 (C2 smooth positive-part), "
             "W = gate_width; g == 1 EXACTLY in the realizable region "
             "mean DP >= 1 and rolls off C2-smoothly below (r74 SS5: no "
             "max(...,0) cliff); the mu1/mu2/gas contributions apply only "
             "when R > 0"),
    "gas_credit": ("dn(gas_species) += R via the small-source -> "
                   "dn_dt*V_poly path, the SAME float as the mu1 drain: "
                   "d(condensed units) + d(gas monomer) = 0 exactly"),
    "moment_law": ("dmu1 = -R; dmu2 = -k*mu0*(g + 2*sp(mu1/mu0 - 1)) -- "
                   "the GATED smooth-pos form of -k*(2*mu1 - mu0), a "
                   "DELIBERATE O(W^2) mu2 under-drain near exhaustion "
                   "(disclosed closure orphan, r78); dmu0 = -k*mu0*p1, "
                   "deliberately UN-gated (applies whenever mu0 > 0, even "
                   "at mean DP < 1: there p1 = 1 while g < 1, so chains "
                   "drain at least as fast as units and mu1 - mu0 is "
                   "pushed back toward the realizable cone mu1 >= mu0 -- "
                   "a cone property ONLY, not a mu1 nonnegativity "
                   "guarantee: an unphysical mu1 = 0, mu0 > 0 state still "
                   "releases a small bounded monomer flux "
                   "(<= k*mu0*W^2) and can drive mu1 slightly negative, "
                   "r78 -- never a permanent dmu0 = 0)"),
    "dp1_closure": ("p1 = min(1, max(0, max(p_gamma, p_floor))), and 0 "
                    "when mu0 <= 1e-30; p_floor = 1 - (3*t^2 - 2*t^3) "
                    "with t = clamp(mu1/mu0 - 1, 0, 1) (C1 smoothstep "
                    "terminal floor over mean DP in [1, 2]); p_gamma = "
                    "max(0, F(1.5) - F(0.5)) / max(0, 1 - F(0.5)) on the "
                    "half-integer-bin gamma CDF F(x) = "
                    "gammainc_reg_lower(k_shape, x/theta) "
                    "(scipy.special.gammainc), zero when the conditioning "
                    "tail 1 - F(0.5) <= 1e-12; k_shape = 1/(PDI - 1), "
                    "theta = (mu1/mu0)/k_shape, PDI = mu2*mu0/mu1^2; the "
                    "gamma leg is 0 (floor-only) when any of mu0/mu1/mu2 "
                    "<= 1e-30 or PDI is nonfinite or PDI <= 1 + 1e-6 or "
                    "k_shape/theta is nonpositive/nonfinite; the scipy "
                    "gammainc leg is NORMATIVE (the scipy-absent discrete "
                    "fallback is not part of this contract)"),
    "out_of_domain": ("a degenerate k(T) -- anything failing "
                      "0 < k(T) < +inf (nonpositive, NaN, overflow) -- "
                      "RAISES (refusing to integrate a poisoned kernel); "
                      "there is no zero-flux fallback"),
    "exclusions": ("mutually exclusive on one pool with legacy "
                   "k_unzip > 0 (the scalar form of the SAME chain-end "
                   "release event), with radical_qssa_unzip (its "
                   "depropagation block IS this lumped channel), and with "
                   "k_homolysis (multi-generation homolysis DEFERRED, r74 "
                   "SS3); k_scission is likewise absent on every "
                   "production carrier (end-radical daughters are born "
                   "with k_scission = 0)"),
    "mass": ("NO condensed-mass formula change: condensed_mass_g stays "
             "mu1*monomer_mw_g_mol (minus mu0*chain_mass_defect_g_mol "
             "only on defect-carrying pools); each unzip event moves "
             "exactly one repeat unit of mass gas_mw_g_mol == the pool's "
             "monomer_mw_g_mol from the condensed basis to gas_species"),
}


def _species_mw_g_mol(spc: Any) -> float:
    """Molar mass [g/mol] of a live core Species, from its first molecule
    graph (kg/mol -> g/mol). Returns 0.0 when no structure is available --
    callers treat that as unresolvable and refuse."""
    mols = getattr(spc, "molecule", None) or []
    getter = getattr(mols[0], "get_molecular_weight", None) if mols else None
    if getter is None:
        return 0.0
    return float(getter()) * 1000.0


def _serialize_end_radical_depropagation_block(
        pool: 'Polymer',
        core_species: Optional[List['Species']] = None
        ) -> Optional[Dict[str, Any]]:
    """Serialize the pool's k_depropagation kernel as the pool-level
    ``end_radical_depropagation`` block (schema 2.8, r74 SS2), or return
    None when the kernel is absent OR the pool is the PARENT declaration
    surface (kernel-free / parent entries emit nothing -- presence-gated
    artifacts stay byte-identical, pinned by test).

    Placement ruling (r78): the block rides the two spawned end-radical
    DAUGHTER pool entries, where the kernel actually integrates. A pool
    carrying both k_depropagation and k_homolysis is the parent's deck
    declaration context (generate_end_radical_daughters copies the triplet
    onto both daughters; validate_configuration excludes the pair on one
    pool and parent configs never carry the triplet), so it emits NO
    block. A pool carrying k_depropagation with NEITHER k_homolysis NOR an
    end-radical daughter identity is a dead-knob shape no deck can produce
    (rmgpy/rmg/input.py requires k_homolysis) -- refuse, the narrowed
    successor of the pre-2.8 blanket producer refusal.

    The triplet is re-normalized through ``validate_k_depropagation`` (the
    shared single source of truth) so a directly-constructed Polymer
    carrying a garbage dict raises here instead of poisoning the artifact.
    The gas monomer routing is resolved from ``monomer_product_species``
    (held BY REFERENCE from the spawn copy) by IDENTITY against the core
    universe, the QSSA routing posture; its molar mass must pin the pool's
    monomer_mw_g_mol exactly (each event moves ONE repeat unit) or the
    artifact would mint/destroy mass on every unzip event."""
    kdep = getattr(pool, "k_depropagation", None)
    if kdep is None:
        return None
    label = getattr(pool, "label", "")
    if getattr(pool, "k_homolysis", None) is not None:
        # Parent declaration surface: the live blocks ride the daughters.
        return None
    if getattr(pool, "end_radical_site", None) not in ("primary",
                                                       "secondary"):
        raise ValueError(
            f"Pool '{label}': carries a k_depropagation kernel but is "
            f"NEITHER a k_homolysis declaration context (parent) NOR an "
            f"end-radical daughter pool (end_radical_site unset). The "
            f"kernel only integrates on spawned end-radical daughters, so "
            f"this is a dead-knob shape no deck can produce "
            f"(rmgpy/rmg/input.py requires k_homolysis beside "
            f"k_depropagation); refusing to serialize an artifact that "
            f"would carry an unreproducible declaration.")
    # Lazy import (the QSSA serializer's cycle-avoidance idiom).
    from rmgpy.solver.polymer import (KDEP_GATE_WIDTH,
                                      validate_k_depropagation)
    trip = validate_k_depropagation(label, kdep)
    mps = getattr(pool, "monomer_product_species", None)
    if mps is None:
        raise ValueError(
            f"Pool '{label}': carries the k_depropagation kernel but no "
            f"monomer_product_species gas routing -- the kernel releases "
            f"ONE monomer volatile per unzip event, and without a "
            f"resolvable emission target the released units would leave "
            f"the condensed phase silently un-conserved. Refusing to "
            f"serialize (consumers hard-reject a routing-free block).")
    if core_species and not any(spc is mps for spc in core_species):
        raise ValueError(
            f"Pool '{label}': k_depropagation gas routing target "
            f"('{getattr(mps, 'label', '')}') is not in the core species "
            f"universe the artifact labels come from (identity check). "
            f"Refusing to serialize a gas routing the consumer cannot "
            f"resolve.")
    gas_mw = _species_mw_g_mol(mps)
    pool_mw = float(getattr(pool, "monomer_mw_g_mol", 0.0) or 0.0)
    if (not (pool_mw > 0.0) or not (gas_mw > 0.0)
            or abs(gas_mw - pool_mw) > 1.0e-6 * pool_mw):
        raise ValueError(
            f"Pool '{label}': k_depropagation gas monomer "
            f"('{getattr(mps, 'label', '')}', M = {gas_mw:g} g/mol) does "
            f"not pin the pool's repeat-unit mass "
            f"(monomer_mw_g_mol = {pool_mw:g}). Each unzip event moves "
            f"exactly ONE repeat unit from the condensed basis into ONE "
            f"mole of the gas monomer, so a diverging molar mass "
            f"mints/destroys mass on every event; refusing to serialize "
            f"(consumers hard-reject it).")
    return {
        "enabled": True,
        "kinetics": {
            "A": float(trip["A"]), "n": float(trip["n"]),
            "Ea": float(trip["Ea"]),
            "units": {"A": "s^-1", "Ea": "J/mol"},
        },
        "gas_species": _artifact_species_label(mps),
        "gas_mw_g_mol": float(gas_mw),
        "gate_width": float(KDEP_GATE_WIDTH),
        "kernel": DEPROPAGATION_KERNEL_NAME,
        "recipe_revision": DEPROPAGATION_BLOCK_RECIPE_REVISION,
        "recipe": dict(DEPROPAGATION_SIDECAR_RECIPE),
    }


def _assert_depropagation_serialization_closure(pools: List[Dict[str, Any]],
                                                carriers: List[Dict[str, Any]],
                                                condensed_labels,
                                                configured_labels,
                                                spawned_labels) -> None:
    """Producer-side closure guard for the schema-2.8 block, the exact
    mirror of the consumer's _check_end_radical_depropagation (the r68
    producer/consumer mirror-property): every carrier must be an
    end-radical DAUGHTER pool entry (ratified suffix + Stage-1 spawn
    provenance pin), solver-configured, never spawned-classified,
    condensed per the artifact's condensed closure; its PARENT entry must
    exist and carry the homolysis_initiation block naming this carrier at
    the matching open-site field (the kernel's initiation feed); its
    SIBLING daughter must carry the block with an IDENTICAL kinetics
    triplet (the producer copies ONE parent-declared triplet onto both);
    the solver's validate_configuration exclusion set must hold on the
    carrier entry (no legacy unzip A > 0, no radical_qssa_unzip channel,
    no homolysis_initiation block -- PLUS the r78-adjudicated k_scission
    rejection: daughters are born with k_scission = 0, so a
    scission-carrying carrier is a direct-config-only shape no generating
    run ever integrated); and the gas routing/mass cross-pins must hold
    (block.gas_species == entry.monomer_routing; block.gas_mw_g_mol ==
    entry.monomer_mw_g_mol)."""
    from rmgpy.solver.polymer import (K_HOMOLYSIS_DAUGHTER_SUFFIXES,
                                      KDEP_GATE_WIDTH)
    by_label = {p.get("label"): p for p in pools}
    configured_labels = set(configured_labels)
    spawned_labels = set(spawned_labels)
    open_site_field = {
        K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]: "open_site_1_radical_pool",
        K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]: "open_site_2_radical_pool",
    }
    for carrier in carriers:
        lab = carrier.get("label", "")
        block = carrier["end_radical_depropagation"]
        suffix = next((s for s in K_HOMOLYSIS_DAUGHTER_SUFFIXES
                       if lab.endswith(s) and len(lab) > len(s)), None)
        if suffix is None:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but is not an end-radical daughter pool (label "
                f"does not follow the ratified "
                f"'<parent><suffix>' convention, suffixes "
                f"{K_HOMOLYSIS_DAUGHTER_SUFFIXES}). The kernel only "
                f"integrates on spawned end-radical daughters; refusing "
                f"to serialize (consumers hard-reject it).")
        # Spawned conjunct FIRST (the 2.6 ordering rationale): the
        # contradiction is named even when the carrier is also missing
        # from the configured set.
        if lab in spawned_labels:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but is classified in the spawned-pool closure "
                f"(conventions.spawned_pools, the configured set's "
                f"complement per schema 2.5) -- contradicts the "
                f"eager-configured daughter design; refusing to serialize "
                f"(consumers hard-reject it).")
        if lab not in configured_labels:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but is not among the solver-configured pools "
                f"(configured_pool_labels). The consumer only builds "
                f"configured pools, so the kernel would be silently "
                f"dropped on load; refusing to serialize (consumers "
                f"hard-reject it).")
        members = carrier.get("phase_species") or []
        missing = sorted(m for m in members if m not in condensed_labels)
        if not members or missing:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but is not condensed per the artifact's condensed "
                f"closure (phase_species={members}, not "
                f"condensed={missing}). The kernel drains this pool's "
                f"condensed moment slots, so an un-condensed carrier is "
                f"the item-16 mass-balance hazard shape; refusing to "
                f"serialize.")
        meta = carrier.get("spawn_event_metadata")
        if not isinstance(meta, dict) or \
                meta.get("source") != HOMOLYSIS_SPAWN_SOURCE:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but lacks the Stage-1 spawn provenance pin "
                f"(expected spawn_event_metadata.source == "
                f"{HOMOLYSIS_SPAWN_SOURCE!r}, got {meta!r}). The kernel's "
                f"home is a producer-spawned end-radical daughter, never "
                f"a user pool; refusing to serialize (consumers "
                f"hard-reject it).")
        # Solver validate_configuration exclusion mirror + r78 k_scission.
        channels = carrier.get("channels") or {}
        unzip_a = float((channels.get("unzip") or {}).get("A", 0.0) or 0.0)
        if unzip_a > 0.0:
            raise ValueError(
                f"Pool '{lab}': serializes an end_radical_depropagation "
                f"block AND legacy unzip A={unzip_a:g} > 0. Legacy "
                f"k_unzip is the scalar form of the SAME chain-end "
                f"release event (validate_configuration excludes the "
                f"pair); refusing to serialize (consumers hard-reject "
                f"it).")
        if "radical_qssa_unzip" in channels:
            raise ValueError(
                f"Pool '{lab}': serializes an end_radical_depropagation "
                f"block AND a radical_qssa_unzip channel. The QSSA "
                f"channel's depropagation block IS this lumped chain-end "
                f"channel (validate_configuration excludes the pair); "
                f"refusing to serialize (consumers hard-reject it).")
        if "homolysis_initiation" in carrier:
            raise ValueError(
                f"Pool '{lab}': serializes an end_radical_depropagation "
                f"block AND a homolysis_initiation block. "
                f"Multi-generation homolysis of radical-ended chains is "
                f"DEFERRED (r74 SS3; validate_configuration excludes the "
                f"pair); refusing to serialize (consumers hard-reject "
                f"it).")
        scission_a = float((channels.get("scission") or {})
                           .get("A", 0.0) or 0.0)
        if scission_a > 0.0:
            raise ValueError(
                f"Pool '{lab}': serializes an end_radical_depropagation "
                f"block AND scission A={scission_a:g} > 0. End-radical "
                f"daughters are born with k_scission = 0 "
                f"(generate_end_radical_daughters), so a "
                f"scission-carrying kernel pool is a direct-config-only "
                f"shape no generating run ever integrated (r78 "
                f"adjudication); refusing to serialize (consumers "
                f"hard-reject it).")
        # Gas routing + mass cross-pins.
        routing = carrier.get("monomer_routing")
        if not routing or routing != block.get("gas_species"):
            raise ValueError(
                f"Pool '{lab}': end_radical_depropagation gas_species="
                f"{block.get('gas_species')!r} does not match the pool "
                f"surface's monomer_routing={routing!r} -- ONE routing, "
                f"cross-pinned (the released monomer would land in an "
                f"unreproducible target); refusing to serialize "
                f"(consumers hard-reject it).")
        carrier_mw = carrier.get("monomer_mw_g_mol")
        gas_mw = block.get("gas_mw_g_mol")
        if (not isinstance(carrier_mw, float) or carrier_mw <= 0.0
                or not isinstance(gas_mw, float)
                or abs(gas_mw - carrier_mw) > 1.0e-6 * carrier_mw):
            raise ValueError(
                f"Pool '{lab}': end_radical_depropagation gas_mw_g_mol="
                f"{gas_mw!r} does not pin the carrier's "
                f"monomer_mw_g_mol={carrier_mw!r}. Each unzip event moves "
                f"exactly ONE repeat unit into ONE mole of gas monomer, "
                f"so a diverging molar mass mints/destroys mass on every "
                f"event; refusing to serialize (consumers hard-reject "
                f"it).")
        # r79 P1: a defect-bearing deprop carrier mints mass. The block's
        # mass contract moves exactly one FULL repeat unit (gas_mw_g_mol
        # == monomer_mw_g_mol) per event, but a defect-bearing pool's
        # condensed mass is mu1*monomer_mw - mu0*defect: a terminal DP1
        # event (dmu0 = dmu1 = -R) drains only R*(monomer_mw - defect) of
        # condensed mass while the gas monomer credits R*monomer_mw --
        # net +R*defect minted per event. The side-group guards
        # deliberately legalize copied defect pools (non-side-group
        # provenance), so THIS is the closing conjunct.
        defect = carrier.get("chain_mass_defect_g_mol")
        if defect is not None and defect != 0:
            raise ValueError(
                f"Pool '{lab}': serializes an end_radical_depropagation "
                f"block AND chain_mass_defect_g_mol={defect!r} (r79 P1). "
                f"The block's mass contract moves one FULL repeat unit "
                f"(gas_mw_g_mol == monomer_mw_g_mol) per unzip event, "
                f"but a defect-bearing pool's condensed mass is "
                f"mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol: a "
                f"terminal DP1 event (dmu0 = dmu1 = -R) drains only "
                f"R*(monomer_mw - defect) of condensed mass while the "
                f"gas monomer credits R*monomer_mw -- minting R*defect "
                f"of mass per event. Depropagation of defect-bearing "
                f"chains (v2) needs a different mass law / gas product; "
                f"refusing to serialize (consumers hard-reject it).")
        if block.get("gate_width") != KDEP_GATE_WIDTH:
            raise ValueError(
                f"Pool '{lab}': end_radical_depropagation gate_width="
                f"{block.get('gate_width')!r} does not pin the solver "
                f"constant KDEP_GATE_WIDTH={KDEP_GATE_WIDTH!r} BITWISE "
                f"(the RHS integrates with exactly that width -- any "
                f"other value is a different law); refusing to serialize "
                f"(consumers hard-reject it).")
        # Parent closure: the initiation feed exists and names this
        # carrier at the matching open-site field.
        parent_lab = lab[:-len(suffix)]
        parent = by_label.get(parent_lab)
        if parent is None or "homolysis_initiation" not in parent:
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but its parent pool '{parent_lab}' is "
                f"{'missing from the serialized registry' if parent is None else 'not a homolysis_initiation carrier'}. "
                f"The kernel's home is a homolysis-spawned daughter with "
                f"a live initiation feed; refusing to serialize "
                f"(consumers hard-reject it).")
        if parent["homolysis_initiation"][open_site_field[suffix]] != lab:
            raise ValueError(
                f"Pool '{lab}': parent pool '{parent_lab}' "
                f"homolysis_initiation does not name this carrier at "
                f"{open_site_field[suffix]!r} -- broken daughter/parent "
                f"closure; refusing to serialize (consumers hard-reject "
                f"it).")
        # Sibling symmetry: ONE parent-declared triplet on BOTH daughters.
        other = next(s for s in K_HOMOLYSIS_DAUGHTER_SUFFIXES
                     if s != suffix)
        sibling_lab = f"{parent_lab}{other}"
        sibling = by_label.get(sibling_lab)
        sib_block = (sibling or {}).get("end_radical_depropagation")
        if not isinstance(sib_block, dict):
            raise ValueError(
                f"Pool '{lab}': carries an end_radical_depropagation "
                f"block but its sibling daughter '{sibling_lab}' does "
                f"not. The producer copies ONE parent-declared triplet "
                f"onto BOTH spawned daughters "
                f"(generate_end_radical_daughters), so a one-sided block "
                f"is a corrupted artifact; refusing to serialize "
                f"(consumers hard-reject it).")
        mine = {k: block["kinetics"][k] for k in ("A", "n", "Ea")}
        theirs = {k: sib_block.get("kinetics", {}).get(k)
                  for k in ("A", "n", "Ea")}
        if mine != theirs:
            raise ValueError(
                f"Pool '{lab}': end_radical_depropagation kinetics "
                f"{mine} diverge from sibling '{sibling_lab}' kinetics "
                f"{theirs}. The producer copies ONE parent-declared "
                f"triplet onto BOTH daughters, so asymmetric siblings are "
                f"a corrupted artifact; refusing to serialize (consumers "
                f"hard-reject it).")


def _serialize_radical_qssa_channel(pool: 'Polymer') -> Optional[Dict[str, Any]]:
    """Serialize the pool's radical_qssa_unzip config for the sidecar, or
    return None when the channel is absent (legacy pools emit nothing).

    The config is re-normalized through validate_radical_qssa_unzip (the
    shared single source of truth for field rules) so a directly-constructed
    Polymer carrying a garbage dict raises here instead of poisoning the
    artifact, and deck-omitted defaults (efficiency/monomer_yield/basis,
    transfer: null) are always filled in the emitted block.

    OLD-CONSUMER LOUD FAILURE (milestone-3 requirement): the block lives
    INSIDE the pool's ``channels`` dict deliberately. The TA loader
    (~/Code/TA/ta/mechanism.py) pins ``SUPPORTED_CHANNELS = frozenset(
    {"scission", "unzip"})`` (mechanism.py:52) and hard-errors on any other
    key inside ``channels`` (mechanism.py:509-517: "unrecognized channel
    key(s) ... A new channel requires a schema minor-version bump and
    consumer support"). RMG stamps schema_version 2.1 on any artifact
    carrying this block (build_polymer_moments_artifact; the minor bump that
    guard demands), and TA's 2.x parser is minor-permissive, so its v2 pool
    parser (the one with that guard) is still the path taken. An old TA
    therefore refuses a QSSA-enabled sidecar loudly instead of silently
    integrating without the channel -- the silent path would produce a
    flat/false TGA.
    """
    qssa = getattr(pool, "radical_qssa_unzip", None)
    if qssa is None:
        return None
    # Lazy import: rmgpy.solver.polymer never imports rmgpy.polymer (cycle
    # note at polymer.pyx:129-131), but keep the compiled-solver import out
    # of this pure-Python module's import time anyway.
    from rmgpy.solver.polymer import QSSA_R_GAS, validate_radical_qssa_unzip
    q = validate_radical_qssa_unzip(getattr(pool, "label", ""), qssa)

    def _triplet(block_name):
        trip = q[block_name]
        if trip is None:
            return None
        return {
            "A": float(trip["A"]), "n": float(trip["n"]), "Ea": float(trip["Ea"]),
            "units": {"A": RADICAL_QSSA_SIDECAR_A_UNITS[block_name],
                      "Ea": "J/mol"},
        }

    # Weak-link allyl/U-state vocabulary (schema 2.2, milestone iv --
    # replaces the milestone-i anti-silent-no-op refusal with the real
    # serialization). validate_radical_qssa_unzip guarantees the group is
    # all-or-nothing and that the legacy SUMMED 'termination' is absent from
    # a weak-link config, so 'initiation_allyl' membership is the complete
    # discriminator. The summed termination slot is emitted explicitly null
    # (split triplets replace it); U0 is STATE, serialized as value + pinned
    # units note, not as a rate triplet. The artifact carrying any of this
    # is stamped schema 2.2 by build_polymer_moments_artifact, so an old
    # consumer/loader fails loudly instead of laundering the channel.
    if "initiation_allyl" in q:
        return {
            "enabled": True,
            "basis": q["basis"],
            "efficiency": float(q["efficiency"]),
            "monomer_yield": float(q["monomer_yield"]),
            "initiation": _triplet("initiation"),
            "depropagation": _triplet("depropagation"),
            "termination": None,
            "transfer": _triplet("transfer"),
            "initiation_allyl": _triplet("initiation_allyl"),
            "termination_recombination":
                _triplet("termination_recombination"),
            "termination_disproportionation":
                _triplet("termination_disproportionation"),
            "unsaturated_tail_ends_initial": {
                "value": float(q["unsaturated_tail_ends_initial"]),
                "units": RADICAL_QSSA_SIDECAR_U0_UNITS,
            },
            "recipe": dict(RADICAL_QSSA_SIDECAR_RECIPE_WEAKLINK),
            "provenance": {
                "radical_balance": (
                    "G_R = 2*f*ki*B + f*ki_allyl*u_active; loss = ktr*R + "
                    "2*kt_total*R^2; Rss no-transfer = "
                    "sqrt(G_R/(2*kt_total))"),
                "moment_closure": "end_shrink_pool_mean/1",
                "R_gas_J_per_mol_K": float(QSSA_R_GAS),
                "concentration_basis": "mol/m^3 condensed volume",
                "transfer_note": ("ktr is pseudo-first-order (s^-1); "
                                  "bimolecular literature k_tr must be "
                                  "premultiplied by substrate concentration "
                                  "before entering this config"),
            },
        }

    return {
        "enabled": True,
        "basis": q["basis"],
        "efficiency": float(q["efficiency"]),
        "monomer_yield": float(q["monomer_yield"]),
        "initiation": _triplet("initiation"),
        "depropagation": _triplet("depropagation"),
        "termination": _triplet("termination"),
        "transfer": _triplet("transfer"),
        "recipe": dict(RADICAL_QSSA_SIDECAR_RECIPE),
        "provenance": {
            "radical_balance": ("G_R = 2*f*ki*B; loss = ktr*R + 2*kt*R^2; "
                                "Rss no-transfer = sqrt(f*ki*B/kt)"),
            "moment_closure": "end_shrink_pool_mean/1",
            "R_gas_J_per_mol_K": float(QSSA_R_GAS),
            "concentration_basis": "mol/m^3 condensed volume",
            "transfer_note": ("ktr is pseudo-first-order (s^-1); bimolecular "
                              "literature k_tr must be premultiplied by "
                              "substrate concentration before entering this "
                              "config"),
        },
    }


def _derive_qssa_monomer_routing(pool: 'Polymer',
                                 core_species: Optional[List['Species']]) -> str:
    """Resolve the monomer-routing artifact label for an enabled-QSSA pool
    that the engine's routing map does not cover (round-25 P1-1).

    The routing target is the pool's own ``monomer_product_species`` -- the
    same live Species the engine's ``monomer_poly_index`` points at for
    configured pools (deck: input.py:432; daughters hold it BY REFERENCE via
    ``_inherit_unzip_channel``). When a core universe is supplied, membership
    is checked by IDENTITY (routing resolution is object-keyed everywhere:
    ``derive_daughter_pool_configs``' spc_map, ``PolymerPool.to_config``).

    Raises ``ValueError`` when the routing cannot be resolved: the
    enabled-QSSA-without-routing artifact shape is defined malformed
    (consumer-side hard reject), so the PRODUCER refuses to emit it.
    """
    label = getattr(pool, "label", "")
    mps = getattr(pool, "monomer_product_species", None)
    if mps is None:
        raise ValueError(
            f"Pool '{label}': radical_qssa_unzip is enabled but no "
            f"monomer_routing could be resolved -- the pool has no engine "
            f"routing entry and no monomer_product_species. Refusing to "
            f"serialize: enabled-QSSA-without-routing is a defined-malformed "
            f"artifact shape (consumers hard-reject it)."
        )
    if core_species and not any(spc is mps for spc in core_species):
        raise ValueError(
            f"Pool '{label}': radical_qssa_unzip is enabled but its "
            f"monomer_routing target ('{getattr(mps, 'label', '')}') is not "
            f"in the core species universe the artifact labels come from "
            f"(identity check). Refusing to serialize enabled QSSA with "
            f"unresolvable monomer_routing."
        )
    return _artifact_species_label(mps)


def _serialize_pool_for_sidecar(pool: 'Polymer',
                                core_species: Optional[List['Species']] = None,
                                monomer_routing: Optional[str] = None,
                                spawned: bool = False,
                                initial_explicit_moles: Optional[Dict[int, float]] = None
                                ) -> Dict[str, Any]:
    """Convert a :class:`Polymer` instance to a JSON-serialisable dict.

    Schema 2.0: 1.0 fields (docs/multi_pool_design.md §6) preserved verbatim;
    additions per docs/polymer_moments_format.md §2.

    ``spawned`` is set by ``build_polymer_moments_artifact``'s
    ``_pool_is_spawned``; direct callers serializing spawned pools must
    pass ``True``.

    ``initial_explicit_moles`` is the pool's ``{dp: moles}`` slice of the
    stage-A solver contract (``HybridPolymerSystem.initial_explicit_species``
    keyed by pool label) — consumed by the schema-2.3 ``explicit_dp`` block
    only; ignored for pools without a handshake species. Absent entries
    default to 0.0 (the species starts empty), never a hole.
    """
    # FR1-K2 fallthrough guard (the successor of the K1 'schema 2.7
    # pending' hard-fail): an X-loss feature-pool identity whose exact
    # mass-defect contract the schema-2.7 vocabulary cannot express
    # (chain_mass_defect_g_mol missing/non-positive/non-finite) must
    # HARD-FAIL here rather than emit an artifact that silently launders
    # the round-70 P1 mass defect. Kernel-carrying pools serialize below
    # (_serialize_side_group_homolysis_block raises its own refusals for
    # any channel shape it cannot serialize).
    if getattr(pool, "side_loss_channel", None):
        _defect = getattr(pool, "chain_mass_defect_g_mol", None)
        try:
            _defect_ok = _defect is not None and math.isfinite(
                float(_defect)) and float(_defect) > 0.0
        except (TypeError, ValueError):
            _defect_ok = False
        if not _defect_ok:
            raise ValueError(
                f"Pool '{getattr(pool, 'label', '')}': carries the X-loss "
                f"feature-pool identity (side_loss_channel="
                f"{pool.side_loss_channel!r}) but no expressible exact "
                f"mass-defect contract (chain_mass_defect_g_mol="
                f"{_defect!r} must be a finite value > 0, the spawning "
                f"channel's gas M_X). Refusing to emit an artifact that "
                f"silently omits the defect (the round-70 P1 mass-minting "
                f"trap on the consumer side).")
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
    # Radical QSSA unzip channel (milestone 3). Absent -> emit nothing (legacy
    # pools stay byte-identical). When present, the enabled block goes INSIDE
    # channels so an old TA fails loudly on it (see
    # _serialize_radical_qssa_channel for the file:line mechanism), and the
    # pool's k_unzip stays 0 (mutual exclusion is an upstream hard error at
    # deck read / to_config / solver validate_configuration). Monomer routing
    # REUSES the existing top-level monomer_routing field: the live main.py
    # hook keys routing on monomer_poly_index (rmgpy/rmg/main.py, the
    # `idx = getattr(p, "monomer_poly_index", None)` loop in save_everything),
    # which is guaranteed non-None for a QSSA pool by three layers of guards
    # (rmgpy/rmg/input.py deck read, PolymerPool.to_config, solver
    # validate_configuration) -- NOT on k_unzip > 0, so a QSSA pool with
    # k_unzip == 0 gets its routing populated the same way an unzip pool does.
    qssa_block = _serialize_radical_qssa_channel(pool)
    if qssa_block is not None:
        d["channels"]["radical_qssa_unzip"] = qssa_block
        # Round-25 P1-1: enabled-QSSA + null routing is a DEFINED-MALFORMED
        # shape (the consumer hard-rejects it; polymer_moments_runner
        # test_rejects_enabled_without_routing). The live save_everything
        # hook keys monomer_routing_by_pool on the ENGINE's configured pools
        # only, so a daughter registered after the last solver rebuild has
        # no entry there -- derive its routing from the pool's own
        # monomer_product_species (held BY REFERENCE from M5 inheritance)
        # against the same core universe the routing labels come from, and
        # HARD-ERROR when that is impossible. The producer never emits the
        # malformed shape.
        if monomer_routing is None:
            monomer_routing = _derive_qssa_monomer_routing(pool, core_species)
    # Explicit-DP handshake block (schema 2.3, stage B). Absent -> emit
    # nothing (flag-OFF pools stay byte-identical, golden-pinned). The block
    # is POOL-LEVEL (state/topology, not kinetics), so unlike the QSSA block
    # it does NOT live inside channels; an old consumer is stopped by the
    # 2.3 schema stamp itself (STRICT-MINOR acceptance, format doc §0/§10)
    # rather than by the channel-key guard.
    explicit_dp_block = _serialize_explicit_dp_block(
        pool, core_species=core_species,
        initial_explicit_moles=initial_explicit_moles)
    if explicit_dp_block is not None:
        d["explicit_dp"] = explicit_dp_block
    # Radical-homolysis initiation block (schema 2.6, Stage 2). Absent ->
    # emit nothing (kernel-free pools stay byte-identical, presence-gated).
    # POOL-LEVEL like explicit_dp (the kernel is solver-owned rate algebra
    # keyed on the pool, not a reactions[] channel); an old consumer is
    # stopped by the 2.6 schema stamp itself (STRICT-MINOR acceptance).
    homolysis_block = _serialize_homolysis_initiation_block(pool)
    if homolysis_block is not None:
        d["homolysis_initiation"] = homolysis_block
    # Side-group homolysis block (schema 2.7, FR1-K2). Absent -> emit
    # nothing (kernel-free pools stay byte-identical, presence-gated).
    # POOL-LEVEL like homolysis_initiation; an old consumer is stopped by
    # the 2.7 schema stamp itself (STRICT-MINOR acceptance).
    side_group_block = _serialize_side_group_homolysis_block(
        pool, core_species=core_species)
    if side_group_block is not None:
        # r75 P1-2 order-stability assertion: the block's serialized
        # site_atom_indices (computed against pool.monomer.atoms) and the
        # entry's monomer_adj_list (the atom order the loader
        # bounds-anchors them in) must be ONE index space -- refuse to
        # emit a kernel whose adj-list write failed or walked a different
        # atoms sequence.
        _assert_side_group_adj_list_alignment(
            d["label"], pool.monomer, monomer_adj_list)
        d["side_group_homolysis"] = side_group_block
    # End-radical depropagation block (schema 2.8, r74 SS2 -- the real
    # serialization that replaced the 'schema 2.8 pending' producer
    # hard-refusal). Absent -> emit nothing (kernel-free pools AND the
    # parent's declaration surface stay byte-identical, presence-gated).
    # POOL-LEVEL like homolysis_initiation; an old consumer is stopped by
    # the 2.8 schema stamp itself (STRICT-MINOR acceptance). Monomer
    # routing REUSES the existing top-level monomer_routing field (the
    # QSSA/unzip idiom): the live save_everything hook keys routing on the
    # ENGINE's monomer_poly_index (non-None for every kernel daughter by
    # three layers of guards), and a daughter serialized without an engine
    # routing entry derives it from its own monomer_product_species (held
    # BY REFERENCE from the spawn copy) -- the block's gas_species and the
    # pool field are cross-pinned by the closure guard.
    deprop_block = _serialize_end_radical_depropagation_block(
        pool, core_species=core_species)
    if deprop_block is not None:
        d["end_radical_depropagation"] = deprop_block
        if monomer_routing is None:
            monomer_routing = deprop_block["gas_species"]
    # X-loss feature-pool mass contract (schema 2.7): the pool-level
    # chain_mass_defect_g_mol POOL PROPERTY (= M_X of the spawning
    # channel's gas_product; producer-pinned). Emitted ONLY when the pool
    # carries a positive defect so every pre-existing pool entry stays
    # byte-identical; consumers apply the NORMATIVE mass formula
    # condensed_mass_g = mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol.
    _pool_defect = getattr(pool, "chain_mass_defect_g_mol", 0.0) or 0.0
    if float(_pool_defect) > 0.0:
        d["chain_mass_defect_g_mol"] = float(_pool_defect)
    phase_species: List[str] = []
    bookkeeping_species: List[str] = []
    if core_species:
        member_bases = {pool.label, f"{pool.label}_mu0", f"{pool.label}_mu1",
                        f"{pool.label}_mu2"}
        explicit_dp_spc = getattr(pool, "explicit_dp_species", None)
        for spc in core_species:
            if _species_base_label(spc) in member_bases:
                artifact_label = _artifact_species_label(spc)
                phase_species.append(artifact_label)
                # Everything THIS branch collects is bookkeeping by
                # construction: the pool's canonical proxy (concentration
                # pinned to 1.0 by the site-scaling rule) and the three
                # µ-dummies (carry ~0 mol; host the moments). Real condensed
                # entries (routed monomer below; explicit-DP chains) never
                # come through this branch — no name/suffix re-classification.
                bookkeeping_species.append(artifact_label)
            elif explicit_dp_spc is not None and spc is explicit_dp_spc:
                # The explicit-DP handshake target is REAL condensed
                # inventory (format doc §2: phase_species lists explicit-DP
                # chains; the bookkeeping subset is exactly the complement's
                # complement — proxy + µ-dummies). Identity match only: no
                # name/suffix classification.
                phase_species.append(_artifact_species_label(spc))
    # monomer_routing is deliberately NOT appended to phase_species (recipe
    # revision 2026-07-03-monomer-gas): the release target is a GAS-phase
    # species -- the unzip/QSSA channel deposits released monomer into the
    # gas amount basis (docs/polymer_moments_format.md §2/§5 revision note).
    d["phase_species"] = phase_species
    d["bookkeeping_species"] = bookkeeping_species
    d["monomer_routing"] = monomer_routing
    d["mu3_closure"] = "log_lagrange/1"

    # --- moments provenance (item #14a, amended 2026-06-12: uniform-t=0) ---
    # pools[].moments are the pool's INITIAL CONDITIONS at t=0 of the
    # simulated experiment, normatively (docs/polymer_moments_format.md s2):
    #   input_declared — declared in the input file; the moments are the
    #                    input-derived t=0 state (the values this serializer
    #                    has always passed through, unchanged);
    #   spawned_empty  — created mid-run (gate-path drain or scission-tail
    #                    registry creation); at t=0 of any consumer
    #                    experiment the pool contains nothing, so [0,0,0] is
    #                    the honest initial condition, not a hole.
    # The moments emission above is deliberately untouched: this field is
    # ADDITIVE (no recipe_revision bump — generation-time semantics, not
    # rate recipe), and no emitted value changes for any pool that exists
    # today.
    d["moments_provenance"] = "spawned_empty" if spawned else "input_declared"
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
    initial_explicit_by_pool=None,
    generation_mass_transfer=None,
    generation_v_poly_m3=None,
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
    initial_explicit_by_pool : dict, optional
        ``{pool_label: {dp: moles}}`` — the stage-A solver contract
        (``HybridPolymerSystem.initial_explicit_species``). Feeds the
        schema-2.3 ``explicit_dp`` block's ``initial_moles``; entries
        default to 0.0 when absent.
    generation_mass_transfer : list of dict, optional
        Deck-declared mass-transfer entries
        (``{"gas_species", "poly_species", "K", "kLa"}`` with artifact
        labels + SI floats). Emitted into the NON-normative
        ``conventions.generation_defaults`` provenance block; omitted deck
        mass_transfer -> key absent.
    generation_v_poly_m3 : float, optional
        The condensed volume the generating run actually used [m^3]
        (``PolymerPhase.calculate_volume``). Same non-normative provenance
        block.

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
        initial_explicit_by_pool=initial_explicit_by_pool,
        generation_mass_transfer=generation_mass_transfer,
        generation_v_poly_m3=generation_v_poly_m3,
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
    int(PolymerFluxArchetype.VOLATILE_EJECTION): "volatile_ejection/1",
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
    VE = int(PolymerFluxArchetype.VOLATILE_EJECTION)

    for rxn in core_reactions:
        arch = int(getattr(rxn, "polymer_flux_archetype", 0))
        src, dst = _resolve_reaction_pools(rxn, pool_set)
        if arch == NONE_ and src is None and dst is None:
            continue  # ordinary chemistry — Cantera handles it untouched

        # Mirror solver demotions (polymer.pyx:557-578) and the r91
        # spawned-pool demotion REFUSAL: a stamped pool-coupled row whose
        # missing required endpoint is attributable to a spawned config-less
        # Polymer participant is refused by the generating solver
        # (conduit-deferred, zero flux), so the artifact must carry the
        # refused marker with the STAMPED archetype -- never a live
        # legacy_mu1/1 row. Computed here (not read off the object): the
        # solver re-derives the refusal per rebuild without stamping
        # polymer_refused, so the emitter mirrors the same predicate against
        # the same configured-pool set.
        unresolved = False
        spawned_refused = False
        if arch == NONE_:
            arch, unresolved = UNR, True
        elif arch == UNR:
            unresolved = True
        elif ((arch in (MIG, SCI, VE) and (src is None or dst is None))
              or (arch == CHIP and src is None)):
            # VOLATILE_EJECTION is cross-pool like MIGRATION: it needs both a
            # configured source and destination pool.
            sides = []
            if src is None:
                sides.append(rxn.reactants or [])
            if dst is None and arch != CHIP:
                sides.append(rxn.products or [])
            attributable = bool(sides)
            for side in sides:
                if not any(isinstance(s, Polymer)
                           and _species_base_label(s) not in pool_set
                           for s in side):
                    attributable = False
            if attributable and (src is not None or dst is not None):
                # Pool-mapped spawned-refused row: stamped archetype kept,
                # refused marker added below.
                spawned_refused = True
            elif attributable:
                # Spawned-refused but NOT pool-mapped (no configured-pool
                # participant): the refused marker is illegal here (loader
                # guard), so keep the legacy demotion emission; the refused
                # block below still warn-louds the over-integration hazard.
                spawned_refused = True
                arch, unresolved = UNR, True
            else:
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
            # Warn for ALL non-Arrhenius entries (retained or dropped).
            # Retained entries still carry kinetics=null; consumers take
            # reversibility from the chem.yaml equation arrow (=> means
            # irreversible — get_reaction_equation mirrors rxn.reversible).
            logging.warning(
                "Polymer artifact: reaction %s has non-Arrhenius kinetics "
                "(%s); entry carries kinetics=null (no A/n/Ea) — consumers "
                "take its rate AND reversibility from chem.yaml (the "
                "equation arrow mirrors rxn.reversible).",
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
        elif arch == VE:
            entry["params"] = {"eject_units": float(getattr(rxn, "polymer_eject_units", 0.0))}
        # Refused-row marker (schema 2.4, format doc §12): the generating
        # solver zeroes this reaction's WHOLE flux (polymer.pyx
        # reaction_refused), so the row must say so or a consumer integrates
        # flux the oracle fabricated nothing for. Stamp-but-keep: the row
        # stays listed (consumers zero its Cantera multiplier from the
        # listing). refused_reason mirrors the solver census reason built
        # from the same stamps (polymer.pyx:1526-1530). Non-refused rows
        # gain NO key (absent, not false — byte-identical artifacts).
        # spawned_refused (r91): the solver refuses the row conduit-deferred
        # at initialize_model (spawned config-less pool endpoint); the
        # emitter mirrors that refusal so consumers zero the row too.
        if getattr(rxn, "polymer_refused", False) or spawned_refused:
            if entry["proxy_reactants"] or entry["proxy_products"]:
                entry["refused"] = True
                entry["refused_reason"] = (
                    "qssa-invalid"
                    if getattr(rxn, "polymer_refused_accumulating", False)
                    else "conduit-deferred")
            else:
                # The marker is legal ONLY on pool-mapped rows (loader
                # guard, both consumers): a refused reaction whose pool is
                # not solver-configured emits unmarked — warn-loud, because
                # a consumer will over-integrate this row (the generating
                # solver still zeroes it).
                logging.warning(
                    "Polymer artifact: refused reaction %s has no "
                    "pool-mapped participant under the configured pools; "
                    "emitting WITHOUT the refused marker (consumers will "
                    "integrate this row while the generating solver zeroed "
                    "it). Configure the pool or expect over-integration.",
                    equation)
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
                                   rmg_commit=None,
                                   initial_explicit_by_pool=None,
                                   generation_mass_transfer=None,
                                   generation_v_poly_m3=None):
    """Assemble the full schema-2.0 polymer moments artifact payload.

    Normative contract: docs/polymer_moments_format.md. The payload mirrors
    the HybridPolymerSystem oracle, including its init-time demotions —
    ``configured_pool_labels`` must be the SOLVER-configured pools (which can
    be a subset of ``pool_registry``: spawned daughters are registry pools
    without solver configs and run as ordinary species).

    ``core_species``/``monomer_routing`` are populated by
    ``build_polymer_moments_artifact`` in the live run (legacy callers omit
    them → ``phase_species: []``, ``monomer_routing: null``). The routing
    label names a GAS-phase species since recipe revision
    2026-07-03-monomer-gas (it is NOT appended to phase_species and NOT in
    ``conventions.condensed_species``; the unzip/QSSA release deposits into
    the gas amount basis).
    """
    if not configured_pool_labels:
        configured_pool_labels = [getattr(p, "label", "") for p in pool_registry]
    monomer_routing_by_pool = monomer_routing_by_pool or {}
    # Stage-A solver contract {pool_label: {dp: moles}} — the same shape
    # HybridPolymerSystem.initial_explicit_species integrates (the live hook
    # passes the engine's attribute verbatim). Feeds the schema-2.3
    # explicit_dp block's initial_moles; pools without an entry default 0.0.
    initial_explicit_by_pool = initial_explicit_by_pool or {}
    # input-vs-spawned for moments_provenance (item #14a amended, spec s4):
    # Path-C scission tails carry NO spawn markers on the pool object (the
    # serializer defaults their spawn_event_metadata to {'source':'input'}),
    # so the PRIMARY signal is the configured set — the live save_everything
    # hook plumbs the engine's polymer_pools, which at HEAD is exactly the
    # input-declared set (compile_polymer_phase is the only config creator).
    # Object spawn-markers are the SECONDARY signal: they keep drained
    # daughters honest in legacy default-label calls, where
    # configured_pool_labels defaults to ALL registry labels. Known
    # limitation of that legacy default: a markerless Path-C scission tail
    # classifies "input_declared" there (its moments are honest [0,0,0]
    # regardless; the live main.py path passes the configured labels and
    # classifies it "spawned_empty"). Item-16
    # coupling (spec s7): when spawned pools gain solver configs this fork
    # must be re-confirmed, not inherited silently.
    configured_set = {str(lbl) for lbl in configured_pool_labels}

    def _pool_is_spawned(p):
        if getattr(p, "label", "") not in configured_set:
            return True
        return bool(getattr(p, "spawn_metadata", None)) or \
            getattr(p, "parent_pool_label", None) is not None

    pools = [
        _serialize_pool_for_sidecar(
            p,
            core_species=core_species,
            monomer_routing=monomer_routing_by_pool.get(getattr(p, "label", "")),
            spawned=_pool_is_spawned(p),
            initial_explicit_moles=initial_explicit_by_pool.get(
                getattr(p, "label", "")),
        )
        for p in pool_registry
    ]
    reactions = compile_polymer_reaction_entries(
        core_reactions or [], core_species or [],
        configured_pool_labels, cantera_index_map)

    # Version contract: the radical_qssa_unzip channel grows the closed
    # channel vocabulary (schema minor bump -> 2.1) and adds new channel/flux
    # algebra (new recipe_revision). Both stamps are conditional on the
    # channel actually appearing in THIS artifact so that artifacts without
    # any QSSA pool keep the legacy 2.0 / 2026-06-10 stamps byte-identically
    # (pinned by test; old consumers keep loading legacy sidecars unchanged).
    qssa_present = any("radical_qssa_unzip" in (p.get("channels") or {})
                       for p in pools)
    # Weak-link vocabulary is a strict superset stamp: the artifact-level
    # schema/recipe stamps are governed by the STRONGEST vocabulary present
    # (a mixed artifact is 2.2). 'initiation_allyl' membership is the
    # complete discriminator (validate_radical_qssa_unzip pins the group
    # all-or-nothing before serialization).
    weaklink_present = any(
        "initiation_allyl" in ((p.get("channels") or {})
                               .get("radical_qssa_unzip") or {})
        for p in pools)
    # Explicit-DP vocabulary (schema 2.3) is pool-level, outside channels;
    # its presence is the strongest SHAPE stamp (2.3 > 2.2 > 2.1 > 2.0), and
    # the recipe token composes with the strongest CHANNEL vocabulary per
    # the monomer-gas 3-way rule (base/qssa/weaklink each have an
    # explicit-dp re-stamp). No explicit-DP anywhere -> the legacy stamps
    # apply byte-identically (golden-pinned).
    explicit_dp_present = any("explicit_dp" in p for p in pools)
    # Refused-row vocabulary (schema 2.4) is ROW-level, on reactions[]; its
    # presence is the strongest SHAPE stamp (2.4 > 2.3 > ...), same
    # presence-based rule as every prior minor: no refused row anywhere ->
    # the older stamps apply byte-identically (golden-pinned). It does NOT
    # touch recipe_revision: zeroing a listed row is consumption semantics
    # carried by new shape vocabulary, not new rate algebra, and the
    # STRICT-MINOR envelope already stops pre-2.4 consumers.
    refused_present = any(e.get("refused") for e in reactions)
    # Spawned-pool closure vocabulary (schema 2.5) is CONVENTIONS-level: the
    # closure complement of configured_pools within the registry, keyed on
    # the PRIMARY signal only (label not in configured set) so it is
    # disjoint from configured_pools by construction. The legacy
    # default-label call (configured defaults to ALL registry labels) has
    # an empty complement -> key absent, stamps untouched — mirroring the
    # documented legacy-default limitation of moments_provenance above.
    spawned_pool_labels = [p["label"] for p in pools
                           if p["label"] not in configured_set]
    if explicit_dp_present:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_EXPLICIT_DP
        recipe_revision = (
            POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP_WEAKLINK
            if weaklink_present
            else (POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP_QSSA
                  if qssa_present
                  else POLYMER_RATE_RECIPE_REVISION_EXPLICIT_DP))
    else:
        schema_version = (
            POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_WEAKLINK if weaklink_present
            else (POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_QSSA if qssa_present
                  else POLYMER_POOLS_SIDECAR_SCHEMA_VERSION))
        recipe_revision = (
            POLYMER_RATE_RECIPE_REVISION_WEAKLINK if weaklink_present
            else (POLYMER_RATE_RECIPE_REVISION_QSSA if qssa_present
                  else POLYMER_RATE_RECIPE_REVISION))
    if refused_present:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_REFUSED
    # Spawned-pool presence is the strongest SHAPE stamp (2.5 > 2.4 > ...);
    # recipe_revision deliberately untouched (2.4 precedent: classification
    # vocabulary, not rate algebra).
    if spawned_pool_labels:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_SPAWNED
    # Homolysis-initiation vocabulary (schema 2.6) is POOL-level, the
    # strongest SHAPE stamp of all (2.6 > 2.5 > ...), presence-gated like
    # every prior minor. It sits BESIDE 2.5: the spawned_pools emission and
    # the condensed closure below are untouched. Artifact-level
    # recipe_revision stays untouched too -- the block carries its own
    # block-local revision token (the explicit-dp 2.3 precedent).
    homolysis_carriers = [p for p in pools if "homolysis_initiation" in p]
    if homolysis_carriers:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_HOMOLYSIS
    # Side-group homolysis vocabulary (schema 2.7, FR1-K2) is POOL-level
    # twice over: the carrier's side_group_homolysis block AND the X-loss
    # feature pools' chain_mass_defect_g_mol mass contract. Either one is
    # the strongest SHAPE stamp of all (2.7 > 2.6 > ...), presence-gated
    # like every prior minor; a homolysis-only artifact stays 2.6
    # byte-identically (negative-control pinned). Artifact-level
    # recipe_revision stays untouched (the block-local token precedent).
    side_group_carriers = [p for p in pools if "side_group_homolysis" in p]
    side_group_present = bool(side_group_carriers) or any(
        "chain_mass_defect_g_mol" in p for p in pools)
    if side_group_present:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_SIDE_GROUP
    # End-radical depropagation vocabulary (schema 2.8, r74 SS2) is
    # POOL-level on the end-radical daughter entries; its presence is the
    # strongest SHAPE stamp of all (2.8 > 2.7 > ...), presence-gated like
    # every prior minor -- a deprop-free artifact keeps its older stamp
    # byte-identically (negative-control pinned). Artifact-level
    # recipe_revision stays untouched (the block-local token precedent).
    deprop_carriers = [p for p in pools
                       if "end_radical_depropagation" in p]
    if deprop_carriers:
        schema_version = POLYMER_POOLS_SIDECAR_SCHEMA_VERSION_DEPROPAGATION

    # conventions.condensed_species closure (schema 2.5): a spawned pool's
    # phase_species (canonical proxy + mu-dummies collected from the same
    # core universe, already declared condensed ROW-side) join the normative
    # condensed list the consumers key on ("Consumers MUST use these lists,
    # not name heuristics", format doc §8) — otherwise TA's normative path
    # classifies a late-spawned pool's bookkeeping species GAS and the
    # condensed mass balance silently breaks. With no spawned pool the set
    # reduces to the caller's condensed_species exactly (byte-identical,
    # golden-pinned); the change in membership semantics ships only inside
    # 2.5-stamped artifacts, which STRICT-MINOR gating hides from every
    # pre-2.5 consumer.
    condensed_labels = {_artifact_species_label(s)
                        for s in (condensed_species or [])}
    if spawned_pool_labels:
        spawned_set = set(spawned_pool_labels)
        condensed_labels.update(
            lbl for p in pools if p["label"] in spawned_set
            for lbl in (p.get("phase_species") or []))

    # Producer-side schema-2.6 closure guard (r68 mirror of the runner's
    # _check_homolysis_initiation): the carrier and both end-radical
    # daughters of every kernel-carrying pool must be solver-configured
    # (and the daughters never spawned-classified), present in pools[],
    # condensed per the FINAL condensed closure (2.5 spawned additions
    # included), and provenance-pinned -- the consumer hard-rejects any
    # other shape, so the producer refuses to emit it.
    if homolysis_carriers:
        _assert_homolysis_serialization_closure(
            pools, homolysis_carriers, condensed_labels,
            configured_set, spawned_pool_labels)
    # Producer-side schema-2.7 closure guard (the r68 mirror-property,
    # side-group edition): carrier configured; every channel's X-loss
    # feature pool present, configured, never spawned-classified,
    # condensed per the FINAL closure, provenance-pinned, monomer_mw ==
    # parent, defect == gas M_X; and no orphan X-loss pool -- the consumer
    # hard-rejects any other shape, so the producer refuses to emit it.
    if side_group_present:
        _assert_side_group_serialization_closure(
            pools, side_group_carriers, condensed_labels,
            configured_set, spawned_pool_labels)
    # Producer-side schema-2.8 closure guard (the r68 mirror-property,
    # depropagation edition): every carrier is a provenance-pinned,
    # configured, condensed end-radical daughter of a serialized
    # homolysis carrier; siblings symmetric; the solver exclusion set
    # (unzip/QSSA/homolysis, plus the r78 k_scission adjudication) holds
    # on the entry; gas routing/MW/gate_width cross-pins hold -- the
    # consumer hard-rejects any other shape, so the producer refuses to
    # emit it.
    if deprop_carriers:
        _assert_depropagation_serialization_closure(
            pools, deprop_carriers, condensed_labels,
            configured_set, spawned_pool_labels)

    conventions = {
        # format_doc mirrors schema_version (same vocabulary fork above):
        # /2.2 when the weak-link vocabulary is present, /2.1 for QSSA-only,
        # /2.0 otherwise so legacy artifacts stay byte-identical
        # (golden-pinned by test). No loader parses this token (grepped
        # rmgpy + TA) — human-facing only.
        "format_doc": ("docs/polymer_moments_format.md "
                       f"(polymer_moments_format/{schema_version})"),
        "recipe_revision": recipe_revision,
        "moment_basis": "extensive mol, DP basis (mu1 = moles of repeat units)",
        "volumes": {
            "V_poly": "constant, consumer-supplied [m^3]",
            "V_gas": "ideal gas, dynamic: V_gas = n_gas*R*T/P (1.0 m^3 floor when n_gas <= 0)",
        },
        "configured_pools": list(configured_pool_labels),
        "condensed_species": sorted(condensed_labels),
        "site_scaling": ("site = max(0, mu_scaling)/V_poly read from the first proxy "
                         "reactant's pool; multiplies ONCE; scales rf AND rr"),
        "chip_site_throttle": ("site = min(max(0,mu0), max(0,mu1)/a)/V_poly when "
                               "archetype=discrete_chip/1 and scaling=mu0 and a>0"),
        "kb_recipe": ("kb = kf/Keq; Keq(T) = (P0/(R*T))^dn * exp(-dG0/(R*T)), "
                      "P0 = 1e5 Pa, dG0 from chem.yaml NASA thermo; dn counts "
                      "ALL species incl. condensed/proxy (format doc s4 step 1)"),
        "mu3_closure": "log_lagrange/1",
        "invariants": {
            "discrete_subset": ("sum_pools(mu1) + sum_chip_species(a_i * n_i) is "
                                "invariant over the discrete-reaction subset only"),
            "with_unzip": ("add + n(monomer_routing) per pool with an active unzip "
                           "channel (unzip moves units from mu1 into that species)"),
        },
    }
    # The configured-pools closure surface (schema 2.5): registry order,
    # disjoint from configured_pools, emitted ONLY when non-empty so
    # spawned-free artifacts stay byte-identical (never an empty list).
    # Rows keep full solver-inertness semantics (format doc §2: no site
    # scaling, no conc:=1.0, no channels integration) — this key carries
    # PHASE-CLASSIFICATION membership only, plus row-side monomer_mw_g_mol /
    # spawn provenance for the consumer's mass accounting.
    if spawned_pool_labels:
        conventions["spawned_pools"] = spawned_pool_labels

    # --- generation_defaults (provenance completeness; NON-normative) ---
    # kLa/K and V_poly are consumer-supplied operating conditions by
    # contract (format doc §7/§8): the artifact never made them normative,
    # which left the generating run's own values unrecorded. Emit them here
    # as explicitly informative provenance — a consumer replaying the
    # generation experiment can start from them, but its OWN experiment
    # config always takes precedence (the "note" says so in-band). Purely
    # ADDITIVE: nothing declared -> key absent, artifacts byte-identical
    # (golden-pinned); no schema/recipe consequences either way (TA-side
    # loaders ignore unknown conventions keys — probed 2026-07-04).
    generation_defaults = {}
    if generation_mass_transfer:
        generation_defaults["mass_transfer"] = [
            {
                "gas_species": str(mt["gas_species"]),
                "poly_species": str(mt["poly_species"]),
                "K": float(mt["K"]),
                "kLa": float(mt["kLa"]),
                "units": {"K": "dimensionless", "kLa": "s^-1"},
            }
            for mt in generation_mass_transfer
        ]
    if generation_v_poly_m3 is not None:
        generation_defaults["V_poly_m3"] = float(generation_v_poly_m3)
    if generation_defaults:
        generation_defaults["note"] = (
            "generation-run values; consumer experiment config takes "
            "precedence")
        conventions["generation_defaults"] = generation_defaults

    return {
        "schema_version": schema_version,
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


# Warn-once registry for withheld/conflicting QSSA channel inheritance
# (round-25 P1-2 / P2-1 / P2-2). Keyed on (parent_label, daughter_label,
# reason) so each pool pair discloses each verdict exactly once instead of
# spamming every enlarge pass. Same idiom as _flux_archetype_warned.
_unzip_inherit_warned = set()


def _same_repeat_chemistry(daughter: 'Polymer', parent: 'Polymer') -> bool:
    """True when the daughter's repeat-unit chemistry is the parent's -- the
    M1 condition under which the elementary initiation/depropagation/
    termination constants transfer (same monomer chemistry AND monomer_mw).

    Round-25 P2-1 probe finding: ``monomer_mw_g_mol`` derives SOLELY from
    ``monomer`` (Polymer.__init__), and every daughter constructor site
    passes the parent's monomer verbatim -- so an mw-only gate would be
    vacuously true for the _mod shape. The truthful discriminator for a
    feature modification is ``feature_monomer``: a daughter whose feature
    unit differs from the parent's carries CHANGED chain chemistry (the
    defect unit participates in initiation/depropagation), so the parent's
    constants do not apply.
    """
    if getattr(daughter, 'monomer_mw_g_mol', None) != \
            getattr(parent, 'monomer_mw_g_mol', None):
        return False
    d_feat = getattr(daughter, 'feature_monomer', None)
    p_feat = getattr(parent, 'feature_monomer', None)
    if (d_feat is None) != (p_feat is None):
        return False
    if d_feat is not None:
        try:
            return bool(d_feat.is_isomorphic(p_feat))
        except Exception:
            return False
    return True


def _inherit_unzip_channel(daughter: 'Polymer', parent: 'Polymer') -> bool:
    """Daughter-pool inheritance of the radical_qssa_unzip channel (M5),
    gated on same repeat chemistry (round-25 P2-1).

    A daughter chain population (scission tail/head, END_MOD fold-back)
    shares the parent's monomer chemistry and monomer_mw, so the same
    elementary initiation/depropagation/termination constants apply -- the
    M1 decision recorded at ``derive_daughter_pool_configs``
    (rmgpy/rmg/polymer_input.py). A feature-modified daughter (changed
    ``feature_monomer``) does NOT share that chemistry: it is created
    channel-free with a once-per-pool WARNING (silent wrong constants are
    worse than a visibly inert pool). Two deliberately different copy
    semantics when inheritance applies:

    - ``radical_qssa_unzip``: DEEP-copied. Post-hoc mutation of the parent's
      channel must never reach the daughter (and vice versa) -- same aliasing
      posture as ``Polymer.copy`` (review round 21, finding 3). Exception:
      ``unsaturated_tail_ends_initial`` (weak-link U-state) is per-pool
      STATE, not a chemistry constant -- it RESETS to 0.0 on the daughter.
    - ``monomer_product_species``: shared BY REFERENCE. Routing resolution
      (the object-keyed ``spc_map`` in ``derive_daughter_pool_configs`` /
      ``PolymerPool.to_config``) needs identity with the live core Species;
      a copy would silently resolve to no core index.

    Channel-free parents bequeath no channel: the daughter stays channel-free
    (no noise). Returns True iff the channel was inherited.
    """
    channel = getattr(parent, 'radical_qssa_unzip', None)
    if not _same_repeat_chemistry(daughter, parent):
        if channel is not None:
            key = (getattr(parent, 'label', ''),
                   getattr(daughter, 'label', ''), 'changed-chemistry')
            if key not in _unzip_inherit_warned:
                _unzip_inherit_warned.add(key)
                logging.warning(
                    "Polymer pool '%s' created channel-free: its repeat-unit "
                    "chemistry (feature/monomer) differs from parent '%s', "
                    "so the parent's radical_qssa_unzip constants were NOT "
                    "inherited.",
                    getattr(daughter, 'label', ''),
                    getattr(parent, 'label', ''))
        return False
    monomer_product = getattr(parent, 'monomer_product_species', None)
    if monomer_product is not None:
        daughter.monomer_product_species = monomer_product
    if channel is None:
        return False
    inherited = deepcopy(channel)
    # CONSTANTS inherit, STATE does not (weak-link milestone i, review P1):
    # unsaturated_tail_ends_initial is the pool's initial U-state amount,
    # not a chemistry constant. State resets on spawn -- copying the
    # parent's would fabricate the same initial U on every daughter pool (a
    # hidden initiation source once the solver U-state lands). A future
    # event-specific spawn law may compute the U transfer explicitly.
    if "unsaturated_tail_ends_initial" in inherited:
        inherited["unsaturated_tail_ends_initial"] = 0.0
    daughter.radical_qssa_unzip = inherited
    return True


def _inherit_spawned_pool_channel(daughter: 'Polymer', parent: 'Polymer',
                                  intent: 'SpawnIntent') -> bool:
    """Certainty-gated QSSA inheritance for spawn-intent daughters
    (round-25 P1-2).

    ``SpawnIntent.parent_pool`` is an attribution SHORTCUT --
    ``process_polymer_candidates_multipool`` queues ``pool_registry[0]``
    verbatim (spec-§3 attribution rules), and the intent carries no record
    of the true source pool (``triggering_reaction_index`` is never
    populated on this path). The QSSA constants bind to CHEMISTRY, so
    inherit only when the intent's TRUE detected motif (``intent.monomer``
    -- the daughter's placeholder ``monomer`` attribute is the parent's and
    proves nothing) matches the attributed parent's own monomer pattern
    (``similarity_merge``, the same predicate the spawn pipeline uses).

    In the live pipeline that match is impossible by construction (Phase E
    is only reached for motifs that failed similarity_merge against EVERY
    pool), so a live-run spawned pool is channel-free: attribution is
    genuinely ambiguous and silent wrong constants are worse than a visibly
    inert pool. A once-per-pool WARNING discloses the withheld channel.
    Returns True iff the channel was inherited.
    """
    motif = getattr(intent, 'monomer', None)
    certain = (parent is not None and motif is not None
               and similarity_merge(motif, [parent]) is parent)
    if certain:
        return _inherit_unzip_channel(daughter, parent)
    if parent is not None and \
            getattr(parent, 'radical_qssa_unzip', None) is not None:
        key = (getattr(parent, 'label', ''),
               getattr(daughter, 'label', ''), 'ambiguous-parentage')
        if key not in _unzip_inherit_warned:
            _unzip_inherit_warned.add(key)
            logging.warning(
                "Spawned polymer pool '%s' created channel-free: parentage "
                "is ambiguous (attributed parent '%s' carries a "
                "radical_qssa_unzip channel, but the spawned motif does not "
                "match its monomer chemistry, so the constants were NOT "
                "inherited).",
                getattr(daughter, 'label', ''),
                getattr(parent, 'label', ''))
    return False


def merge_unzip_channel_on_dedup(existing: 'Polymer',
                                 incoming: 'Polymer') -> None:
    """Fingerprint-dedup channel merge (round-25 P2-2).

    ``CoreEdgeReactionModel._register_polymer`` returns the EXISTING Polymer
    on a fingerprint match and discards the incoming object -- first-writer
    -wins on every attribute. A daughter registered channel-free first must
    not silently swallow a later channel-bearing equivalent, so the verdict
    is transferred onto the canonical object (same posture as the durable
    gas-veto props transfer in ``make_new_species``, commit c133b34e1):

    - existing lacks the channel, incoming has it: channel DEEP-copied onto
      the existing object (the incoming is discarded; aliasing posture as
      ``_inherit_unzip_channel``).
    - ``monomer_product_species``: transferred BY REFERENCE whenever the
      existing lacks it (identity is load-bearing for routing resolution).
    - BOTH carry channels and they differ: keep the existing one and WARN
      (disclosed first-writer-wins) once per pool.
    """
    inc_mps = getattr(incoming, 'monomer_product_species', None)
    if inc_mps is not None and \
            getattr(existing, 'monomer_product_species', None) is None:
        existing.monomer_product_species = inc_mps
    inc = getattr(incoming, 'radical_qssa_unzip', None)
    if inc is None:
        return
    ex = getattr(existing, 'radical_qssa_unzip', None)
    if ex is None:
        existing.radical_qssa_unzip = deepcopy(inc)
        return
    if ex != inc:
        key = (getattr(existing, 'label', ''), 'dedup-channel-conflict')
        if key not in _unzip_inherit_warned:
            _unzip_inherit_warned.add(key)
            logging.warning(
                "Polymer '%s': fingerprint-deduped equivalent arrived with a "
                "DIFFERENT radical_qssa_unzip channel; keeping the existing "
                "channel (first-writer-wins, disclosed).",
                getattr(existing, 'label', ''))


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
        # Honest-empty seeding (item #14a, amended 2026-06-12 uniform-t=0):
        # a just-spawned daughter genuinely contains nothing, and the
        # artifact's pools[].moments are INITIAL CONDITIONS at t=0 of the
        # simulated experiment — so mu = [0, 0, 0] is the physically correct
        # seed, exactly unifying with the Path-C scission-tail convention
        # (initial_mass=0 -> zero moments). Nothing recomputes moments at
        # emission: the sidecar passes pool.moments through verbatim. NOTE:
        # the Polymer(...) constructor above derives interim moments from
        # Mn/Mw/initial_mass=0.001 (spec section-7 ruling: KEPT — the
        # constant parameterizes only that derivation and never reaches the
        # artifact); this assignment overrides them. Mn/Mw stay the
        # PARENT's, as lineage metadata only (no DP is derived from them).
        new_pool.moments = np.zeros(3, dtype=np.float64)
        new_pool.parent_pool_label = parent.label
        # Daughter-pool channel inheritance (radical_qssa_unzip M5), gated on
        # parentage certainty (round-25 P1-2): intent.parent_pool is the
        # pool_registry[0] attribution shortcut, so the channel transfers
        # only when the intent's TRUE motif is the parent's own monomer
        # chemistry; otherwise the pool spawns channel-free with a
        # once-per-pool WARNING.
        _inherit_spawned_pool_channel(new_pool, parent, intent)
        new_pool.spawn_iteration = iteration
        new_pool.end_groups_str = list(intent.end_groups)
        new_pool.mu_indices = (next_idx, next_idx + 1, next_idx + 2)
        next_idx += 3
        new_pool.spawn_metadata = {
            # triggering_dp is motif METADATA (repeats per proxy-scale
            # triggering molecule — repeats-per-chain is genuinely useful);
            # it is NEVER a moment multiplier (item #14a section 2: a
            # quantity measured on a proxy-scale object is representation,
            # not chemistry — nothing multiplies it anywhere).
            # triggering_moles is deleted outright: never emitted by any
            # real deck (clean-delete compat ruling, spec section 2).
            "triggering_dp": int(intent.triggering_dp),
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


def clear_polymer_proxy(obj: Union['Species', Molecule]) -> None:
    """Clear the ``is_polymer_proxy`` flag on a Molecule or Species.

    Inverse of :func:`_tag_polymer_proxy`. The reaction-generation handshake
    (``family.py:1665``) blanket-stamps every product ``is_polymer_proxy=True``
    when any reactant is a proxy (the polymer pool always is). For a product
    that the polymer machinery does NOT convert to a :class:`Polymer` -- i.e. a
    genuine discrete gas-phase species or a demoted discrete chip -- that stamp
    is stale and would wrongly count the species as a melt reference-state
    participant in the solver. This helper scrubs it: the attribute, the
    ``props`` entry (if a ``props`` dict exists), and the same on every
    constituent :class:`Molecule` when ``obj`` is a Species.

    Order matters: clear the constituent molecules FIRST, then the
    object-level flag LAST. :attr:`Species.is_polymer_proxy` is a sticky
    lazy-cache property (``species.py``): its setter calls
    ``propagate_polymer_proxy_to_molecules``, whose getter re-derives ``True``
    from ANY still-proxy molecule and re-caches it. Clearing the species flag
    before the molecules would let that getter re-stamp the cache ``True`` off
    the not-yet-cleared molecules, leaving the species ``_is_polymer_proxy``
    stuck ``True`` (the molecules then go ``False`` but the cached species flag
    never re-clears) -- which ``make_new_species`` (model.py) then ORs onto the
    solver-visible Species. Molecules-first lets the final setter settle the
    cache to ``False``.
    """
    if getattr(obj, "molecule", None):
        for m in obj.molecule:
            m.is_polymer_proxy = False
            if isinstance(getattr(m, "props", None), dict):
                m.props["is_polymer_proxy"] = False
    obj.is_polymer_proxy = False
    if isinstance(getattr(obj, "props", None), dict):
        obj.props["is_polymer_proxy"] = False


#: props key carrying the DURABLE "this is a genuine discrete gas volatile"
#: verdict. Unlike :attr:`is_polymer_proxy` -- a monotonic multi-writer sticky
#: cache re-stamped by ``family.py:1665`` and the ``species.py`` sticky getter,
#: with no authoritative "gas" clear point -- this is a POSITIVE veto set ONCE
#: at the discrete-product creation point (the polymer handshake / chip
#: demotion) and never touched by the proxy stamping machinery. It lives in
#: ``props`` because ``Species.copy`` deep-copies ``props`` (species.py) while
#: ad-hoc attributes and ``Molecule.props`` are NOT preserved across copies.
#: The solver reference-state melt gate (``polymer.pyx``) reads it: a species
#: is a melt tag-branch participant only if proxy AND chain-scale MW AND NOT
#: this veto -- so a genuine gas volatile that got proxy-contaminated is
#: correctly excluded from the melt reference-state sum.
POLYMER_REFERENCE_STATE_GAS_VETO_KEY = "polymer_reference_state_gas_veto"

#: Species-level heavy-atom (non-H) count carried in ``props`` for species
#: whose STRUCTURE does not cross a data boundary (r89 dual-axis melt gate).
#: Consumer-world species (rmgpy/tools/polymer_moments_runner.py
#: ``_species_from_yaml``) are label + thermo only (``molecule == []``); their
#: MW is reconstructed from the chem.yaml elemental composition, and the
#: solver reference-state melt gate's HEAVY-ATOM axis
#: (``polymer.pyx _dual_axis_polymer_sized``, mirroring
#: :func:`_discrete_is_polymer_sized`) needs the heavy count from the same
#: composition -- the r89 adjudication forbids approximating the heavy axis
#: from mass. Lives in ``props`` for the same reason as the gas veto above:
#: ``Species.copy`` deep-copies ``props`` while ad-hoc attributes are lost.
#: The solver reads the literal (pinned alongside the veto-key literal); a
#: species with structure never needs it (the gate computes the count from
#: ``molecule[0]`` first). Absent + structureless => the heavy axis is
#: UNDECIDABLE and the gate answers conservative-gas with a census warning.
POLYMER_HEAVY_ATOM_COUNT_KEY = "polymer_heavy_atom_count"


def set_polymer_gas_veto(obj: Union['Species', Molecule]) -> None:
    """Stamp the durable gas-volatile veto on a Molecule or Species.

    See :data:`POLYMER_REFERENCE_STATE_GAS_VETO_KEY`. Sets the veto in ``props``
    on ``obj`` and -- when ``obj`` is a Species -- on every constituent
    :class:`Molecule` too (belt-and-suspenders for the make_new_species path,
    which may read the verdict off either the incoming object or its molecule).
    The Species-level ``props`` entry is the load-bearing one: it is the copy of
    the verdict that survives ``Species.copy``.
    """
    if not isinstance(getattr(obj, "props", None), dict):
        obj.props = {}
    obj.props[POLYMER_REFERENCE_STATE_GAS_VETO_KEY] = True
    if getattr(obj, "molecule", None):
        for m in obj.molecule:
            if not isinstance(getattr(m, "props", None), dict):
                m.props = {}
            m.props[POLYMER_REFERENCE_STATE_GAS_VETO_KEY] = True


def strip_rmg_index_suffix(label: str) -> str:
    """Canonical base-label convention (PP run-5 erratum, 2026-07-05): strip
    ONLY a trailing RMG index suffix ``(<int>)`` from a species label
    (``'PS(2)' -> 'PS'``, ``'C[CH]CC(C)C(6)' -> 'C[CH]CC(C)C'``); any other
    label -- including SMILES-derived labels whose parentheses are BRANCHING
    syntax, e.g. ``'C[CH]CC(C)C'`` -- is returned unchanged.

    This replaces the former ``label.partition('(')[0]`` convention
    everywhere base labels are produced or consumed (polymer_input's
    ``_base_label``, the solver's pool binding and edge-daughter condensed
    application, ``_species_base_label`` below, and
    ``derive_daughter_pool_configs``). Truncating at the FIRST ``'('``
    aliased structurally different SMILES-labelled species onto one base:
    a genuine C9H19 H-loss edge daughter ``'C[CH]CC(C)CC(C)C'`` registered
    base ``'C[CH]CC'``, which the get_gas_mask section-D application then
    matched against the unrelated CORE C6H13 tail ``'C[CH]CC(C)C'`` -- but
    only in the combined chain(core, edge) call, so the prospective mask
    core prefix diverged from gas_species_mask and RIDER R1 raised (the PP
    run-5 PROSPECTIVE-MASK TRIPWIRE crash). A parenthesised group of bare
    digits is never valid SMILES, so trailing-index stripping cannot
    truncate a structural label.

    ONE convention, never duplicated: every producer and consumer of base
    labels must call this (directly or via a thin delegate) so set
    membership stays lockstep across module boundaries.
    """
    if label and label.endswith(')'):
        head, sep, tail = label[:-1].rpartition('(')
        if sep and tail.isdigit():
            return head
    return label


def is_h_loss_radical_daughter(molecule: Molecule, proxy_element_counts) -> bool:
    """Structural core of the H-loss radical-daughter qualifier (ratified
    2026-07-03; veto-scoping extraction 2026-07-04). Returns True iff
    ``molecule`` is an H-loss radical daughter of a proxy whose element
    counts appear in ``proxy_element_counts`` (an iterable of
    ``Molecule.get_element_count()`` dicts):

      (ii)  radical-bearing and neutral;
      (iii) same non-H element composition (heavy-atom skeleton) as a proxy;
      (iv)  differs from that proxy by H-loss only
            (H_proxy - H_daughter == its radical count >= 1,
            e.g. C9H20 -> C9H19*).

    ONE predicate, never duplicated -- shared by BOTH consumers:

    * ``PolymerPhase.get_h_loss_radical_daughter_bases``
      (rmgpy/rmg/polymer_input.py), which adds condition (v)
      (not durably gas-vetoed) and the base-label bookkeeping;
    * the handshake veto scoping in ``_handshake_structures``
      (rmgpy/data/kinetics/family.py), which adds the chain-window MW
      conjunct (MW >= monomer_mw + slack) before EXEMPTING a refused
      product from the durable gas veto -- without that exemption the veto
      defeats condition (v) of the very qualifier above (the PP run-2
      self-defeat loop, gate/conduit diagnosis Phenomenon 1).

    Conditions (ii)-(iv) fail closed on any structure query error.
    """
    try:
        n_rad = molecule.get_radical_count()
        net_charge = molecule.get_net_charge()
        comp = molecule.get_element_count()
    except Exception:
        return False
    if n_rad < 1 or net_charge != 0:                # (ii)
        return False
    heavy = {el: n for el, n in comp.items() if el != 'H'}
    n_h = comp.get('H', 0)
    for pcomp in proxy_element_counts:
        p_heavy = {el: n for el, n in pcomp.items() if el != 'H'}
        # (iii) same heavy skeleton + (iv) H-loss only, one H per radical site
        if p_heavy == heavy and pcomp.get('H', 0) - n_h == n_rad:
            return True
    return False


def compute_h_loss_feature_verdicts(reactants, products, polymer_reactants):
    """Per-product H-loss conduit verdicts (stage S2, feature-pool conduit
    arc, adjudicated design): a list of bools parallel to ``products``,
    ``True`` iff that product may route through the radical-feature producer
    path (``create_reacted_copy(..., h_loss_feature=True)``).

    Computed at the ONE place where resolved reactants and raw products are
    both visible (the ``make_new_reaction`` handshake call site in
    rmgpy/rmg/model.py) -- NEVER inferred inside :class:`Polymer` from a
    single molecule. A product's verdict is True iff ALL of:

    * exactly ONE polymer reactant source (``len(polymer_reactants) == 1``);
    * the product is a same-heavy-skeleton SINGLE-H-loss radical daughter of
      that source's reacting proxy (:func:`is_h_loss_radical_daughter` with
      radical count exactly 1);
    * abstraction co-product EVIDENCE: the non-polymer side of the reaction
      (all non-polymer participants except this daughter) nets a gain of
      exactly ONE hydrogen and no heavy atoms -- the H atom / H2 / RH
      co-product actually carries the abstracted H, not just a formula
      coincidence. Fails closed if any non-polymer participant cannot be
      weighed;
    * the daughter radical is QSSA-eliminating
      (:func:`is_qssa_eliminating_radical`); accumulating
      (resonance-stabilized) daughters stay refused and pick up the
      ``qssa-invalid`` refuse stamp downstream.
    """
    verdicts = [False] * len(products)
    if len(polymer_reactants) != 1:
        return verdicts
    proxy_mols = getattr(polymer_reactants[0], 'molecule', None) or []
    if not proxy_mols or proxy_mols[0] is None:
        return verdicts
    try:
        proxy_comp = proxy_mols[0].get_element_count()
    except Exception:
        return verdicts

    def _participant_mol(item):
        """Molecule of a non-polymer participant; None for Polymers (pool
        chains do not count) and a ValueError for unweighable non-polymers
        (fail closed)."""
        if isinstance(item, Polymer) or getattr(item, 'is_polymer', False):
            return None
        if isinstance(item, Molecule):
            return item
        mols = getattr(item, 'molecule', None)
        if mols and mols[0] is not None:
            return mols[0]
        raise ValueError("unweighable non-polymer participant")

    for i, item in enumerate(products):
        try:
            mol = _participant_mol(item)
        except ValueError:
            continue
        if mol is None:
            continue
        try:
            if mol.get_radical_count() != 1:
                continue
        except Exception:
            continue
        if not is_h_loss_radical_daughter(mol, [proxy_comp]):
            continue
        # Abstraction co-product evidence: net element delta of the
        # NON-POLYMER side excluding this daughter must be exactly +1 H.
        delta: Dict[str, int] = {}
        try:
            for j, other in enumerate(products):
                if j == i:
                    continue
                omol = _participant_mol(other)
                if omol is None:
                    continue
                for el, n in omol.get_element_count().items():
                    delta[el] = delta.get(el, 0) + n
            for r in reactants:
                rmol = _participant_mol(r)
                if rmol is None:
                    continue
                for el, n in rmol.get_element_count().items():
                    delta[el] = delta.get(el, 0) - n
        except Exception:
            continue
        if any(n != 0 for el, n in delta.items() if el != 'H'):
            continue
        if delta.get('H', 0) != 1:
            continue
        if not is_qssa_eliminating_radical(mol):
            continue
        verdicts[i] = True
    return verdicts


def has_polymer_gas_veto(obj) -> bool:
    """Return True if ``obj`` (a Molecule or Species) carries the durable
    gas-volatile veto in ``props`` (or, for a Species, on any of its
    molecules). See :data:`POLYMER_REFERENCE_STATE_GAS_VETO_KEY`."""
    props = getattr(obj, "props", None)
    if isinstance(props, dict) and props.get(POLYMER_REFERENCE_STATE_GAS_VETO_KEY):
        return True
    for m in getattr(obj, "molecule", None) or []:
        mprops = getattr(m, "props", None)
        if isinstance(mprops, dict) and mprops.get(POLYMER_REFERENCE_STATE_GAS_VETO_KEY):
            return True
    return False


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
