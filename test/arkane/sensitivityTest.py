#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the :mod:`arkane.sensitivity` module.

`PDepSensitivity` normally requires a fully executed :class:`PressureDependenceJob`
to instantiate, so these tests exercise its `perturb` method as an unbound method
against a lightweight stub object that only exposes the attributes `perturb` reads:
`self.job.network.path_reactions` and `self.perturbation`.
"""

import logging

import pytest

import rmgpy.quantity as quantity
from rmgpy.kinetics import Arrhenius, MultiArrhenius, Chebyshev
from rmgpy.reaction import Reaction
from rmgpy.species import Species, TransitionState
from rmgpy.statmech import Conformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator

from arkane.sensitivity import PDepSensitivity, SYNTHETIC_TS_E0_TOLERANCE, _ts_e0_is_synthetic


class _Network(object):
    def __init__(self, path_reactions):
        self.path_reactions = path_reactions


class _Job(object):
    def __init__(self, path_reactions):
        self.network = _Network(path_reactions)


class _StubPDepSensitivity(object):
    """A minimal stand-in exposing only what `PDepSensitivity.perturb` reads."""

    def __init__(self, path_reactions, perturbation):
        self.job = _Job(path_reactions)
        self.perturbation = perturbation


def _make_modeless_ts(e0_value_si):
    """A TS with E0 only, and no statmech modes -- i.e. `can_tst()` is False."""
    conformer = Conformer(E0=(e0_value_si, 'J/mol'), modes=[])
    return TransitionState(label='TS', conformer=conformer)


def _make_full_ts(e0_value_si):
    """A TS with real statmech modes -- i.e. `can_tst()` is True."""
    modes = [
        IdealGasTranslation(mass=(28.0, 'amu')),
        NonlinearRotor(inertia=([1.0, 2.0, 3.0], 'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([500.0, 1000.0], 'cm^-1')),
    ]
    conformer = Conformer(E0=(e0_value_si, 'J/mol'), modes=modes)
    return TransitionState(label='TS', conformer=conformer)


def _make_species_with_e0(label, e0_value_si):
    """A Species with only a conformer E0 set (as `_ts_e0_is_synthetic` needs)."""
    return Species(label=label, conformer=Conformer(E0=(e0_value_si, 'J/mol'), modes=[]))


class TestPDepSensitivityPerturb:
    """
    Unit tests for :meth:`PDepSensitivity.perturb`, focused on the ILT-based path
    reaction case where perturbing the TS's E0 alone is a structural no-op for the
    resulting microcanonical rate, and the reaction's Arrhenius `Ea` must also move.
    """

    def setup_method(self):
        self.perturbation = quantity.Quantity(1, 'kcal/mol')

    def test_ilt_reaction_perturbs_e0_and_ea(self):
        """For an ILT-based path reaction (no TS modes, n >= 0.25), both E0 and Ea should move."""
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        # E0(TS) = sum(E0(reactants)) + Ea, so the synthetic relation holds and Ea is perturbed too.
        reactant = _make_species_with_e0('A', e0_0 - ea_0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        rxn = Reaction(reactants=[reactant], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()
        assert kinetics.n.value_si >= 0.25

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + p)
        assert kinetics.Ea.value_si == pytest.approx(ea_0 + p)

        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert ts.conformer.E0.value_si == pytest.approx(e0_0)
        assert kinetics.Ea.value_si == pytest.approx(ea_0)

    def test_rrkm_reaction_does_not_perturb_ea(self):
        """For an RRKM-based path reaction (TS has modes), Ea must not be touched, only E0."""
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_full_ts(e0_0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        assert rxn.can_tst()

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + p)
        assert kinetics.Ea.value_si == pytest.approx(ea_0)  # unchanged

        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert ts.conformer.E0.value_si == pytest.approx(e0_0)

    def test_network_kinetics_is_perturbed_when_present(self):
        """When `network_kinetics` is set, it (not `kinetics`) must be the one perturbed."""
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        network_ea_0 = 30000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        # Chosen so E0(TS) = sum(E0(reactants)) + network_ea_0 (the kinetics actually used for k(E)).
        reactant = _make_species_with_e0('A', e0_0 - network_ea_0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        network_kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(network_ea_0, 'J/mol'))
        rxn = Reaction(reactants=[reactant], products=[Species(label='B')], transition_state=ts,
                        kinetics=kinetics, network_kinetics=network_kinetics)
        assert not rxn.can_tst()

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert network_kinetics.Ea.value_si == pytest.approx(network_ea_0 + p)
        assert kinetics.Ea.value_si == pytest.approx(ea_0)  # the non-selected kinetics is untouched

    def test_multi_arrhenius_falls_through_to_unsupported_kinetics_warning(self, caplog):
        """
        `MultiArrhenius`-kinetics path reactions cannot genuinely reach the ILT sensitivity path:
        `apply_inverse_laplace_transform_method` (rmgpy/pdep/reaction.pyx) takes a Cython-typed
        `Arrhenius kinetics` argument, and `MultiArrhenius` is a distinct `KineticsModel` subclass
        (not an `Arrhenius` subclass), so passing a `MultiArrhenius` reaction's kinetics into it
        raises a `TypeError` before Arkane's sensitivity code is ever reached (verified directly
        against the compiled extension). Since a real `MultiArrhenius` ILT path reaction can never
        reach `perturb`, it is intentionally NOT given a dedicated perturbation branch here (that
        would have been an unverified vector-perturbation claim); it must fall through to the
        same unsupported-kinetics warning as any other unhandled kinetics type, with only E0
        perturbed.
        """
        e0_0 = 10000.0  # J/mol
        ea_0a, ea_0b = 15000.0, 25000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        arrh_a = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0a, 'J/mol'))
        arrh_b = Arrhenius(A=(1e10, 's^-1'), n=2.0, Ea=(ea_0b, 'J/mol'))
        kinetics = MultiArrhenius(arrhenius=[arrh_a, arrh_b])
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        with caplog.at_level(logging.WARNING):
            PDepSensitivity.perturb(stub, ts)

        assert any('MultiArrhenius' in record.message for record in caplog.records)
        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + p)
        assert arrh_a.Ea.value_si == pytest.approx(ea_0a)  # untouched
        assert arrh_b.Ea.value_si == pytest.approx(ea_0b)  # untouched

    def test_shared_transition_state_perturbs_ea_for_every_owning_reaction(self):
        """
        Arkane input files may point more than one `reaction` block at the same `transitionState`
        label, so two path reactions can share a single `TransitionState` object. The E0
        perturbation above already reaches every reaction using that shared object, so the Ea
        perturbation must match that same blast radius: previously a `break` after the first
        matching reaction meant only one of two owning reactions' Ea moved, desynchronizing it
        from the (already shared) E0 perturbation.
        """
        e0_0 = 10000.0  # J/mol
        ea_0_1 = 20000.0  # J/mol
        ea_0_2 = 5000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        reactant_1 = _make_species_with_e0('A', e0_0 - ea_0_1)
        reactant_2 = _make_species_with_e0('C', e0_0 - ea_0_2)
        kinetics_1 = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0_1, 'J/mol'))
        kinetics_2 = Arrhenius(A=(1e12, 's^-1'), n=1.2, Ea=(ea_0_2, 'J/mol'))
        rxn_1 = Reaction(reactants=[reactant_1], products=[Species(label='B')],
                          transition_state=ts, kinetics=kinetics_1)
        rxn_2 = Reaction(reactants=[reactant_2], products=[Species(label='D')],
                          transition_state=ts, kinetics=kinetics_2)
        assert not rxn_1.can_tst() and not rxn_2.can_tst()
        assert rxn_1.transition_state is rxn_2.transition_state is ts

        stub = _StubPDepSensitivity([rxn_1, rxn_2], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + p)
        assert kinetics_1.Ea.value_si == pytest.approx(ea_0_1 + p)
        assert kinetics_2.Ea.value_si == pytest.approx(ea_0_2 + p)

        # Round-trip: unperturb must restore both Ea values (and E0) exactly.
        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert ts.conformer.E0.value_si == pytest.approx(e0_0, abs=1e-9)
        assert kinetics_1.Ea.value_si == pytest.approx(ea_0_1, abs=1e-9)
        assert kinetics_2.Ea.value_si == pytest.approx(ea_0_2, abs=1e-9)

    def test_modeless_ts_with_independent_e0_is_not_perturbed_for_ea(self, caplog):
        """
        A hand-authored, modeless TS whose E0 is an independent physical value (i.e. it does NOT
        satisfy E0(TS) = sum(E0(reactants)) + Ea, nor the products-side equivalent) must have only
        its E0 perturbed; moving Ea as well would not be physically justified for such a TS.
        """
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        # Reactant/product E0s chosen so neither side comes anywhere near satisfying the relation.
        reactant = _make_species_with_e0('A', 500000.0)
        product = _make_species_with_e0('B', 600000.0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        rxn = Reaction(reactants=[reactant], products=[product], transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()
        assert not _ts_e0_is_synthetic(rxn, ts.conformer.E0.value_si, kinetics.Ea.value_si)

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        with caplog.at_level(logging.INFO):
            PDepSensitivity.perturb(stub, ts)

        assert any('independent' in record.message for record in caplog.records)
        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + p)
        assert kinetics.Ea.value_si == pytest.approx(ea_0)  # untouched

        # Round-trip: unperturb must restore E0 exactly, and Ea must remain untouched throughout,
        # bit-for-bit (the same verdict -- "not synthetic" -- must be reached on both calls).
        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert ts.conformer.E0.value_si == e0_0
        assert kinetics.Ea.value_si == ea_0

    def test_ts_e0_is_synthetic_tries_products_side_when_reactants_side_fails(self):
        """
        A path reaction may be stored in either direction relative to which side the TS sits on;
        `_ts_e0_is_synthetic` must try the products-side sum when the reactants-side sum does not
        satisfy the relation.
        """
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        reactant = _make_species_with_e0('A', 999999.0)  # reactants-side check fails
        product = _make_species_with_e0('B', e0_0 - ea_0)  # products-side check matches
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        rxn = Reaction(reactants=[reactant], products=[product], transition_state=ts, kinetics=kinetics)

        assert _ts_e0_is_synthetic(rxn, ts.conformer.E0.value_si, kinetics.Ea.value_si)

    def test_ts_e0_is_synthetic_tolerance_boundary(self):
        """Sanity-check the module-level tolerance constant is applied as documented."""
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        just_inside = _make_species_with_e0('A', e0_0 - ea_0 + 0.9 * SYNTHETIC_TS_E0_TOLERANCE)
        just_outside = _make_species_with_e0('A', e0_0 - ea_0 + 1.1 * SYNTHETIC_TS_E0_TOLERANCE)
        rxn_inside = Reaction(reactants=[just_inside], products=[Species(label='B')],
                               transition_state=ts, kinetics=kinetics)
        rxn_outside = Reaction(reactants=[just_outside], products=[Species(label='B')],
                                transition_state=ts, kinetics=kinetics)

        assert _ts_e0_is_synthetic(rxn_inside, ts.conformer.E0.value_si, kinetics.Ea.value_si)
        assert not _ts_e0_is_synthetic(rxn_outside, ts.conformer.E0.value_si, kinetics.Ea.value_si)

    def test_unsupported_kinetics_type_logs_warning(self, caplog):
        """An ILT-based reaction with an unsupported kinetics type should log a warning, not fail silently."""
        e0_0 = 10000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        kinetics = Chebyshev()
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        with caplog.at_level(logging.WARNING):
            PDepSensitivity.perturb(stub, ts)

        assert any('Chebyshev' in record.message for record in caplog.records)
        # E0 is still perturbed even though Ea perturbation was skipped
        assert ts.conformer.E0.value_si == pytest.approx(e0_0 + self.perturbation.value_si)
