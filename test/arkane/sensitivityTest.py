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

from arkane.sensitivity import PDepSensitivity


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
    """A Species with only a conformer E0 set."""
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
        reactant = _make_species_with_e0('A', 5000.0)
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
        reactant = _make_species_with_e0('A', 5000.0)
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

    def test_ilt_reaction_with_unrelated_ts_e0_still_perturbs_ea(self):
        """
        Regression test built from a real RMG-written network (T3 PDep network13_2,
        CH2CNH <=> CH3CN): the written network file carries a TS E0 that does NOT satisfy
        E0(TS) = sum(E0(reactants)) + Ea on either side (the relation misses by 137 and
        23 kJ/mol respectively -- RMG mutates Ea via `fix_barrier_height` after synthesizing
        the TS E0, and Arkane strips `energy_correction` when writing the network file, so
        the relation is generally unrecoverable from the file). An earlier revision gated the
        Ea perturbation on that arithmetic relation, which only 24 of 6,434 resolvable path
        reactions across 400 real RMG networks satisfied (0.4%); real networks therefore got a
        structurally zero TS sensitivity coefficient. The Ea perturbation must apply whenever
        the reaction is ILT-based with Arrhenius kinetics, regardless of any relation between
        the TS E0 and the species energies.
        """
        e0_ts = 331959.0  # J/mol, from the real network file
        ea_0 = 294135.2  # J/mol
        ts = _make_modeless_ts(e0_ts)
        reactant = _make_species_with_e0('CH2CNH', 174694.0)
        product = _make_species_with_e0('CH3CN', 61223.1)
        kinetics = Arrhenius(A=(2.5e13, 's^-1'), n=0.0, Ea=(ea_0, 'J/mol'), T0=(1, 'K'))
        rxn = Reaction(reactants=[reactant], products=[product], transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()
        # The defining property of this regression case: the old synthetic-E0 relation fails
        # on both sides, by far more than any plausible tolerance.
        for side in (rxn.reactants, rxn.products):
            e0_sum = sum(spc.conformer.E0.value_si for spc in side)
            assert abs((e0_ts - ea_0) - e0_sum) > 20000.0

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert ts.conformer.E0.value_si == pytest.approx(e0_ts + p)
        assert kinetics.Ea.value_si == pytest.approx(ea_0 + p)  # perturbed despite the broken relation

        # Round-trip: unperturb must restore both exactly (the condition gating the Ea
        # perturbation is structural -- can_tst() and the kinetics type -- so it is invariant
        # across the perturb/unperturb pair by construction).
        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert ts.conformer.E0.value_si == pytest.approx(e0_ts, abs=1e-9)
        assert kinetics.Ea.value_si == pytest.approx(ea_0, abs=1e-9)

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

    def test_perturbing_one_ts_does_not_perturb_other_ts_reactions(self):
        """
        A network with two DISTINCT transition states, each owning its own ILT-based reaction:
        perturbing TS A must move only TS A's reaction Ea, leaving TS B's reaction Ea and E0 fully
        untouched. This pins the `rxn.transition_state is entry` scoping in the Ea loop -- without
        it, perturbing one TS would shift Ea for every ILT reaction anywhere in the network and
        silently corrupt every other TS's sensitivity coefficient computed in the same job.
        """
        e0_a, e0_b = 10000.0, 12000.0  # J/mol
        ea_a, ea_b = 20000.0, 30000.0  # J/mol
        ts_a = _make_modeless_ts(e0_a)
        ts_b = _make_modeless_ts(e0_b)
        kinetics_a = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_a, 'J/mol'))
        kinetics_b = Arrhenius(A=(1e12, 's^-1'), n=1.2, Ea=(ea_b, 'J/mol'))
        rxn_a = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                          transition_state=ts_a, kinetics=kinetics_a)
        rxn_b = Reaction(reactants=[Species(label='C')], products=[Species(label='D')],
                          transition_state=ts_b, kinetics=kinetics_b)
        assert ts_a is not ts_b
        assert not rxn_a.can_tst() and not rxn_b.can_tst()

        stub = _StubPDepSensitivity([rxn_a, rxn_b], self.perturbation)
        PDepSensitivity.perturb(stub, ts_a)

        p = self.perturbation.value_si
        assert ts_a.conformer.E0.value_si == pytest.approx(e0_a + p)
        assert kinetics_a.Ea.value_si == pytest.approx(ea_a + p)
        # TS B's reaction is entirely untouched.
        assert ts_b.conformer.E0.value_si == pytest.approx(e0_b)
        assert kinetics_b.Ea.value_si == pytest.approx(ea_b)

    def test_shared_kinetics_object_is_perturbed_only_once(self):
        """
        If two path reactions of the same TS share a single Arrhenius object, it represents one
        barrier: the E0 shift already touches it once, so its Ea must move by exactly the
        perturbation, not by twice it (once per owning reaction). The Ea loop keys on kinetics-object
        identity to guard this.
        """
        e0_0 = 10000.0  # J/mol
        ea_0 = 20000.0  # J/mol
        ts = _make_modeless_ts(e0_0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(ea_0, 'J/mol'))
        rxn_1 = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                          transition_state=ts, kinetics=kinetics)
        rxn_2 = Reaction(reactants=[Species(label='C')], products=[Species(label='D')],
                          transition_state=ts, kinetics=kinetics)
        assert rxn_1.kinetics is rxn_2.kinetics is kinetics
        assert not rxn_1.can_tst() and not rxn_2.can_tst()

        stub = _StubPDepSensitivity([rxn_1, rxn_2], self.perturbation)
        PDepSensitivity.perturb(stub, ts)

        p = self.perturbation.value_si
        assert kinetics.Ea.value_si == pytest.approx(ea_0 + p)  # shifted once, not 2*p

        PDepSensitivity.perturb(stub, ts, unperturb=True)
        assert kinetics.Ea.value_si == pytest.approx(ea_0, abs=1e-9)

    def test_kinetics_shared_across_transition_states_warns(self, caplog):
        """
        Two ILT reactions with DIFFERENT transition states that share one Arrhenius object is an
        unsupported input: perturbing TS A shifts the shared Ea, collaterally perturbing TS B's
        reaction rate and contaminating its sensitivity row. perturb() cannot undo the aliasing, but
        it must surface it loudly rather than report a silently corrupted coefficient.
        """
        ts_a = _make_modeless_ts(10000.0)
        ts_b = _make_modeless_ts(12000.0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(20000.0, 'J/mol'))
        rxn_a = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                          transition_state=ts_a, kinetics=kinetics)
        rxn_b = Reaction(reactants=[Species(label='C')], products=[Species(label='D')],
                          transition_state=ts_b, kinetics=kinetics)
        assert rxn_a.kinetics is rxn_b.kinetics is kinetics
        assert ts_a is not ts_b

        stub = _StubPDepSensitivity([rxn_a, rxn_b], self.perturbation)
        with caplog.at_level(logging.WARNING):
            PDepSensitivity.perturb(stub, ts_a)

        assert any('shares its Arrhenius kinetics object' in record.message for record in caplog.records)

    def test_unsupported_kinetics_warning_not_duplicated_on_unperturb(self, caplog):
        """
        The unsupported-kinetics warning must fire once per perturbation, not once on the perturb
        call and again on the matching unperturb call (which re-runs the same loop). It is guarded by
        `not unperturb` for exactly this reason.
        """
        ts = _make_modeless_ts(10000.0)
        kinetics = Chebyshev()
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        assert not rxn.can_tst()

        stub = _StubPDepSensitivity([rxn], self.perturbation)
        with caplog.at_level(logging.WARNING):
            PDepSensitivity.perturb(stub, ts)
            PDepSensitivity.perturb(stub, ts, unperturb=True)

        n_warnings = sum('Not perturbing Ea' in record.message for record in caplog.records)
        assert n_warnings == 1


class TestClassifyIltTransitionStates:
    """
    Unit tests for :meth:`PDepSensitivity._classify_ilt_transition_states`, the single source that
    save() and plot() use to decide which TS rows carry the combined-E0+Ea (ILT) semantics and which
    are contaminated by an Arrhenius kinetics object shared across transition states.
    """

    def setup_method(self):
        self.perturbation = quantity.Quantity(1, 'kcal/mol')

    def test_arrhenius_ilt_ts_is_classified_ilt_not_contaminated(self):
        ts = _make_modeless_ts(10000.0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(20000.0, 'J/mol'))
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        stub = _StubPDepSensitivity([rxn], self.perturbation)
        ilt, contaminated = PDepSensitivity._classify_ilt_transition_states(stub)
        assert ilt == ['TS']
        assert contaminated == []

    def test_rrkm_ts_is_not_classified_ilt(self):
        ts = _make_full_ts(10000.0)  # has statmech modes -> can_tst() is True
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(20000.0, 'J/mol'))
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=kinetics)
        stub = _StubPDepSensitivity([rxn], self.perturbation)
        ilt, contaminated = PDepSensitivity._classify_ilt_transition_states(stub)
        assert ilt == []
        assert contaminated == []

    def test_non_arrhenius_ilt_ts_is_not_classified_ilt(self):
        """
        The classification gate must match perturb()'s Ea shift: a non-Arrhenius ILT reaction has
        its Ea left untouched, so its row is a plain E0-only derivative and must NOT be marked ILT --
        otherwise the YAML/header/plot would over-promise the combined-E0+Ea semantics for a row that
        never got the Ea perturbation.
        """
        ts = _make_modeless_ts(10000.0)
        rxn = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                        transition_state=ts, kinetics=Chebyshev())
        stub = _StubPDepSensitivity([rxn], self.perturbation)
        ilt, contaminated = PDepSensitivity._classify_ilt_transition_states(stub)
        assert ilt == []
        assert contaminated == []

    def test_kinetics_shared_across_distinct_ts_is_contaminated(self):
        """One Arrhenius object shared across two DIFFERENT transition states contaminates both."""
        ts_a = _make_modeless_ts(10000.0)
        ts_b = _make_modeless_ts(12000.0)
        ts_a.label, ts_b.label = 'TS_A', 'TS_B'
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(20000.0, 'J/mol'))
        rxn_a = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                          transition_state=ts_a, kinetics=kinetics)
        rxn_b = Reaction(reactants=[Species(label='C')], products=[Species(label='D')],
                          transition_state=ts_b, kinetics=kinetics)
        stub = _StubPDepSensitivity([rxn_a, rxn_b], self.perturbation)
        ilt, contaminated = PDepSensitivity._classify_ilt_transition_states(stub)
        assert sorted(ilt) == ['TS_A', 'TS_B']
        assert sorted(contaminated) == ['TS_A', 'TS_B']

    def test_kinetics_shared_within_one_ts_is_not_contaminated(self):
        """An Arrhenius object shared across reactions of the SAME TS is one barrier, not contamination."""
        ts = _make_modeless_ts(10000.0)
        kinetics = Arrhenius(A=(1e13, 's^-1'), n=1.5, Ea=(20000.0, 'J/mol'))
        rxn_1 = Reaction(reactants=[Species(label='A')], products=[Species(label='B')],
                          transition_state=ts, kinetics=kinetics)
        rxn_2 = Reaction(reactants=[Species(label='C')], products=[Species(label='D')],
                          transition_state=ts, kinetics=kinetics)
        stub = _StubPDepSensitivity([rxn_1, rxn_2], self.perturbation)
        ilt, contaminated = PDepSensitivity._classify_ilt_transition_states(stub)
        assert ilt == ['TS']
        assert contaminated == []
