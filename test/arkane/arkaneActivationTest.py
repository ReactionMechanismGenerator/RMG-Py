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
Unit tests of the optional Arkane activation block (:mod:`arkane.activation`).

Hand values are from a reviewed prototype of this block, computed for the shipped examples
``examples/arkane/reactions/H+C2H4=C2H5`` (dn = -1) and ``NH2+N2H3=NH+N2H4`` (dn = 0) at 298.15 K.
"""

import math
import os
import re
import shutil
import sys

import pytest

import rmgpy.constants as constants
from rmgpy.exceptions import InputError, ReactionError
from rmgpy.statmech.translation import IdealGasTranslation

import arkane.activation
from arkane.activation import (KC_BASIS_FACTOR, P_STATMECH, STANDARD_STATES, compute_activation,
                               format_activation_block, k_in_state_units, standard_concentration,
                               validate_standard_state)
from arkane.kinetics import KineticsJob
from arkane.main import Arkane

T298 = 298.15
STATES = ('1 atm', '1 bar', '1 M')
EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                        'examples', 'arkane')
GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'activation',
                      'H_C2H4_C2H5_output.py.golden')
H_C2H4 = 'H+C2H4=C2H5'
NH2_N2H3 = 'NH2+N2H3=NH+N2H4'
EXAMPLE_SPECIES = {H_C2H4: ('H', 'C2H4', 'C2H5'), NH2_N2H3: ()}


def _copy_example(name, root):
    """Copy a shipped reaction example and the species it reads to `root`, keeping relative paths."""
    shutil.copytree(os.path.join(EXAMPLES, 'reactions', name), os.path.join(root, 'reactions', name))
    for spc in EXAMPLE_SPECIES[name]:
        shutil.copytree(os.path.join(EXAMPLES, 'species', spc), os.path.join(root, 'species', spc))
    return os.path.join(root, 'reactions', name, 'input.py')


def _add_kinetics_kwargs(input_file, kwargs):
    """Insert keyword arguments into the last ``kinetics(`` call of an Arkane input file."""
    with open(input_file) as f:
        text = f.read()
    i = text.rindex('\nkinetics(') + len('\nkinetics(')
    with open(input_file, 'w') as f:
        f.write(text[:i] + kwargs + ', ' + text[i:])


def _run(input_file, output_directory):
    os.makedirs(output_directory)
    arkane = Arkane(input_file=input_file, output_directory=output_directory)
    arkane.plot = False
    arkane.execute()
    job = [j for j in arkane.job_list if isinstance(j, KineticsJob)][0]
    return job, output_directory


def _read(path):
    with open(path) as f:
        return f.read()


_NUMBER = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
_NASA_BLOCK = re.compile(r'^    thermo = NASA\(\n.*?^    \),\n', re.S | re.M)


def _text_and_numbers(output):
    """Split output.py text, without its fitted NASA blocks, into the text with numbers masked and the numbers."""
    output = _NASA_BLOCK.sub('', output)
    return _NUMBER.sub('#N', output), [float(x) for x in _NUMBER.findall(output)]


_INPUT_GLOBALS = ('species_dict', 'transition_state_dict', 'reaction_dict', 'network_dict', 'job_list')


@pytest.fixture(scope='module', autouse=True)
def _restore_arkane_input_globals():
    """Running Arkane in-process rebinds arkane.input's module-level dicts (load_input_file); restore them after
    this module, so later modules that call species()/transition_state() directly do not see its labels."""
    import arkane.input
    saved = {name: getattr(arkane.input, name) for name in _INPUT_GLOBALS}
    yield
    for name, value in saved.items():
        setattr(arkane.input, name, value)


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    """Run each example once with the default input and once with the block on (1 M)."""
    out = {}
    for name in (H_C2H4, NH2_N2H3):
        root = str(tmp_path_factory.mktemp(name.replace('+', '_').replace('=', '_')))
        input_file = _copy_example(name, os.path.join(root, 'default'))
        out[(name, 'default')] = _run(input_file, os.path.join(root, 'out_default'))
        input_file = _copy_example(name, os.path.join(root, 'block'))
        _add_kinetics_kwargs(input_file, "activation_block=True, standard_state='1 M'")
        out[(name, 'block')] = _run(input_file, os.path.join(root, 'out_block'))
    return out


def _block(runs, name, state, T=T298):
    job = runs[(name, 'default')][0]
    return compute_activation(job.reaction, T, state)


class TestDefaultOutputUnchanged:
    def test_default_output_matches_main(self, runs):
        """Without the option, output.py of the shipped example equals the one made by official/main
        (golden generated at 032e9d20e): all text exactly, every number to rel 1e-4. The fitted NASA polynomials
        are left out; their breakpoint and coefficients depend on the platform's optimizer, not on this PR."""
        out_dir = runs[(H_C2H4, 'default')][1]
        out_text, out_numbers = _text_and_numbers(_read(os.path.join(out_dir, 'output.py')))
        golden_text, golden_numbers = _text_and_numbers(_read(GOLDEN))
        assert out_text == golden_text
        assert len(out_numbers) == len(golden_numbers)
        for a, b in zip(out_numbers, golden_numbers):
            assert a == pytest.approx(b, rel=1e-4, abs=1e-12)

    def test_default_output_compare_is_not_vacuous(self):
        """The comparison above sees an extra line and a changed number."""
        golden_text, golden_numbers = _text_and_numbers(_read(GOLDEN))
        assert len(golden_numbers) > 400
        assert _text_and_numbers(_read(GOLDEN) + '# x\n')[0] != golden_text
        changed = _read(GOLDEN).replace('A = (4.00615e+08', 'A = (4.01e+08', 1)
        assert _text_and_numbers(changed)[1] != golden_numbers

    def test_block_only_appends(self, runs):
        """With the block on, output.py is the default output plus the block; chem.inp is unchanged."""
        default_dir = runs[(H_C2H4, 'default')][1]
        job, block_dir = runs[(H_C2H4, 'block')]
        default = _read(os.path.join(default_dir, 'output.py'))
        with_block = _read(os.path.join(block_dir, 'output.py'))
        assert job.activation_text
        assert job.activation_text in with_block
        assert with_block.replace(job.activation_text, '', 1) == default
        assert _read(os.path.join(block_dir, 'chem.inp')) == _read(os.path.join(default_dir, 'chem.inp'))
        assert repr(job.reaction.kinetics) == repr(runs[(H_C2H4, 'default')][0].reaction.kinetics)

    def test_block_text_states_state_and_units(self, runs):
        job = runs[(H_C2H4, 'block')][0]
        text = job.activation_text
        # default temperatures: 298.15 K and the ends of the fit range (400-1200 K)
        assert job.get_activation_temperatures() == [298.15, 400.0, 1200.0]
        assert text.count('# Standard state: 1 M (ideal gas, c0 = 1 mol/L)') == 3
        for T in ('298.15 K', '400 K', '1200 K'):
            assert 'H + C2H4 <=> C2H5 at {0}'.format(T) in text
        assert 'relative to the reactants of each direction, each at the standard state' in text
        assert 'M^-1 s^-1' in text
        assert 'kinetics(' not in text and 'Arrhenius(' not in text  # ARC parses those from output.py
        assert all(line.startswith('#') for line in text.splitlines() if line)


class TestHandValues:
    """H + C2H4 <=> C2H5 at 298.15 K, three states (reviewed prototype values)."""
    # state: (dG++ fwd kJ/mol, dS++ fwd J/(mol K), dG0 kJ/mol, dS0 J/(mol K), Keq, k_fwd in state units)
    STATE_VALUES = {
        '1 atm': (32.670, -84.778, -126.860, -76.539, 1.678084e22, 1.987189e7),
        '1 bar': (32.702, -84.887, -126.827, -76.649, 1.656141e22, 1.961203e7),
        '1 M': (24.744, -58.194, -134.785, -49.956, 4.105506e23, 4.861743e8),
    }

    @pytest.mark.parametrize('state', STATES)
    def test_hand_values(self, runs, state):
        r = _block(runs, H_C2H4, state)
        fwd, rev, rxn = r['directions']['forward'], r['directions']['reverse'], r['thermo']
        dG, dS, dG0, dS0, keq, k_state = self.STATE_VALUES[state]
        assert fwd['dG'] / 1000 == pytest.approx(dG, abs=1e-3)
        assert fwd['dS'] == pytest.approx(dS, abs=1e-3)
        assert rxn['dG0'] / 1000 == pytest.approx(dG0, abs=1e-3)
        assert rxn['dS0'] == pytest.approx(dS0, abs=1e-3)
        assert rxn['Keq'] == pytest.approx(keq, rel=1e-6)
        assert fwd['k_state'] == pytest.approx(k_state, rel=1e-6)
        # the same in every state
        assert r['dn'] == -1 and fwd['m'] == 2 and rev['m'] == 1
        assert fwd['dH'] / 1000 == pytest.approx(7.393, abs=1e-3)
        assert rev['dH'] / 1000 == pytest.approx(157.073, abs=1e-3)
        assert rev['dG'] / 1000 == pytest.approx(159.529, abs=1e-3)
        assert rev['dS'] == pytest.approx(-8.239, abs=1e-3)
        assert rxn['dH0'] / 1000 == pytest.approx(-149.680, abs=1e-3)
        assert r['kappa'] == pytest.approx(1.69221, rel=1e-5)
        assert r['tunneling'] == 'Eckart'
        assert fwd['k_cgs'] == pytest.approx(4.861743e11, rel=1e-6)
        assert fwd['k_unit'] == {'1 atm': 'atm^-1 s^-1', '1 bar': 'bar^-1 s^-1', '1 M': 'M^-1 s^-1'}[state]
        assert rev['k_state'] == pytest.approx(1.184201e-15, rel=1e-6)
        assert rev['k_unit'] == 's^-1'


class TestEyring:
    @pytest.mark.parametrize('name', (H_C2H4, NH2_N2H3))
    @pytest.mark.parametrize('state', STATES)
    @pytest.mark.parametrize('direction', ('forward', 'reverse'))
    def test_eyring_matches_tst(self, runs, name, state, direction):
        """kappa (kB T/h) C0^(1-m) exp(-dG++/RT) equals the TST k. On the reverse of the dn != 0 reaction,
        this catches a missing (P/RT)^dn factor in Ka -> Kc."""
        r = _block(runs, name, state)
        d = r['directions'][direction]
        assert abs(d['eyring_SI'] / d['k_SI'] - 1.0) < 1e-10
        assert abs(d['eyring_rel_diff']) < 1e-10

    def test_reverse_k_is_independent_reverse_tst(self, runs):
        """k_r = kf/Kc equals the reverse TST from the TS and the products directly (dn = -1)."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        T = T298
        r = compute_activation(rxn, T, '1 atm')
        v = constants.R * T / P_STATMECH
        q_prod = 1.0
        e0 = rxn.transition_state.conformer.E0.value_si
        for spc in rxn.products:
            q_prod *= spc.get_partition_function(T) / v
            e0 -= spc.conformer.E0.value_si
        q_ts = rxn.transition_state.get_partition_function(T) / v
        k_r = (constants.kB * T / constants.h * q_ts / q_prod * math.exp(-e0 / (constants.R * T))
               * rxn.transition_state.calculate_tunneling_factor(T))
        assert r['directions']['reverse']['k_SI'] == pytest.approx(k_r, rel=1e-10)


class TestDeltaNZero:
    def test_state_free_quantities(self, runs):
        """NH2 + N2H3 <=> NH + N2H4 (dn = 0): k_c, Keq and dG0 are the same in all three states."""
        blocks = [_block(runs, NH2_N2H3, state) for state in STATES]
        assert all(b['dn'] == 0 for b in blocks)
        ref = blocks[0]
        for b in blocks[1:]:
            for direction in ('forward', 'reverse'):
                assert b['directions'][direction]['k_SI'] == pytest.approx(ref['directions'][direction]['k_SI'],
                                                                           rel=1e-12)
            assert b['thermo']['Keq'] == pytest.approx(ref['thermo']['Keq'], rel=1e-12)
            assert b['thermo']['dG0'] == pytest.approx(ref['thermo']['dG0'], rel=1e-12)
        assert ref['thermo']['Keq'] == pytest.approx(7.186301e-09, rel=1e-6)
        assert ref['thermo']['dG0'] / 1000 == pytest.approx(46.483, abs=1e-3)


class TestEquilibriumConstantBasis:
    def test_block_kc_vs_get_equilibrium_constant(self, runs):
        """The block's Kc is from Q. get_equilibrium_constant uses P0 = 1e5 Pa on a 1 atm dG0 (a separate
        issue, not fixed here), so for dn = -1 the ratio is the named factor KC_BASIS_FACTOR. The 1e-6
        tolerance covers the NASA-fit residual of the example's thermo (about 3e-7), not the factor."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        r = compute_activation(rxn, T298, '1 atm')
        ratio = rxn.get_equilibrium_constant(T298, type='Kc') / r['Kc']
        assert r['Kc_ratio'] == pytest.approx(ratio, rel=1e-12)
        assert KC_BASIS_FACTOR == 101325. / 1e5
        assert ratio == pytest.approx(KC_BASIS_FACTOR ** (-r['dn']), rel=1e-6)
        assert abs(ratio - 1.0) > 1e-2  # the factor is really there; a tolerance does not hide it

    def test_statmech_basis_is_one_atm(self):
        """P_STATMECH must match the reference pressure of the translational partition function."""
        T = 500.
        tr = IdealGasTranslation(mass=(28.0313, 'amu'))
        mass = tr.mass.value_si
        expected = (2 * math.pi * mass * constants.kB * T / constants.h ** 2) ** 1.5 * constants.kB * T / P_STATMECH
        assert tr.get_partition_function(T) == pytest.approx(expected, rel=1e-12)


class TestInputValidation:
    @pytest.mark.parametrize('bad', ('1 Bar', '1M', '1 mol/L', '', 1, 1.0, None))
    def test_bad_state_raises(self, bad):
        with pytest.raises(InputError):
            validate_standard_state(bad)

    @pytest.mark.parametrize('state', STATES)
    def test_good_state(self, state):
        assert validate_standard_state(state) == state

    def test_bad_state_raises_from_input_file(self, tmp_path):
        input_file = _copy_example(NH2_N2H3, str(tmp_path / 'ex'))
        _add_kinetics_kwargs(input_file, "activation_block=True, standard_state='1 mbar'")
        os.makedirs(str(tmp_path / 'out'))
        arkane = Arkane(input_file=input_file, output_directory=str(tmp_path / 'out'))
        arkane.plot = False
        with pytest.raises(InputError, match='Unknown standard state'):
            arkane.execute()

    def test_kineticsjob_validation(self, runs):
        rxn = runs[(NH2_N2H3, 'default')][0].reaction
        with pytest.raises(InputError):
            KineticsJob(rxn, activation_block='yes')
        with pytest.raises(InputError):  # state given while the block is off would print nothing
            KineticsJob(rxn, standard_state='1 M')
        with pytest.raises(InputError):
            KineticsJob(rxn, activation_Tlist=([300], 'K'))
        with pytest.raises(InputError):
            KineticsJob(rxn, activation_block=True, standard_state='STP')
        with pytest.raises(InputError):
            KineticsJob(rxn, activation_block=True, activation_Tlist=([0, 300], 'K'))
        for bad_t in (float('nan'), float('inf')):
            with pytest.raises(InputError):
                KineticsJob(rxn, activation_block=True, activation_Tlist=([300, bad_t], 'K'))
            with pytest.raises(InputError):
                compute_activation(rxn, bad_t, '1 atm')
        job = KineticsJob(rxn, activation_block=True, activation_Tlist=(500, 'K'))  # a scalar is one temperature
        assert job.get_activation_temperatures() == [500.0]
        job = KineticsJob(rxn, activation_block=True)
        assert job.standard_state == '1 atm'
        job = KineticsJob(rxn)
        assert job.activation_block is False and job.activation_text is None

    def test_default_temperatures_drop_ends_near_298(self, runs):
        """The default fit range (298-2500 K) gives 298.15 and 2500 K, not a 298 K section 0.15 K away."""
        rxn = runs[(NH2_N2H3, 'default')][0].reaction
        job = KineticsJob(rxn, activation_block=True)
        assert job.get_activation_temperatures() == [298.15, 2500.0]
        job = KineticsJob(rxn, activation_block=True, Tmin=(299.2, 'K'), Tmax=(1000, 'K'))
        assert job.get_activation_temperatures() == [298.15, 299.2, 1000.0]  # 1.05 K away: kept

    def test_write_output_failure_raises_with_block(self, runs, tmp_path, monkeypatch):
        """With the block on, a write_output failure raises from execute (the block would be lost); with it off,
        execute keeps the existing warn-only behavior."""
        import copy
        rxn = copy.deepcopy(runs[(NH2_N2H3, 'default')][0].reaction)
        rxn.kinetics = None  # so generate_kinetics runs TST

        def fail(self, output_directory):
            raise ValueError('No valid NASA polynomial at temperature 3500 K.')

        monkeypatch.setattr(KineticsJob, 'write_output', fail)
        monkeypatch.setattr(KineticsJob, 'write_chemkin', lambda self, output_directory: None)
        job = KineticsJob(rxn, Tlist=([400, 600, 800, 1000], 'K'), activation_block=True, activation_Tlist=([500], 'K'))
        with pytest.raises(ValueError, match='NASA'):
            job.execute(output_directory=str(tmp_path), plot=False)
        rxn.kinetics = None
        job = KineticsJob(rxn, Tlist=([400, 600, 800, 1000], 'K'))
        job.execute(output_directory=str(tmp_path), plot=False)  # warns only
        assert job.usedTST

    def test_state_revalidated_at_use(self, runs):
        """A state changed on the job after construction is checked again when the block is built."""
        job = runs[(NH2_N2H3, 'default')][0]
        rxn = job.reaction
        new_job = KineticsJob(rxn, activation_block=True)
        new_job.usedTST = True
        new_job.standard_state = '2 atm'
        with pytest.raises(InputError):
            new_job.generate_activation_block()

    def test_block_without_tst_raises_from_execute(self, runs):
        """A block on a reaction with given kinetics raises from execute; it is not logged and dropped."""
        from rmgpy.kinetics import Arrhenius
        from rmgpy.reaction import Reaction
        rxn = Reaction(reactants=[], products=[],
                       kinetics=Arrhenius(A=(1e10, 's^-1'), n=0, Ea=(10, 'kJ/mol'), T0=(1, 'K')))
        rxn.transition_state = runs[(NH2_N2H3, 'default')][0].reaction.transition_state
        job = KineticsJob(rxn, activation_block=True)
        with pytest.raises(InputError, match='needs TST kinetics'):
            job.execute(output_directory=None, plot=False)


class TestExtremeTemperature:
    def test_overflowing_kc_does_not_abort(self, runs):
        """At 10 K, Kc of H + C2H4 <=> C2H5 overflows a float and k_r underflows. The block is still built
        from logarithms: dG0 stays finite, and the Eyring check holds for both directions."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        r = compute_activation(rxn, 10., '1 M')
        assert math.isinf(r['Kc'])
        assert math.isfinite(r['thermo']['dG0'])
        assert r['thermo']['dG0'] == pytest.approx(r['directions']['forward']['dG'] - r['directions']['reverse']['dG'],
                                                   rel=1e-12)
        for direction in ('forward', 'reverse'):
            assert abs(r['directions'][direction]['eyring_rel_diff']) < 1e-10
        assert 'Standard state: 1 M' in format_activation_block([r])


class _ReactionProxy:
    """Delegate to a real reaction, but override some methods (Reaction is a compiled class)."""

    def __init__(self, reaction, **overrides):
        self._reaction = reaction
        self._overrides = overrides

    def __getattr__(self, name):
        if name in self._overrides:
            return self._overrides[name]
        return getattr(self._reaction, name)

    def __str__(self):
        return str(self._reaction)


class TestUnderflow:
    def test_forward_tst_underflow_keeps_reverse_k(self, runs):
        """If Arkane's TST kf underflows to 0, kf and kr come from the same TST expression in log space."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        ref = compute_activation(rxn, T298, '1 atm')
        proxy = _ReactionProxy(rxn, calculate_tst_rate_coefficient=lambda T: 0.0)
        r = compute_activation(proxy, T298, '1 atm')
        for direction in ('forward', 'reverse'):
            assert r['directions'][direction]['k_SI'] == pytest.approx(ref['directions'][direction]['k_SI'],
                                                                       rel=1e-10)
            assert abs(r['directions'][direction]['eyring_rel_diff']) < 1e-10

    def test_k_underflowing_only_in_si_is_kept(self, runs, monkeypatch):
        """A forward k below the smallest float in m^3/(mol s) but not in L/(mol s) or cm^3/(mol s) is printed
        nonzero in those units: the conversions are done on ln k, not on the underflowed SI value."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        ts = rxn.transition_state
        ln_kf_true = math.log(compute_activation(rxn, T298, '1 M')['directions']['forward']['k_SI'])
        delta = (ln_kf_true + 748.) * constants.R * T298  # raise the TS energy so ln k_SI = -748
        real_terms = arkane.activation._species_terms

        def shifted_terms(obj, T):
            d = real_terms(obj, T)
            if obj is ts:
                d['E0'] += delta
                d['H'] += delta
            return d

        monkeypatch.setattr(arkane.activation, '_species_terms', shifted_terms)
        proxy = _ReactionProxy(rxn, calculate_tst_rate_coefficient=lambda T: 0.0)  # TST kf from the shifted Q
        fwd = compute_activation(proxy, T298, '1 M')['directions']['forward']
        assert fwd['m'] == 2
        assert fwd['k_SI'] == 0.0
        assert fwd['k_state'] > 0.0 and fwd['k_cgs'] > 0.0
        assert fwd['k_state'] == pytest.approx(fwd['eyring_state'], rel=0.05)  # subnormal floats: few digits
        # Arkane's TST kf a subnormal float (few digits): the log-space TST value is used instead
        proxy = _ReactionProxy(rxn, calculate_tst_rate_coefficient=lambda T: 5e-324)
        fwd = compute_activation(proxy, T298, '1 M')['directions']['forward']
        assert abs(fwd['eyring_rel_diff']) < 1e-10

    def test_forward_tst_overflow_keeps_reverse_k(self, runs):
        """A submerged barrier at 10 K makes Reaction.calculate_tst_rate_coefficient raise OverflowError (math.exp).
        The block is still built: kf comes from the log-space TST expression and prints as inf, and both directions
        pass the Eyring check, which compares logarithms."""
        import copy
        rxn = copy.deepcopy(runs[(H_C2H4, 'default')][0].reaction)
        ts = rxn.transition_state
        ts.tunneling = None
        e0_reac = sum(spc.conformer.E0.value_si for spc in rxn.reactants)
        ts.conformer.E0 = (e0_reac - 60000., 'J/mol')  # 60 kJ/mol below the reactants
        with pytest.raises(OverflowError):
            rxn.calculate_tst_rate_coefficient(10.)
        r = compute_activation(rxn, 10., '1 atm')
        assert math.isinf(r['directions']['forward']['k_SI'])
        for direction in ('forward', 'reverse'):
            assert abs(r['directions'][direction]['eyring_rel_diff']) < 1e-10
        assert math.isfinite(r['thermo']['dG0'])
        assert 'Standard state: 1 atm' in format_activation_block([r])

    def test_get_equilibrium_constant_outside_nasa_range_is_not_fatal(self, runs):
        """At 3500 K, above the 3000 K end of Arkane's NASA fits, get_equilibrium_constant raises ValueError.
        The comparison is then nan and the block, computed from Q, is still built."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        with pytest.raises(ValueError):
            rxn.get_equilibrium_constant(3500., type='Kc')
        r = compute_activation(rxn, 3500., '1 bar')
        assert math.isnan(r['Kc_ratio'])
        for direction in ('forward', 'reverse'):
            assert abs(r['directions'][direction]['eyring_rel_diff']) < 1e-10
        assert 'ratio = nan' in format_activation_block([r])

    def test_subnormal_boltzmann_factor_in_normal_kf(self, runs):
        """C2H5 => H + C2H4 with dE0 = 735 RT at 50 K: Arkane's kf is a normal float, but its exp(-dE0/RT) was
        subnormal and lost digits (3e-5 relative). The log-space TST value is used, so the Eyring check holds."""
        import copy
        from rmgpy.reaction import Reaction
        fwd = copy.deepcopy(runs[(H_C2H4, 'default')][0].reaction)
        ts = fwd.transition_state
        ts.tunneling = None
        rxn = Reaction(reactants=fwd.products, products=fwd.reactants, transition_state=ts)
        T = 50.
        e0_reac = sum(spc.conformer.E0.value_si for spc in rxn.reactants)
        ts.conformer.E0 = (e0_reac + 735. * constants.R * T, 'J/mol')
        kf = rxn.calculate_tst_rate_coefficient(T)
        assert kf >= sys.float_info.min  # the case: a normal kf
        r = compute_activation(rxn, T, '1 atm')
        for direction in ('forward', 'reverse'):
            assert abs(r['directions'][direction]['eyring_rel_diff']) < 1e-10

    @pytest.mark.parametrize('kappa', (0.0, -1.0, float('inf'), float('nan')))
    def test_bad_tunneling_factor_raises(self, runs, kappa):
        """An invalid tunneling factor raises an InputError instead of being printed as an underflowed k = 0."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        ts = _ReactionProxy(rxn.transition_state, calculate_tunneling_factor=lambda T: kappa)
        proxy = _ReactionProxy(rxn, transition_state=ts)
        with pytest.raises(InputError, match='tunneling factor'):
            compute_activation(proxy, T298, '1 atm')

    def test_rel_diff_does_not_overflow(self):
        """A discrepancy too large for a float is reported as inf, not raised."""
        assert arkane.activation._rel_diff(800., 0.) == float('inf')
        assert arkane.activation._rel_diff(0., 800.) == pytest.approx(-1.0)

    @pytest.mark.parametrize('q', (0.0, float('inf'), float('nan')))
    def test_bad_partition_function_raises(self, runs, q):
        """A partition function that is 0, inf or nan as a float raises an InputError naming the species,
        instead of printing nan."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        ts = _ReactionProxy(rxn.transition_state, get_partition_function=lambda T: q)
        proxy = _ReactionProxy(rxn, transition_state=ts)
        with pytest.raises(InputError, match='partition function of'):
            compute_activation(proxy, T298, '1 atm')

    def test_get_equilibrium_constant_error_is_not_fatal(self, runs):
        """get_equilibrium_constant raises ReactionError when its K underflows to 0 (C2H5 => H + C2H4 at
        10 K). The comparison is then nan and the block is still built."""
        from rmgpy.reaction import Reaction
        fwd = runs[(H_C2H4, 'default')][0].reaction
        rxn = Reaction(reactants=fwd.products, products=fwd.reactants, transition_state=fwd.transition_state)
        with pytest.raises(ReactionError):
            rxn.get_equilibrium_constant(10., type='Kc')
        r = compute_activation(rxn, 10., '1 M')
        assert math.isnan(r['Kc_get_equilibrium_constant']) and math.isnan(r['Kc_ratio'])
        assert math.isfinite(r['thermo']['dG0'])
        assert 'ratio = nan' in format_activation_block([r])

    @pytest.mark.parametrize('kc_gec', (0.0, float('inf')))
    def test_degenerate_get_equilibrium_constant_formats(self, runs, kc_gec):
        """A get_equilibrium_constant Kc of 0 or inf makes the comparison ratio nan; formatting does not raise."""
        rxn = runs[(H_C2H4, 'default')][0].reaction
        proxy = _ReactionProxy(rxn, get_equilibrium_constant=lambda T, type='Kc': kc_gec)
        r = compute_activation(proxy, T298, '1 atm')
        assert math.isnan(r['Kc_ratio'])
        assert 'ratio = nan' in format_activation_block([r])


class TestUnits:
    def test_k_in_state_units(self):
        T = 1000.
        k_si = 1.0e-6  # m^3/(mol s) = 1 cm^3/(mol s)
        assert k_in_state_units(k_si, 2, '1 M', T) == pytest.approx(1.0e-3)
        assert k_in_state_units(k_si, 2, '1 bar', T) == pytest.approx(k_si / (constants.R * T) * 1e5)
        assert k_in_state_units(k_si, 2, '1 atm', T) == pytest.approx(k_si / (constants.R * T) * 101325.)
        assert k_in_state_units(3.0, 1, '1 atm', T) == 3.0
        assert standard_concentration('1 M', T) == 1000.
        assert standard_concentration('1 bar', T) == pytest.approx(1e5 / (constants.R * T))
        assert set(STANDARD_STATES) == set(STATES)

    def test_format_has_ratio_line_only_for_nonzero_dn(self, runs):
        text_h = format_activation_block([_block(runs, H_C2H4, '1 bar')])
        text_n = format_activation_block([_block(runs, NH2_N2H3, '1 bar')])
        assert 'ratio = 1.013250' in text_h
        assert 'expected from the reference-pressure basis: 1.013250' in text_h
        assert 'remainder ratio/expected - 1' in text_h
        assert 'get_equilibrium_constant' not in text_n
