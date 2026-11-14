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
Activation and reaction thermochemistry of an Arkane kinetics job in a
user-chosen standard state.

This is an optional report. It changes nothing that Arkane computes or writes
as data (fitted Arrhenius, libraries, Chemkin, YAML).

Convention
----------
All Arkane partition functions Q are on a 1 atm basis: the translational Q
uses the molecular volume kB*T/101325 Pa (``rmgpy/statmech/translation.pyx``).
For a standard state with standard concentration C0 (C0 = P0/RT for a
pressure state, C0 = c0 for 1 M) the standard-state partition function is
``Q0 = Q * P_STATMECH / (R T C0)``. Activation quantities are relative to the
reactants of each direction, each at the standard state::

    K++   = (Q0_TS / prod Q0_r) exp(-(E0_TS - sum E0_r)/RT),   dG++ = -RT ln K++
    dH++  = H_TS - sum H_r,   H_i = E0_i + H_thermal,i(T)      (no state)
    dS++  = (dH++ - dG++) / T

so that ``k_c = kappa (kB T/h) C0^(1-m) K++`` in every state. The tunneling
factor kappa is reported apart and is not folded into dG++.

Reaction quantities use the same Q::

    Kc  = prod (Q_p/V) e^(-E0_p/RT) / prod (Q_r/V) e^(-E0_r/RT),   V = RT/P_STATMECH
    Keq = Kc / C0^dn   (= Kp / P0^dn for a pressure state),        dG0 = -RT ln Keq
    dH0 = sum H_p - sum H_r,                                       dS0 = (dH0 - dG0)/T

The reverse rate coefficient is ``k_r = k_f / Kc`` with the same kappa.
Kc is deliberately *not* taken from ``Reaction.get_equilibrium_constant``,
which combines a 1 atm dG0 with P0 = 1e5 Pa (see ``KC_BASIS_FACTOR``).
"""

import math
import sys

import rmgpy.constants as constants
from rmgpy.exceptions import InputError, ReactionError

# The reference pressure of Arkane's translational partition function, in Pa.
# Must match ``rmgpy/statmech/translation.pyx`` (IdealGasTranslation) and the
# ``R*T/101325.`` divisor in ``Reaction.calculate_tst_rate_coefficient``.
P_STATMECH = 101325.

# Known factor, not fixed here: ``Reaction.get_equilibrium_constant`` converts a
# dG0 on the 1 atm basis of Arkane's partition functions to Kc with P0 = 1e5 Pa,
# so its Kc differs from the Kc of this module by (101325/1e5)^(-dn).
# This is a separate issue in RMG-Py; when it is fixed, this factor becomes 1.
KC_BASIS_FACTOR = 101325. / 1.0e5

# label -> (kind, value in SI): a pressure in Pa or a concentration in mol/m^3
STANDARD_STATES = {
    '1 atm': ('pressure', 101325.),
    '1 bar': ('pressure', 1.0e5),
    '1 M': ('concentration', 1000.),
}

DEFAULT_STANDARD_STATE = '1 atm'

_STATE_DESCRIPTIONS = {
    '1 atm': '1 atm (ideal gas, P0 = 101325 Pa)',
    '1 bar': '1 bar (ideal gas, P0 = 100000 Pa)',
    '1 M': '1 M (ideal gas, c0 = 1 mol/L)',
}

_STATE_UNITS = {'1 atm': 'atm', '1 bar': 'bar', '1 M': 'M'}


def validate_standard_state(state):
    """
    Return `state` if it is an accepted standard state label, else raise an InputError.
    """
    if not isinstance(state, str) or state not in STANDARD_STATES:
        raise InputError('Unknown standard state {0!r}. Accepted values are: {1}.'.format(
            state, ', '.join(repr(s) for s in STANDARD_STATES)))
    return state


def standard_concentration(state, T):
    """
    Return the standard concentration C0 in mol/m^3 of `state` at temperature `T` in K.
    """
    kind, value = STANDARD_STATES[validate_standard_state(state)]
    return value / (constants.R * T) if kind == 'pressure' else value


def k_in_state_units(k_si, m, state, T):
    """
    Convert a rate coefficient `k_si` of molecularity `m` from SI concentration units
    ((m^3/mol)^(m-1) s^-1) to the units of `state`: M^(1-m) s^-1 for a concentration state,
    or bar^(1-m) s^-1 / atm^(1-m) s^-1 for a pressure state (k_p = k_c (RT)^(1-m)).
    """
    kind, value = STANDARD_STATES[validate_standard_state(state)]
    if kind == 'concentration':
        return k_si * value ** (m - 1)
    k_pa = k_si * (constants.R * T) ** (1 - m)
    return k_pa * value ** (m - 1)


def k_state_unit_label(m, state):
    """
    Return the unit label of a rate coefficient of molecularity `m` in the units of `state`.
    """
    if m == 1:
        return 's^-1'
    return '{0}^{1:d} s^-1'.format(_STATE_UNITS[validate_standard_state(state)], 1 - m)


def _exp(x):
    """exp(x), returning inf instead of raising OverflowError (a report must not abort on a huge K)."""
    try:
        return math.exp(x)
    except OverflowError:
        return float('inf')


def _rel_diff(ln_a, ln_b):
    """a/b - 1 from logarithms; nan when either value is zero (log of 0 is -inf) or undefined."""
    if not (math.isfinite(ln_a) and math.isfinite(ln_b)):
        return float('nan')
    try:
        return math.expm1(ln_a - ln_b)
    except OverflowError:
        return float('inf')


def _ln(x):
    return math.log(x) if x > 0 else float('-inf')


def check_temperature(T):
    """Return `T` as a float if it is a finite positive temperature in K, else raise an InputError."""
    try:
        T = float(T)
    except (TypeError, ValueError):
        raise InputError('The activation block temperature must be a number, got {0!r}.'.format(T))
    if not math.isfinite(T) or T <= 0:
        raise InputError('The activation block temperature must be finite and positive, got {0} K.'.format(T))
    return T


def _require_conformer(obj, role, reaction):
    if getattr(obj, 'conformer', None) is None:
        raise InputError('Cannot compute the activation block of reaction {0}: the {1} {2!r} has no '
                         'conformer (statmech) data.'.format(reaction, role, obj.label))


def _species_terms(obj, T):
    """E0, H = E0 + thermal H, and ln Q (1 atm basis) of a species or transition state, SI."""
    conformer = obj.conformer
    q = float(obj.get_partition_function(T))
    if not (0 < q < float('inf')):
        raise InputError('Cannot compute the activation block: the partition function of {0!r} at {1} K is {2}, '
                         'not a finite positive number.'.format(obj.label, T, q))
    return {'E0': conformer.E0.value_si,
            'H': conformer.E0.value_si + conformer.get_enthalpy(T),
            'lnQ': math.log(q)}


def _ln_k_dimensionless(final, initial, T, c0):
    """
    ln of the dimensionless equilibrium constant between two sets of species,
    each at standard concentration `c0` (mol/m^3): ln Q0_i = ln Q_i + ln(P_STATMECH/(R T c0)).
    """
    shift = math.log(P_STATMECH / (constants.R * T * c0))
    total = 0.0
    for d in final:
        total += d['lnQ'] + shift - d['E0'] / (constants.R * T)
    for d in initial:
        total -= d['lnQ'] + shift - d['E0'] / (constants.R * T)
    return total


def compute_activation(reaction, T, state):
    """
    Compute the activation block of `reaction` at temperature `T` (K) in standard state `state`.

    The reaction needs conformers for all reactants, products and the transition state, and the
    tunneling model of the transition state must already be set up (as done in
    ``KineticsJob.generate_kinetics``). All returned values are in SI units.
    """
    validate_standard_state(state)
    T = check_temperature(T)
    ts = reaction.transition_state
    if ts is None:
        raise InputError('Cannot compute the activation block of reaction {0}: no transition state.'.format(
            reaction))
    for spc in reaction.reactants:
        _require_conformer(spc, 'reactant', reaction)
    for spc in reaction.products:
        _require_conformer(spc, 'product', reaction)
    _require_conformer(ts, 'transition state', reaction)

    R, kB, h = constants.R, constants.kB, constants.h
    reac = [_species_terms(s, T) for s in reaction.reactants]
    prod = [_species_terms(s, T) for s in reaction.products]
    tsd = _species_terms(ts, T)
    m_f, m_r = len(reac), len(prod)
    dn = m_r - m_f
    kappa = float(ts.calculate_tunneling_factor(T))
    if not (0 < kappa < float('inf')):
        raise InputError('Cannot compute the activation block of reaction {0}: the tunneling factor at {1} K is {2}, '
                         'not a finite positive number.'.format(reaction, T, kappa))
    tunneling = type(ts.tunneling).__name__ if ts.tunneling is not None else 'None'

    # Arkane's TST k, the concentration-basis k_c; it includes kappa and does not depend on the state
    kf = reaction.calculate_tst_rate_coefficient(T)
    # State-free Kc (mol/m^3)^dn from Q: Ka on the P_STATMECH basis times (P_STATMECH/RT)^dn
    c_statmech = P_STATMECH / (R * T)
    # Kept in log space, so an extreme Kc at low T cannot overflow before k_r is formed
    ln_kc = _ln_k_dimensionless(prod, reac, T, c_statmech) + dn * math.log(c_statmech)
    # ln kf from Arkane's TST k; if that underflowed at low T (to 0, or to a subnormal float with few digits left),
    # the same TST expression from Q in log space, so a representable reverse k is not lost
    ln_kf = _ln(kf)
    if not sys.float_info.min <= kf < float('inf'):
        ln_kf = (_ln(kappa) + math.log(kB * T / h) + (1 - m_f) * math.log(c_statmech)
                 + _ln_k_dimensionless([tsd], reac, T, c_statmech))
    kc = _exp(ln_kc)
    kr = _exp(ln_kf - ln_kc)

    c0 = standard_concentration(state, T)
    sum_h_reac = sum(d['H'] for d in reac)
    sum_h_prod = sum(d['H'] for d in prod)
    directions = {}
    for name, m, initial, ln_k_si, dH in (('forward', m_f, reac, ln_kf, tsd['H'] - sum_h_reac),
                                          ('reverse', m_r, prod, ln_kf - ln_kc, tsd['H'] - sum_h_prod)):
        ln_k_dagger = _ln_k_dimensionless([tsd], initial, T, c0)
        dG = -R * T * ln_k_dagger
        # In the state's own unit the prefactor C0^(1-m) is 1^(1-m)
        ln_eyring_state = _ln(kappa) + math.log(kB * T / h) + ln_k_dagger
        ln_eyring_si = ln_eyring_state + (1 - m) * math.log(c0)
        # Unit conversions in log space, so a k that underflows in SI but not in the state's or cgs units is kept
        ln_k_state = ln_k_si + math.log(k_in_state_units(1.0, m, state, T))
        ln_k_cgs = ln_k_si + (m - 1) * math.log(1.0e6)
        directions[name] = {
            'm': m, 'dH': dH, 'dG': dG, 'dS': (dH - dG) / T,
            'k_SI': _exp(ln_k_si), 'k_state': _exp(ln_k_state), 'k_unit': k_state_unit_label(m, state),
            'k_cgs': _exp(ln_k_cgs),
            'eyring_state': _exp(ln_eyring_state),
            'eyring_SI': _exp(ln_eyring_si),
            'eyring_rel_diff': _rel_diff(ln_eyring_state, ln_k_state),
            'eyring_rel_diff_SI': _rel_diff(ln_eyring_si, ln_k_si),
        }
    ln_keq = ln_kc - dn * math.log(c0)
    keq = _exp(ln_keq)
    dH0 = sum_h_prod - sum_h_reac
    dG0 = -R * T * ln_keq
    try:
        kc_gec = float(reaction.get_equilibrium_constant(T, type='Kc'))
    except ReactionError:
        # get_equilibrium_constant raises when its K underflows to 0; the comparison is then unavailable
        kc_gec = float('nan')
    ln_kc_gec = _ln(kc_gec)
    return {
        'reaction': str(reaction), 'T': T, 'state': state, 'C0': c0,
        'kappa': kappa, 'tunneling': tunneling, 'dn': dn,
        'Kc': kc, 'Kc_get_equilibrium_constant': kc_gec,
        # ratio get_equilibrium_constant / this block, from logs; nan if either Kc is 0, inf or unavailable
        'Kc_ratio': (_exp(ln_kc_gec - ln_kc) if math.isfinite(ln_kc_gec) and math.isfinite(ln_kc)
                     else float('nan')),
        'directions': directions,
        'thermo': {'dH0': dH0, 'dG0': dG0, 'dS0': (dH0 - dG0) / T, 'Keq': keq},
    }


def format_activation_block(results):
    """
    Return the activation block text (``#``-commented lines) for a list of `results`
    from :func:`compute_activation`, one section per temperature.
    """
    sections = [_format_one(r) for r in results]
    return '\n\n'.join(sections) + '\n\n'


def _format_one(r):
    kJ, kcal, cal = 1000.0, 4184.0, 4.184
    state = r['state']
    rows = r['directions']
    rr = r['thermo']
    lines = [
        '# Standard state: {0}'.format(_STATE_DESCRIPTIONS[state]),
        '# Activation and reaction thermochemistry for {0} at {1:g} K'.format(r['reaction'], r['T']),
        '# Activation quantities are relative to the reactants of each direction, each at the standard state;',
        '# dS++ and dG++ depend on the standard state unless m = 1. kappa is not folded into dG++.',
        '#   C0 = {0:.6g} mol/m^3   kappa ({1}) = {2:.6g}'.format(r['C0'], r['tunneling'], r['kappa']),
        '#',
        '#   Direction  m   dH++          dS++               dG++                 '
        'k ({0} units)              k'.format(state),
        '#                  kJ/mol        J/(mol*K)          kJ/mol    kcal/mol   '
        '                              cm3, mol, s',
    ]
    for name in ('forward', 'reverse'):
        d = rows[name]
        lines.append('#   {0:9s} {1:d}  {2:9.3f}  {3:9.3f}  {4:9.3f}  {5:9.3f}  {6:12.6e} {7:12s}  '
                     '{8:12.6e}'.format(name, d['m'], d['dH'] / kJ, d['dS'], d['dG'] / kJ,
                                        d['dG'] / kcal, d['k_state'], d['k_unit'], d['k_cgs']))
    lines.append('#   (dH++ in kcal/mol: forward {0:.3f}, reverse {1:.3f}; dS++ in cal/(mol*K): '
                 'forward {2:.3f}, reverse {3:.3f})'.format(
                     rows['forward']['dH'] / kcal, rows['reverse']['dH'] / kcal,
                     rows['forward']['dS'] / cal, rows['reverse']['dS'] / cal))
    keq_def = 'Kc/c0^dn' if STANDARD_STATES[state][0] == 'concentration' else 'Kp/P0^dn'
    lines.extend([
        '#',
        '#   Reaction   dn = {0:+d}'.format(r['dn']),
        '#     dH0 = {0:10.3f} kJ/mol  ({1:9.3f} kcal/mol)'.format(rr['dH0'] / kJ, rr['dH0'] / kcal),
        '#     dS0 = {0:10.3f} J/(mol*K) ({1:9.3f} cal/(mol*K))'.format(rr['dS0'], rr['dS0'] / cal),
        '#     dG0 = {0:10.3f} kJ/mol  ({1:9.3f} kcal/mol)'.format(rr['dG0'] / kJ, rr['dG0'] / kcal),
        '#     Keq = {0} = {1:.6e} (dimensionless)'.format(keq_def, rr['Keq']),
        '#',
        '#   Eyring check: kappa*(kB*T/h)*C0^(1-m)*exp(-dG++/RT) vs TST k, relative difference',
    ])
    for name in ('forward', 'reverse'):
        d = rows[name]
        lines.append('#     {0:8s} {1:.12e} {2:12s} vs {3:.12e}   {4:+.2e}'.format(
            name, d['eyring_state'], d['k_unit'], d['k_state'], d['eyring_rel_diff']))
    if r['dn'] != 0:
        lines.append('#   Kc in SI units (this block, from Q) = {0:.6e}; '
                     'Reaction.get_equilibrium_constant Kc = {1:.6e}; ratio = {2:.6f}'.format(
                         r['Kc'], r['Kc_get_equilibrium_constant'], r['Kc_ratio']))
        expected = KC_BASIS_FACTOR ** (-r['dn'])
        lines.append('#   (expected from the reference-pressure basis: {0:.6f} = {1:.5f}^(-dn), because '
                     'get_equilibrium_constant uses P0 = 1e5 Pa on a 1 atm dG0;'.format(expected, KC_BASIS_FACTOR))
        lines.append('#    remainder ratio/expected - 1 = {0:+.2e}, from its thermo (e.g. NASA fit) vs the partition '
                     'functions. The reverse k above is kf/Kc from Q.)'.format(
                         r['Kc_ratio'] / expected - 1.0))
    return '\n'.join(lines)
