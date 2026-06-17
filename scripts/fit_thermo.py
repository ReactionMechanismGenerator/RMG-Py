# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: rmg_env18
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Fitting a 14-parameter NASA polynomial one coefficient at a time
#
# We import a Chemkin-format thermo entry for `C4H6-2` (1,2-butadiene) as an RMG
# [`NASA`](../rmgpy/thermo/nasa.pyx) object, then ask a *sensitivity* question:
#
# > Starting from the published polynomial, if we are allowed to adjust **only one**
# > of the 14 coefficients, how much closer can we get to a set of target
# > thermochemistry data?
#
# So we perform **14 independent single-parameter optimizations**. Each one starts
# from the original polynomial, optimizes exactly one coefficient, records the
# result, and then *restores* that coefficient before moving to the next — every
# fitted polynomial therefore differs from the original in exactly one number.
#
# The data we fit to:
#
# | quantity   | target            |
# |------------|-------------------|
# | ΔHf(298 K) | 34.75 kcal/mol    |
# | S(298 K)   | 67.66 cal/mol/K   |
# | Cp(298 K)  | 18.63 cal/mol/K   |
# | Cp(500 K)  | 26.36 cal/mol/K   |
# | Cp(800 K)  | 35.13 cal/mol/K   |
# | Cp(1000 K) | 39.20 cal/mol/K   |
# | Cp(1500 K) | 45.77 cal/mol/K   |
#
# plus the *additional goal* that Cp, H and S are continuous where the two
# polynomials meet (the common temperature `Tint = 1000 K`), included as a soft
# penalty in the objective.
#
# We reuse RMG's own machinery throughout: `rmgpy.chemkin.read_thermo_entry` to
# parse the Chemkin block, and `NASAPolynomial.get_heat_capacity / get_enthalpy /
# get_entropy` (which return SI units, J/mol·K and J/mol) to evaluate the curves.

# %%
import copy

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from rmgpy.chemkin import read_thermo_entry

# Importing rmgpy forces the non-interactive "Agg" matplotlib backend, so request
# the inline backend *after* that import to get figures embedded in the notebook.
# %matplotlib inline

# Unit conversions: RMG thermo getters return SI (J/mol, J/mol/K).
CAL = 4.184           # J per cal      -> divide J/mol/K  by this for cal/mol/K
KCAL = 4184.0         # J per kcal     -> divide J/mol     by this for kcal/mol

# %% [markdown]
# ## 1. Import the Chemkin thermo entry
#
# `read_thermo_entry` returns the species label, a `NASA` object (two
# `NASAPolynomial` pieces — note the *high*-T polynomial is written first in the
# Chemkin format), and the elemental formula.

# %%
entry = """C4H6-2            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
  9.0338133E+00  8.2124510E-03  7.1753952E-06 -5.8834334E-09  1.0343915E-12    2
  1.4335068E+04 -2.0985762E+01  2.1373338E+00  2.6486229E-02 -9.0568711E-06    3
 -5.5386397E-19  2.1281884E-22  1.5710902E+04  1.3529426E+01  1.7488676E+04    4"""

species, thermo0, formula = read_thermo_entry(entry)

# The data we are fitting to is reported at 298 K, just below the published lower
# limit of 300 K. Per the exercise, we lower the low-T validity limit to 298 K so
# the polynomial can be evaluated there.
thermo0.Tmin = (298.0, "K")
thermo0.polynomials[0].Tmin = (298.0, "K")

Tint = thermo0.polynomials[0].Tmax.value_si  # common temperature where the two polys meet (1000 K)

print(f"species : {species}")
print(f"formula : {formula}")
print(f"valid   : {thermo0.Tmin.value_si:.0f}-{thermo0.Tmax.value_si:.0f} K, joined at Tint = {Tint:.0f} K")
print(f"low-T  coeffs (298-{Tint:.0f} K): {thermo0.polynomials[0].coeffs}")
print(f"high-T coeffs ({Tint:.0f}-{thermo0.Tmax.value_si:.0f} K): {thermo0.polynomials[1].coeffs}")

# %% [markdown]
# ## 2. Target data
#
# Heat capacities at five temperatures, plus the 298 K enthalpy of formation and
# entropy. Each Cp point is associated with the polynomial that owns that
# temperature (≤ `Tint` → low-T piece, otherwise high-T piece); `Cp(1000 K)` sits
# exactly on the boundary and is assigned to the low-T piece.

# %%
T_REF = 298.15  # standard reference temperature for the "298 K" quantities

Cp_T = np.array([298.15, 500.0, 800.0, 1000.0, 1500.0])   # K
Cp_target = np.array([18.63, 26.36, 35.13, 39.20, 45.77])  # cal/mol/K
H_target = 34.75   # kcal/mol  (ΔHf at 298 K)
S_target = 67.66   # cal/mol/K (S at 298 K)


def piece(thermo, T):
    """Return the NASAPolynomial that owns temperature `T` (low-T owns Tint)."""
    return thermo.polynomials[0] if T <= Tint else thermo.polynomials[1]


# %% [markdown]
# ## 3. Objective function: data residuals + soft continuity penalty
#
# Every residual is made dimensionless by dividing by a fixed reference scale, so
# the terms are comparable and the objective is a clean sum of squares:
#
# * **Data residuals** — relative error of each Cp point and of H(298), S(298).
# * **Continuity residuals** — the Cp, H and S gaps between the two polynomials at
#   `Tint`, normalised by fixed reference magnitudes (the target Cp(1000) and the
#   original H, S at `Tint`). A weight `w_cont` controls their importance.
#
# Because Cp, H and S are all *linear* in the NASA coefficients and every
# denominator here is a constant, the objective is exactly quadratic in any single
# coefficient — which makes each 1-D optimization well-behaved.

# %%
# Fixed reference magnitudes for normalising the continuity gaps (computed once
# from the original polynomial so they never change during fitting).
H_REF = abs(thermo0.polynomials[0].get_enthalpy(Tint) / KCAL)   # ~ kcal/mol at Tint
S_REF = abs(thermo0.polynomials[0].get_entropy(Tint) / CAL)     # ~ cal/mol/K at Tint

RESIDUAL_LABELS = (
    [f"Cp({T:.0f})" for T in Cp_T]
    + ["H(298)", "S(298)"]
    + ["cont:Cp", "cont:H", "cont:S"]
)


def residuals(thermo, w_cont=1.0):
    """Vector of dimensionless residuals (data points then continuity gaps)."""
    w_cont = 10.0 # over-riding it here!
    lo, hi = thermo.polynomials
    r = [piece(thermo, T).get_heat_capacity(T) / CAL / c - 1.0
         for T, c in zip(Cp_T, Cp_target)]
    r.append(lo.get_enthalpy(T_REF) / KCAL / H_target - 1.0)
    r.append(lo.get_entropy(T_REF) / CAL / S_target - 1.0)
    # continuity gaps at Tint, normalised by fixed reference scales
    r.append(w_cont * (lo.get_heat_capacity(Tint) - hi.get_heat_capacity(Tint)) / CAL / Cp_target[3])
    r.append(w_cont * (lo.get_enthalpy(Tint) - hi.get_enthalpy(Tint)) / KCAL / H_REF)
    r.append(w_cont * (lo.get_entropy(Tint) - hi.get_entropy(Tint)) / CAL / S_REF)
    return np.array(r)


def objective(thermo, w_cont=1.0):
    """Scalar least-squares objective: sum of squared residuals."""
    return float(np.sum(residuals(thermo, w_cont) ** 2))


# %% [markdown]
# How good is the published polynomial before we touch anything?

# %%
def report(thermo, w_cont=1.0):
    print(f"{'quantity':>10s}  {'model':>10s}  {'target':>10s}  {'rel.err':>9s}")
    for T, c in zip(Cp_T, Cp_target):
        m = piece(thermo, T).get_heat_capacity(T) / CAL
        print(f"  Cp({T:5.0f})  {m:10.3f}  {c:10.3f}  {m / c - 1:+9.2%}")
    h = thermo.polynomials[0].get_enthalpy(T_REF) / KCAL
    s = thermo.polynomials[0].get_entropy(T_REF) / CAL
    print(f"   H(298)   {h:10.3f}  {H_target:10.3f}  {h / H_target - 1:+9.2%}")
    print(f"   S(298)   {s:10.3f}  {S_target:10.3f}  {s / S_target - 1:+9.2%}")
    lo, hi = thermo.polynomials
    print(f"\n  continuity gaps at {Tint:.0f} K (low - high):")
    print(f"    Cp: {(lo.get_heat_capacity(Tint) - hi.get_heat_capacity(Tint)) / CAL:+.4f} cal/mol/K")
    print(f"     H: {(lo.get_enthalpy(Tint) - hi.get_enthalpy(Tint)) / KCAL:+.4f} kcal/mol")
    print(f"     S: {(lo.get_entropy(Tint) - hi.get_entropy(Tint)) / CAL:+.4f} cal/mol/K")
    print(f"\n  objective = {objective(thermo, w_cont):.6f}")


report(thermo0)
OBJ0 = objective(thermo0)

# %% [markdown]
# The published polynomial is already a good fit; the most visible flaws are a
# ~1.7% low S(298) and clear H and S discontinuities at 1000 K.

# %% [markdown]
# ## 4. Optimize each of the 14 coefficients independently
#
# The 14 parameters are the 7 coefficients of each polynomial. In RMG's ordering
# `coeffs = [c0, c1, c2, c3, c4, c5, c6]`:
#
# * `c0…c4` set the Cp(T) = R·(c0 + c1T + c2T² + c3T³ + c4T⁴) shape,
# * `c5` is the enthalpy integration constant (shifts H),
# * `c6` is the entropy integration constant (shifts S).
#
# Each coefficient is allowed to vary **freely** — far from its starting value and
# across zero if that helps — so a coefficient that starts at ~1e-19 is not
# confined to stay tiny. To keep all 14 searches numerically well-conditioned
# despite the coefficients spanning ~22 orders of magnitude, we optimize a
# perturbation `delta` measured in units of each coefficient's *natural* scale
# (the value that makes a unit change produce an O(1) change in Cp/R, H/R or S/R),
# which is set by the coefficient's role in the polynomial — **not** by its
# starting magnitude. The optimizer (`scipy.optimize.minimize_scalar`, unbounded
# Brent) is then free to send `delta` anywhere.

# %%
# Enumerate the 14 parameters as (polynomial index, coefficient index).
PARAMS = [(p, i) for p in (0, 1) for i in range(7)]
PARAM_NAMES = [f"{'low' if p == 0 else 'high'} c{i}" for p, i in PARAMS]

# Natural scale of each coefficient, independent of its starting value. In
# Cp/R = c0 + c1·T + c2·T² + c3·T³ + c4·T⁴, coefficient cᵢ multiplies Tⁱ, so its
# natural size is Tref^(-i); c5 is the H/R offset (units of K ~ Tref) and c6 is
# the dimensionless S/R offset (~1).
TREF = 1000.0


def natural_scale(i):
    if i <= 4:
        return TREF ** (-i)
    return TREF if i == 5 else 1.0


def optimize_single(p, i, w_cont=1.0):
    """Optimize only coefficient `i` of polynomial `p`, starting from the original.

    Returns a fresh NASA object that differs from `thermo0` in exactly one
    coefficient (the original is never mutated).
    """
    thermo = copy.deepcopy(thermo0)
    v0 = thermo0.polynomials[p].coeffs[i]
    scale = natural_scale(i)  # fixed scale, independent of the starting value

    def set_value(value):
        coeffs = thermo.polynomials[p].coeffs.copy()
        coeffs[i] = value
        thermo.polynomials[p].coeffs = coeffs

    def f(delta):
        set_value(v0 + delta * scale)
        return objective(thermo, w_cont)

    # Unbounded Brent: delta may grow large and cross zero, so the coefficient is
    # free to change sign and move orders of magnitude from where it started.
    res = minimize_scalar(f, method="brent")
    set_value(v0 + res.x * scale)  # leave `thermo` holding the optimum
    return thermo, v0, thermo.polynomials[p].coeffs[i], res.fun


results = []
for (p, i), name in zip(PARAMS, PARAM_NAMES):
    thermo_opt, v_old, v_new, obj = optimize_single(p, i)
    results.append({"name": name, "p": p, "i": i, "thermo": thermo_opt,
                    "v_old": v_old, "v_new": v_new, "obj": obj})

# Sanity check: the original object was never modified.
assert objective(thermo0) == OBJ0

# %% [markdown]
# ## 5. Summary: which single coefficient helps the most?

# %%
order = sorted(range(len(results)), key=lambda k: results[k]["obj"])
print(f"original objective = {OBJ0:.6f}\n")
print(f"{'parameter':>10s}  {'old value':>13s}  {'new value':>13s}  {'objective':>10s}  {'reduction':>9s}")
for k in order:
    r = results[k]
    print(f"{r['name']:>10s}  {r['v_old']:13.4e}  {r['v_new']:13.4e}  "
          f"{r['obj']:10.6f}  {1 - r['obj'] / OBJ0:+9.1%}")

# %%
fig, ax = plt.subplots(figsize=(9, 4.5))
names = [results[k]["name"] for k in order]
objs = [results[k]["obj"] for k in order]
ax.bar(range(len(objs)), objs, color="steelblue")
ax.axhline(OBJ0, color="crimson", ls="--", label=f"original = {OBJ0:.4f}")
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("objective (sum of squared residuals)")
ax.set_title("Best achievable objective when varying one coefficient at a time")
ax.legend()
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Cp, H and S curves for each single-parameter fit
#
# For every parameter we plot the resulting Cp, H and S curves (solid) against the
# original polynomial (grey dashed) and the target data (markers). Each polynomial
# piece is drawn over its own temperature range, so any discontinuity at
# `Tint = 1000 K` shows up as a visible jump.

# %%
T_lo = np.linspace(thermo0.Tmin.value_si, Tint, 200)
T_hi = np.linspace(Tint, thermo0.Tmax.value_si, 200)


def curves(thermo, prop):
    """Evaluate a property ('Cp', 'H', 'S') over each polynomial's own range."""
    lo, hi = thermo.polynomials
    if prop == "Cp":
        f, sc = "get_heat_capacity", CAL
    elif prop == "H":
        f, sc = "get_enthalpy", KCAL
    else:
        f, sc = "get_entropy", CAL
    y_lo = np.array([getattr(lo, f)(T) / sc for T in T_lo])
    y_hi = np.array([getattr(hi, f)(T) / sc for T in T_hi])
    return y_lo, y_hi


def plot_fit(result):
    """Three-panel Cp/H/S comparison for one single-parameter fit."""
    thermo = result["thermo"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    panels = [
        ("Cp", "Cp [cal/mol/K]", Cp_T, Cp_target),
        ("H", "H [kcal/mol]", np.array([T_REF]), np.array([H_target])),
        ("S", "S [cal/mol/K]", np.array([T_REF]), np.array([S_target])),
    ]
    for ax, (prop, ylabel, T_pts, y_pts) in zip(axes, panels):
        o_lo, o_hi = curves(thermo0, prop)
        n_lo, n_hi = curves(thermo, prop)
        ax.plot(T_lo, o_lo, color="grey", ls="--", lw=1)
        ax.plot(T_hi, o_hi, color="grey", ls="--", lw=1, label="original")
        ax.plot(T_lo, n_lo, color="C0", lw=1.8)
        ax.plot(T_hi, n_hi, color="C0", lw=1.8, label="optimized")
        ax.plot(T_pts, y_pts, "rs", ms=7, label="target", zorder=5)
        ax.axvline(Tint, color="k", lw=0.6, alpha=0.4)
        ax.set_xlabel("T [K]")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
    fig.suptitle(f"Vary '{result['name']}':  {result['v_old']:.4e} → "
                 f"{result['v_new']:.4e}   |   objective {OBJ0:.5f} → {result['obj']:.5f}",
                 fontsize=11)
    fig.tight_layout()
    plt.show()


for result in results:
    plot_fit(result)
