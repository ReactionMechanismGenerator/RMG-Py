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
# # C4H6-2 (but-2-yne) Thermochemistry Fitting
#
# # Fitting a NASA polynomial one parameter at a time
#
# We import a Chemkin-format thermo entry for `C4H6-2` as an RMG
# [`NASA`](../rmgpy/thermo/nasa.pyx) object then ask:
#
# > Starting from the published polynomial, if we are allowed to adjust **only one**
# > of the 14 coefficients, how much closer can we get to a set of target
# > thermochemistry data?
#
# So we perform **independent single-parameter optimizations**. Each one starts
# from the original polynomial, optimizes exactly one parameter, records the
# result, and then *restores* it before moving to the next — every fitted
# polynomial therefore differs from the original in exactly one number.
#
# There are 15 parameters: the 7 coefficients of each polynomial (14 total) plus a
# **15th**, the *join temperature* `Tint` where the two pieces meet. Sections 4-6
# treat all 15 on the same footing; section 7 looks at the join temperature on its
# own. Holding `Tint` fixed at 1000 K (sections 1-3) is just the conventional
# default that the 15th parameter relaxes.
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
# Came from a Hidakwa J. Phys. Chem. 1993, 97, 10977-10983 10.1021/j100144a014 which attributed to 
# Handbook of Chemistry, 3rd ed.; Ed. by the Chemical Society of Japan: Maruzen: Tokyo, 1984; Chapter 9.
#
#
# plus the *additional goal* that Cp, H and S are continuous where the two
# polynomials meet (the common temperature `Tint = 1000 K`), included as a soft
# penalty in the objective.
#
# We reuse RMG's own machinery throughout: `rmgpy.chemkin.read_thermo_entry` to
# parse the Chemkin block, and `NASAPolynomial.get_heat_capacity / get_enthalpy /
# get_entropy` (which return SI units, J/mol·K and J/mol) to evaluate the curves.

# %% [markdown]
# "Detailed kinetic modeling of 1,3-butadiene oxidation at high temperatures" International Journal of Chemical Kinetics, 32, pp. 589-614 (2000).
# https://ignis.usc.edu:80/Mechanisms/C4H6/c4h6.html
#
# ```
# C4H612            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
#   0.1781557E 02 -0.4257502E-02  0.1051185E-04 -0.4473844E-08  0.5848138E-12    2
#   0.1267342E 05 -0.6982662E 02  0.1023467E 01  0.3495919E-01 -0.2200905E-04    3
#   0.6942272E-08 -0.7879187E-12  0.1811799E 05  0.1975066E 02  0.1950807E+05    4
# C4H6-2            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
#   9.0338133E+00 8.2124510E-03  7.1753952E-06 -5.8834334E-09  1.0343915E-12    2
#   1.4335068E+04 -2.0985762E+01  2.1373338E+00  2.6486229E-02 -9.0568711E-06    3
#  -5.5386397E-19 2.1281884E-22  1.5710902E+04  1.3529426E+01  1.7488676E+04    4
#  ```
#
#  note the line lengths are weird.
#
#
#  The thermodynamic properties of the species were taken mostly
# from previous compilations [4,46–48] and in part
# from individual works [2,31,36,49]. The heats of formation of several species were estimated using the
# NIST Structures and Properties code [50].
#
#
# 21. Hidaka, Y.; Higashihara, T.; Ninomiya, N.; Masaoka,
# H.; Nakamura, T.; Kawano, H. Int J Chem Kinet 1996,
# 28, 137.
# 31. Hidaka, Y.; Higashihara, T.; Ninomiya, N.; Oshita, H.;
# Kawano, H. J Phys Chem 1993, 97, 10977

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


def piece(thermo, T, tint=None):
    """Return the NASAPolynomial that owns temperature `T` (low-T owns the join).

    `tint` is the join temperature; it defaults to the published `Tint` (1000 K)
    but section 7 passes a different value to let the join float.
    """
    if tint is None:
        tint = Tint
    return thermo.polynomials[0] if T <= tint else thermo.polynomials[1]


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


def residuals(thermo, w_cont=1.0, tint=None):
    """Vector of dimensionless residuals (data points then continuity gaps).

    `tint` is the join temperature where the two pieces meet and where the
    continuity gaps are measured; it defaults to the published `Tint` (1000 K).
    """
    if tint is None:
        tint = Tint
    w_cont = 10.0 # over-riding it here!
    lo, hi = thermo.polynomials
    r = [piece(thermo, T, tint).get_heat_capacity(T) / CAL / c - 1.0
         for T, c in zip(Cp_T, Cp_target)]
    r.append(lo.get_enthalpy(T_REF) / KCAL / H_target - 1.0)
    r.append(lo.get_entropy(T_REF) / CAL / S_target - 1.0)
    # continuity gaps at the join, normalised by fixed reference scales
    r.append(w_cont * (lo.get_heat_capacity(tint) - hi.get_heat_capacity(tint)) / CAL / Cp_target[3])
    r.append(w_cont * (lo.get_enthalpy(tint) - hi.get_enthalpy(tint)) / KCAL / H_REF)
    r.append(w_cont * (lo.get_entropy(tint) - hi.get_entropy(tint)) / CAL / S_REF)
    return np.array(r)


def objective(thermo, w_cont=1.0, tint=None):
    """Scalar least-squares objective: sum of squared residuals."""
    return float(np.sum(residuals(thermo, w_cont, tint) ** 2))


# %% [markdown]
# How good is the published polynomial before we touch anything?

# %%
def report(thermo, w_cont=1.0, tint=None):
    if tint is None:
        tint = Tint
    print(f"{'quantity':>10s}  {'model':>10s}  {'target':>10s}  {'rel.err':>9s}")
    for T, c in zip(Cp_T, Cp_target):
        m = piece(thermo, T, tint).get_heat_capacity(T) / CAL
        print(f"  Cp({T:5.0f})  {m:10.3f}  {c:10.3f}  {m / c - 1:+9.2%}")
    h = thermo.polynomials[0].get_enthalpy(T_REF) / KCAL
    s = thermo.polynomials[0].get_entropy(T_REF) / CAL
    print(f"   H(298)   {h:10.3f}  {H_target:10.3f}  {h / H_target - 1:+9.2%}")
    print(f"   S(298)   {s:10.3f}  {S_target:10.3f}  {s / S_target - 1:+9.2%}")
    lo, hi = thermo.polynomials
    print(f"\n  continuity gaps at {tint:.0f} K (low - high):")
    print(f"    Cp: {(lo.get_heat_capacity(tint) - hi.get_heat_capacity(tint)) / CAL:+.4f} cal/mol/K")
    print(f"     H: {(lo.get_enthalpy(tint) - hi.get_enthalpy(tint)) / KCAL:+.4f} kcal/mol")
    print(f"     S: {(lo.get_entropy(tint) - hi.get_entropy(tint)) / CAL:+.4f} cal/mol/K")
    print(f"\n  objective = {objective(thermo, w_cont, tint):.6f}")


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


def objective_vs_tint(tint, w_cont=1.0):
    """Objective of the *original* polynomial as a function of the join `tint`."""
    return objective(thermo0, w_cont, tint)


def optimize_tint(w_cont=1.0, n_grid=2000):
    """Optimize only the join temperature, the 15th parameter, with all 14
    coefficients held at their published values.

    Moving the join changes two things: (a) which polynomial piece owns each Cp
    data point, and (b) the temperature at which the continuity gaps are
    measured. Effect (a) makes the objective piecewise — it jumps each time
    `tint` crosses one of the Cp data temperatures — so a single Brent search
    could land in the wrong piece. We therefore scan a dense grid across the
    full validity range for the global region, then refine locally with a
    bounded search inside the smooth interval around the grid winner.

    Returns a fresh NASA object whose join boundary has been moved to the
    optimum (the original `thermo0` is never mutated).
    """
    thermo = copy.deepcopy(thermo0)
    lo, hi = thermo.polynomials
    t_min = thermo0.Tmin.value_si
    t_max = thermo0.Tmax.value_si

    grid = np.linspace(t_min + 1.0, t_max - 1.0, n_grid)
    vals = np.array([objective(thermo, w_cont, t) for t in grid])
    t_seed = grid[int(np.argmin(vals))]

    # Refine within one grid step of the winner, staying inside one smooth piece.
    step = grid[1] - grid[0]
    res = minimize_scalar(lambda t: objective(thermo, w_cont, t),
                          method="bounded", bounds=(t_seed - step, t_seed + step))
    t_opt = res.x if res.fun <= objective(thermo, w_cont, t_seed) else t_seed

    # Leave `thermo` holding the optimum: move the shared join boundary so the
    # returned object is self-consistent (its own select_polynomial uses t_opt).
    lo.Tmax = (t_opt, "K")
    hi.Tmin = (t_opt, "K")
    return thermo, Tint, t_opt, objective(thermo, w_cont, t_opt)


results = []
for (p, i), name in zip(PARAMS, PARAM_NAMES):
    thermo_opt, v_old, v_new, obj = optimize_single(p, i)
    results.append({"name": name, "p": p, "i": i, "thermo": thermo_opt,
                    "v_old": v_old, "v_new": v_new, "obj": obj, "tint": Tint})

# Sanity check: the original object was never modified.
assert objective(thermo0) == OBJ0

# The 15th parameter: the join temperature itself (see section 7 for the full
# landscape). We append it to `results` so the summary table, bar chart and
# per-fit plots below treat it on the same footing as the 14 coefficients.
thermo_tint, t_old, t_new, obj_tint = optimize_tint()
results.append({"name": "Tint", "p": None, "i": None, "thermo": thermo_tint,
                "v_old": t_old, "v_new": t_new, "obj": obj_tint, "tint": t_new})

# %% [markdown]
# ## 5. Summary: which single parameter helps the most?
#
# The table and bar chart below rank all **15** single-parameter fits — the 14
# coefficients plus the join temperature `Tint`. For `Tint` the "old/new value"
# columns are temperatures in K rather than coefficient values.

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
ax.set_title("Best achievable objective when varying one parameter at a time")
ax.legend()
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Cp, H and S curves for each single-parameter fit
#
# For every parameter we plot the resulting Cp, H and S curves (solid) against the
# original polynomial (grey dashed) and the target data (markers). Each polynomial
# piece is drawn over its own temperature range, so any discontinuity at the join
# shows up as a visible jump. The last figure is the `Tint` fit: the two original
# pieces are unchanged but the join (black line) has moved off 1000 K.

# %%
def curves(thermo, prop, tint=None):
    """Evaluate a property ('Cp', 'H', 'S') over each polynomial's own range.

    The low piece is drawn over [Tmin, tint] and the high piece over [tint, Tmax],
    so moving the join (the 15th parameter) is reflected directly in the split.
    Returns the four arrays (T_lo, y_lo, T_hi, y_hi).
    """
    if tint is None:
        tint = Tint
    lo, hi = thermo.polynomials
    T_lo = np.linspace(thermo.Tmin.value_si, tint, 200)
    T_hi = np.linspace(tint, thermo.Tmax.value_si, 200)
    if prop == "Cp":
        f, sc = "get_heat_capacity", CAL
    elif prop == "H":
        f, sc = "get_enthalpy", KCAL
    else:
        f, sc = "get_entropy", CAL
    y_lo = np.array([getattr(lo, f)(T) / sc for T in T_lo])
    y_hi = np.array([getattr(hi, f)(T) / sc for T in T_hi])
    return T_lo, y_lo, T_hi, y_hi


def plot_fit(result):
    """Three-panel Cp/H/S comparison for one single-parameter fit."""
    thermo = result["thermo"]
    tint = result.get("tint", Tint)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    panels = [
        ("Cp", "Cp [cal/mol/K]", Cp_T, Cp_target),
        ("H", "H [kcal/mol]", np.array([T_REF]), np.array([H_target])),
        ("S", "S [cal/mol/K]", np.array([T_REF]), np.array([S_target])),
    ]
    for ax, (prop, ylabel, T_pts, y_pts) in zip(axes, panels):
        oT_lo, o_lo, oT_hi, o_hi = curves(thermo0, prop)          # original join (1000 K)
        nT_lo, n_lo, nT_hi, n_hi = curves(thermo, prop, tint)     # this fit's join
        ax.plot(oT_lo, o_lo, color="grey", ls="--", lw=1)
        ax.plot(oT_hi, o_hi, color="grey", ls="--", lw=1, label="original")
        ax.plot(nT_lo, n_lo, color="C0", lw=1.8)
        ax.plot(nT_hi, n_hi, color="C0", lw=1.8, label="optimized")
        ax.plot(T_pts, y_pts, "rs", ms=7, label="target", zorder=5)
        ax.axvline(tint, color="k", lw=0.6, alpha=0.4)
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

# %% [markdown]
# ## 7. The 15th parameter: a floating join temperature
#
# So far the join between the two polynomial pieces has been fixed at the
# published `Tint = 1000 K`. Here we treat that join as the **15th optimizable
# parameter** and ask the same sensitivity question: *holding all 14 coefficients
# at their published values, where should the two pieces meet to fit the data
# best?*
#
# This is qualitatively different from moving a coefficient. The two polynomials
# are unchanged; we only slide the temperature at which we switch from one to the
# other. That changes the objective in two ways:
#
# * **which piece owns a Cp point** — as `Tint` sweeps past one of the Cp data
#   temperatures, that point is suddenly evaluated on the *other* polynomial,
#   producing a step in the objective;
# * **where continuity is enforced** — the Cp/H/S gaps in the penalty are measured
#   at the join, and the two pieces drift apart away from 1000 K.
#
# The result is a piecewise objective with kinks/jumps at each Cp data
# temperature, which is why `optimize_tint` scans a dense grid before refining.

# %%
report(thermo_tint, tint=t_new)
print(f"\njoin temperature: {t_old:.1f} K -> {t_new:.1f} K")
print(f"objective: {OBJ0:.6f} -> {obj_tint:.6f}  ({1 - obj_tint / OBJ0:+.1%})")

# %% [markdown]
# ### The objective as a function of the join temperature
#
# Scanning `Tint` across the full validity range shows the landscape the optimizer
# searches: a mostly smooth curve (the continuity penalty growing as the pieces
# drift apart) punctuated by jumps where the join crosses a Cp data point.

# %%
t_scan = np.linspace(thermo0.Tmin.value_si + 1.0, thermo0.Tmax.value_si - 1.0, 2000)
obj_scan = np.array([objective_vs_tint(t) for t in t_scan])

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(t_scan, obj_scan, color="steelblue", lw=1.5)
ax.axhline(OBJ0, color="grey", ls=":", lw=1, label=f"published Tint = {OBJ0:.4f}")
ax.axvline(Tint, color="grey", ls=":", lw=1)
ax.plot([t_new], [obj_tint], "r*", ms=14, zorder=5,
        label=f"optimum: {t_new:.0f} K, obj = {obj_tint:.4f}")
for T in Cp_T:
    ax.axvline(T, color="k", lw=0.5, alpha=0.15)
ax.set_xlabel("join temperature Tint [K]")
ax.set_ylabel("objective (sum of squared residuals)")
ax.set_title("Objective vs. the 15th parameter (faint lines mark Cp data temperatures)")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Conclusion
# There is no single parameter typo that can fix the inconsistencies in the thermodynamic data.
# Let's instead replace the Cp(T) curves with more reliable ones, but preserve the H and S values that these models might rely on

# %% [markdown]
# ## Burcat 2005
#
# Burcat, A. & Ruscic, B. Third Millenium Ideal Gas and Condensed Phase Thermochemical Database for Combustion with Updates from Active Thermochemical Tables. (2005).
# https://publications.anl.gov/anlpubs/2005/07/53802.pdf
# https://respecth.elte.hu/burcat/BURCAT.THR.txt
#
#   
#

# %%
Burcat = """
503-17-3
C4H6  2-BUTAYN (DIMETHYLACETYLENE) STATWT=1  SIGMA=2  IA=1.0605  IB=IC=25.0029
Ir=0.2620  ROSYM=3 [V(3)=25. cm-1 REF=B3LYP]  (ONE ROTATION ONLY)   NU=213(2),
371(2),725,1125,1029(2),1050(2),1380(2),1448(2),1468(2),2270,2916,2966(2),
2976(3)  REF=YOST OSBORNE GARNER  JACS 63,(1941),3492   HF298=34.97 kcal
REF= Stull, Westrum & Sinke 1969  {HF298=146.314+/-4. kJ REF=Burcat G3B3 calc;
HF298=145.767+/-0.769 kJ  REF=ATcT A}   MAX LST SQ ERROR Cp @ 1300 K  0.60 %.
C4H6 Dimethyl Ac  A 1/05C  4.H  6.   0.   0.G   200.000  6000.000  B  54.09044 1
 7.26055302E+00 1.80160845E-02-6.47062409E-06 1.04411453E-09-6.24741250E-14    2
 1.39644246E+04-1.29484347E+01 5.39211846E+00 2.98346178E-03 5.22542032E-05    3
-6.64726627E-08 2.56305331E-11 1.55148209E+04 1.71080366E+00 1.75974868E+04    4
"""


# %%
from rmgpy.constants import R
R
H = 1.75974868E+04 * R / KCAL
print(f"Burcat has a enthalpy of formation of {H:.2f} kcal/mol")
# %%
Burcat_CK = """
C4H6-2            A 1/05C   4H   6    0    0G   200.00   6000.00   1000.0      1
 7.26055302E+00 1.80160845E-02-6.47062409E-06 1.04411453E-09-6.24741250E-14    2
 1.39644246E+04-1.29484347E+01 5.39211846E+00 2.98346178E-03 5.22542032E-05    3
-6.64726627E-08 2.56305331E-11 1.55148209E+04 1.71080366E+00                   4
""".strip()
read_thermo_entry(Burcat_CK)

# %%
# Import the Burcat 2005 Chemkin entry the same way as the original (section 1).
burcat_species, thermo_burcat, burcat_formula = read_thermo_entry(Burcat_CK)
burcat_tint = thermo_burcat.polynomials[0].Tmax.value_si
print(f"species : {burcat_species}")
print(f"formula : {burcat_formula}")
print(f"valid   : {thermo_burcat.Tmin.value_si:.0f}-{thermo_burcat.Tmax.value_si:.0f} K, "
      f"joined at {burcat_tint:.0f} K")

# Plot the original (USC-Mech) and Burcat 2005 polynomials together, each piece
# drawn over its own range, with the target data as markers.
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
panels = [
    ("Cp", "Cp [cal/mol/K]", Cp_T, Cp_target),
    ("H", "H [kcal/mol]", np.array([T_REF]), np.array([H_target])),
    ("S", "S [cal/mol/K]", np.array([T_REF]), np.array([S_target])),
]
curve_sources = [
    (thermo0, Tint, "C0", "original (USC)"),
    (thermo_burcat, burcat_tint, "C1", "Burcat 2005"),
]
for ax, (prop, ylabel, T_pts, y_pts) in zip(axes, panels):
    for thermo, tint, color, label in curve_sources:
        T_lo, y_lo, T_hi, y_hi = curves(thermo, prop, tint)
        mask = T_hi < 2000
        T_hi = T_hi[mask]
        y_hi = y_hi[mask]
        ax.plot(T_lo, y_lo, color=color, lw=1.8)
        ax.plot(T_hi, y_hi, color=color, lw=1.8, label=label)
    ax.plot(T_pts, y_pts, "rs", ms=7, label="target", zorder=5)
    ax.set_xlabel("T [K]")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)
fig.suptitle("C4H6-2 thermochemistry: original (USC-Mech) vs. Burcat 2005", fontsize=11)
fig.tight_layout()
plt.show()



# %%
# Gibbs free energy G = H - T*S over each piece's own range, both curves.
def gibbs_curves(thermo, tint):
    """Return (T_lo, G_lo, T_hi, G_hi) in kcal/mol for each polynomial piece."""
    lo, hi = thermo.polynomials
    T_lo = np.linspace(thermo.Tmin.value_si, tint, 200)
    T_hi = np.linspace(tint, thermo.Tmax.value_si, 200)
    G_lo = np.array([(lo.get_enthalpy(T) - T * lo.get_entropy(T)) / KCAL for T in T_lo])
    G_hi = np.array([(hi.get_enthalpy(T) - T * hi.get_entropy(T)) / KCAL for T in T_hi])
    return T_lo, G_lo, T_hi, G_hi


fig, ax = plt.subplots(figsize=(7, 4.5))
for thermo, tint, color, label in curve_sources:
    T_lo, G_lo, T_hi, G_hi = gibbs_curves(thermo, tint)
    mask = T_hi < 1500
    T_hi = T_hi[mask]
    G_hi = G_hi[mask]
    ax.plot(T_lo, G_lo, color=color, lw=1.8)
    ax.plot(T_hi, G_hi, color=color, lw=1.8, label=label)
ax.set_xlabel("T [K]")
ax.set_ylabel("G [kcal/mol]")
ax.set_title("C4H6-2 Gibbs free energy: original (USC) vs. Burcat 2005")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()


# %% [markdown]
# # What to do
# Use the original USC-mech entry (from the top of the notebook) to evaluate the H and S values at 298 K.
# Report these in kcal/mol and cal/mol/K.
# Derive a new set of polynomials that use the Burcat 2005 Cp data, but shifted so that the H and S values at 298 K match those from the original USC-mech entry. The joins at 100K should be continuous.
#

# %% [markdown]
# ## A "NEW" polynomial: Burcat Cp shape anchored to USC-mech H/S at 298 K
#
# We keep the Burcat 2005 **Cp shape** (the `c0…c4` of each piece) untouched and
# move only the two integration constants — `c5` (the enthalpy offset) and `c6`
# (the entropy offset). H and S are *linear* in those constants
# (`H = R·(… + c5)`, `S = R·(… + c6)`, so `dH/dc5 = dS/dc6 = R`), so shifting
# them is exact, reproduces any target H/S exactly, and leaves Cp unchanged.
#
# 1. Read H(298) and S(298) off the original USC-mech low-T piece.
# 2. Shift the Burcat **low-T** `c5`, `c6` so the new low-T piece reproduces those
#    two numbers at 298 K.
# 3. Shift the Burcat **high-T** `c5`, `c6` so H and S are continuous with the
#    (already shifted) low-T piece at the 1000 K join. (Cp is continuous there
#    already because we never touch `c0…c4`.)

# %%
from rmgpy.constants import R

# 1. Target H and S at 298 K, taken from the original USC-mech low-T polynomial.
H298_usc_si = thermo0.polynomials[0].get_enthalpy(T_REF)   # J/mol
S298_usc_si = thermo0.polynomials[0].get_entropy(T_REF)    # J/mol/K
print(f"USC-mech at 298 K:  H = {H298_usc_si / KCAL:.3f} kcal/mol   "
      f"S = {S298_usc_si / CAL:.3f} cal/mol/K")

# 2. Start from a fresh copy of Burcat and shift its integration constants.
thermo_new = copy.deepcopy(thermo_burcat)
lo_new, hi_new = thermo_new.polynomials


def shift_offsets(poly, T, H_target_si, S_target_si):
    """Move c5/c6 of `poly` so its H and S equal the SI targets at temperature T.

    Cp (set by c0…c4) is untouched, so only the H and S offsets move.
    """
    coeffs = poly.coeffs.copy()
    coeffs[5] += (H_target_si - poly.get_enthalpy(T)) / R   # dH/dc5 = R
    coeffs[6] += (S_target_si - poly.get_entropy(T)) / R    # dS/dc6 = R
    poly.coeffs = coeffs


# Low-T piece: match USC-mech H, S at 298 K.
shift_offsets(lo_new, T_REF, H298_usc_si, S298_usc_si)

# 3. High-T piece: match the (now shifted) low-T piece at the 1000 K join so H
#    and S are continuous across the boundary.
shift_offsets(hi_new, burcat_tint,
              lo_new.get_enthalpy(burcat_tint), lo_new.get_entropy(burcat_tint))

print(f"NEW       at 298 K:  H = {lo_new.get_enthalpy(T_REF) / KCAL:.3f} kcal/mol   "
      f"S = {lo_new.get_entropy(T_REF) / CAL:.3f} cal/mol/K")
print(f"\ncontinuity gaps at the {burcat_tint:.0f} K join (low - high):")
print(f"  Cp: {(lo_new.get_heat_capacity(burcat_tint) - hi_new.get_heat_capacity(burcat_tint)) / CAL:+.2e} cal/mol/K")
print(f"   H: {(lo_new.get_enthalpy(burcat_tint) - hi_new.get_enthalpy(burcat_tint)) / KCAL:+.2e} kcal/mol")
print(f"   S: {(lo_new.get_entropy(burcat_tint) - hi_new.get_entropy(burcat_tint)) / CAL:+.2e} cal/mol/K")

# %% [markdown]
# ### Three-panel comparison: USC-mech, Burcat 2005 and NEW
#
# Same Cp/H/S panels as before, now with all three curves and the temperature
# range limited to 2000 K. By construction the NEW curve shares Burcat's Cp shape
# while passing through the USC-mech H and S at 298 K.

# %%
TMAX_PLOT = 2500.0
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
panels = [
    ("Cp", "Cp [cal/mol/K]", Cp_T, Cp_target),
    ("H", "H [kcal/mol]", np.array([T_REF]), np.array([H_target])),
    ("S", "S [cal/mol/K]", np.array([T_REF]), np.array([S_target])),
]
curve_sources = [
    (thermo0, Tint, "C0", "original (USC)"),
    (thermo_burcat, burcat_tint, "C1", "Burcat 2005"),
    (thermo_new, burcat_tint, "C2", "NEW (Burcat Cp, USC H/S)"),
]
for ax, (prop, ylabel, T_pts, y_pts) in zip(axes, panels):
    for thermo, tint, color, label in curve_sources:
        T_lo, y_lo, T_hi, y_hi = curves(thermo, prop, tint)
        # Draw each piece over its own range, clipped to the 2000 K plot limit.
        # Only the high piece carries the legend label (one entry per curve).
        for T_arr, y_arr, lab in ((T_lo, y_lo, None), (T_hi, y_hi, label)):
            mask = T_arr <= TMAX_PLOT
            ax.plot(T_arr[mask], y_arr[mask], color=color, lw=1.8, label=lab)
    ax.plot(T_pts, y_pts, "rs", ms=7, label="target", zorder=5)
    ax.set_xlim(right=TMAX_PLOT)
    ax.set_xlabel("T [K]")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)
fig.suptitle("C4H6-2: original (USC) vs. Burcat 2005 vs. NEW "
             "(Burcat Cp shifted to USC H/S at 298 K)", fontsize=11)
fig.tight_layout()
plt.show()

# %%
thermo_new.comment="H298 and S298 values from USC-Mech-ii,JetSurF,AramcoMech1.3,etc.;\nCp from Burcat 2005, because USC-Mech inconsistent and discontinuous."


# %%
thermo_new

# %%
try:
    import black
    print(black.format_str(repr(thermo_new), mode=black.Mode(line_length=110)))
except ImportError:
    print("Black is not installed. Skipping formatting.")

# %% [markdown]
# ## Differences between the NEW and original entries
#
# Plot the change (NEW − original) in H, S and G as a function of temperature, up
# to 2500 K. The `NASA` object selects the correct polynomial piece for each T, so
# these curves span both pieces of each entry. By construction H and S match at
# 298 K (so the H and S differences start near zero there); they drift apart at
# higher T because the two entries use different Cp shapes.

# %%
T_diff = np.linspace(T_REF, TMAX_PLOT, 300)
dH = np.array([(thermo_new.get_enthalpy(T) - thermo0.get_enthalpy(T)) / KCAL for T in T_diff])
dS = np.array([(thermo_new.get_entropy(T) - thermo0.get_entropy(T)) / CAL for T in T_diff])
dG = np.array([((thermo_new.get_enthalpy(T) - T * thermo_new.get_entropy(T))
                - (thermo0.get_enthalpy(T) - T * thermo0.get_entropy(T))) / KCAL for T in T_diff])

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, y, ylabel in zip(
    axes,
    (dH, dS, dG),
    ("ΔH [kcal/mol]", "ΔS [cal/mol/K]", "ΔG [kcal/mol]"),
):
    ax.plot(T_diff, y, color="C3", lw=1.8)
    ax.axhline(0.0, color="grey", ls=":", lw=1)
    ax.set_xlabel("T [K]")
    ax.set_ylabel(ylabel)
fig.suptitle("NEW − original (USC-mech) differences vs. T", fontsize=11)
fig.tight_layout()
plt.show()

# %%
