"""Numpy-only reference consumer for the polymer moments artifact.

Implements docs/polymer_moments_format.md §4-§7 with numpy + stdlib ONLY.
THIS MODULE MUST NOT IMPORT rmgpy (that is the artifact's entire point);
test_consumer_module_is_rmgpy_free enforces it.

Rate constants come from each entry's own `kinetics` (every entry in the
test decks is emitted cantera-null), so no Cantera is needed here either;
Keq for reversible entries is computed from caller-supplied NASA7 data via
the documented recipe.
"""

import numpy as np

# MUST equal the oracle's rmgpy.constants.R (CODATA-2006), NOT the 2018-SI
# exact value 8.31446261815324: a mismatched R shifts every gas-phase / Keq /
# mass-transfer term by a constant 1.13e-6 (format doc §4 definitions).
R = 8.314472          # J/(mol K)
P0 = 1.0e5            # Pa, Keq reference pressure
SMALL_EPS = 1e-30
LN_EXP_OVERFLOW_GUARD = 700.0


def safe_mu3(mu0, mu1, mu2):
    """log_lagrange/1 closure with realizability guard (format doc §6)."""
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    if mu1 < mu0:
        return 0.0
    ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
    if ln_mu3 > LN_EXP_OVERFLOW_GUARD:
        return float("inf")
    return float(np.exp(ln_mu3))


def nasa_g_over_rt(coeffs, T):
    """G/(R*T) from one NASA7 row [a0..a4, a5, a6]."""
    a0, a1, a2, a3, a4, a5, a6 = coeffs
    h_rt = (a0 + a1 * T / 2.0 + a2 * T**2 / 3.0 + a3 * T**3 / 4.0
            + a4 * T**4 / 5.0 + a5 / T)
    s_r = (a0 * np.log(T) + a1 * T + a2 * T**2 / 2.0 + a3 * T**3 / 3.0
           + a4 * T**4 / 4.0 + a6)
    return h_rt - s_r


def _base(label):
    return label.partition("(")[0]


class ArtifactConsumer:
    """Integrates the artifact's pool moments + species ODEs.

    Parameters
    ----------
    artifact : dict           parsed polymer_pools.json (schema 2.0)
    species_order : [str]     chem.yaml labels defining the y-vector layout
                              (same order as the oracle's core species)
    P : float                 pressure [Pa]
    V_poly : float            condensed-phase volume [m^3]
    mass_transfer : [dict]    [{"gas": label, "poly": label, "K": f, "kLa": f}]
    nasa : {label: {"Tmid": f, "low": [7], "high": [7]}}, optional
                              NASA7 data for Keq of reversible entries
    """

    def __init__(self, artifact, species_order, P, V_poly,
                 mass_transfer=None, nasa=None):
        self.P = float(P)
        self.V_poly = float(V_poly)
        self.nasa = nasa or {}
        self.idx = {lab: i for i, lab in enumerate(species_order)}
        n = len(species_order)

        conv = artifact["conventions"]
        condensed = set(conv["condensed_species"])
        self.gas_mask = np.array([lab not in condensed for lab in species_order],
                                 dtype=bool)
        self.configured = list(conv["configured_pools"])

        # pools: label -> dict(mu indices, channels, monomer routing index)
        self.pools = {}
        for p in artifact["pools"]:
            lab = p["label"]
            if lab not in self.configured:
                continue
            mu_idx = tuple(self.idx[f"{lab}_mu{k}"] for k in range(3))
            routing = p.get("monomer_routing")
            self.pools[lab] = {
                "mu": mu_idx,
                "k_s": p["channels"]["scission"]["A"],
                "k_u": p["channels"]["unzip"]["A"],
                "routing": self.idx[routing] if routing else None,
            }

        # reactions: precompute index forms
        self.entries = []
        for e in artifact["reactions"]:
            kin = e["kinetics"]
            assert kin is not None, f"entry {e['id']} has no kinetics"
            ridx = [self.idx[s] for s in e["reactants"]]
            pidx = [self.idx[s] for s in e["products"]]
            pool_mapped = ({self.idx[s] for s in e["proxy_reactants"]}
                           | {self.idx[s] for s in e["proxy_products"]})
            self.entries.append({
                "A": kin["A"], "n": kin["n"], "Ea": kin["Ea"],
                "reversible": kin["reversible"],
                "ridx": ridx, "pidx": pidx, "pool_mapped": pool_mapped,
                "r_labels": list(e["reactants"]), "p_labels": list(e["products"]),
                "proxy_r_pools": [_base(s) for s in e["proxy_reactants"]],
                "proxy_p_pools": [_base(s) for s in e["proxy_products"]],
                "scaling": e["scaling"],
                "src": e["src_pool"], "dst": e["dst_pool"],
                "arch": e["archetype"],
                "a": int(e.get("params", {}).get("a", 0)),
            })

        self.mass_transfer = []
        for mt in (mass_transfer or []):
            self.mass_transfer.append((self.idx[mt["gas"]], self.idx[mt["poly"]],
                                       float(mt["K"]), float(mt["kLa"])))

    # ----- helpers ------------------------------------------------------

    def _keq(self, entry, T):
        def g(label):
            d = self.nasa[_label_lookup(self.nasa, label)]
            row = d["low"] if T <= d["Tmid"] else d["high"]
            return nasa_g_over_rt(row, T)  # G/(R T)
        dg_rt = sum(g(s) for s in entry["p_labels"]) - sum(g(s) for s in entry["r_labels"])
        # dn counts ALL species as written (condensed/proxy included) —
        # format doc §4 step 1.
        dn = len(entry["p_labels"]) - len(entry["r_labels"])
        return (P0 / (R * T)) ** dn * np.exp(-dg_rt)

    def _chain_bundle(self, pool, y, end_group):
        i0, i1, i2 = self.pools[pool]["mu"]
        mu0 = max(0.0, y[i0]) / self.V_poly
        mu1 = max(0.0, y[i1]) / self.V_poly
        mu2 = max(0.0, y[i2]) / self.V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = safe_mu3(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False

    # ----- RHS (format doc §4-§7) ----------------------------------------

    def rhs(self, y, T):
        dn = np.zeros_like(y)
        n_gas = float(np.sum(np.clip(y[self.gas_mask], 0.0, None)))
        V_gas = n_gas * R * T / self.P if n_gas > 0 else 1.0
        Vp = self.V_poly
        C = np.where(self.gas_mask,
                     np.clip(y, 0.0, None) / V_gas,
                     np.clip(y, 0.0, None) / Vp)

        for e in self.entries:
            kf = e["A"] * T ** e["n"] * np.exp(-e["Ea"] / (R * T))
            kb = kf / self._keq(e, T) if e["reversible"] else 0.0

            # step 2: phase + gate
            is_poly_event = any(not self.gas_mask[i] for i in e["ridx"])
            V_rxn = Vp if is_poly_event else V_gas
            has_poly_prod = any(not self.gas_mask[i] for i in e["pidx"])
            if is_poly_event and not has_poly_prod:
                continue
            if (not is_poly_event) and has_poly_prod:
                continue

            # step 3: concentration products
            rf = kf
            for i in e["ridx"]:
                rf *= 1.0 if i in e["pool_mapped"] else C[i]
            rr = kb
            for i in e["pidx"]:
                rr *= 1.0 if i in e["pool_mapped"] else C[i]

            # step 4: site scaling
            if e["src"] is not None:
                i0, i1, _ = self.pools[e["src"]]["mu"]
                if e["arch"] == "discrete_chip/1" and e["scaling"] == "mu0" and e["a"] > 0:
                    site = min(max(0.0, y[i0]), max(0.0, y[i1]) / e["a"]) / Vp
                else:
                    mi = i0 if e["scaling"] == "mu0" else i1
                    site = max(0.0, y[mi]) / Vp
                rf *= site
                rr *= site

            r_mol = (rf - rr) * V_rxn

            # step 5: stoichiometric flux for non-pool-mapped species
            for i in e["ridx"]:
                if i not in e["pool_mapped"]:
                    dn[i] -= r_mol
            for i in e["pidx"]:
                if i not in e["pool_mapped"]:
                    dn[i] += r_mol

            # step 6: archetype bundles
            arch = e["arch"]
            if arch == "migration/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    for ev, frm, to in ((rf, src, dst), (rr, dst, src)):
                        if ev <= 0.0:
                            continue
                        ev_mol = ev * V_rxn
                        b0, b1, b2, ok = self._chain_bundle(frm, y, e["scaling"] == "mu0")
                        if b0 == 0.0:
                            continue
                        f = self.pools[frm]["mu"]
                        t = self.pools[to]["mu"]
                        dn[f[0]] -= ev_mol * b0
                        dn[f[1]] -= ev_mol * b1
                        dn[t[0]] += ev_mol * b0
                        dn[t[1]] += ev_mol * b1
                        if ok:
                            dn[f[2]] -= ev_mol * b2
                            dn[t[2]] += ev_mol * b2
            elif arch == "scission_fragment/1":
                src, dst = e["src"], e["dst"]
                if src and dst and src != dst:
                    s = self.pools[src]["mu"]
                    d = self.pools[dst]["mu"]
                    mu0p = max(0.0, y[s[0]]) / Vp
                    mu1p = max(0.0, y[s[1]]) / Vp
                    mu2p = max(0.0, y[s[2]]) / Vp
                    ok = mu1p > SMALL_EPS
                    if ok and r_mol < 0.0:
                        if (max(0.0, y[d[0]]) / Vp <= SMALL_EPS
                                or max(0.0, y[d[1]]) / Vp <= SMALL_EPS):
                            ok = False
                    if ok:
                        e_n = mu2p / mu1p
                        dn[s[1]] -= r_mol * e_n / 2.0
                        dn[d[0]] += r_mol
                        dn[d[1]] += r_mol * e_n / 2.0
                        mu3p = safe_mu3(mu0p, mu1p, mu2p)
                        if np.isfinite(mu3p):
                            e_n2 = mu3p / mu1p
                            dn[s[2]] -= r_mol * (2.0 / 3.0) * e_n2
                            dn[d[2]] += r_mol * e_n2 / 3.0
            elif arch == "discrete_chip/1":
                src = e["src"]
                if src:
                    a = float(e["a"])
                    b0, b1, _b2, _ok = self._chain_bundle(src, y, e["scaling"] == "mu0")
                    if b0 != 0.0:
                        s = self.pools[src]["mu"]
                        e_n = b1
                        if rf > 0.0:
                            rf_mol = rf * V_rxn
                            dn[s[1]] -= rf_mol * a
                            dmu2 = 2.0 * a * e_n - a * a
                            if dmu2 > 0.0:
                                dn[s[2]] -= rf_mol * dmu2
                        if rr > 0.0:
                            rr_mol = rr * V_rxn
                            dn[s[1]] += rr_mol * a
                            dn[s[2]] += rr_mol * (2.0 * a * e_n + a * a)
            elif arch == "legacy_mu1/1":
                for pool in e["proxy_r_pools"]:
                    dn[self.pools[pool]["mu"][1]] -= r_mol
                for pool in e["proxy_p_pools"]:
                    dn[self.pools[pool]["mu"][1]] += r_mol
            # same_pool/1: no moment flux

        # step 7a: channels (format doc §5)
        for pool in self.pools.values():
            i0, i1, i2 = pool["mu"]
            mu0 = max(0.0, y[i0]) / Vp
            mu1 = max(0.0, y[i1]) / Vp
            mu2 = max(0.0, y[i2]) / Vp
            dmu0 = dmu1 = dmu2 = 0.0
            if pool["k_s"] > 0.0:
                mu3 = safe_mu3(mu0, mu1, mu2)
                dmu0 += pool["k_s"] * (mu1 - mu0)
                if np.isfinite(mu3):
                    dmu2 += pool["k_s"] * (mu1 - mu3) / 3.0
            if pool["k_u"] > 0.0:
                r_ev = pool["k_u"] * mu0
                dmu1 -= r_ev
                dmu2 -= pool["k_u"] * (2.0 * mu1 - mu0)
                if pool["routing"] is not None:
                    dn[pool["routing"]] += r_ev * Vp
            dn[i0] += dmu0 * Vp
            dn[i1] += dmu1 * Vp
            dn[i2] += dmu2 * Vp

        # step 7b: mass transfer (format doc §7)
        for ig, ip, K, kLa in self.mass_transfer:
            J = kLa * (C[ip] - K * C[ig])
            dq = J * Vp
            dn[ig] += dq
            dn[ip] -= dq

        return dn

    def integrate_euler(self, y0, T, dt, n_steps, record_every=1):
        """Fixed-step forward Euler; returns (times, trajectory[n_rec, n_spc])."""
        y = np.array(y0, dtype=float)
        times, traj = [0.0], [y.copy()]
        for k in range(n_steps):
            y = y + dt * self.rhs(y, T)
            if (k + 1) % record_every == 0:
                times.append((k + 1) * dt)
                traj.append(y.copy())
        return np.array(times), np.array(traj)


def _label_lookup(nasa, label):
    """NASA data may be keyed by full chem.yaml label or its base form."""
    if label in nasa:
        return label
    base = _base(label)
    if base in nasa:
        return base
    raise KeyError(f"no NASA thermo for {label!r}")
