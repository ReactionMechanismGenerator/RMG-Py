#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Proton_Electron_Reduction_Alpha/training"
shortDesc = u"Reaction kinetics used to generate rate rules"
longDesc = u"""
Put kinetic parameters for specific reactions in this file to use as a
training set for generating rate rules to populate this kinetics family.
"""

# entry(
#     index = 1,
#     label = "CX + H <=> CHX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.20, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.24, 'V'), # reference potential
#         Ea = (0.61, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 2,
    label = "CX + H <=> CHX_p",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.20, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.44, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 1,
#     label = "CX + H <=> CHX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.06, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.29, 'V'), # reference potential
#         Ea = (0.19, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 1,
#     label = "CX + H <=> CHX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.06, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.5, 'V'), # reference potential
#         Ea = (0.13, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 3,
#     label = "CHX + H <=> CH2X_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.31, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.32, 'V'), # reference potential
#         Ea = (0.77, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 3,
    label = "CHX + H <=> CH2X_p",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.31, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.44, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 3,
#     label = "CHX + H <=> CH2X_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.05, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.30, 'V'), # reference potential
#         Ea = (0.59, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 4,
#     label = "CHX + H <=> CH2X_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.05, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.5, 'V'), # reference potential
#         Ea = (0.53, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 6,
#     label = "CH2X + H <=> CH3X",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.23, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.38, 'V'), # reference potential
#         Ea = (0.62, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 4,
    label = "CH2X + H <=> CH3X",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.23, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.37, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 8,
#     label = "NX + H <=> HNX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.07, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.12, 'V'), # reference potential
#         Ea = (0.15, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2017.01.050""",
#     longDesc = u"""
# """,
#     metal = "Cu",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 9,
#     label = "NX + H <=> HNX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.23, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.24, 'V'), # reference potential
#         Ea = (0.78, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 5,
    label = "NX + H <=> HNX_p",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.23, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.59, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 10,
#     label = "NX + H <=> HNX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.07, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.3, 'V'), # reference potential
#         Ea = (0.17, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 4,
#     label = "NX + H <=> HNX_p",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.07, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.5, 'V'), # reference potential
#         Ea = (0.09, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 11,
#     label = "HNX + H <=> H2NX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.27, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.28, 'V'), # reference potential
#         Ea = (1.20, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 6,
    label = "HNX + H <=> H2NX",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.27, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.97, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 11,
#     label = "HNX + H <=> H2NX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.64, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.31, 'V'), # reference potential
#         Ea = (0.99, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 5,
#     label = "HNX + H <=> H2NX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.64, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.5, 'V'), # reference potential
#         Ea = (0.43, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 12,
#     label = "OX + H <=> HOX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.41, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.31, 'V'), # reference potential
#         Ea = (0.87, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Tafel""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

entry(
    index = 7,
    label = "OX + H <=> HOX",
    degeneracy = 1,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.42, # charge transfer coeff
        A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
        n = 0, # temperature coeff
        V0 = (-0.5, 'V'), # reference potential
        Ea = (0.48, 'eV/molecule'), # activation energy
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
        electrons = -1, # electron stochiometric coeff 
    ),
    shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
    longDesc = u"""Tafel""",
    metal = "Pt",
    facet = "111",
    site = "",
    rank = 5,
)

# entry(
#     index = 12,
#     label = "OX + H <=> HOX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.02, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (0.57, 'V'), # reference potential
#         Ea = (0.06, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )

# entry(
#     index = 6,
#     label = "OX + H <=> HOX",
#     degeneracy = 1,
#     kinetics = SurfaceChargeTransfer(
#         alpha = 0.02, # charge transfer coeff
#         A = (2.5e14, 'm^3/(mol*s)'), # pre-exponential factor estimate 10^11 s^-1 * 2.5e5 m^2/mol / 1000 m^3/mol H+
#         n = 0, # temperature coeff
#         V0 = (-0.5, 'V'), # reference potential
#         Ea = (0.03, 'eV/molecule'), # activation energy
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#         electrons = -1, # electron stochiometric coeff 
#     ),
#     shortDesc = u"""https://doi.org/10.1016/j.cattod.2018.03.048""",
#     longDesc = u"""Heyrovsky""",
#     metal = "Pt",
#     facet = "111",
#     site = "",
#     rank = 5,
# )
