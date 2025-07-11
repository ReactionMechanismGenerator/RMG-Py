#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Abstraction_vdW/training"
shortDesc = "Reaction kinetics used to generate rate rules"
longDesc = """
Put kinetic parameters for specific reactions in this file to use as a
training set for generating rate rules to populate this kinetics family.
"""

# reverse of 16, below
# entry(
#     index = 34,
#     label = "H2O* + O* <=> OH_2* + OH_4*",
#     degeneracy = 2,
#     kinetics = SurfaceArrhenius(
#         A = (8.14E20, 'm^2/(mol*s)'),
#         n = -0.274,
#         Ea = (218400, 'J/mol'),
#         Tmin = (200, 'K'),
#         Tmax = (3000, 'K'),
#     ),
#     rank = 10,
#     shortDesc = u"""Default""",
#     longDesc = u"""R34
#     test surface mechanism: based upon Olaf Deutschmann's work:
#     "Surface Reaction Kinetics of Steam- and CO2-Reforming as well as Oxidation of Methane over Nickel-Based Catalysts"
#     Delgado et al
#     Catalysts, 2015, 5, 871-904""",
# 	  metal = "Ni",
# )

# reverse of 34, above
entry(
    index=16,
    label="OH_2* + OH_4* <=> H2O* + O*",
    degeneracy=2,
    kinetics=SurfaceArrhenius(
        A=(5.691e16, "m^2/(mol*s)"),
        n=0.0,
        Ea=(14.0669343, "kcal/mol"),
        Tmin=(298, "K"),
        Tmax=(2000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""
Reaction 16 from table 2 in "Mechanism of Methanol Synthesis on Cu through CO2
and CO Hydrogenation", Grabow and Mavrikakis.  doi:10.1021/cs200055d

A factor from paper / surface site density of Cu
1.675e12 m^4/(mol^2 * s) / 2.943e‐5 mol/m^2 = 5.691e16 m^2/(mol*s)
""",
    metal="Cu",
)

entry(
    index=21,
    label="CH4* + O* <=> CH3* + OH_4*",
    degeneracy=4,
    kinetics=SurfaceArrhenius(
        A=(5.62e20, "m^2/(mol*s)"),
        n=-0.101,
        Ea=(92700, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""R21
    test surface mechanism: based upon Olaf Deutschmann's work:
    "Surface Reaction Kinetics of Steam- and CO2-Reforming as well as Oxidation of Methane over Nickel-Based Catalysts"
    Delgado et al
    Catalysts, 2015, 5, 871-904""",
    metal="Ni",
)

entry(
    index=40,
    label="OH_2* + HCO* <=> H2O* + CO*",
    degeneracy=1,
    kinetics=SurfaceArrhenius(
        A=(3.261e17, "m^2/(mol*s)"),
        n=0.0,
        Ea=(6.9181644, "kcal/mol"),
        Tmin=(298, "K"),
        Tmax=(2000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""
Reaction 40 from table 2 in "Mechanism of Methanol Synthesis on Cu through CO2
and CO Hydrogenation", Grabow and Mavrikakis.  doi:10.1021/cs200055d

A factor from paper / surface site density of Cu
9.597e12 m^4/(mol^2 * s) / 2.943e‐5 mol/m^2 = 3.261e17 m^2/(mol*s)
""",
    metal="Cu",
)

entry(
    index=41,
    label="HCOO_1* + HCO* <=> HCOOH* + CO*",
    degeneracy=1,
    kinetics=SurfaceArrhenius(
        A=(7.475e18, "m^2/(mol*s)"),
        n=0.0,
        Ea=(13.8363288, "kcal/mol"),
        Tmin=(298, "K"),
        Tmax=(2000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""
Reaction 41 from table 2 in "Mechanism of Methanol Synthesis on Cu through CO2
and CO Hydrogenation", Grabow and Mavrikakis.  doi:10.1021/cs200055d

A factor from paper / surface site density of Cu
2.2e14 m^4/(mol^2 * s) / 2.943e‐5 mol/m^2 = 7.475e18 m^2/(mol*s)
""",
    metal="Cu",
)

entry(
    index=45,
    label="CH3O* + HCO* <=> CH3OH* + CO*",
    degeneracy=1,
    kinetics=SurfaceArrhenius(
        A=(6.572e16, "m^2/(mol*s)"),
        n=0.0,
        Ea=(8.76300824, "kcal/mol"),
        Tmin=(298, "K"),
        Tmax=(2000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""
Reaction 45 from table 2 in "Mechanism of Methanol Synthesis on Cu through CO2
and CO Hydrogenation", Grabow and Mavrikakis.  doi:10.1021/cs200055d

A factor from paper / surface site density of Cu
1.934e12 m^4/(mol^2 * s) / 2.943e‐5 mol/m^2 = 6.572e16 m^2/(mol*s)
""",
    metal="Cu",
)

entry(
    index=46,
    label="CH3O* + HCOO_5* <=> HCOOCH3* + O*",
    degeneracy=1,
    kinetics=SurfaceArrhenius(
        A=(2.356e16, "m^2/(mol*s)"),
        n=0.0,
        Ea=(28.5950795, "kcal/mol"),
        Tmin=(298, "K"),
        Tmax=(2000, "K"),
    ),
    rank=10,
    shortDesc="""Default""",
    longDesc="""
Reaction 46 from table 2 in "Mechanism of Methanol Synthesis on Cu through CO2
and CO Hydrogenation", Grabow and Mavrikakis.  doi:10.1021/cs200055d

A factor from paper / surface site density of Cu
6.934e11 m^4/(mol^2 * s) / 2.943e‐5 mol/m^2 = 2.356e16 m^2/(mol*s)
""",
    metal="Cu",
)
