#!/usr/bin/env python
# encoding: utf-8

name = "Surface"
shortDesc = ""
longDesc = """
test surface mechanism: based upon Olaf Deutschmann's work:
"Surface Reaction Kinetics of Steam- and CO2-Reforming as well as Oxidation of Methane over Nickel-Based Catalysts"
Delgado et al
Catalysts, 2015, 5, 871-904
"""

entry(
    index=1,
    label="O2 + Ni + Ni <=> OX + OX",
    kinetics=StickingCoefficient(
        A=4.36e-2,
        n=-0.206,
        Ea=(1.5e3, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
        coverage_dependence={
            "OX": {"E": (0.0, "J/mol"), "m": 0.0, "a": 0.0},
        },
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)


entry(
    index=2,
    label="O2 + Ni <=> O2X",
    kinetics=StickingCoefficient(
        A=0.0,
        n=0,
        Ea=(0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)

entry(
    index=3,
    label="CH4 + Ni + Ni <=> CH3X + HX",
    kinetics=StickingCoefficient(
        A=8.0e-3,
        n=0,
        Ea=(0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)

entry(
    index=4,
    label="H2 + Ni + Ni <=> HX + HX",
    kinetics=StickingCoefficient(
        A=3.2e-2,
        n=0,
        Ea=(0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)


entry(
    index=5,
    label="HX + HOX <=> H2O + Ni + Ni",
    kinetics=SurfaceArrhenius(
        A=(1.85e16, "m^2/(mol*s)"),
        n=0.086,
        Ea=(41500.0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)


entry(
    index=6,
    label="CO + Ni <=> OCX",
    kinetics=StickingCoefficient(
        A=5.0e-1,
        n=0,
        Ea=(0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)


entry(
    index=7,
    label="OCX + OX <=> CO2 + Ni + Ni",
    kinetics=SurfaceArrhenius(
        A=(2.00e15, "m^2/(mol*s)"),
        n=0.0,
        Ea=(123600.0, "J/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="""Default""",
    longDesc="""Made up""",
    metal="Ni",
)

entry(
    index=8,
    label="CH4 + OX + Ni <=> CH3X + HOX",
    kinetics = StickingCoefficientBEP(
        A = 5e-6,
        n = 0,
        alpha = 0,
        E0 = (68.66, 'kJ/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc="StickingCoefficientBEP Eley-Rideal example",
    longDesc="""Because a gas phase molecule (CH4) reactss with an adsorbate,
                we use StickingCoefficientBEP.
                Copied from Surface_EleyRideal_Addition_Multiple_Bond.
                Note the kinetics are therefore for a reaction with one surface site
                but being used for one with two. They might be very inappropriate.""",
    metal="Ni",
)

entry(
    index=10,
    label="H2 + Ni + Ni <=> HX + HX",
    kinetics=StickingCoefficientBEP(
        A=0.01,
        n=0,
        alpha=0,
        E0=(10, "kcal/mol"),
        Tmin=(200, "K"),
        Tmax=(3000, "K"),
    ),
    shortDesc="BEP sticking coefficient example",
    longDesc="""This is an example reaction using the BEP kinetics model
                for a sticking coefficient.
                Copied from the Surface_Adsorption_Dissociative family.""",
    metal="Ni",
)
