#!/usr/bin/env python
# encoding: utf-8

name = "Baeyer-Villiger_step1_cat/training"
shortDesc = "Kinetics used to train group additivity values"
longDesc = """
Put kinetic parameters for reactions to use as a training set for fitting
group additivity values in this file.
"""

entry(
    index=1,
    label="acetone + peracetic_acid + acetic_acid1 <=> acetone_peracetic_criegee + acetic_acid2",
    degeneracy=1.0,
    kinetics=Arrhenius(
        A=(1.07543e-11, "cm^6/(mol^2*s)"),
        n=5.47295,
        Ea=(-38.5379, "kJ/mol"),
        T0=(1, "K"),
    ),
    rank=6,
    shortDesc="""CBS-QB3 calculation without HR""",
    longDesc="""
CBS-QB3 calculation without HR fitted over range from 300-600 K
""",
)

entry(
    index=2,
    label="cyclohexanone + peracetic_acid + acetic_acid1 <=> cyclohexanone_peracetic_criegee + acetic_acid2",
    degeneracy=1.0,
    kinetics=Arrhenius(
        A=(1.32822e-11, "cm^6/(mol^2*s)"),
        n=5.49341,
        Ea=(-44.5298, "kJ/mol"),
        T0=(1, "K"),
    ),
    rank=6,
    shortDesc="""CBS-QB3 calculation without HR""",
    longDesc="""
CBS-QB3 calculation without HR fitted over range from 300-600 K
""",
)

entry(
    index=3,
    label="acetone + methylhydroperoxide + acetic_acid1 <=> acetone_methyl_criegee + acetic_acid2",
    degeneracy=1.0,
    kinetics=Arrhenius(
        A=(2.6104e-09, "cm^6/(mol^2*s)"),
        n=4.30497,
        Ea=(-30.1492, "kJ/mol"),
        T0=(1, "K"),
    ),
    rank=6,
    shortDesc="""CBS-QB3 calculation without HR""",
    longDesc="""
CBS-QB3 calculation without HR fitted over range from 300-600 K
""",
)

entry(
    index=4,
    label="cyclohexanone + methylhydroperoxide + acetic_acid1 <=> cyclohexanone_methyl_criegee + acetic_acid2",
    degeneracy=1.0,
    kinetics=Arrhenius(
        A=(5.58493e-09, "cm^6/(mol^2*s)"),
        n=4.34471,
        Ea=(-35.857, "kJ/mol"),
        T0=(1, "K"),
    ),
    rank=6,
    shortDesc="""CBS-QB3 calculation without HR""",
    longDesc="""
CBS-QB3 calculation without HR fitted over range from 300-600 K
""",
)
