#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_Carbene_Intra_Disproportionation/training"
shortDesc = u"Kinetics used to train group additivity values"
longDesc = u"""
Put kinetic parameters for reactions to use as a training set for fitting
group additivity values in this file.
"""



entry(
    index = 1,
    label = "C6H6 <=> C6H6-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.067e+10, 's^-1'), n=0.649, Ea=(8.03, 'kcal/mol'), T0=(1, 'K')),
    rank = 3,
    shortDesc = u"""Training reaction from kinetics library: 2003_Miller_Propargyl_Recomb_High_P""",
    longDesc = 
u"""
Taken from entry: A <=> IV
""",
)

entry(
    index = 2,
    label = "C6H6-3 <=> C6H6-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.454e+12, 's^-1'), n=0.178, Ea=(0.205, 'kcal/mol'), T0=(1, 'K')),
    rank = 3,
    shortDesc = u"""Training reaction from kinetics library: 2003_Miller_Propargyl_Recomb_High_P""",
    longDesc = 
u"""
Taken from entry: IX <=> VII
""",
)

entry(
    index = 3,
    label = "C6H6-5 <=> C6H6-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.865e+11, 's^-1'), n=0.577, Ea=(29.169, 'kcal/mol'), T0=(1, 'K')),
    rank = 3,
    shortDesc = u"""Training reaction from kinetics library: 2003_Miller_Propargyl_Recomb_High_P""",
    longDesc = 
u"""
Taken from entry: X <=> IX
""",
)

entry(
    index = 4,
    label = "C6H6-7 <=> C6H6-8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.355e+12, 's^-1'), n=0.294, Ea=(35.954, 'kcal/mol'), T0=(1, 'K')),
    rank = 3,
    shortDesc = u"""Training reaction from kinetics library: 2003_Miller_Propargyl_Recomb_High_P""",
    longDesc = 
u"""
Taken from entry: X <=> XI
""",
)

