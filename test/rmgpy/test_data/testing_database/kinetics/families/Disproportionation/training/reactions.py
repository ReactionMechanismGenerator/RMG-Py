#!/usr/bin/env python
# encoding: utf-8

name = "Disproportionation/training"
shortDesc = u"Reaction kinetics used to generate rate rules"
longDesc = u"""
Put kinetic parameters for specific reactions in this file to use as a
training set for generating rate rules to populate this kinetics family.
"""
entry(
    index = 0,
    label = "C2H + CH3O <=> C2H2 + CH2O",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (3.61e+13, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kcal/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 9,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: C2H + CH2OH --> C2H2 + CH2O

pg. 504: Discussion on evaluated data

Entry 39,21 (a): CH2OH + C2H --> C2H2 + CH2O

Author suggest a disproportionation rate coefficient of 6.0x10^-11 cm3/molecule/s, due

to very exothermic rxn.  No data available at the time.
MRH 30-Aug-2009
""",
)

entry(
    index = 1,
    label = "C2H3 + O2 = C2H2_1 + HO2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (1.04e+16, 'cm^3/(mol*s)', '*|/', 5),
        n = -1.26,
        Ea = (3.31, 'kcal/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 10,
    shortDesc = u"""S.S. Merchant estimate""",
    longDesc = 
u"""
This rate rule is a estimate taken from NIST, ref: Aromatic and Polycyclic Aromatic
Hydrocarbon Formation in a Laminar Premixed n-butane Flame
Derived from fitting to a complex mechanism for C2H3 + O2 = C2H2 + HO2
""",
)

entry(
    index = 2,
    label = "C2H5 + O2 <=> HO2 + C2H4",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (4.338e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (66.9022, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (700, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""[AJ] Miyoshi 2011 (Table 4, Node 'sp') dx.doi.org/10.1021/jp112152n""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: i-C3H7 + O2 = HO2 + C3H6

pg. 931-932: Discussion on evaluated data

Entry 42,3 (a): Author appears to be skeptical of the only experimentally reported

value.  Author notes that more recent work on C2H5+O2 suggested that the
addition and disproportionation rxns may be coupled through a common intermediate.
For the time being, the author decided to recommend the only experimentally
reported rate coefficient, only for temperatures above 700K, as they note the
addition rxn should be the predominant rxn at lower temperatures.
MRH 30-Aug-2009

Divide the rate constant by 12 to account for symmetry of 2 (O2) and 6 (i-C3H7, carbons #1 and 3).  The final result is 1.05833e+10 cm3/mol/s.
JDM 31-Mar-2010

Converted to training reaction from rate rule: O2b;Cmethyl_Csrad
""",
)

entry(
    index = 3,
    label = "CH2 + C2H5 <=> CH3 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (9.03e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2(triplet) + i-C3H7 --> C3H6 + CH3

pg. 944: Discussion on evaluated data

Entry 42,26: No data available at the time.  Author suggests this is a minor channel,

stating the main process should be combination, leading to chemically activated
i-butyl radical.  Rate coefficient is estimate.
MRH 30-Aug-2009

Converted to training reaction from rate rule: CH2_triplet;Cmethyl_Csrad
""",
)

entry(
    index = 4,
    label = "H + C2H5 <=> H2 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.083e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + i-C3H7 --> C3H6 + H2

pg. 932: Discussion on evaluated data

Entry 42,4 (a): No data available at the time.  Author recommends a rate coefficient

expression equal to double the rate expression of H+C2H5=H2+C2H4.
MRH 30-Aug-2009

Converted to training reaction from rate rule: H_rad;Cmethyl_Csrad
""",
)

entry(
    index = 5,
    label = "CH3_r1 + C2H5 <=> CH4 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (6.57e+14, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.68,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH3 + i-C3H7 --> C3H6 + CH4

pg. 937: Discussion on evaluated data

Entry 42,16 (b): Author notes that measurements performed by Arthur and Anastasi on

the rate coefficient of total CH3+i-C3H7 decomposition matches very well with
the coefficient derived from the recommended rate of CH3+CH3 decomposition, the 
recommended rate of i-C3H7+i-C3H7 decomposition, and the geometric rule.  The author
recommends a high-pressure rate expression of 4.7x10^-11*(300/T)^0.68 cm3/molecule/s
for the addition rexn (based on the geometric mean rule???) and recommends the 
branching ratio of 0.16 reported by Gibian and Corley (1973).
NOTE: Previous data entry appeared to compute A and n as such:

A = 0.16 * 4.7x10^-11 * (1/300)^0.68
n = 0.68
However, MRH would compute A and n:

A = 0.16 * 4.7x10^-11 * (300)^0.68
n = -0.68
These are the values that now reside in the database.  The online NIST database

(kinetics.nist.gov) agree with what I have calculated.
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_methyl;Cmethyl_Csrad
""",
)

entry(
    index = 6,
    label = "C2H5 + C2H5-2 <=> C2H6 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (6.9e+13, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.35,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H5 + i-C3H7 --> C3H6 + C2H6

pg. 937-938: Discussion on evaluated data

Entry 42,17 (c): No data available at the time.  The rate coefficient expression for

the combination rxn is computed using the geometric mean rule and is reported as
2.6x10^-11 * (300/T)^0.35 cm3/molecule/s.  The recommended branching ratio for 
disproportionation to addition is that reported by Gibian and Corley (1973).
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cs;Cmethyl_Csrad
""",
)

entry(
    index = 7,
    label = "C3H5 + C2H5 <=> C3H6 + C2H4",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (1.374e+14, 'cm^3/(mol*s)', '*|/', 3),
        n = -0.35,
        Ea = (-0.54392, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [93] Literature review.""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + iC3H7 --> C3H6 + C3H6

pg.268: Discussion on evaluated data

Entry 47,42(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for "alkyl radicals" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.2.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC3H7+iC3H7-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cd;Cmethyl_Csrad
""",
)

entry(
    index = 8,
    label = "CH3O-2 + C2H5 <=> CH4O + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (8.67e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + i-C3H7 --> C3H6 + CH3OH

pg. 945: Discussion on evaluated data

Entry 42,39 (c): No data available at the time.  Author recommends a rate coefficient

of 4.8x10^-12 based on the rate expression of i-C3H7+C2H5=C2H6+C3H6
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/O;Cmethyl_Csrad
""",
)

entry(
    index = 9,
    label = "C2H5 + C3H7 <=> C3H8 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (6.33e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.7,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  i-C3H7 + i-C3H7 --> C3H6 + C3H8

pg. 946-947: Discussion on evaluated data

Entry 42,42 (b): No high-Temperature data available.  Author has fit rate coefficient

expression for addition rxn to 4 sets of experimental data.  Recommended branching
ratio agrees well with most of the experimental data.
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H/NonDeC;Cmethyl_Csrad
""",
)

entry(
    index = 10,
    label = "C2H5 + C4H9 <=> C4H10 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (8.58e+15, 'cm^3/(mol*s)', '*|/', 1.7),
        n = -1.1,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + i-C3H7 --> C3H6 + i-C4H10

pg. 46: Discussion on evaluated data

Entry 44,42 (a): The author computes the combination rate expression using the geometric

mean rule (of the rxns t-C4H9+t-C4H9-->adduct and i-C3H7+i-C3H7-->adduct).  The
disproportionation rate coefficient expression was then computed using the
reported branching ratio.
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/Cs3;Cmethyl_Csrad
""",
)

entry(
    index = 11,
    label = "C2H3-2 + C2H5 <=> C2H4-2 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (4.56e+14, 'cm^3/(mol*s)', '*|/', 1.5),
        n = -0.7,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H3 + i-C3H7 --> C3H6 + C2H4

pg. 939-940: Discussion on evaluated data

Entry 42,19 (a): No data available at the time.  Author recommends the rate coefficient

expression of C2H5+i-C3H7 for the rate expression for C2H3+i-C3H7.  Author also
recommends the branching ratio of disproportionation to addition of the 
C2H5+i-C3H7 system for the C2H3+i-C3H7 system.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Cd_pri_rad;Cmethyl_Csrad
""",
)

entry(
    index = 12,
    label = "C2H + C2H5 <=> C2H2 + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.083e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H + i-C3H7 --> C3H6 + C2H2

pg. 941-942: Discussion on evaluated data

Entry 42,21 (a): No data available at the time.  Author recommends a rate coefficient

of 6x10^-12 cm3/molecule/s, a "typical" disproportionation rate.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Ct_rad/Ct;Cmethyl_Csrad
""",
)

entry(
    index = 13,
    label = "HO + C2H5 <=> H2O + C2H4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (7.23e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + i-C3H7 --> C3H6 + H2O

pg. 934: Discussion on evaluated data

Entry 42,6: No data available at the time.  Author notes that both a H-atom abstraction

rxn and an addition + hot adduct decomposition rxn will result in the same products.
MRH 30-Aug-2009

Converted to training reaction from rate rule: O_pri_rad;Cmethyl_Csrad
""",
)

entry(
    index = 14,
    label = "H + CH3O-3 <=> H2 + CH2O-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (5.43e+13, 'cm^3/(mol*s)', '*|/', 3.16),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1000, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Baulch et al [95] literature review.""",
    longDesc = 
u"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; Troe, J.; Walker, R.W.; Warnatz, J.; Journal of Physical and Chemical Reference Data (1992), 21(3), 411-734.
pg.523: Discussion of evaluated data

H+CH3O --> H2+CH2O: Authors state that no new data have been reported for this reaction.

MRH assumes the recommended value comes from a previous review article published
by authors.  In any case, recommended data fits the reported data well.
MRH 31-Aug-2009

Converted to training reaction from rate rule: H_rad;Cmethyl_Orad
""",
)

entry(
    index = 15,
    label = "C3H7-2 + O2 <=> HO2 + C3H6-2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (1.833e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (62.1324, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (500, 'K'),
        Tmax = (900, 'K'),
    ),
    rank = 10,
    shortDesc = u"""[AJ] Miyoshi 2011 (Table 4, Node 'ss') dx.doi.org/10.1021/jp112152n""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: n-C3H7 + O2 = HO2 + C3H6

pg. 914-915: Discussion on evaluated data

Entry 41,3 (a): The author suggests a rate coefficient based on those reported in the

literature.  The author notes that the data reported in the literature suggests
the formation of C3H6 is controlled by the addition rxn.  The author further
notes that it is surprising that p-dependence effects are not observed for
C3H6 formation.
MRH 30-Aug-2009

Divide the rate constant by 4 to account for symmetry of 2 (O2) and 2 (n-C3H7, carbon #2).  The final result is 2.25825e+10 cm3/mol/s.
JDM 31-Mar-2010

Converted to training reaction from rate rule: O2b;C/H2/Nd_Csrad
""",
)

entry(
    index = 16,
    label = "CH2 + C3H7-2 <=> CH3 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.62e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2_triplet + n-C3H7 --> C3H6 + CH3

pg. 925: Discussion on evaluated data

Entry 41,26: No data available at the time.  Author estimates the rate coefficient

expression of the addition rxn.  The author then recommends that the disproportionation
rate coefficient not exceed 10% of the combination rate.  Thus, the rate coefficient
is an upper limit.
MRH 30-Aug-2009

Converted to training reaction from rate rule: CH2_triplet;C/H2/Nd_Csrad
""",
)

entry(
    index = 17,
    label = "C2H3S + C3H7-2 <=> C2H4S + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.438, 'cm^3/(mol*s)'),
        n = 3.13,
        Ea = (-15.2716, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 6,
    shortDesc = u"""CAC calc CBS-QB3, 1dhr""",
    longDesc = 
u"""
Converted to training reaction from rate rule: S_rad/OneDe;C/H2/Nd_Csrad
""",
)

entry(
    index = 18,
    label = "H + C3H7-2 <=> H2 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.62e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + n-C3H7 --> C3H6 + H2

pg. 915-916: Discussion on evaluated data

Entry 41,4 (a): No data available at the time.  Author recommends the rate coefficient

of the H+C2H5=C2H4+H2 rxn for the H+n-C3H7=C3H6+H2 rxn.
MRH 30-Aug-2009

Converted to training reaction from rate rule: H_rad;C/H2/Nd_Csrad
""",
)

entry(
    index = 19,
    label = "C4H7 + C2H5S <=> C4H8 + C2H4S-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.526e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (-2.3012, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Rough estimate based on 1/10 of #3026 in R_Recombination""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH3 + n-C3H7 --> C3H6 + CH4

pg. 920: Discussion on evaluated data

Entry 41,16 (b): No direct measurements for either the addition or disproportionation

rxns.  Author recommends a rate coefficient expression for the addition rxn, based
on the geometric mean rule of the rxns CH3+CH3=>adduct and n-C3H7+n-C3H7=>adduct.
Furthermore, author recommends a branching ratio for disproportionation to
addition of 0.06 (which appears to MRH to be consistent with the experimentally
measured branching ratios)
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H/OneDeC;C/H2/Nd_Srad
""",
)

entry(
    index = 20,
    label = "CH3_r1 + C3H7-2 <=> CH4 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.3e+13, 'cm^3/(mol*s)', '*|/', 1.7),
        n = -0.32,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH3 + n-C3H7 --> C3H6 + CH4

pg. 920: Discussion on evaluated data

Entry 41,16 (b): No direct measurements for either the addition or disproportionation

rxns.  Author recommends a rate coefficient expression for the addition rxn, based
on the geometric mean rule of the rxns CH3+CH3=>adduct and n-C3H7+n-C3H7=>adduct.
Furthermore, author recommends a branching ratio for disproportionation to
addition of 0.06 (which appears to MRH to be consistent with the experimentally
measured branching ratios)
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_methyl;C/H2/Nd_Csrad
""",
)

entry(
    index = 21,
    label = "C3H7-2 + C2H5-2 <=> C2H6 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.9e+12, 'cm^3/(mol*s)', '*|/', 1.4),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H5 + n-C3H7 --> C3H6 + C2H6

pg. 937-938: Discussion on evaluated data

Entry 42,17 (b): No direct measurements for either the addition or disproportionation

rxns.  Author recommends a rate coefficient expression for the addition rxn, based
on the geometric mean rule of the rxns C2H5+C2H5=>adduct and n-C3H7+n-C3H7=>adduct.
Furthermore, author recommends a branching ratio for disproportionation to
addition of 0.073 (which is an average of the 2 experimentally determined
branching ratios)
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cs;C/H2/Nd_Csrad
""",
)

entry(
    index = 22,
    label = "C5H7 + C2H5S <=> C5H8 + C2H4S-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.88e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (1.50624, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Rough estimate based on 1/10 of #3027 in R_Recombination""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H5 + n-C3H7 --> C3H6 + C2H6

pg. 937-938: Discussion on evaluated data

Entry 42,17 (b): No direct measurements for either the addition or disproportionation

rxns.  Author recommends a rate coefficient expression for the addition rxn, based
on the geometric mean rule of the rxns C2H5+C2H5=>adduct and n-C3H7+n-C3H7=>adduct.
Furthermore, author recommends a branching ratio for disproportionation to
addition of 0.073 (which is an average of the 2 experimentally determined
branching ratios)
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H/TwoDe;C/H2/Nd_Srad
""",
)

entry(
    index = 23,
    label = "C3H5 + C3H7-2 <=> C3H6 + C3H6-2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (5.8e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (-0.54392, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [93] Literature review.""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + nC3H7 --> C3H6 + C3H6

pg.268: Discussion on evaluated data

Entry 47,41(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for "alkyl radicals" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.07.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and nC3H7+nC3H7-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cd;C/H2/Nd_Csrad
""",
)

entry(
    index = 24,
    label = "C2H3S + C4H9-2 <=> C2H4S + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (7.63e+11, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (-2.3012, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""VERY Rough estimate based on 1/10 of #3026 in R_Recombination""",
    longDesc = 
u"""
Converted to training reaction from rate rule: S_rad/OneDe;C/H/NdNd_Csrad
""",
)

entry(
    index = 25,
    label = "CH3O-2 + C3H7-2 <=> CH4O + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (9.64e+11, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + n-C3H7 --> C3H6 + CH3OH

pg. 926: Discussion on evaluated data

Entry 41,39 (c): No data available at the time.  Author estimates the rate coefficient

for the addition rxn to be similar to the rate for n-C3H7+n-C3H7=>adduct.  Author
also estimates the branching ratio of disproportionation to addition as 0.051
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/O;C/H2/Nd_Csrad
""",
)

entry(
    index = 26,
    label = "C4H9-2 + CH3S <=> CH4S + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (7.63e+11, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (-2.3012, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Rough estimate based on 1/10 of #3026 in R_Recombination""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + n-C3H7 --> C3H6 + CH3OH

pg. 926: Discussion on evaluated data

Entry 41,39 (c): No data available at the time.  Author estimates the rate coefficient

for the addition rxn to be similar to the rate for n-C3H7+n-C3H7=>adduct.  Author
also estimates the branching ratio of disproportionation to addition as 0.051
MRH 30-Aug-2009

Converted to training reaction from rate rule: S_rad/NonDeC;C/H/NdNd_Csrad
""",
)

entry(
    index = 27,
    label = "C3H7-2 + C3H7 <=> C3H8 + C3H6-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (5.13e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.35,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  i-C3H7 + n-C3H7 --> C3H6 + C3H8

pg. 945-946: Discussion on evaluated data

Entry 42,41 (b): No data available at the time.  Author estimates the rate coefficient

expression of the addition rxn using the rate for i-C3H7+i-C3H7=>adduct, the rate
for n-C3H7+n-C3H7=>adduct, and the geometric mean rule.  The author recommends
the branching ratio of disproportionation to addition reported by Gibian and
Corley (1973).
MRH 30-Aug-2009

Degeneracy not recalculated

Converted to training reaction from rate rule: C_rad/H/NonDeC;C/H2/Nd_Csrad
""",
)

entry(
    index = 28,
    label = "C3H7-2 + C4H9 <=> C4H10 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.32e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.75,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + n-C3H7 --> C3H6 + i-C4H10

pg. 45: Discussion on evaluated data

Entry 44,41 (a): No data available at the time.  Author estimates the rate expression

for the combination rxn using the geometric mean rule (of the rxns t-C4H9+t-C4H9-->adduct
and n-C3H7+n-C3H7-->adduct).  The author then estimates the disproportionation
rate expression using the branching ratio; the branching ratio is from "analogous
processes".
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/Cs3;C/H2/Nd_Csrad
""",
)

entry(
    index = 29,
    label = "C2H3-2 + C3H7-2 <=> C2H4-2 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.42e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H3 + n-C3H7 --> C3H6 + C2H4

pg. 922: Discussion on evaluated data

Entry 41,19 (a): No data available at the time.  Author estimates the rate coefficient

based on the rxn C2H5+n-C3H7=C3H6=C2H6.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Cd_pri_rad;C/H2/Nd_Csrad
""",
)

entry(
    index = 30,
    label = "HS2 + C5H9 <=> H2S2 + C5H8-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.288e+09, 'cm^3/(mol*s)'),
        n = 1.19,
        Ea = (2.13384, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Very rough based on R_Recomb #491""",
    longDesc = 
u"""
Converted to training reaction from rate rule: S_rad/NonDeS;C/H2/Nd_Csrad/H/Cd
""",
)

entry(
    index = 31,
    label = "C2H + C3H7-2 <=> C2H2 + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.206e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H + n-C3H7 --> C3H6 + C2H2

pg. 923: Discussion on evaluated data

Entry 41,21 (a): No data available at the time.  Author notes that the rxn is more exothermic

than the rxn CH3+n-C3H7=C3H6+CH4 and suggests a rate coefficient 3x larger,
namely 1.0x10^-11 cm3/molecule/s.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Ct_rad/Ct;C/H2/Nd_Csrad
""",
)

entry(
    index = 32,
    label = "HS2 + CH3S-2 <=> H2S2 + CH2S",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (6.44e+08, 'cm^3/(mol*s)'),
        n = 1.19,
        Ea = (2.13384, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Very rough based on R_Recomb #491""",
    longDesc = 
u"""
Converted to training reaction from rate rule: S_rad/NonDeS;S_Csrad
""",
)

entry(
    index = 33,
    label = "HO + C3H7-2 <=> H2O + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.82e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + n-C3H7 --> C3H6 + H2O

pg. 917: Discussion on evaluated data

Entry 41,6 (a): No data available at the time.  Author estimates rate coefficient based

on the rate coefficient for OH+C2H5=C2H4+H2O, namely 4.0x10^-11 cm3/molecule/s.
MRH 30-Aug-2009

Converted to training reaction from rate rule: O_pri_rad;C/H2/Nd_Csrad
""",
)

entry(
    index = 34,
    label = "C4H9-2 + O2 <=> HO2 + C4H8-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.4088e+10, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (600, 'K'),
        Tmax = (1000, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: O2 + iC4H9 --> iC4H8 + HO2

pg. 52-53: Discussion on evaluated data

Entry 45,3 (a): The author recommends a rate coefficient based on the experiments performed

by Baker et al. (yielding a disproportionation-to-decomposition ratio) and the
current (Tsang) study's recommended iC4H9 unimolecular decomposition rate.
MRH 31-Aug-2009

Divide the rate constant by 2 to account for symmetry of 2 (O2) and 1 (i-C4H9, carbon #2).  The final result is 1.2044e+10 cm3/mol/s.
JDM 31-Mar-2010

Converted to training reaction from rate rule: O2b;C/H/NdNd_Csrad
""",
)

entry(
    index = 35,
    label = "C2H + C4H9-2 <=> C2H2 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H + i-C4H9 --> i-C4H8 + C2H2

pg. 61: Discussion on evaluated data

Entry 45,21: No data available at the time.  The author estimates the rate of 

disproportionation to be 1x10^-11 cm3/molecule/s.
*** NOTE: RMG_database previously had CH2_triplet as Y_rad_birad node, not Ct_rad ***

MRH 30-Aug-2009

Converted to training reaction from rate rule: Ct_rad/Ct;C/H/NdNd_Csrad
""",
)

entry(
    index = 36,
    label = "H + C4H9-2 <=> H2 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (9.04e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: H + i-C4H9 --> i-C4H8 + H2

pg. 53: Discussion on evaluated data

Entry 45,4 (c): No data available at the time.  The author estimates the disproportionation

rate coefficent as half the rate of H+n-C3H7=C3H6+H2 (due to the presence of 2
H-atoms on the alpha-carbon in n-C3H7 and only 1 on the alpha-carbon of i-C4H9).
The author also states that the branching ratio is pressure-dependent and supplies
fall-off tables and collisional efficiencies.
MRH 30-Aug-2009

Converted to training reaction from rate rule: H_rad;C/H/NdNd_Csrad
""",
)

entry(
    index = 37,
    label = "CH3_r1 + C4H9-2 <=> CH4 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (6.02e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.32,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: CH3 + i-C4H9 --> i-C4H8 + CH4

pg. 58: Discussion on evaluated data

Entry 45,16 (b): No data available at the time.  The author estimates the disproportionation

rate coefficient as half the rate of CH3+n-C3H7=C3H6+H2 (due to half as many H-atoms
on the alpha-carbon).
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_methyl;C/H/NdNd_Csrad
""",
)

entry(
    index = 38,
    label = "C4H9-2 + C2H5-2 <=> C2H6 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (8.43e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H5 + i-C4H9 --> i-C4H8 + C2H6

pg. 59: Discussion on evaluated data

Entry 45,17 (a): No direct measurements of either the addition or disproportionation rxns.

The combination rate coefficient was computed using the geometric mean rule (of the
rxns C2H5+C2H5-->adduct and i-C4H9+i-C4H9-->adduct).  The disproportionation rate
coefficient was computed using the disproportionation-to-combination ratio reported
by Gibian and Corley (1973).
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cs;C/H/NdNd_Csrad
""",
)

entry(
    index = 39,
    label = "CH3O-2 + C4H9-2 <=> CH4O + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.41e+11, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: CH2OH + i-C4H9 --> i-C4H8 + CH3OH

pg. 64: Discussion on evaluated data

Entry 45,39 (c): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of CH2OH+n-C3H7=C3H6+CH3OH (due to half as many H-atoms
on the alpha-carbon).
*** NOTE: Although author states the the rate coefficient of CH2OH+i-C4H9=i-C4H8+CH3OH

is half that of CH2OH+n-C3H7=C3H6+CH3OH, MRH finds them to be equal, both in the electronic
references and the online NIST database (kinetics.nist.gov).  I am therefore
cutting the A in the RMG_database in two. ***
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/O;C/H/NdNd_Csrad
""",
)

entry(
    index = 40,
    label = "C3H5 + C4H9-2 <=> C3H6 + C4H8-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.566e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (-0.54392, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [93] Literature review.""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + iC4H9 --> iC4H8 + C3H6

pg.270: Discussion on evaluated data

Entry 47,45(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC4H9+iC4H9-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cd;C/H/NdNd_Csrad
""",
)

entry(
    index = 41,
    label = "C4H9-2 + C3H7 <=> C3H8 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.56e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.35,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: i-C3H7 + i-C4H9 --> i-C4H8 + C3H8

pg. 65: Discussion on evaluated data

Entry 45,42 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of i-C3H7+n-C3H7=C3H6+C3H8 (due to half as many H-atoms
on the alpha-carbon).
*** NOTE: MRH computes half the rate of i-C3H7+n-C3H7=C3H6+C3H8 as 0.52x10^-11 * (300/T)^0.35,

not 0.58x10^-11 * (300/T)^0.35.  However, there may be a reason for the relatively
small discrepancy between the author's stated and implemented calculation. ***
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H/NonDeC;C/H/NdNd_Csrad
""",
)

entry(
    index = 42,
    label = "C4H9-2 + C4H9 <=> C4H10 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.08e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.75,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + i-C4H9 --> i-C4H8 + i-C4H10

pg. 66: Discussion on evaluated data

Entry 45,44 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of t-C4H9+n-C3H7=C3H6+i-C4H10 (due to half as many H-atoms
on the alpha-carbon).
*** NOTE: Although author states the the rate coefficient of t-C4H9+i-C4H9=i-C4H8+i-C4H10

is half that of t-C4H9+n-C3H7=C3H6+i-C4H10, MRH finds them to be equal, both in the electronic
references and the online NIST database (kinetics.nist.gov).  I am therefore
cutting the A in the RMG_database in two. ***
MRH 30-Aug-2009

Degeneracy not recalculated

Converted to training reaction from rate rule: C_rad/Cs3;C/H/NdNd_Csrad
""",
)

entry(
    index = 43,
    label = "C2H3-2 + C4H9-2 <=> C2H4-2 + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (8.43e+11, 'cm^3/(mol*s)', '*|/', 4),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H3 + i-C4H9 --> i-C4H8 + C2H4

pg. 60: Discussion on evaluated data

Entry 45,19 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient based on the rate of C2H5+i-C4H9=i-C4H8+C2H6.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Cd_pri_rad;C/H/NdNd_Csrad
""",
)

entry(
    index = 44,
    label = "HO + C4H9-2 <=> H2O + C4H8-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.21e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: OH + i-C4H9 --> i-C4H8 + H2O

pg. 55: Discussion on evaluated data

Entry 45,6 (a): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of OH+n-C3H7=C3H6+H2O (due to half as many H-atoms
on the alpha-carbon).
MRH 30-Aug-2009

Converted to training reaction from rate rule: O_pri_rad;C/H/NdNd_Csrad
""",
)

entry(
    index = 45,
    label = "C3H5-2 + O2 <=> HO2 + C3H4",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.2044e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (56.6932, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [93] Literature review.""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: O2 + C3H5 --> H2C=C=CH2 + HO2

pg.251: Discussion on evaluated data

*** UPPER LIMIT ***

Entry 47,3(b): The author states that there is uncertainty whether this rxn is appreciable

at high temperatures.  There were conflicting results published regarding the
significance above 461K (Morgan et al. and Slagle and Gutman).  The author thus
decides to place an upper limit on the rate coefficient of 2x10^-12 * exp(-6820/T)
cm3/molecule/s.  The author further notes that this upper limit assumes no
contribution from a complex rearrangement of the adduct.  Finally, the author
notes that this rxn should not be significant in combustion situations.
MRH 31-Aug-2009

Divide the rate constant by 2 to account for symmetry of 2 (O2) and 1 (allyl, carbon #2). The final result is 6.022e+11 cm3/mol/s, Ea = 13.55 kcal/mol.
JDM 31-Mar-2010

Converted to training reaction from rate rule: O2b;Cdpri_Csrad
""",
)

entry(
    index = 46,
    label = "CH3_r1 + C3H5-2 <=> CH4 + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (3.01e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: CH3 + C3H5 --> H2C=C=CH2 + CH4

pg.257: Discussion on evaluated data

Entry 47,16(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.); this leads to disproportionation-
to-addition ratio of 0.03.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and CH3+CH3-->adduct.
NOTE: The Ea reported in the discussion is Ea/R=-132 Kelvin.  However, in the table near

the beginning of the review article (summarizing all reported data) and in the NIST
online database (kinetics.nist.gov), the reported Ea/R=-66 Kelvin.  MRH took the
geometric mean of the allyl combination rxn (1.70x10^-11 * exp(132/T)) and methyl
combination rxn (1.68x10^-9 * T^-0.64) to obtain 1.69x10^-11 * T^-0.32 * exp(66/T).
Multiplying by 0.03 results in the recommended rate coefficient expression.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_methyl;Cdpri_Csrad
""",
)

entry(
    index = 47,
    label = "C3H5-2 + C2H5-2 <=> C2H6 + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (9.64e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C2H5 + C3H5 --> H2C=C=CH2 + C2H6

pg.259: Discussion on evaluated data

Entry 47,17(a): The recommended rate expression is derived from the experimentally-

determined disproportionation-to-addition ratio of 0.047 (James and Troughton)
and the addition rate rule (C2H5+C3H5-->adduct) calculated using the geometric
mean rule of the rxns C2H5+C2H5-->adduct and C3H5+C3H5-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cs;Cdpri_Csrad
""",
)

entry(
    index = 48,
    label = "C3H5 + C3H5-2 <=> C3H6 + C3H4",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.686e+11, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + C3H5 --> H2C=C=CH2 + C3H6

pg.271-272: Discussion on evaluated data

Entry 47,47(b): The recommended rate expression is derived from the experimentally-

determined disproportionation-to-addition ratio of 0.008 (James and Kambanis)
and the addition rate rule (C3H5+C3H5-->adduct) calculated based on the results
of Tulloch et al.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cd;Cdpri_Csrad
""",
)

entry(
    index = 49,
    label = "C3H5-2 + C3H7 <=> C3H8 + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (4.58e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: iC3H7 + C3H5 --> H2C=C=CH2 + C3H8

pg.268: Discussion on evaluated data

Entry 47,42(b): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for "alkyl radicals" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC3H7+iC3H7-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H/NonDeC;Cdpri_Csrad
""",
)

entry(
    index = 50,
    label = "C3H5-2 + C4H9 <=> C4H10 + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.89e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: tC4H9 + C3H5 --> H2C=C=CH2 + iC4H10

pg.269: Discussion on evaluated data

Entry 47,44(b): No data available at the time.  Recommended rate coefficient expression

based on "allyl and alkyl radicals behaving in similar fashion" (possibly referencing
Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and tC4H9+tC4H9-->adduct.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/Cs3;Cdpri_Csrad
""",
)

entry(
    index = 51,
    label = "C2H3-2 + C3H5-2 <=> C2H4-2 + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.41e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C2H3 + C3H5 --> H2C=C=CH2 + C2H4

pg.261-262: Discussion on evaluated data

Entry 47,19(d): No data available at the time.  Author recommends a rate coefficient

of 4x10^-12 cm3/molecule/s for the disproportionation rxn.
MRH 31-Aug-2009

Converted to training reaction from rate rule: Cd_pri_rad;Cdpri_Csrad
""",
)

entry(
    index = 52,
    label = "HO + C3H5-2 <=> H2O + C3H4",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (25.104, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: OH + C3H5 --> H2C=C=CH2 + H2O

pg.253: Discussion on evaluated data

Entry 47,6(a): No data available at the time.  Author recommends a rate coefficient

of 1x10^-11 cm3/molecule/s, based on "comparable rxns".
MRH 31-Aug-2009

Converted to training reaction from rate rule: O_pri_rad;Cdpri_Csrad
""",
)

entry(
    index = 53,
    label = "CH3O + O2 <=> HO2 + CH2O",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.14418e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (298, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Atkinson et al [98] literature review.""",
    longDesc = 
u"""
[98] Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F., Jr.; Kerr, J.A.; Rossi, M.J.; Troe, J. "Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry,", 2001.
Literature review: CH3CHOH + O2 --> CH3CHO + HO2

Recommended value is k298.  This reference just gives a table of results,

with no discussion on how the preferred numbers were arrived at.
MRH 31-Aug-2009

Divide the rate constant by 2 to account for symmetry of 2 (O2) and 1 (CH3CHOH, oxygen atom). The final result is 5.7209e+12 cm3/mol/s.
JDM 31-Mar-2010

Converted to training reaction from rate rule: O2b;O_Csrad
""",
)

entry(
    index = 54,
    label = "CH3O + O <=> HO-2 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (9.04e+13, 'cm^3/(mol*s)', '+|-', 3.01e+13),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (298, 'K'),
    ),
    rank = 6,
    shortDesc = u"""Grotheer et al [189].""",
    longDesc = 
u"""
[189] Grotheer, H.; Riekert, G.; Walter, D.; Just, T. Symp. Int. Combust. Proc. 1989, 22, 963.
Absolute value measured directly. Excitation: discharge, analysis: mass spectroscopy. Original uncertainty 3.0E+13

O + CH2OC --> OH + CH2O, O + CH3CHOH --> OH + CH3CHO

O+CH2OH --> OH+CH2O && O+CH3CHOH --> OH+CH3CHO

pg.963: Measured rate coefficients mentioned in abstract as k_2M and k_2E.

pg.965-967: Discussion on measured rate coefficients.

MRH 1-Sept-2009

Converted to training reaction from rate rule: O_atom_triplet;O_Csrad
""",
)

entry(
    index = 55,
    label = "CH2 + CH3O <=> CH3 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.21e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2 + CH2OH --> CH3 + CH2O

pg. 505: Discussion on evaluated data

Entry 39,26 (b): CH2OH + CH2(triplet) --> CH3 + CH2O

Author estimates the rate of disproportionation as 2.0x10^-12 cm3/molecule/s.  No data at the time.

MRH 30-Aug-2009

Converted to training reaction from rate rule: CH2_triplet;O_Csrad
""",
)

entry(
    index = 56,
    label = "H + CH3O <=> H2 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2e+13, 'cm^3/(mol*s)', '+|-', 1e+13),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (295, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Edelbuttel-Einhaus et al [190].""",
    longDesc = 
u"""
[190] Edelbuttel-Einhaus, J.; Hoyermann, K.; Rohde, G.; Seeba, J. Symp. Int. Combust. Proc. 1992, 22, 661.
Data derived from fitting to a complex mechanism. Excitation: discharge, analysis: mass spectroscopy. Original uncertainty 1.0E+13

H + CH3CHOH --> H2 + CH3CHO

H+CH3CHOH --> H2+CH3CHO

pg.661: Measured rate coefficient mentioned in abstract as k6.

pg.665-666: Discussion on measured rate coefficient.  The reported rate coefficient is

for H+CH3CHOH --> products, making this an UPPER LIMIT.  The rate coefficient
was calculated based on the rate coefficient of the rxn C2H5+H --> CH3+CH3; the
value the authors used was 3.6x10^13 cm3/mol/s.
MRH 1-Sept-2009

Converted to training reaction from rate rule: H_rad;O_Csrad
""",
)

entry(
    index = 57,
    label = "CH3O + CH3_r1 <=> CH4 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (8.49e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (298, 'K'),
    ),
    rank = 6,
    shortDesc = u"""Pagsberg et al [191].""",
    longDesc = 
u"""
[191] Pagsberg, P.; Munk, J.; Sillesen, A.; Anastasi, C. Chem. Phys. Lett. 1988, 146, 375.
Absolute value measured directly. Excitatio: electron beam, analysis: Vis-UV absorption.

CH2OH + CH3 --> CH2O + CH4

pg.378 Table 2: Formation and decay rates of CH2OH, CH3, and OH observed by pulse radiolysis of

gas mixtures of varying composition.  Chemical composition of systems A-E as in Table 1.
The authors note below Table 2 that the reported rate coefficient for CH3+CH2OH is an

"adjustment of model to reproduce the observed decay rates of CH3 and CH2OH".
MRH is skeptical of data, as this specific rxn is not directly referenced in the article,

nor do the authors address whether other channels besides -->CH4+CH2O exist / are significant.
The value of A in the database is consistent with that reported in Table 2.
MRH 1-Sept-2009

Converted to training reaction from rate rule: C_methyl;O_Csrad
""",
)

entry(
    index = 58,
    label = "CH3O + C2H5-2 <=> C2H6 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.41e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: C2H5 + CH2OH --> C2H6 + CH2O

pg. 502: Discussion on evaluated data

Entry 39,17 (b): C2H5 + CH2OH --> C2H6 + CH2O

Author estimates the disproportionation rate coefficient as 4x10^-12 cm3/molecule/s.

No data at the time.
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cs;O_Csrad
""",
)

entry(
    index = 59,
    label = "CH3O + C3H5 <=> C3H6 + CH2O",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.62e+13, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [93] Literature review.""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + CH2OH --> CH2O + C3H6

pg.267: Discussion on evaluated data

Entry 47,39: No data available at the time.  Author notes that combination of these two

reactants will form 3-butene-1-ol which should decompose under combustion conditions
to form C3H6 + CH2O (same products).  The author therefore recommends a rate
coefficient of 3x10^-11 cm3/molecule/s.
MRH 31-Aug-2009

Converted to training reaction from rate rule: C_rad/H2/Cd;O_Csrad
""",
)

entry(
    index = 60,
    label = "CH3O-2 + CH3O <=> CH4O + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (4.82e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2OH + CH2OH --> CH3OH + CH2O

pg. 506: Discussion on evaluated data

Entry 39,39 (b): CH2OH + CH2OH --> CH3OH + CH2O

Meier, et al. (1985) measured the rate of addition + disproportionation.  Tsang estimates

a disproportionation to combination ratio of 0.5
NOTE: Rate coefficient given in table at beginning of reference (summarizing all data

presented) gives k_a+b = 2.4x10^-11, leading to k_b = 8x10^-12.  NIST's online
database (kinetics.nist.gov) reports this number as well.  However, the discussion
on pg. 506 suggests k_a+b = 1.5x10^-11, leading to k_b = 5x10^-12.
MRH 30-Aug-2009

*** NEED TO INVESTIGATE ***

Converted to training reaction from rate rule: C_rad/H2/O;O_Csrad
""",
)

entry(
    index = 61,
    label = "CH3O + C3H7 <=> C3H8 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.35e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: CH2OH + i-C3H7 = C3H8 + CH2O

pg. 945: Discussion on evaluated data

Entry 42,39 (b): No data available at the time.  Author suggests rate coefficient based

on rxn C2H5+i-C3H7=C3H8+C2H4, namely 3.9x10^-12 cm3/molecule/s
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/H/NonDeC;O_Csrad
""",
)

entry(
    index = 62,
    label = "CH3O + C4H9 <=> C4H10 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (3.47e+14, 'cm^3/(mol*s)', '*|/', 3),
        n = -0.75,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + CH2OH = CH2O + i-C4H10

pg. 44: Discussion on evaluated data

Entry 44,39 (a): No data available at the time.  Author estimates the addition rxn rate

coefficient based on the rate for t-C4H9+C2H5-->adduct.  The author uses a
disproportionation-to-addition ratio of 0.52 to obtain the reported rate coefficient
expression.
*** NOTE: Previous value in RMG was for k_c (the addition rxn).  I have changed it to match

the rate for the disproportionation rxn. ***
MRH 30-Aug-2009

Converted to training reaction from rate rule: C_rad/Cs3;O_Csrad
""",
)

entry(
    index = 63,
    label = "CH3O + C2H3-2 <=> C2H4-2 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (3.01e+13, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2OH + C2H3 --> C2H4 + CH2O

pg. 503: Discussion on evaluated data

Entry 39,19 (a): CH2OH + C2H3 --> C2H4 + CH2O

Author suggests a disproportionation rate coefficient near the collision limit, due

to rxn's exothermicity.  No data available at the time.
MRH 30-Aug-2009

Converted to training reaction from rate rule: Cd_pri_rad;O_Csrad
""",
)

entry(
    index = 64,
    label = "CHO + CH3O <=> CH2O-3 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.81e+14, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: HCO + CH2OH --> CH2O + CH2O

pg. 500: Discussion on evaluated data

Entry 39,15 (b): CH2OH + HCO --> 2 CH2O

Author estimates a disproportionation rate coefficient of 3x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009

Converted to training reaction from rate rule: CO_pri_rad;O_Csrad
""",
)

entry(
    index = 65,
    label = "HO + CH3O <=> H2O + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: OH + CH2OH --> H2O + CH2O

pg. 497: Discussion on evaluated data

Entry 39,6: CH2OH + OH --> H2O + CH2O

Author estimates a disproportionation rate coefficient of 4x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009

Converted to training reaction from rate rule: O_pri_rad;O_Csrad
""",
)

entry(
    index = 66,
    label = "CH3O + CH3O-4 <=> CH4O-2 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH3O + CH2OH --> CH3OH + CH2O

pg. 505: Discussion on evaluated data

Entry 39,24: CH2OH + CH3O --> CH3OH + CH2O

Author estimates a disproportionation rate coefficient of 4x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009

Degeneracy not recalculated

Converted to training reaction from rate rule: O_rad/NonDeC;O_Csrad
""",
)

entry(
    index = 67,
    label = "HO2-2 + CH3O <=> H2O2 + CH2O",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.21e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: HO2 + CH2OH --> CH3OH + H2O2

pg. 498: Discussion on evaluated data

Entry 39,7: CH2OH + HO2 --> H2O2 + CH2O

Author recommends a disproportionation rate coefficient of 2x10^-11 cm3/molecules/s.

No data available at the time.
MRH 30-Aug-2009

Converted to training reaction from rate rule: O_rad/NonDeO;O_Csrad
""",
)

entry(
    index = 68,
    label = "CH3S + CH3S-3 <=> CH4S + CH2S-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (2.937e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (298, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Tycholiz et al [A].""",
    longDesc = 
u"""
Converted to training reaction from rate rule: S_rad/NonDeC;Cmethyl_Srad
""",
)

entry(
    index = 69,
    label = "C2H5S-2 + C3H7-2 <=> C2H6S + C3H6-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (6.74e-06, 'cm^3/(mol*s)'),
        n = 4.35,
        Ea = (4.76976, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 6,
    shortDesc = u"""CAC calc CBS-QB3 1dhr""",
    longDesc = 
u"""
Converted to training reaction from rate rule: C_rad/H/CsS;C/H2/Nd_Csrad
""",
)

entry(
    index = 70,
    label = "C4H7-2 + O2 <=> HO2 + C4H6",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (4.338e+13, 'cm^3/(mol*s)', '*|/', 10),
        n = 0,
        Ea = (92.048, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 11,
    shortDesc = u"""S.S. Merchant estimate""",
    longDesc = 
u"""
SSM estimate based on Miyoshi rate rule for secondary carbon in dx.doi.org/10.1021/jp112152n, 
modified to account for allylic stability (+7 kcal)

Converted to training reaction from rate rule: O2b;Cmethyl_Csrad/H/Cd
""",
)

entry(
    index = 71,
    label = "C4H7-3 + O2 <=> HO2 + C4H6-2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (4e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H2/De_Csrad
""",
)

entry(
    index = 72,
    label = "C3H7-2 + O2 <=> HO2 + C3H6-2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (4e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H2/Nd_Rrad
""",
)

entry(
    index = 73,
    label = "C4H7-3 + O2 <=> HO2 + C4H6-2",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (4e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H2/De_Rrad
""",
)

entry(
    index = 74,
    label = "C4H9-2 + O2 <=> HO2 + C4H8-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H/NdNd_Rrad
""",
)

entry(
    index = 75,
    label = "C5H9-2 + O2 <=> HO2 + C5H8-3",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H/NdDe_Rrad
""",
)

entry(
    index = 76,
    label = "C6H9 + O2 <=> HO2 + C6H8",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2e+10, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 10,
    shortDesc = u"""Estimated value, AG Vandeputte""",
    longDesc = 
u"""
Converted to training reaction from rate rule: O2b;C/H/DeDe_Rrad
""",
)

entry(
    index = 77,
    label = "HO2-3 + H2N <=> H3N + O2-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2 + HO2 = NH3 + O2 (B&D #14d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Degeneracy not recalculated

Converted to training reaction from rate rule: NH2_rad;O_Orad
""",
)

entry(
    index = 78,
    label = "HN2 + O2 <=> HO2 + N2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.4e+12, 'cm^3/(mol*s)'),
        n = -0.34,
        Ea = (0.6276, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + O2 = N2 + HO2 (B&D #28b1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O2b;N3d/H_d_Nrad
""",
)

entry(
    index = 79,
    label = "H + HN2 <=> H2 + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + H = N2 + H2 (B&D #28c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3d/H_d_Nrad
""",
)

entry(
    index = 80,
    label = "HO + HN2 <=> H2O + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + OH = N2 + H2O (B&D #28d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3d/H_d_Nrad
""",
)

entry(
    index = 81,
    label = "HN2 + O <=> HO-2 + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + O = N2 + OH (B&D #28e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3d/H_d_Nrad
""",
)

entry(
    index = 82,
    label = "H2N + HN2 <=> H3N + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + NH2 = N2 + NH3 (B&D #28f) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;N3d/H_d_Nrad
""",
)

entry(
    index = 83,
    label = "HO2-2 + HN2 <=> H2O2 + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (14000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + HO2 = N2 + H2O2 (B&D #28g1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;N3d/H_d_Nrad
""",
)

entry(
    index = 84,
    label = "HN2 + NO <=> HNO + N2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + NO = N2 + HNO (B&D #28h) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: N3d_rad/O;N3d/H_d_Nrad
""",
)

entry(
    index = 85,
    label = "H + H2N2 <=> H2 + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (9.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + H = NNH + H2 (B&D #30c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3s/H2_s_Nbirad
""",
)

entry(
    index = 86,
    label = "H2N2 + O <=> HO-2 + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (6.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + O = NNH + OH (B&D #30d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3s/H2_s_Nbirad
""",
)

entry(
    index = 87,
    label = "HO + H2N2 <=> H2O + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + OH = NNH + H2O (B&D #30e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3s/H2_s_Nbirad
""",
)

entry(
    index = 88,
    label = "H2N2 + CH3_r1 <=> CH4 + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.2e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (0.54392, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + CH3 = NNH + CH4 (B&D #30f3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;N3s/H2_s_Nbirad
""",
)

entry(
    index = 89,
    label = "H2N + H2N2 <=> H3N + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.6e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + NH2 = NNH + NH3 (B&D #30g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;N3s/H2_s_Nbirad
""",
)

entry(
    index = 90,
    label = "HO2-2 + H2N2 <=> H2O2 + HN2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (58000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + HO2 = NNH + H2O2 (B&D #30h2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;N3s/H2_s_Nbirad
""",
)

entry(
    index = 91,
    label = "H + H3N2 <=> H2 + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + H = N2H2 + H2 (B&D #31b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3s/H2_s_Nrad
""",
)

entry(
    index = 92,
    label = "H3N2 + O <=> HO-2 + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-2.7196, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + O = N2H2 + OH (B&D #31c3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3s/H2_s_Nrad
""",
)

entry(
    index = 93,
    label = "HO + H3N2 <=> H2O + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + OH = N2H2 + H2O (B&D #31d1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3s/H2_s_Nrad
""",
)

entry(
    index = 94,
    label = "H3N2 + CH3_r1 <=> CH4 + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.64e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (7.61488, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + CH3 = N2H2 + CH4 (B&D #31e1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;N3s/H2_s_Nrad
""",
)

entry(
    index = 95,
    label = "H2N + H3N2 <=> H3N + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (9.2e+05, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + NH2 = N2H2 + NH3 (B&D #31f1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;N3s/H2_s_Nrad
""",
)

entry(
    index = 96,
    label = "HO2-2 + H3N2 <=> H2O2 + H2N2-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (58000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + HO2 = N2H2 + H2O2 (B&D #31g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;N3s/H2_s_Nrad
""",
)

entry(
    index = 97,
    label = "HO2-3 + H3N2-2 <=> H4N2 + O2-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (8.91192, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + HO2 = N2H4 + O2 (B&D #31g3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Degeneracy not recalculated

Converted to training reaction from rate rule: N3s_rad/H/NonDeN;O_Orad
""",
)

entry(
    index = 98,
    label = "H + H2NO <=> H2 + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (9.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (6.52704, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + H = HNO + H2 (B&D #37c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3s/H2_s_Orad
""",
)

entry(
    index = 99,
    label = "H2NO + O <=> HO-2 + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (6.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (2.05016, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + O = HNO + OH (B&D #37d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3s/H2_s_Orad
""",
)

entry(
    index = 100,
    label = "HO + H2NO <=> H2O + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + OH = HNO + H2O (B&D #37e) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3s/H2_s_Orad
""",
)

entry(
    index = 101,
    label = "H2NO + CH3_r1 <=> CH4 + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.2e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (12.3846, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + CH3 = CH4 + HNO (B&D #37f2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;N3s/H2_s_Orad
""",
)

entry(
    index = 102,
    label = "H2N + H2NO <=> H3N + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.6e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + NH2 = HNO + NH3 (B&D #37g) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;N3s/H2_s_Orad
""",
)

entry(
    index = 103,
    label = "HO2-2 + H2NO <=> H2O2 + HNO-2",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (58000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + HO2 = HNO + H2O2 (B&D #37h1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;N3s/H2_s_Orad
""",
)

entry(
    index = 104,
    label = "HO2-3 + H2NO-2 <=> H3NO + O2-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + HO2 = NH2OH + O2 (B&D #37h2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Degeneracy not recalculated

Converted to training reaction from rate rule: O_rad/NonDeN;O_Orad
""",
)

entry(
    index = 105,
    label = "H + H2NO-3 <=> H2 + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (1.58992, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + H = HNO + H2 (B&D #38b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;O_Nrad
""",
)

entry(
    index = 106,
    label = "H2NO-3 + O <=> HO-2 + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-1.50624, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + O = HNO + OH (B&D #38c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;O_Nrad
""",
)

entry(
    index = 107,
    label = "HO + H2NO-3 <=> H2O + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + OH = HNO + H2O (B&D #38d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;O_Nrad
""",
)

entry(
    index = 108,
    label = "H2NO-3 + CH3_r1 <=> CH4 + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (8.7864, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + CH3 = CH4 + HNO (B&D #38e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;O_Nrad
""",
)

entry(
    index = 109,
    label = "H2N + H2NO-3 <=> H3N + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + NH2 = HNO + NH3 (B&D #38f3)  in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;O_Nrad
""",
)

entry(
    index = 110,
    label = "HO2-2 + H2NO-3 <=> H2O2 + HNO-3",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (29400, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + HO2 = HNO + H2O2 (B&D #38g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;O_Nrad
""",
)

entry(
    index = 111,
    label = "HO2-3 + H2NO-4 <=> H3NO-2 + O2-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.6944, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + HO2 = NH2OH + O2 (B&D #38g3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Degeneracy not recalculated

Converted to training reaction from rate rule: N3s_rad/H/NonDeO;O_Orad
""",
)

entry(
    index = 112,
    label = "HO2-2 + CH2N <=> H2O2 + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (28000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.73624, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + HO2 = HCN + H2O2 (B&D #45b1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_rad/NonDeO;Cds/H2_d_N3rad
""",
)

entry(
    index = 113,
    label = "HO2-3 + CH2N-2 <=> CH3N + O2-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (14000, 'cm^3/(mol*s)'),
        n = 2.69,
        Ea = (-6.73624, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + HO2 = H2CNH + O2 (B&D #45b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Degeneracy not recalculated

Converted to training reaction from rate rule: N3d_rad/C;O_Orad
""",
)

entry(
    index = 114,
    label = "CH3_r1 + CH2N <=> CH4 + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.62e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-4.64424, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + CH3 = HCN + CH4 (B&D #45d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;Cds/H2_d_N3rad
""",
)

entry(
    index = 115,
    label = "HO + CH2N <=> H2O + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + OH = HCN + H2O (B&D #45e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;Cds/H2_d_N3rad
""",
)

entry(
    index = 116,
    label = "H + CH2N <=> H2 + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + H = HCN + H2 (B&D #45g) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;Cds/H2_d_N3rad
""",
)

entry(
    index = 117,
    label = "H2N + CH2N <=> H3N + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (1.84e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + NH2 = HCN + NH3 (B&D #45h) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;Cds/H2_d_N3rad
""",
)

entry(
    index = 118,
    label = "CH2N + O <=> HO-2 + CHN",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + O = HCN + OH (B&D #45i1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;Cds/H2_d_N3rad
""",
)

entry(
    index = 119,
    label = "H + CH2N-3 <=> H2 + CHN-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + H = HCN + H2 (B&D #46a2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3d/H_d_Crad
""",
)

entry(
    index = 120,
    label = "CH2N-3 + O <=> HO-2 + CHN-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + O = HCN + OH (B&D #46b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3d/H_d_Crad
""",
)

entry(
    index = 121,
    label = "HO + CH2N-3 <=> H2O + CHN-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + OH = HCN + H2O (B&D #46c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3d/H_d_Crad
""",
)

entry(
    index = 122,
    label = "CH2N-3 + CH3_r1 <=> CH4 + CHN-2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (820000, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-4.64424, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + CH3 = HCN + CH4 (B&D #46d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;N3d/H_d_Crad
""",
)

entry(
    index = 123,
    label = "H + CH4N <=> H2 + CH3N-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (2.16e+09, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + H = H2CNH + H2 (B&D #49b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;Cmethyl_Nrad
""",
)

entry(
    index = 124,
    label = "CH4N + O <=> HO-2 + CH3N-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.5e+09, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + O = H2CNH + OH (B&D #49c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;Cmethyl_Nrad
""",
)

entry(
    index = 125,
    label = "HO + CH4N <=> H2O + CH3N-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.08e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + OH = H2CNH + H2O (B&D #49d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;Cmethyl_Nrad
""",
)

entry(
    index = 126,
    label = "CH4N + CH3_r1 <=> CH4 + CH3N-2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (7.2e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-4.64424, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + CH3 = H2CNH + CH4 (B&D #49e) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;Cmethyl_Nrad
""",
)

entry(
    index = 127,
    label = "H + CH4N-2 <=> H2 + CH3N-3",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + H = H2CNH + H2 (B&D #50b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;N3s/H2_s_Cssrad
""",
)

entry(
    index = 128,
    label = "CH4N-2 + O <=> HO-2 + CH3N-3",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (6.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + O = H2CNH + OH (B&D #50c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_atom_triplet;N3s/H2_s_Cssrad
""",
)

entry(
    index = 129,
    label = "HO + CH4N-2 <=> H2O + CH3N-3",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + OH = H2CNH + H2O (B&D #50d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;N3s/H2_s_Cssrad
""",
)

entry(
    index = 130,
    label = "CH3_r1 + CH4N-2 <=> CH4 + CH3N-3",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.2e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-2.63592, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + CH3 = H2CNH + CH4 (B&D #50e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;N3s/H2_s_Cssrad
""",
)

entry(
    index = 131,
    label = "H + CH2NO <=> H2 + CHNO",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (9.6e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-3.72376, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + H = HCNO + H2
The reacting structures are CH2=[N.+][O-] + R = [CH]#[N+][O-] + RH
(D&B #57c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: H_rad;Cds/H2_d_N5dcrad/O
""",
)

entry(
    index = 132,
    label = "HO + CH2NO <=> H2O + CHNO",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (4.8e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-4.97896, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + OH = HCNO + H2O
(D&B #57e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: O_pri_rad;Cds/H2_d_N5dcrad/O
""",
)

entry(
    index = 133,
    label = "CH3_r1 + CH2NO <=> CH4 + CHNO",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.2e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        Ea = (-4.64424, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + CH3 = HCNO + CH4
(D&B #57f2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: C_methyl;Cds/H2_d_N5dcrad/O
""",
)

entry(
    index = 134,
    label = "H2N + CH2NO <=> H3N + CHNO",
    degeneracy = 2.0,
    kinetics = Arrhenius(
        A = (3.6e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        Ea = (-4.8116, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 2,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + NH2 = HCNO + NH3
(D&B #57g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",

Converted to training reaction from rate rule: NH2_rad;Cds/H2_d_N5dcrad/O
""",
)

entry(
    index = 135,
    label = "C5H7 + C4H7-2 <=> C5H8 + C4H6",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.5e+11, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 11,
    shortDesc = u"""Estimated by S.S. Merchant""",
    longDesc = 
u"""
Estimating rate coefficient for cyclopentadienyl radical + butadieneyl radical
NIST estimate for allyl + iso-butyl is 8E+11 at 1000 K, however in our system the butadieneyl radical is also resonance stabilized
and it will be harder to break the bond to give butadiene + cyclopentadiene. Currently estimate it to be a factor of 5 slower.

Converted to training reaction from rate rule: C_rad/H/TwoDe;Cmethyl_Csrad/H/Cd
""",
)

entry(
    index = 136,
    label = "C6H5 + C6H5-2 <=> C6H6 + C6H4",
    degeneracy = 2.0,
    kinetics = Arrhenius(A=(1350, 'cm^3/(mol*s)'), n=2.7, Ea=(-4.403, 'kcal/mol'), T0=(1, 'K')),
    reference = Article(
        authors = ['Tranter, R. S.', 'Klippenstein, S. J.', 'Harding, L. B.', 'Giri, B. R.', 'Yang, X.', 'Kiefer, J. H.'],
        title = 'Experimental and Theoretical Investigation of the Self-Reaction of Phenyl Radicals',
        journal = 'The Journal of Physical Chemistry A',
        volume = '114 (32)',
        pages = '8240-8261',
        year = '2010',
    ),
    referenceType = "theory",
    rank = 3,
    longDesc = 
u"""
CASPT2(2e,2o)/cc-pvdz (VRC-TST)
""",
)

