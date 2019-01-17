#!/usr/bin/env python
# encoding: utf-8

name = "Disproportionation/rules"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 485,
    label = "Y_rad_birad_trirad_quadrad;XH_Rrad_birad",
    kinetics = ArrheniusEP(
        A = (3e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 0,
    shortDesc = u"""Default""",
)

entry(
    index = 487,
    label = "O2b;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (7.23e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (15.99, 'kcal/mol'),
        Tmin = (700, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 488,
    label = "CH2_triplet;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (3.01e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 489,
    label = "H_rad;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (3.61e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + i-C3H7 --> C3H6 + H2

pg. 932: Discussion on evaluated data

Entry 42,4 (a): No data available at the time.  Author recommends a rate coefficient

expression equal to double the rate expression of H+C2H5=H2+C2H4.
MRH 30-Aug-2009
""",
)

entry(
    index = 490,
    label = "H_rad;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (1.81e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [89] Literature review.""",
    longDesc = 
u"""
[89] Tsang, W.; Hampson, R.F.; Journal of Physical and Chemical Reference Data (1986) 15(3), 1087-1279.
Literature review.  H + C2H5 --> C2H4 + H2

pg. 1174: Discussion on evaluated data

Entry 17,4 (c): Author recommends rate coefficient from study performed by 

Camilleri, et al. (1974)
MRH 30-Aug-2009
""",
)

entry(
    index = 491,
    label = "C_methyl;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.19e+14, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.68,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 492,
    label = "C_rad/H2/Cs;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.3e+13, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.35,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 493,
    label = "C_rad/H2/Cd;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.29e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = -0.35,
        alpha = 0,
        E0 = (-0.13, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 494,
    label = "C_rad/H2/O;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.89e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + i-C3H7 --> C3H6 + CH3OH

pg. 945: Discussion on evaluated data

Entry 42,39 (c): No data available at the time.  Author recommends a rate coefficient

of 4.8x10^-12 based on the rate expression of i-C3H7+C2H5=C2H6+C3H6
MRH 30-Aug-2009
""",
)

entry(
    index = 495,
    label = "C_rad/H/NonDeC;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.11e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.7,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 496,
    label = "C_rad/Cs3;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.86e+15, 'cm^3/(mol*s)', '*|/', 1.7),
        n = -1.1,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 497,
    label = "Cd_pri_rad;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (1.52e+14, 'cm^3/(mol*s)', '*|/', 1.5),
        n = -0.7,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 498,
    label = "Ct_rad/Ct;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (3.61e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H + i-C3H7 --> C3H6 + C2H2

pg. 941-942: Discussion on evaluated data

Entry 42,21 (a): No data available at the time.  Author recommends a rate coefficient

of 6x10^-12 cm3/molecule/s, a "typical" disproportionation rate.
MRH 30-Aug-2009
""",
)

entry(
    index = 499,
    label = "O_pri_rad;Cmethyl_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + i-C3H7 --> C3H6 + H2O

pg. 934: Discussion on evaluated data

Entry 42,6: No data available at the time.  Author notes that both a H-atom abstraction

rxn and an addition + hot adduct decomposition rxn will result in the same products.
MRH 30-Aug-2009
""",
)

entry(
    index = 500,
    label = "H_rad;Cmethyl_Orad",
    kinetics = ArrheniusEP(
        A = (1.81e+13, 'cm^3/(mol*s)', '*|/', 3.16),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1000, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Baulch et al [95] literature review.""",
    longDesc = 
u"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; Troe, J.; Walker, R.W.; Warnatz, J.; Journal of Physical and Chemical Reference Data (1992), 21(3), 411-734.
pg.523: Discussion of evaluated data

H+CH3O --> H2+CH2O: Authors state that no new data have been reported for this reaction.

MRH assumes the recommended value comes from a previous review article published
by authors.  In any case, recommended data fits the reported data well.
MRH 31-Aug-2009
""",
)

entry(
    index = 501,
    label = "O2b;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (4.5825e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (14.85, 'kcal/mol'),
        Tmin = (500, 'K'),
        Tmax = (900, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 502,
    label = "CH2_triplet;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.81e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 502,
    label = "S_rad/OneDe;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (0.719, 'cm^3/(mol*s)'),
        n = 3.13,
        alpha = 0,
        E0 = (-3.65, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 3,
    shortDesc = u"""CAC calc CBS-QB3, 1dhr""",
)

entry(
    index = 503,
    label = "H_rad;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.81e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + n-C3H7 --> C3H6 + H2

pg. 915-916: Discussion on evaluated data

Entry 41,4 (a): No data available at the time.  Author recommends the rate coefficient

of the H+C2H5=C2H4+H2 rxn for the H+n-C3H7=C3H6+H2 rxn.
MRH 30-Aug-2009
""",
)

entry(
    index = 504,
    label = "C_rad/H/OneDeC;C/H2/Nd_Srad",
    kinetics = ArrheniusEP(
        A = (7.63e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (-0.55, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 504,
    label = "C_methyl;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.15e+13, 'cm^3/(mol*s)', '*|/', 1.7),
        n = -0.32,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 505,
    label = "C_rad/H2/Cs;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.45e+12, 'cm^3/(mol*s)', '*|/', 1.4),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 505,
    label = "C_rad/H/TwoDe;C/H2/Nd_Srad",
    kinetics = ArrheniusEP(
        A = (1.94e+12, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0.36, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 506,
    label = "C_rad/H2/Cd;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.45e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (-0.13, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 506,
    label = "S_rad/OneDe;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (7.63e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (-0.55, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""VERY Rough estimate based on 1/10 of #3026 in R_Recombination""",
)

entry(
    index = 507,
    label = "C_rad/H2/O;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (4.82e+11, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 507,
    label = "S_rad/NonDeC;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (7.63e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (-0.55, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 508,
    label = "C_rad/H/NonDeC;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (5.13e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.35,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 509,
    label = "C_rad/Cs3;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (2.16e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.75,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 510,
    label = "Cd_pri_rad;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.21e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H3 + n-C3H7 --> C3H6 + C2H4

pg. 922: Discussion on evaluated data

Entry 41,19 (a): No data available at the time.  Author estimates the rate coefficient

based on the rxn C2H5+n-C3H7=C3H6=C2H6.
MRH 30-Aug-2009
""",
)

entry(
    index = 510,
    label = "S_rad/NonDeS;C/H2/Nd_Csrad/H/Cd",
    kinetics = ArrheniusEP(
        A = (6.44e+08, 'cm^3/(mol*s)'),
        n = 1.19,
        alpha = 0,
        E0 = (0.51, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""Very rough based on R_Recomb #491""",
)

entry(
    index = 511,
    label = "Ct_rad/Ct;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 511,
    label = "S_rad/NonDeS;S_Csrad",
    kinetics = ArrheniusEP(
        A = (6.44e+08, 'cm^3/(mol*s)'),
        n = 1.19,
        alpha = 0,
        E0 = (0.51, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""Very rough based on R_Recomb #491""",
)

entry(
    index = 512,
    label = "O_pri_rad;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + n-C3H7 --> C3H6 + H2O

pg. 917: Discussion on evaluated data

Entry 41,6 (a): No data available at the time.  Author estimates rate coefficient based

on the rate coefficient for OH+C2H5=C2H4+H2O, namely 4.0x10^-11 cm3/molecule/s.
MRH 30-Aug-2009
""",
)

entry(
    index = 513,
    label = "O2b;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.2044e+10, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (600, 'K'),
        Tmax = (1000, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 514,
    label = "Ct_rad/Ct;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 515,
    label = "H_rad;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (9.04e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 516,
    label = "C_methyl;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (6.02e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.32,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 517,
    label = "C_rad/H2/Cs;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (8.43e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 518,
    label = "C_rad/H2/O;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+11, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 519,
    label = "C_rad/H2/Cd;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (7.83e+11, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (-0.13, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 520,
    label = "C_rad/H/NonDeC;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (2.56e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.35,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 521,
    label = "C_rad/Cs3;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.08e+14, 'cm^3/(mol*s)', '*|/', 2),
        n = -0.75,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 522,
    label = "Cd_pri_rad;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (8.43e+11, 'cm^3/(mol*s)', '*|/', 4),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [92] Literature review.""",
    longDesc = 
u"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H3 + i-C4H9 --> i-C4H8 + C2H4

pg. 60: Discussion on evaluated data

Entry 45,19 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient based on the rate of C2H5+i-C4H9=i-C4H8+C2H6.
MRH 30-Aug-2009
""",
)

entry(
    index = 523,
    label = "O_pri_rad;C/H/NdNd_Csrad",
    kinetics = ArrheniusEP(
        A = (1.21e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 524,
    label = "O2b;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (6.022e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (13.55, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 526,
    label = "C_methyl;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (3.01e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 527,
    label = "C_rad/H2/Cs;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (9.64e+11, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 528,
    label = "C_rad/H2/Cd;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (8.43e+10, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 529,
    label = "C_rad/H/NonDeC;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (4.58e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 530,
    label = "C_rad/Cs3;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (2.89e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
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
""",
)

entry(
    index = 531,
    label = "Cd_pri_rad;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C2H3 + C3H5 --> H2C=C=CH2 + C2H4

pg.261-262: Discussion on evaluated data

Entry 47,19(d): No data available at the time.  Author recommends a rate coefficient

of 4x10^-12 cm3/molecule/s for the disproportionation rxn.
MRH 31-Aug-2009
""",
)

entry(
    index = 532,
    label = "O_pri_rad;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""SSM estimate. Original value with 6 kcal barrier""",
    longDesc = 
u"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: OH + C3H5 --> H2C=C=CH2 + H2O

pg.253: Discussion on evaluated data

Entry 47,6(a): No data available at the time.  Author recommends a rate coefficient

of 1x10^-11 cm3/molecule/s, based on "comparable rxns".
MRH 31-Aug-2009
""",
)

entry(
    index = 533,
    label = "O2b;O_Csrad",
    kinetics = ArrheniusEP(
        A = (5.7209e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (298, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 534,
    label = "O2b;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.92067e+12, 'cm^3/(mol*s)', '*|/', 1.3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (298, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Atkinson et al [98] literature review.""",
    longDesc = 
u"""
[98] Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F., Jr.; Kerr, J.A.; Rossi, M.J.; Troe, J. "Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry,", 2001.
Literature review: CH2OH + O2 --> CH2O + HO2

Recommended value is k298.  This reference just gives a table of results,

with no discussion on how the preferred numbers were arrived at.
MRH 31-Aug-2009

Divide the rate constant by 2 to account for symmetry of 2 (O2) and 1 (CH2OH, oxygen atom). The final result is 2.92067e+12 cm3/mol/s.
JDM 31-Mar-2010
""",
)

entry(
    index = 535,
    label = "O2b;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.74001e+12, 'cm^3/(mol*s)', '*|/', 1.3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol', '+|-', 0.4),
        Tmin = (200, 'K'),
        Tmax = (300, 'K'),
    ),
    rank = 4,
    shortDesc = u"""DeMore et al [183] literature review.""",
    longDesc = 
u"""
[183] DeMore, W.B.; Sander, S.P.; Golden, D.M.; Hampson, R.F.; Kurylo, M.J.; Howard, C.J.; Ravishankara, A.R.; Kolb, C.E.; Molina, M.J.; JPL Publication 97-4
Literature review: CH2OH + O2 --> CH2O + HO2

pg.62 D38: Discussion on evaluated data

pg.22: Recommended A-factor and E/R parameter values

MRH 1-Sept-2009

Divide the rate constant by 2 to account for symmetry of 2 (O2) and 1 (CH2OH, oxygen atom). The final result is 2.74001e+12 cm3/mol/s.
JDM 31-Mar-2010
""",
)

entry(
    index = 536,
    label = "O_atom_triplet;O_Csrad",
    kinetics = ArrheniusEP(
        A = (9.04e+13, 'cm^3/(mol*s)', '+|-', 3.01e+13),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (298, 'K'),
    ),
    rank = 3,
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
""",
)

entry(
    index = 537,
    label = "CH2_triplet;O_Csrad",
    kinetics = ArrheniusEP(
        A = (1.21e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2 + CH2OH --> CH3 + CH2O

pg. 505: Discussion on evaluated data

Entry 39,26 (b): CH2OH + CH2(triplet) --> CH3 + CH2O

Author estimates the rate of disproportionation as 2.0x10^-12 cm3/molecule/s.  No data at the time.

MRH 30-Aug-2009
""",
)

entry(
    index = 538,
    label = "H_rad;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2e+13, 'cm^3/(mol*s)', '+|-', 1e+13),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (295, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 539,
    label = "H_rad;O_Csrad",
    kinetics = ArrheniusEP(
        A = (6.03e+12, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: H + CH2OH --> H2 + CH2O

pg. 496-497: Discussion on evaluated data

Entry 39,4 (a): CH2OH + H --> H2 + CH2O

Author estimates disproportionation rate will be faster than the H+C2H5=H2+C2H4 reaction

and reports rate coefficient as 1.0x10^-11 cm3/molecule/s.  No data at the time.
MRH 30-Aug-2009
""",
)

entry(
    index = 540,
    label = "C_methyl;O_Csrad",
    kinetics = ArrheniusEP(
        A = (8.49e+13, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (298, 'K'),
    ),
    rank = 3,
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
""",
)

entry(
    index = 541,
    label = "C_methyl;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [90] Literature review.""",
    longDesc = 
u"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH3 + CH2OH --> CH4 + CH2O

pg. 500-501: Discussion on evaluated data

Entry 39,16 (b): CH2OH + CH3 --> CH4 + CH2O

Author estimates ratio of disproportionation rate to addition rate to be 0.2,

namely 4x10^-12 cm3/molecule/s.  No data at the time.
MRH 30-Aug-2009
""",
)

entry(
    index = 542,
    label = "C_rad/H2/Cs;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 543,
    label = "C_rad/H2/Cd;O_Csrad",
    kinetics = ArrheniusEP(
        A = (1.81e+13, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 544,
    label = "C_rad/H2/O;O_Csrad",
    kinetics = ArrheniusEP(
        A = (4.82e+12, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 545,
    label = "C_rad/H/NonDeC;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.35e+12, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tsang [91] Literature review.""",
    longDesc = 
u"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: CH2OH + i-C3H7 = C3H8 + CH2O

pg. 945: Discussion on evaluated data

Entry 42,39 (b): No data available at the time.  Author suggests rate coefficient based

on rxn C2H5+i-C3H7=C3H8+C2H4, namely 3.9x10^-12 cm3/molecule/s
MRH 30-Aug-2009
""",
)

entry(
    index = 546,
    label = "C_rad/Cs3;O_Csrad",
    kinetics = ArrheniusEP(
        A = (3.47e+14, 'cm^3/(mol*s)', '*|/', 3),
        n = -0.75,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 547,
    label = "Cd_pri_rad;O_Csrad",
    kinetics = ArrheniusEP(
        A = (3.01e+13, 'cm^3/(mol*s)', '*|/', 2.5),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 549,
    label = "CO_pri_rad;O_Csrad",
    kinetics = ArrheniusEP(
        A = (1.81e+14, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 550,
    label = "O_pri_rad;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 551,
    label = "O_rad/NonDeC;O_Csrad",
    kinetics = ArrheniusEP(
        A = (2.41e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 552,
    label = "O_rad/NonDeO;O_Csrad",
    kinetics = ArrheniusEP(
        A = (1.21e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
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
""",
)

entry(
    index = 553,
    label = "S_rad/NonDeC;Cmethyl_Srad",
    kinetics = ArrheniusEP(
        A = (9.79e+11, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (298, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Tycholiz et al [A].""",
)

entry(
    index = 554,
    label = "C_rad/H/CsS;C/H2/Nd_Csrad",
    kinetics = ArrheniusEP(
        A = (3.37e-06, 'cm^3/(mol*s)'),
        n = 4.35,
        alpha = 0,
        E0 = (1.14, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 3,
    shortDesc = u"""CAC calc CBS-QB3 1dhr""",
)

entry(
    index = 555,
    label = "O2b;O_Csrad",
    kinetics = ArrheniusEP(
        A = (3.7e+16, 'cm^3/(mol*s)'),
        n = -1.63,
        alpha = 0,
        E0 = (3.4, 'kcal/mol'),
        Tmin = (250, 'K'),
        Tmax = (1000, 'K'),
    ),
    rank = 4,
    shortDesc = u"""S.S. Merchant estimate""",
    longDesc = 
u"""
Estimate on basis of C3H7 + O2 rate from NIST kinetic datbase, Measurements, Theory, and Modeling of OH Formation in Ethyl + O2 and Propyl + O2 Reactions
ref: DOI: 10.1021/jp0221946
""",
)

entry(
    index = 555,
    label = "O2b;Cd_Cdrad",
    kinetics = ArrheniusEP(
        A = (1.3e+15, 'cm^3/(mol*s)', '*|/', 5),
        n = -1.26,
        alpha = 0,
        E0 = (3.31, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 5,
    shortDesc = u"""S.S. Merchant estimate""",
    longDesc = 
u"""
This rate rule is a estimate taken from NIST, ref: Aromatic and Polycyclic Aromatic
Hydrocarbon Formation in a Laminar Premixed n-butane Flame
Derived from fitting to a complex mechanism for C2H3 + O2 = C2H2 + HO2
""",
)

entry(
    index = 555,
    label = "O2b;Cmethyl_Csrad/H/Cd",
    kinetics = ArrheniusEP(
        A = (7.23e+12, 'cm^3/(mol*s)', '*|/', 10),
        n = 0,
        alpha = 0,
        E0 = (22, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
    ),
    rank = 5,
    shortDesc = u"""S.S. Merchant estimate""",
    longDesc = 
u"""
SSM estimate based on Miyoshi rate rule for secondary carbon in dx.doi.org/10.1021/jp112152n, 
modified to account for allylic stability (+7 kcal)
""",
)

entry(
    index = 556,
    label = "O2b;XH_Rrad_birad",
    kinetics = ArrheniusEP(
        A = (5e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""A.G. Vandeputte estimated value""",
)

entry(
    index = 556,
    label = "Y_rad_birad_trirad_quadrad;Cdpri_Csrad",
    kinetics = ArrheniusEP(
        A = (1e+09, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 557,
    label = "O2b;C/H2/De_Csrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 558,
    label = "O2b;C/H2/Nd_Rrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 559,
    label = "O2b;C/H2/De_Rrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 560,
    label = "O2b;C/H/NdNd_Rrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 561,
    label = "O2b;C/H/NdDe_Rrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 562,
    label = "O2b;C/H/DeDe_Rrad",
    kinetics = ArrheniusEP(
        A = (1e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 4,
    shortDesc = u"""Estimated value, AG Vandeputte""",
)

entry(
    index = 600,
    label = "NH2_rad;O_Orad",
    kinetics = ArrheniusEP(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2 + HO2 = NH3 + O2 (B&D #14d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 601,
    label = "O2b;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (1.2e+12, 'cm^3/(mol*s)'),
        n = -0.34,
        alpha = 0,
        E0 = (0.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + O2 = N2 + HO2 (B&D #28b1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 602,
    label = "H_rad;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + H = N2 + H2 (B&D #28c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 603,
    label = "O_pri_rad;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + OH = N2 + H2O (B&D #28d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 604,
    label = "O_atom_triplet;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + O = N2 + OH (B&D #28e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 605,
    label = "NH2_rad;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + NH2 = N2 + NH3 (B&D #28f) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 606,
    label = "O_rad/NonDeO;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (14000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + HO2 = N2 + H2O2 (B&D #28g1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 607,
    label = "N3d_rad/O;N3d/H_d_Nrad",
    kinetics = ArrheniusEP(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NNH + NO = N2 + HNO (B&D #28h) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 608,
    label = "H_rad;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + H = NNH + H2 (B&D #30c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 609,
    label = "O_atom_triplet;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + O = NNH + OH (B&D #30d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 610,
    label = "O_pri_rad;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + OH = NNH + H2O (B&D #30e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 611,
    label = "C_methyl;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (0.13, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + CH3 = NNH + CH4 (B&D #30f3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 612,
    label = "NH2_rad;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + NH2 = NNH + NH3 (B&D #30g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 613,
    label = "O_rad/NonDeO;N3s/H2_s_Nbirad",
    kinetics = ArrheniusEP(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2NN + HO2 = NNH + H2O2 (B&D #30h2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 614,
    label = "H_rad;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + H = N2H2 + H2 (B&D #31b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 615,
    label = "O_atom_triplet;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.65, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + O = N2H2 + OH (B&D #31c3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 616,
    label = "O_pri_rad;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + OH = N2H2 + H2O (B&D #31d1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 617,
    label = "C_methyl;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (820000, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (1.82, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + CH3 = N2H2 + CH4 (B&D #31e1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 618,
    label = "NH2_rad;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + NH2 = N2H2 + NH3 (B&D #31f1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 619,
    label = "O_rad/NonDeO;N3s/H2_s_Nrad",
    kinetics = ArrheniusEP(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + HO2 = N2H2 + H2O2 (B&D #31g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 620,
    label = "N3s_rad/H/NonDeN;O_Orad",
    kinetics = ArrheniusEP(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (2.13, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: N2H3 + HO2 = N2H4 + O2 (B&D #31g3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 621,
    label = "H_rad;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (1.56, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + H = HNO + H2 (B&D #37c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 622,
    label = "O_atom_triplet;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (0.49, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + O = HNO + OH (B&D #37d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 623,
    label = "O_pri_rad;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + OH = HNO + H2O (B&D #37e) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 624,
    label = "C_methyl;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (2.96, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + CH3 = CH4 + HNO (B&D #37f2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 625,
    label = "NH2_rad;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + NH2 = HNO + NH3 (B&D #37g) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 626,
    label = "O_rad/NonDeO;N3s/H2_s_Orad",
    kinetics = ArrheniusEP(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + HO2 = HNO + H2O2 (B&D #37h1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 627,
    label = "O_rad/NonDeN;O_Orad",
    kinetics = ArrheniusEP(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: NH2O + HO2 = NH2OH + O2 (B&D #37h2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 628,
    label = "H_rad;O_Nrad",
    kinetics = ArrheniusEP(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (0.38, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + H = HNO + H2 (B&D #38b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 629,
    label = "O_atom_triplet;O_Nrad",
    kinetics = ArrheniusEP(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.36, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + O = HNO + OH (B&D #38c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 630,
    label = "O_pri_rad;O_Nrad",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + OH = HNO + H2O (B&D #38d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 631,
    label = "C_methyl;O_Nrad",
    kinetics = ArrheniusEP(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (2.1, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + CH3 = CH4 + HNO (B&D #38e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 632,
    label = "NH2_rad;O_Nrad",
    kinetics = ArrheniusEP(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + NH2 = HNO + NH3 (B&D #38f3)  in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 633,
    label = "O_rad/NonDeO;O_Nrad",
    kinetics = ArrheniusEP(
        A = (29400, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + HO2 = HNO + H2O2 (B&D #38g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 634,
    label = "N3s_rad/H/NonDeO;O_Orad",
    kinetics = ArrheniusEP(
        A = (29000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.6, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HNOH + HO2 = NH2OH + O2 (B&D #38g3) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 635,
    label = "O_rad/NonDeO;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (14000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.61, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + HO2 = HCN + H2O2 (B&D #45b1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 636,
    label = "N3d_rad/C;O_Orad",
    kinetics = ArrheniusEP(
        A = (14000, 'cm^3/(mol*s)'),
        n = 2.69,
        alpha = 0,
        E0 = (-1.61, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + HO2 = H2CNH + O2 (B&D #45b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 637,
    label = "C_methyl;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (810000, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (-1.11, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + CH3 = HCN + CH4 (B&D #45d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 638,
    label = "O_pri_rad;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + OH = HCN + H2O (B&D #45e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 639,
    label = "H_rad;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + H = HCN + H2 (B&D #45g) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 640,
    label = "NH2_rad;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (920000, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + NH2 = HCN + NH3 (B&D #45h) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 641,
    label = "O_atom_triplet;Cds/H2_d_N3rad",
    kinetics = ArrheniusEP(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: H2CN + O = HCN + OH (B&D #45i1) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 642,
    label = "H_rad;N3d/H_d_Crad",
    kinetics = ArrheniusEP(
        A = (2.4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + H = HCN + H2 (B&D #46a2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 643,
    label = "O_atom_triplet;N3d/H_d_Crad",
    kinetics = ArrheniusEP(
        A = (1.7e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + O = HCN + OH (B&D #46b2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 644,
    label = "O_pri_rad;N3d/H_d_Crad",
    kinetics = ArrheniusEP(
        A = (1.2e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + OH = HCN + H2O (B&D #46c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 645,
    label = "C_methyl;N3d/H_d_Crad",
    kinetics = ArrheniusEP(
        A = (820000, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (-1.11, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: HCNH + CH3 = HCN + CH4 (B&D #46d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 646,
    label = "H_rad;Cmethyl_Nrad",
    kinetics = ArrheniusEP(
        A = (7.2e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + H = H2CNH + H2 (B&D #49b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 647,
    label = "O_atom_triplet;Cmethyl_Nrad",
    kinetics = ArrheniusEP(
        A = (5e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + O = H2CNH + OH (B&D #49c) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 648,
    label = "O_pri_rad;Cmethyl_Nrad",
    kinetics = ArrheniusEP(
        A = (3.6e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + OH = H2CNH + H2O (B&D #49d) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 649,
    label = "C_methyl;Cmethyl_Nrad",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (-1.11, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH3NH + CH3 = H2CNH + CH4 (B&D #49e) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 650,
    label = "H_rad;N3s/H2_s_Cssrad",
    kinetics = ArrheniusEP(
        A = (4e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + H = H2CNH + H2 (B&D #50b) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 651,
    label = "O_atom_triplet;N3s/H2_s_Cssrad",
    kinetics = ArrheniusEP(
        A = (3.3e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + O = H2CNH + OH (B&D #50c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 652,
    label = "O_pri_rad;N3s/H2_s_Cssrad",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + OH = H2CNH + H2O (B&D #50d2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 653,
    label = "C_methyl;N3s/H2_s_Cssrad",
    kinetics = ArrheniusEP(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (-0.63, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NH2 + CH3 = H2CNH + CH4 (B&D #50e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 654,
    label = "H_rad;Cds/H2_d_N5ddrad/O",
    kinetics = ArrheniusEP(
        A = (4.8e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        alpha = 0,
        E0 = (-0.89, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + H = HCNO + H2 (B&D #57c2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 655,
    label = "O_pri_rad;Cds/H2_d_N5ddrad/O",
    kinetics = ArrheniusEP(
        A = (2.4e+06, 'cm^3/(mol*s)'),
        n = 2,
        alpha = 0,
        E0 = (-1.19, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + OH = HCNO + H2O (B&D #57e2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 656,
    label = "C_methyl;Cds/H2_d_N5ddrad/O",
    kinetics = ArrheniusEP(
        A = (1.6e+06, 'cm^3/(mol*s)'),
        n = 1.87,
        alpha = 0,
        E0 = (-1.11, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + CH3 = HCNO + CH4 (B&D #57f2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 657,
    label = "NH2_rad;Cds/H2_d_N5ddrad/O",
    kinetics = ArrheniusEP(
        A = (1.8e+06, 'cm^3/(mol*s)'),
        n = 1.94,
        alpha = 0,
        E0 = (-1.15, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 1,
    shortDesc = u"""Added by Beat Buesser""",
    longDesc = 
u"""
Added by Beat Buesser, value for reaction: CH2NO + NH2 = HCNO + NH3 (B&D #57g2) in 'Gas-Phase Combustion Chemistry' (ISBN: 978-1-4612-7088-1), chapter 2, 'Combustion Chemistry of Nitrogen', Anthony M. Dean, Joseph W. Bozzelli",
""",
)

entry(
    index = 658,
    label = "C_rad/H/TwoDe;Cmethyl_Csrad/H/Cd",
    kinetics = ArrheniusEP(
        A = (5e+10, 'cm^3/(mol*s)'),
        n = 0,
        alpha = 0,
        E0 = (0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
    ),
    rank = 5,
    shortDesc = u"""Estimated by S.S. Merchant""",
    longDesc = 
u"""
Estimating rate coefficient for cyclopentadienyl radical + butadieneyl radical
NIST estimate for allyl + iso-butyl is 8E+11 at 1000 K, however in our system the butadieneyl radical is also resonance stabilized
and it will be harder to break the bond to give butadiene + cyclopentadiene. Currently estimate it to be a factor of 5 slower.
""",
)

