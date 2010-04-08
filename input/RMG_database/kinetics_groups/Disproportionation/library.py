# encoding: utf-8
header = """
Disproportionation

Y_rad_birad + XH_Rrad -> Y_H + X_R


Reverse name: Molecular_Addition

(1) FORM_BOND		{*1,S,*4}
(2) BREAK_BOND		{*2,S,*4}
(3) CHANGE_BOND		{*2,1,*3}
(4) LOSE_RADICAL 	{*1,1}
(5) LOSE_RADICAL	{*3,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Disproportionation"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f09a: Disproportionation reaction
// original from rate_library_4.txt, Cath, 03/07/28

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f09a_Disproportionation
// changing root number into 3e11, 0, 0 according to Wing Tsang . J. Phys. Chem. Ref. Data 1991, 20, 221-273, JS, july 27, 2003
// add 10001 to give HO2. + R. -> O2 + HR, JS, july 27, 2003

// JS, Aug, 5, 2003
// move the following rate constants of O2 + XH_Rrad into the new reaction family: Disproportionation_O2d
// move #: 487, 501, 513, 514, 533, 534, 535

// Catherina Wijaya thesis, pg 157

// [90] Tsang, W. J. Phys. Chem. Ref. Data 1987, 16, 471.
// [91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
// [92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
// [93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
// [95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; 
//		Troe, J.; Walker, R.W.; Warnatz, J. J. Phys. Chem. Ref. Data 1992, 21, 411.
// [189] Grotheer, H.; Riekert, G.; Walter, D.; Just, T. Symp. Int. Combust. Proc. 1989, 22, 963.
// [190] Edelbuttel-Einhaus, J.; Hoyermann, K.; Rohde, G.; Seeba, J. Symp. Int. Combust. Proc. 1992, 22, 661.
// [191] Pagsberg, P.; Munk, J.; Sillesen, A.; Anastasi, C. Chem. Phys. Lett. 1988, 146, 375.

 
//No.		Y_rad_birad		XH_Rrad				Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments
//491.		C_methyl		Cmethyl_Csrad		300-2500	9.41E+10	0.68	0		0		*1.5	0		0		0		4		Tsang [91] Literature review.
//514.		CH2_triplet		C/H/NdNd_Csrad		300-2500	6.03E+12	0		0		0		*3.0	0		0		0		4		Tsang [92] Literature review. C2H + iso-C4H9--> C2H2 + iso-C4H8


"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
Y_rad_birad
Union {Y_2centeradjbirad, Y_1centerbirad, Y_rad}
""",
  group2 = 
"""
XH_Rrad
1. *2 {R!H} 0 {2,S}, {3,S}
2. *3 {R!H} 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "485",
  short_comment = "Default",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
O_sec_rad
1 *1 O 1 {2,S}
2 {R!H} 0 {1,S}
""",
  group2 = 
"""
O_Orad
1. *2 O 0 {2,S}, {3,S}
2. *3 O 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.75E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-3.275,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "486",
  short_comment = "[8] Curran\'s estimation in reaction type 13, RO2 + HO2 --> RO2H + O2",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
CH2_triplet
1 *1 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(3.01E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "488",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2(triplet) + i-C3H7 --> C3H6 + CH3

pg. 944: Discussion on evaluated data

Entry 42,26: No data available at the time.  Author suggests this is a minor channel,

stating the main process should be combination, leading to chemically activated
i-butyl radical.  Rate coefficient is estimate.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(3.61E+12,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "489",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + i-C3H7 --> C3H6 + H2

pg. 932: Discussion on evaluated data

Entry 42,4 (a): No data available at the time.  Author recommends a rate coefficient

expression equal to double the rate expression of H+C2H5=H2+C2H4.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "490",
  short_comment = "Tsang [89] Literature review.",
  long_comment = 
"""
[89] Tsang, W.; Hampson, R.F.; Journal of Physical and Chemical Reference Data (1986) 15(3), 1087-1279.
Literature review.  H + C2H5 --> C2H4 + H2

pg. 1174: Discussion on evaluated data

Entry 17,4 (c): Author recommends rate coefficient from study performed by 

Camilleri, et al. (1974)
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.19E+14,A_UNITS,"*/",1.1),
                 n=(-0.68,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "491",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
C_rad/H2/Cs
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.3E+13,A_UNITS,"*/",1.1),
                 n=(-0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "492",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H5 + i-C3H7 --> C3H6 + C2H6

pg. 937-938: Discussion on evaluated data

Entry 42,17 (c): No data available at the time.  The rate coefficient expression for

the combination rxn is computed using the geometric mean rule and is reported as
2.6x10^-11 * (300/T)^0.35 cm3/molecule/s.  The recommended branching ratio for 
disproportionation to addition is that reported by Gibian and Corley (1973).
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
C_rad/H2/Cd
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.29E+13,A_UNITS,"*/",3.0),
                 n=(-0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "493",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + iC3H7 --> C3H6 + C3H6

pg.268: Discussion on evaluated data

Entry 47,42(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for \"alkyl radicals\" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.2.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC3H7+iC3H7-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
C_rad/H2/O
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.89E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "494",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + i-C3H7 --> C3H6 + CH3OH

pg. 945: Discussion on evaluated data

Entry 42,39 (c): No data available at the time.  Author recommends a rate coefficient

of 4.8x10^-12 based on the rate expression of i-C3H7+C2H5=C2H6+C3H6
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
C_rad/H/NonDeC
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.11E+14,A_UNITS,"*/",2.0),
                 n=(-0.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "495",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  i-C3H7 + i-C3H7 --> C3H6 + C3H8

pg. 946-947: Discussion on evaluated data

Entry 42,42 (b): No high-Temperature data available.  Author has fit rate coefficient

expression for addition rxn to 4 sets of experimental data.  Recommended branching
ratio agrees well with most of the experimental data.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
C_rad/Cs3
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.86E+15,A_UNITS,"*/",1.7),
                 n=(-1.10,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "496",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + i-C3H7 --> C3H6 + i-C4H10

pg. 46: Discussion on evaluated data

Entry 44,42 (a): The author computes the combination rate expression using the geometric

mean rule (of the rxns t-C4H9+t-C4H9-->adduct and i-C3H7+i-C3H7-->adduct).  The
disproportionation rate coefficient expression was then computed using the
reported branching ratio.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
rate(
  group1 = 
"""
Cd_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(1.52E+14,A_UNITS,"*/",1.5),
                 n=(-0.70,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "497",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H3 + i-C3H7 --> C3H6 + C2H4

pg. 939-940: Discussion on evaluated data

Entry 42,19 (a): No data available at the time.  Author recommends the rate coefficient

expression of C2H5+i-C3H7 for the rate expression for C2H3+i-C3H7.  Author also
recommends the branching ratio of disproportionation to addition of the 
C2H5+i-C3H7 system for the C2H3+i-C3H7 system.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
rate(
  group1 = 
"""
Ct_rad
1 *1 C 1 {2,T}
2 C 0 {1,T}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(3.61E+12,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "498",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H + i-C3H7 --> C3H6 + C2H2

pg. 941-942: Discussion on evaluated data

Entry 42,21 (a): No data available at the time.  Author recommends a rate coefficient

of 6x10^-12 cm3/molecule/s, a \"typical\" disproportionation rate.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
O_pri_rad
1 *1 O 1 {2,S}
2 H 0 {1,S}
""",
  group2 = 
"""
Cmethyl_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+13,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "499",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + i-C3H7 --> C3H6 + H2O

pg. 934: Discussion on evaluated data

Entry 42,6: No data available at the time.  Author notes that both a H-atom abstraction

rxn and an addition + hot adduct decomposition rxn will result in the same products.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 15
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
Cmethyl_Orad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 O 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. H 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+13,A_UNITS,"*/",3.16),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1000),
  rank = 4,
  old_id = "500",
  short_comment = "Baulch et al [95] literature review.",
  long_comment = 
"""
[95] Baulch, D.L.; Cobos, C.J.; Cox, R.A.; Esser, C.; Frank, P.; Just, T.; Kerr, J.A.; Pilling, M.J.; Troe, J.; Walker, R.W.; Warnatz, J.; Journal of Physical and Chemical Reference Data (1992), 21(3), 411-734.
pg.523: Discussion of evaluated data

H+CH3O --> H2+CH2O: Authors state that no new data have been reported for this reaction.

MRH assumes the recommended value comes from a previous review article published
by authors.  In any case, recommended data fits the reported data well.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 16
rate(
  group1 = 
"""
CH2_triplet
1 *1 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "502",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2_triplet + n-C3H7 --> C3H6 + CH3

pg. 925: Discussion on evaluated data

Entry 41,26: No data available at the time.  Author estimates the rate coefficient

expression of the addition rxn.  The author then recommends that the disproportionation
rate coefficient not exceed 10% of the combination rate.  Thus, the rate coefficient
is an upper limit.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 17
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+12,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "503",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  H + n-C3H7 --> C3H6 + H2

pg. 915-916: Discussion on evaluated data

Entry 41,4 (a): No data available at the time.  Author recommends the rate coefficient

of the H+C2H5=C2H4+H2 rxn for the H+n-C3H7=C3H6+H2 rxn.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 18
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.15E+13,A_UNITS,"*/",1.7),
                 n=(-0.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "504",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 19
rate(
  group1 = 
"""
C_rad/H2/Cs
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.45E+12,A_UNITS,"*/",1.4),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "505",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 20
rate(
  group1 = 
"""
C_rad/H2/Cd
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.45E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "506",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + nC3H7 --> C3H6 + C3H6

pg.268: Discussion on evaluated data

Entry 47,41(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for \"alkyl radicals\" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.07.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and nC3H7+nC3H7-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 21
rate(
  group1 = 
"""
C_rad/H2/O
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(4.82E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "507",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  CH2OH + n-C3H7 --> C3H6 + CH3OH

pg. 926: Discussion on evaluated data

Entry 41,39 (c): No data available at the time.  Author estimates the rate coefficient

for the addition rxn to be similar to the rate for n-C3H7+n-C3H7=>adduct.  Author
also estimates the branching ratio of disproportionation to addition as 0.051
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 22
rate(
  group1 = 
"""
C_rad/H/NonDeC
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(5.13E+13,A_UNITS,"*/",2.0),
                 n=(-0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "508",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 23
rate(
  group1 = 
"""
C_rad/Cs3
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.16E+14,A_UNITS,"*/",2.0),
                 n=(-0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "509",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: t-C4H9 + n-C3H7 --> C3H6 + i-C4H10

pg. 45: Discussion on evaluated data

Entry 44,41 (a): No data available at the time.  Author estimates the rate expression

for the combination rxn using the geometric mean rule (of the rxns t-C4H9+t-C4H9-->adduct
and n-C3H7+n-C3H7-->adduct).  The author then estimates the disproportionation
rate expression using the branching ratio; the branching ratio is from \"analogous
processes\".
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 24
rate(
  group1 = 
"""
Cd_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "510",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H3 + n-C3H7 --> C3H6 + C2H4

pg. 922: Discussion on evaluated data

Entry 41,19 (a): No data available at the time.  Author estimates the rate coefficient

based on the rxn C2H5+n-C3H7=C3H6=C2H6.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 25
rate(
  group1 = 
"""
Ct_rad
1 *1 C 1 {2,T}
2 C 0 {1,T}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "511",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  C2H + n-C3H7 --> C3H6 + C2H2

pg. 923: Discussion on evaluated data

Entry 41,21 (a): No data available at the time.  Author notes that the rxn is more exothermic

than the rxn CH3+n-C3H7=C3H6+CH4 and suggests a rate coefficient 3x larger,
namely 1.0x10^-11 cm3/molecule/s.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 26
rate(
  group1 = 
"""
O_pri_rad
1 *1 O 1 {2,S}
2 H 0 {1,S}
""",
  group2 = 
"""
C/H2/Nd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. H 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+13,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "512",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review.  OH + n-C3H7 --> C3H6 + H2O

pg. 917: Discussion on evaluated data

Entry 41,6 (a): No data available at the time.  Author estimates rate coefficient based

on the rate coefficient for OH+C2H5=C2H4+H2O, namely 4.0x10^-11 cm3/molecule/s.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 27
rate(
  group1 = 
"""
Ct_rad
1 *1 C 1 {2,T}
2 C 0 {1,T}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "514",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H + i-C4H9 --> i-C4H8 + C2H2

pg. 61: Discussion on evaluated data

Entry 45,21: No data available at the time.  The author estimates the rate of 

disproportionation to be 1x10^-11 cm3/molecule/s.
*** NOTE: RMG_database previously had CH2_triplet as Y_rad_birad node, not Ct_rad ***

MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 28
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(9.04E+11,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "515",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 29
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(6.02E+12,A_UNITS,"*/",2.0),
                 n=(-0.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "516",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: CH3 + i-C4H9 --> i-C4H8 + CH4

pg. 58: Discussion on evaluated data

Entry 45,16 (b): No data available at the time.  The author estimates the disproportionation

rate coefficient as half the rate of CH3+n-C3H7=C3H6+H2 (due to half as many H-atoms
on the alpha-carbon).
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 30
rate(
  group1 = 
"""
C_rad/H2/Cs
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(8.43E+11,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "517",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 31
rate(
  group1 = 
"""
C_rad/H2/O
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "518",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 32
rate(
  group1 = 
"""
C_rad/H2/Cd
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(7.83E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "519",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + iC4H9 --> iC4H8 + C3H6

pg.270: Discussion on evaluated data

Entry 47,45(a): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC4H9+iC4H9-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 33
rate(
  group1 = 
"""
C_rad/H/NonDeC
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(2.56E+13,A_UNITS,"*/",2.0),
                 n=(-0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "520",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: i-C3H7 + i-C4H9 --> i-C4H8 + C3H8

pg. 65: Discussion on evaluated data

Entry 45,42 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of i-C3H7+n-C3H7=C3H6+C3H8 (due to half as many H-atoms
on the alpha-carbon).
*** NOTE: MRH computes half the rate of i-C3H7+n-C3H7=C3H6+C3H8 as 0.52x10^-11 * (300/T)^0.35,

not 0.58x10^-11 * (300/T)^0.35.  However, there may be a reason for the relatively
small discrepancy between the author\'s stated and implemented calculation. ***
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 34
rate(
  group1 = 
"""
C_rad/Cs3
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.08E+14,A_UNITS,"*/",2.0),
                 n=(-0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "521",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 35
rate(
  group1 = 
"""
Cd_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(8.43E+11,A_UNITS,"*/",4.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "522",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: C2H3 + i-C4H9 --> i-C4H8 + C2H4

pg. 60: Discussion on evaluated data

Entry 45,19 (b): No data available at the time.  Author estimates the disproportionation rate

coefficient based on the rate of C2H5+i-C4H9=i-C4H8+C2H6.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 36
rate(
  group1 = 
"""
O_pri_rad
1 *1 O 1 {2,S}
2 H 0 {1,S}
""",
  group2 = 
"""
C/H/NdNd_Csrad
1. *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
4. {Cs,O} 0 {1,S}
5. {Cs,O} 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+13,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "523",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: OH + i-C4H9 --> i-C4H8 + H2O

pg. 55: Discussion on evaluated data

Entry 45,6 (a): No data available at the time.  Author estimates the disproportionation rate

coefficient as half the rate of OH+n-C3H7=C3H6+H2O (due to half as many H-atoms
on the alpha-carbon).
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 37
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+13,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "525",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: H + C3H5 --> H2C=C=CH2 + H2

pg.252: Discussion on evaluated data

Entry 47,4(c): No data available at the time.  Author assigns a rate coefficient of 

3x10^-11 cm3/molecule/s for the disproportionation rxn.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 38
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.01E+12,A_UNITS,"*/",3.0),
                 n=(-0.32,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "526",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 39
rate(
  group1 = 
"""
C_rad/H2/Cs
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.64E+11,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "527",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C2H5 + C3H5 --> H2C=C=CH2 + C2H6

pg.259: Discussion on evaluated data

Entry 47,17(a): The recommended rate expression is derived from the experimentally-

determined disproportionation-to-addition ratio of 0.047 (James and Troughton)
and the addition rate rule (C2H5+C3H5-->adduct) calculated using the geometric
mean rule of the rxns C2H5+C2H5-->adduct and C3H5+C3H5-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 40
rate(
  group1 = 
"""
C_rad/H2/Cd
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.43E+10,A_UNITS,"*/",2.5),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.26,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "528",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + C3H5 --> H2C=C=CH2 + C3H6

pg.271-272: Discussion on evaluated data

Entry 47,47(b): The recommended rate expression is derived from the experimentally-

determined disproportionation-to-addition ratio of 0.008 (James and Kambanis)
and the addition rate rule (C3H5+C3H5-->adduct) calculated based on the results
of Tulloch et al.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 41
rate(
  group1 = 
"""
C_rad/H/NonDeC
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.58E+12,A_UNITS,"*/",3.0),
                 n=(-0.35,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "529",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: iC3H7 + C3H5 --> H2C=C=CH2 + C3H8

pg.268: Discussion on evaluated data

Entry 47,42(b): No data available at the time.  Recommended rate coefficient expression

based on rxn C3H5+C2H5=C2H4+C3H6 (James, D.G.L. and Troughton, G.E.) and values
for \"alkyl radicals\" (Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and iC3H7+iC3H7-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 42
rate(
  group1 = 
"""
C_rad/Cs3
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.89E+13,A_UNITS,"*/",3.0),
                 n=(-0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(-0.13,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "530",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: tC4H9 + C3H5 --> H2C=C=CH2 + iC4H10

pg.269: Discussion on evaluated data

Entry 47,44(b): No data available at the time.  Recommended rate coefficient expression

based on \"allyl and alkyl radicals behaving in similar fashion\" (possibly referencing
Gibian M.J. and Corley R.C.); this leads to disproportionation-
to-addition ratio of 0.04.  The addition rate expression was derived using the geometric
mean rule for the rxns C3H5+C3H5-->adduct and tC4H9+tC4H9-->adduct.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 43
rate(
  group1 = 
"""
Cd_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "531",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C2H3 + C3H5 --> H2C=C=CH2 + C2H4

pg.261-262: Discussion on evaluated data

Entry 47,19(d): No data available at the time.  Author recommends a rate coefficient

of 4x10^-12 cm3/molecule/s for the disproportionation rxn.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 44
rate(
  group1 = 
"""
O_pri_rad
1 *1 O 1 {2,S}
2 H 0 {1,S}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "532",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: OH + C3H5 --> H2C=C=CH2 + H2O

pg.253: Discussion on evaluated data

Entry 47,6(a): No data available at the time.  Author recommends a rate coefficient

of 1x10^-11 cm3/molecule/s, based on \"comparable rxns\".
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 45
rate(
  group1 = 
"""
O_atom_triplet
1 *1 O 2T
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(9.04E+13,A_UNITS,"+-",3.01e+13),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 3,
  old_id = "536",
  short_comment = "Grotheer et al [189].",
  long_comment = 
"""
[189] Grotheer, H.; Riekert, G.; Walter, D.; Just, T. Symp. Int. Combust. Proc. 1989, 22, 963.
Absolute value measured directly. Excitation: discharge, analysis: mass spectroscopy. Original uncertainty 3.0E+13

O + CH2OC --> OH + CH2O, O + CH3CHOH --> OH + CH3CHO

O+CH2OH --> OH+CH2O && O+CH3CHOH --> OH+CH3CHO

pg.963: Measured rate coefficients mentioned in abstract as k_2M and k_2E.

pg.965-967: Discussion on measured rate coefficients.

MRH 1-Sept-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 46
rate(
  group1 = 
"""
CH2_triplet
1 *1 C 2T {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "537",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2 + CH2OH --> CH3 + CH2O

pg. 505: Discussion on evaluated data

Entry 39,26 (b): CH2OH + CH2(triplet) --> CH3 + CH2O

Author estimates the rate of disproportionation as 2.0x10^-12 cm3/molecule/s.  No data at the time.

MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 47
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.0E+13,A_UNITS,"+-",1e+13),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (295,None),
  rank = 4,
  old_id = "538",
  short_comment = "Edelbuttel-Einhaus et al [190].",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 48
rate(
  group1 = 
"""
H_rad
1 *1 H 1
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(6.03E+12,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "539",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: H + CH2OH --> H2 + CH2O

pg. 496-497: Discussion on evaluated data

Entry 39,4 (a): CH2OH + H --> H2 + CH2O

Author estimates disproportionation rate will be faster than the H+C2H5=H2+C2H4 reaction

and reports rate coefficient as 1.0x10^-11 cm3/molecule/s.  No data at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 49
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(8.49E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 3,
  old_id = "540",
  short_comment = "Pagsberg et al [191].",
  long_comment = 
"""
[191] Pagsberg, P.; Munk, J.; Sillesen, A.; Anastasi, C. Chem. Phys. Lett. 1988, 146, 375.
Absolute value measured directly. Excitatio: electron beam, analysis: Vis-UV absorption.

CH2OH + CH3 --> CH2O + CH4

pg.378 Table 2: Formation and decay rates of CH2OH, CH3, and OH observed by pulse radiolysis of

gas mixtures of varying composition.  Chemical composition of systems A-E as in Table 1.
The authors note below Table 2 that the reported rate coefficient for CH3+CH2OH is an

\"adjustment of model to reproduce the observed decay rates of CH3 and CH2OH\".
MRH is skeptical of data, as this specific rxn is not directly referenced in the article,

nor do the authors address whether other channels besides -->CH4+CH2O exist / are significant.
The value of A in the database is consistent with that reported in Table 2.
MRH 1-Sept-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 50
rate(
  group1 = 
"""
C_methyl
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "541",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH3 + CH2OH --> CH4 + CH2O

pg. 500-501: Discussion on evaluated data

Entry 39,16 (b): CH2OH + CH3 --> CH4 + CH2O

Author estimates ratio of disproportionation rate to addition rate to be 0.2,

namely 4x10^-12 cm3/molecule/s.  No data at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 51
rate(
  group1 = 
"""
C_rad/H2/Cs
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "542",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: C2H5 + CH2OH --> C2H6 + CH2O

pg. 502: Discussion on evaluated data

Entry 39,17 (b): C2H5 + CH2OH --> C2H6 + CH2O

Author estimates the disproportionation rate coefficient as 4x10^-12 cm3/molecule/s.

No data at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 52
rate(
  group1 = 
"""
C_rad/H2/Cd
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Cd 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+13,A_UNITS,"*/",2.5),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "543",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
[93] Tsang, W.; Journal of Physical and Chemical Reference Data (1991), 20(2), 221-273.
Literature review: C3H5 + CH2OH --> CH2O + C3H6

pg.267: Discussion on evaluated data

Entry 47,39: No data available at the time.  Author notes that combination of these two

reactants will form 3-butene-1-ol which should decompose under combustion conditions
to form C3H6 + CH2O (same products).  The author therefore recommends a rate
coefficient of 3x10^-11 cm3/molecule/s.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 53
rate(
  group1 = 
"""
C_rad/H2/O
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 O 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.82E+12,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "544",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2OH + CH2OH --> CH3OH + CH2O

pg. 506: Discussion on evaluated data

Entry 39,39 (b): CH2OH + CH2OH --> CH3OH + CH2O

Meier, et al. (1985) measured the rate of addition + disproportionation.  Tsang estimates

a disproportionation to combination ratio of 0.5
NOTE: Rate coefficient given in table at beginning of reference (summarizing all data

presented) gives k_a+b = 2.4x10^-11, leading to k_b = 8x10^-12.  NIST\'s online
database (kinetics.nist.gov) reports this number as well.  However, the discussion
on pg. 506 suggests k_a+b = 1.5x10^-11, leading to k_b = 5x10^-12.
MRH 30-Aug-2009

*** NEED TO INVESTIGATE ***
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 54
rate(
  group1 = 
"""
C_rad/H/NonDeC
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.35E+12,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "545",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: CH2OH + i-C3H7 = C3H8 + CH2O

pg. 945: Discussion on evaluated data

Entry 42,39 (b): No data available at the time.  Author suggests rate coefficient based

on rxn C2H5+i-C3H7=C3H8+C2H4, namely 3.9x10^-12 cm3/molecule/s
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 55
rate(
  group1 = 
"""
C_rad/Cs3
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.47E+14,A_UNITS,"*/",3.0),
                 n=(-0.75,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "546",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
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
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 56
rate(
  group1 = 
"""
Cd_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 C 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.01E+13,A_UNITS,"*/",2.5),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "547",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH2OH + C2H3 --> C2H4 + CH2O

pg. 503: Discussion on evaluated data

Entry 39,19 (a): CH2OH + C2H3 --> C2H4 + CH2O

Author suggests a disproportionation rate coefficient near the collision limit, due

to rxn\'s exothermicity.  No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 57
rate(
  group1 = 
"""
Ct_rad
1 *1 C 1 {2,T}
2 C 0 {1,T}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.61E+13,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "548",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: C2H + CH2OH --> C2H2 + CH2O

pg. 504: Discussion on evaluated data

Entry 39,21 (a): CH2OH + C2H --> C2H2 + CH2O

Author suggest a disproportionation rate coefficient of 6.0x10^-11 cm3/molecule/s, due

to very exothermic rxn.  No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 58
rate(
  group1 = 
"""
CO_pri_rad
1 *1 C 1 {2,D}, {3,S}
2 O 0 {1,D}
3 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.81E+14,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "549",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: HCO + CH2OH --> CH2O + CH2O

pg. 500: Discussion on evaluated data

Entry 39,15 (b): CH2OH + HCO --> 2 CH2O

Author estimates a disproportionation rate coefficient of 3x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 59
rate(
  group1 = 
"""
O_pri_rad
1 *1 O 1 {2,S}
2 H 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "550",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: OH + CH2OH --> H2O + CH2O

pg. 497: Discussion on evaluated data

Entry 39,6: CH2OH + OH --> H2O + CH2O

Author estimates a disproportionation rate coefficient of 4x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 60
rate(
  group1 = 
"""
O_rad/NonDeC
1 *1 O 1 {2,S}
2 Cs 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.41E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "551",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: CH3O + CH2OH --> CH3OH + CH2O

pg. 505: Discussion on evaluated data

Entry 39,24: CH2OH + CH3O --> CH3OH + CH2O

Author estimates a disproportionation rate coefficient of 4x10^-11 cm3/molecule/s.

No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 61
rate(
  group1 = 
"""
O_rad/NonDeO
1 *1 O 1 {2,S}
2 O 0 {1,S}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "552",
  short_comment = "Tsang [90] Literature review.",
  long_comment = 
"""
[90] Tsang, W.; Journal of Physical and Chemical Reference Data (1987), 16(3), 471-508.
Literature review: HO2 + CH2OH --> CH3OH + H2O2

pg. 498: Discussion on evaluated data

Entry 39,7: CH2OH + HO2 --> H2O2 + CH2O

Author recommends a disproportionation rate coefficient of 2x10^-11 cm3/molecules/s.

No data available at the time.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 62
rate(
  group1 = 
"""
Y_rad
1 *1 R 1
""",
  group2 = 
"""
O_Orad
1. *2 O 0 {2,S}, {3,S}
2. *3 O 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(3E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "10001",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


