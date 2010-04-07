# encoding: utf-8
header = """
Disproportionation_O2d

O2d + XH_Rrad -> HO2 + X_R


Reverse name: HO2_Addition

(1) CHANGE_BOND		{*1,-1,*5}
(2) GAIN_RADICAL	{*5,1}
(3) FORM_BOND		{*1,S,*4}
(4) BREAK_BOND		{*2,S,*4}
(5) CHANGE_BOND		{*2,1,*3}
(6) LOSE_RADICAL	{*3,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Disproportionation_O2d"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f09b: Disproportionation_O2d reaction
// move rate constants from original disproportionation rate library, JS, Aug, 05, 2003

// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f09b_Disproportionation_O2d
// 20001: usinging root number of all disproporionation reacitons

// Catherina Wijaya Thesis pg 157 - 158
// [91] Tsang, W. J. Phys. Chem. Ref. Data 1988, 17, 887.
// [92] Tsang, W. J. Phys. Chem. Ref. Data 1990, 19, 1.
// [93] Tsang, W. J. Phys. Chem. Ref. Data 1991, 20, 221.
// [98] Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F., Jr.; Kerr, J.A.; Rossi, M.J.; Troe, J. 
// \"summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry,\", 2001.
// [183] DeMore, W.B.; Sander, S.P.; Golden, D.M.; Hampson, R.F.; Kurylo, M.J.; Howard, C.J.; Ravishankara, A.R.; Kolb, C.E.; Molina, M.J. JPL Publication 97-4 1997, 1. 

//No.		Y_rad_birad		XH_Rrad			Temp.		A			N		a		E0		DA		Dn		Da		DE0		Rank	Comments


"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
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
  old_id = "20001",
  short_comment = "",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 2
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
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
  kf = Arrhenius(A=(1.26E+11,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (700,2500),
  rank = 4,
  old_id = "487",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
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
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
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
  kf = Arrhenius(A=(9.04E+10,A_UNITS,"*/",3.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (500,900),
  rank = 4,
  old_id = "501",
  short_comment = "Tsang [91] Literature review.",
  long_comment = 
"""
[91] Tsang, W.; Journal of Physical and Chemical Reference Data (1988), 17(2), 887-951.
Literature review: n-C3H7 + O2 = HO2 + C3H6

pg. 914-915: Discussion on evaluated data

Entry 41,3 (a): The author suggests a rate coefficient based on those reported in the

literature.  The author notes that the data reported in the literature suggests
the formation of C3H6 is controlled by the addition rxn.  The author further
notes that it is surprising that p-dependence effects are not observed for
C3H6 formation.
MRH 30-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
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
  kf = Arrhenius(A=(2.41E+10,A_UNITS,"*/",5.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (600,1000),
  rank = 4,
  old_id = "513",
  short_comment = "Tsang [92] Literature review.",
  long_comment = 
"""
[92] Tsang, W.; Journal of Physical and Chemical Reference Data (1990), 19(1), 1-68.
Literature review: O2 + iC4H9 --> iC4H8 + HO2

pg. 52-53: Discussion on evaluated data

Entry 45,3 (a): The author recommends a rate coefficient based on the experiments performed

by Baker et al. (yielding a disproportionation-to-decomposition ratio) and the
current (Tsang) study\'s recommended iC4H9 unimolecular decomposition rate.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
""",
  group2 = 
"""
Cdpri_Csrad
1. *2 Cd 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.21E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.55,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2500),
  rank = 4,
  old_id = "524",
  short_comment = "Tsang [93] Literature review.",
  long_comment = 
"""
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
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.14E+13,A_UNITS,"*/",2.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 4,
  old_id = "533",
  short_comment = "Atkinson et al [98] literature review.",
  long_comment = 
"""
[98] Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F., Jr.; Kerr, J.A.; Rossi, M.J.; Troe, J. \"Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry,\", 2001.
Literature review: CH3CHOH + O2 --> CH3CHO + HO2

Recommended value is k298.  This reference just gives a table of results,

with no discussion on how the preferred numbers were arrived at.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.78E+12,A_UNITS,"*/",1.3),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 4,
  old_id = "534",
  short_comment = "Atkinson et al [98] literature review.",
  long_comment = 
"""
[98] Atkinson, R.; Baulch, D.L.; Cox, R.A.; Crowley, J.N.; Hampson, R.F., Jr.; Kerr, J.A.; Rossi, M.J.; Troe, J. \"Summary of Evaluated Kinetic and Photochemical Data for Atmospheric Chemistry,\", 2001.
Literature review: CH2OH + O2 --> CH2O + HO2

Recommended value is k298.  This reference just gives a table of results,

with no discussion on how the preferred numbers were arrived at.
MRH 31-Aug-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
O2d
1 *1 O 0 {2,D}
2 *5 O 0 {1,D}
""",
  group2 = 
"""
O_Csrad
1. *2 O 0 {2,S}, {3,S}
2. *3 Cs 1 {1,S}
3. *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(5.48E+12,A_UNITS,"*/",1.3),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.4)
                 ),
  temperature_range = (200,300),
  rank = 4,
  old_id = "535",
  short_comment = "DeMore et al [183] literature review.",
  long_comment = 
"""
[183] DeMore, W.B.; Golden, D.M.; Hampson, R.F.; Howard, C.J.; Kolb, C.E.; Molina, M.J.; JPL Publication 97-4
Literature review: CH2OH + O2 --> CH2O + HO2

pg.62 D38: Discussion on evaluated data

pg.22: Recommended A-factor and E/R parameter values

MRH 1-Sept-2009
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


