# encoding: utf-8
header = """
1+2_Cycloaddition

multiplebond + elec_def -> cycle


Reverse name: Three_Ring_Cleavage

(1) CHANGE_BOND		{*1,-1,*2}
(2) FORM_BOND		{*1,S,*3}
(3) FORM_BOND		{*2,S,*3}
(4) LOSE_RADICAL 	{*3,2}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "1+2_Cycloaddition"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f15: 1+2 cycloaddition

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f15_1+2_cycloaddition

// Catherina Wijaya Thesis, pg 131
// [106] Cobos, C. Chem. - Phys. 1985, 83, 1010 pages 3-4
// [192] Frey, H. M. - J. Am. Chem. Soc. 1957, 79, 6373.
// [194] Gaedtke, H. Symp. Int. Combust. Proc. 1973, 14, 295.
// [195] Herbrechtsmeier, P. Reactions of O(3P) Atoms with Unsaturated C3 Hydrocarbons. In Combust. Inst. European Symp., 1973; pp13.
// [196] Smith, I.W.M. Trans. Faraday Soc. 1968, 64, 378.
// [197] Cvetanovic, R. J. Chem. Phys. 1959, 30, 19.

// rate constants from rate_library_4.txt, Cath, 03/07/28

//No.		elec_def	multiplebond		Temp.		A			n		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
elec_def
Union {carbene, me_carbene, dime_carbene, ph_carbene, o_atom}
""",
  group2 = 
"""
multiplebond
1 *1 {Cd,CO,O} 0 {2,D}
2 *2 {Cd,O} 0 {1,D}
""",
  kf = Arrhenius(A=(1E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "576",
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
carbene
1 *3 C 2 {2,S} {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  group2 = 
"""
mb_db_unsub
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.98E+12,A_UNITS,"*/",3.2),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(5.29,E_UNITS,"+-",0.26)
                 ),
  temperature_range = (296,728),
  rank = 3,
  old_id = "577",
  short_comment = "Frey et al [192]",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_o2_doublebond
1 *1 O 0 {2,D}
2 *2 O 0 {1,D}
""",
  kf = Arrhenius(A=(3.49E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.46,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,2100),
  rank = 3,
  old_id = "578",
  short_comment = "Cobos et al [106] Transition state theory. (pages 3-4)",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_unsub
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  kf = Arrhenius(A=(7.0E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,None),
  rank = 4,
  old_id = "579",
  short_comment = "Gaedtke et al [194]",
  long_comment = 
"""
[194] Gaedtke, H. Symp. Int. Combust. Proc. 1973, 14, 295. 
Excitation: direct photolysis, analysis: UV-Vis absorption, Pressure 0.1 - 1000 atm. O + C2H4 --> Oxirane
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_monosub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(2.9E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,None),
  rank = 4,
  old_id = "580",
  short_comment = "Gaedtke et al [194]",
  long_comment = 
"""
[194] Gaedtke, H. Symp. Int. Combust. Proc. 1973, 14, 295. 
Excitation: direct photolysis, analysis: UV-Vis absorption, Pressure 0.1 - 1000 atm. O + CH3CH=CH2 --> methyloxirane
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_monosub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(4.2E+12,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.50,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (275,360),
  rank = 4,
  old_id = "581",
  short_comment = "Herbrechtsmeier et al [195]",
  long_comment = 
"""
[195] Herbrechtsmeier, P. Reactions of O(3P) Atoms with Unsaturated C3 Hydrocarbons. In Combust. Inst. European Symp., 1973; pp13.
Absolute values measured directly. Excitation: discharge, analysis :GC, Pressure 0.01 atm. O + CH3CH=CH2 --> methyloxirane
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_monosub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.9E+12,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.80,E_UNITS,"+-",0.4)
                 ),
  temperature_range = (298,410),
  rank = 4,
  old_id = "582",
  short_comment = "Smith [196]",
  long_comment = 
"""
[196] Smith, I.W.M. Trans. Faraday Soc. 1968, 64, 378.
Data derived from fitting to a complex mechanism. Excitation: flash photolysis, analysis : UV-Vis absorption. Pressure 0.13 atm

O + 1-C4H8 --> ethyloxirane. Original uncertainty 3.0E+11
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_onecdisub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 {Cs,O} 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(7.6E+12,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0.10,E_UNITS,"+-",0.4)
                 ),
  temperature_range = (298,410),
  rank = 4,
  old_id = "583",
  short_comment = "Smith [196]",
  long_comment = 
"""
[196] Smith, I.W.M. Trans. Faraday Soc. 1968, 64, 378.
Data derived from fitting to a complex mechanism. Excitation: flash photolysis, analysis : UV-Vis absorption. Pressure 0.13 atm

O + iso-C4H8 --> 2,2- dimethyloxirane. Original uncertainty 1.2E+12
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_twocdisub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.54E+13,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 4,
  old_id = "584",
  short_comment = "Cvetanovic [197]",
  long_comment = 
"""
[197] Cvetanovic, R. J. Chem. Phys. 1959, 30, 19.
Relative value measured (O + (Z)-2-C4H8 --> cis-2,3-dimethyloxirane/O + C2H4 = Oxirane --> 2.2E+01) 

Pressure 0.39 atm. Excitation : sensitized photolysis, analysis :GC. 
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
o_atom
1 *3 O 2
""",
  group2 = 
"""
mb_db_tetrasub_Nd
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 {Cs,O} 0 {1,S}
4 {Cs,O} 0 {1,S}
5 {Cs,O} 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(3.18E+13,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (298,None),
  rank = 4,
  old_id = "585",
  short_comment = "Cvetanovic [197]",
  long_comment = 
"""
[197] Cvetanovic, R. J. Chem. Phys. 1959, 30, 19.
Relative value measured (O + (CH3)2C=C(CH3)2 --> tetramethyl-oxirane/O + iso-C4H8 --> 2,2-Dimethyloxirane = 4.18)  

Pressure 0.39 atm. Excitation : sensitized photolysis, analysis :GC.
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


