# encoding: utf-8
header = """
1,2_Insertion

CO_birad + RR' -> R_CO_R'


Reverse name: 1,1_Elimination

(1) BREAK_BOND		{*2,S,*3}
(2) FORM_BOND		{*1,S,*2}
(3) FORM_BOND		{*1,S,*3}
(4) LOSE_RADICAL 	{*1,2}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "1,2_Insertion"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f11: 1,2 insertion
// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP

Arrhenius_EP

// f11_1,2_insertion
// Catherina Wijaya Thesis pg 130
// [87] Sumathi et al (2003) - CBS-QB3 calculations. 
// rate constants from rate_library_4.txt, Cath, 03/07/28

//No.	CO_birad	RR\'					Temp.		A			N		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
RR'
Union {R_H,R_R'}
""",
  kf = Arrhenius(A=(1E+05,A_UNITS,"+-",0.0),
                 n=(2,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(80,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "553",
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
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
C_methyl_C_methyl
1 *2 Cs 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *3 Cs 0 {1,S}, {6,S}, {7,S}, {8,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
6 H 0 {2,S}
7 H 0 {2,S}
8 H 0 {2,S}
""",
  kf = Arrhenius(A=(5.38E+02,A_UNITS,"+-",0.0),
                 n=(3.29,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(104.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "554",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
H2
1  *2 H 0 {2,S}
2  *3 H 0 {1,S}
""",
  kf = Arrhenius(A=(2.89E+09,A_UNITS,"+-",0.0),
                 n=(1.16,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(82.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "555",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
C_methane
1 *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *3 H 0 {1,S}
3    H 0 {1,S}
4    H 0 {1,S}
5    H 0 {1,S}
""",
  kf = Arrhenius(A=(1.64E+04,A_UNITS,"+-",0.0),
                 n=(2.86,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(86.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "556",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
C_pri/NonDeC
1 *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *3 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(9.14E+04,A_UNITS,"+-",0.0),
                 n=(2.53,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(85.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "557",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
C/H2/NonDeC
1 *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *3 H 0 {1,S}
3    H 0 {1,S}
4   Cs 0 {1,S}
5   Cs 0 {1,S}
""",
  kf = Arrhenius(A=(7.66E+05,A_UNITS,"+-",0.0),
                 n=(2.07,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(82.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "558",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
CO_birad
1 *1 C 2 {2,D}
2    O 0 {1,D}
""",
  group2 = 
"""
C/H/Cs3
1 *2 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *3 H 0 {1,S}
3 Cs 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(8.89E+07,A_UNITS,"+-",0.0),
                 n=(1.51,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(79.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "559",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


