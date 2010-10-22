# encoding: utf-8
header = """
1,3_Insertion_CO2

CO2 + RR' -> R_(CO2)_R'


Reverse name: 1,2_Elimination_CO2

(1) BREAK_BOND		{*3,S,*4}
(2) CHANGE_BOND		{*1,-1,*2}
(3) FORM_BOND		{*1,S,*3}
(4) FORM_BOND		{*2,S,*4}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "1,3_Insertion_CO2"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f13: 1,3 insertion for CO2



// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP

Arrhenius_EP



//f13_1,3_insertion_CO2
// Catherina Wijaya Thesis pg 130

// [87] Sumathi et al (2003) - CBS-QB3 calculations. 

// rate constants from rate_library_4.txt, Cath, 03/07/28

//No.		CO2			RR\'				Temp.		A			n	a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
CO2
Union {CO2_Od,CO2_Cdd}
""",
  group2 = 
"""
RR'
Union {R_H,R_R'}
""",
  kf = Arrhenius(A=(1E+05,A_UNITS,"+-",0.0),
                 n=(2,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(75,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "571",
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
CO2_Cdd
1 *1 Cdd 0 {2,D} {3,D}
2 *2 Od  0 {1,D}
3    Od  0 {1,D}
""",
  group2 = 
"""
H2
1  *3 H 0 {2,S}
2  *4 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.51E+09,A_UNITS,"+-",0.0),
                 n=(1.23,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(73.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "572",
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
CO2_Cdd
1 *1 Cdd 0 {2,D} {3,D}
2 *2 Od  0 {1,D}
3    Od  0 {1,D}
""",
  group2 = 
"""
C_methane
1 *3 Cs 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *4 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.53E+03,A_UNITS,"+-",0.0),
                 n=(2.83,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(79.2,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "573",
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
CO2_Cdd
1 *1 Cdd 0 {2,D} {3,D}
2 *2 Od  0 {1,D}
3    Od  0 {1,D}
""",
  group2 = 
"""
C_pri/NonDeC
1 *3 Cs 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *4 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.09E+04,A_UNITS,"+-",0.0),
                 n=(2.56,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(76.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "574",
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
CO2_Cdd
1 *1 Cdd 0 {2,D} {3,D}
2 *2 Od  0 {1,D}
3    Od  0 {1,D}
""",
  group2 = 
"""
C/H2/NonDeC
1 *3 Cs 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *4 H 0 {1,S}
3 H 0 {1,S}
4 Cs 0 {1,S}
5 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.06E+05,A_UNITS,"+-",0.0),
                 n=(2.13,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(77.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "575",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


