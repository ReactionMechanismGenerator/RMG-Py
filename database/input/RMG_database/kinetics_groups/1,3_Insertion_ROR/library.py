# encoding: utf-8
header = """
1,3_Insertion_ROR

doublebond + R_OR -> R_doublebond_OR


Reverse name: 1,2_Elimination_ROR

(1) BREAK_BOND		{*3,S,*4}
(2) CHANGE_BOND		{*1,-1,*2}
(3) FORM_BOND		{*1,S,*3}
(4) FORM_BOND		{*2,S,*4}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "1,3_Insertion_ROR"

# determines permitted units for rate expression:
reaction_order = 2

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f13: 1,3 insertion for ROR



// JS, define key word for format of the rate: either Arrhenius or Arrhenius_EP

Arrhenius_EP



//f13_1,3_insertion_ROR
// Catherina Wijaya Thesis pg 130

// [87] Sumathi et al (2003) - CBS-QB3 calculations. 

// rate constants from rate_library_4.txt, Cath, 03/07/28

//No.	doublebond			R_OR			Temp.		A			N		a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
doublebond
Union {Cd_Cdd, Cdd_Cd, Cd_Cd}
""",
  group2 = 
"""
R_OR
Union {H_OR,R_OH}
""",
  kf = Arrhenius(A=(1E+02,A_UNITS,"+-",0.0),
                 n=(3,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(45,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "560",
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
Cd/H2_Cd/Nd2
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 *2 C 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 {Cs,O} 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_OCmethyl
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}, {5,S}, {6,S}
4 H 0 {3,S}
5 H 0 {3,S}
6 H 0 {3,S}
""",
  kf = Arrhenius(A=(9.36E+01,A_UNITS,"+-",0.0),
                 n=(2.85,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "561",
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
Cd/H2_Cd/H/Nd
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 *2 C 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_OCmethyl
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}, {5,S}, {6,S}
4 H 0 {3,S}
5 H 0 {3,S}
6 H 0 {3,S}
""",
  kf = Arrhenius(A=(2.48E+01,A_UNITS,"+-",0.0),
                 n=(2.93,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(43.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "562",
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
Cd/unsub_Cd/unsub
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_OCmethyl
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 Cs 0 {2,S}, {4,S}, {5,S}, {6,S}
4 H 0 {3,S}
5 H 0 {3,S}
6 H 0 {3,S}
""",
  kf = Arrhenius(A=(4.73E+01,A_UNITS,"+-",0.0),
                 n=(3.00,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(47.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "563",
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
cco_2H
1 *1 Cd 0 {2,D}, {4,S}, {5,S}
2 *2 Cdd 0 {1,D}, {3,D}
3 Od 0 {2,D}
4 H 0 {1,S}
5 H 0 {1,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.57E+02,A_UNITS,"+-",0.0),
                 n=(3.04,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(39.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "564",
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
cco_HNd
1 *1 Cd 0 {2,D}, {4,S}, {5,S}
2 *2 Cdd 0 {1,D}, {3,D}
3 Od 0 {2,D}
4 H 0 {1,S}
5 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(5.15E+01,A_UNITS,"+-",0.0),
                 n=(3.05,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(41.0,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "565",
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
cco_Nd2
1 *1 Cd 0 {2,D}, {4,S}, {5,S}
2 *2 Cdd 0 {1,D}, {3,D}
3 Od 0 {2,D}
4 {Cs,O} 0 {1,S}
5 {Cs,O} 0 {1,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.45E+03,A_UNITS,"+-",0.0),
                 n=(2.80,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(40.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "566",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
Cd/unsub_Cd/unsub
1 *1 Cd 0 {2,D} {3,S} {4,S}
2 *2 Cd 0 {1,D} {5,S} {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.47E+02,A_UNITS,"+-",0.0),
                 n=(2.94,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(53.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "567",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
Cd/H/Nd_Cd/H2
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 *2 C 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 {Cs,O} 0 {1,S}
5 H 0 {2,S}
6 H 0 {2,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(2.27E+02,A_UNITS,"+-",0.0),
                 n=(2.74,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(56.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "568",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
Cd/H2_Cd/H/Nd
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 *2 C 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 H 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(6.52E+01,A_UNITS,"+-",0.0),
                 n=(2.92,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(50.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "569",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
rate(
  group1 = 
"""
Cd/H2_Cd/Nd2
1 *1 C 0 {2,D}, {3,S}, {4,S}
2 *2 C 0 {1,D}, {5,S}, {6,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 {Cs,O} 0 {2,S}
6 {Cs,O} 0 {2,S}
""",
  group2 = 
"""
H_OH
1 *3 H 0 {2,S}
2 *4 Os 0 {1,S}, {3,S}
3 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.25E+03,A_UNITS,"+-",0.0),
                 n=(2.76,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(48.5,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "570",
  short_comment = "[87]CBS-QB3 calculations from Sumathi 2003.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


