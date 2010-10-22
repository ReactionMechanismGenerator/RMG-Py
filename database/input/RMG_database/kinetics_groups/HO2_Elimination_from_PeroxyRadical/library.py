# encoding: utf-8
header = """
HO2_Elimination_from_PeroxyRadical

R2OO -> R=R + OOH


Reverse name: none

(1) BREAK_BOND		{*1,S,*5}
(2) BREAK_BOND		{*2,S,*3}
(3) CHANGE_BOND		{*1,1,*2}
(4) FORM_BOND		{*4,S,*5}
(5) GAIN_RADICAL	{*3,1}
(6) LOSE_RADICAL 	{*4,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "HO2_Elimination_from_PeroxyRadical"

# determines permitted units for rate expression:
reaction_order = 1

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f34: HO2 elimination from peroxy radical

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

// f34 HO2 elimination from peroxy radical
// rate constants from rate_library_4.txt, cath, 03/07/28

//		R2OO			Temp.		A			n	a		E0		DA		Dn		Da		DE0		Rank	Comments

"""

# Set some units for all the rates in this file

A_UNITS = "1/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
R2OO
1 *1 {C,Si} 0 {2,S}, {5,S}
2 *2 {C,Si} 0 {1,S}, {3,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
""",
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(1,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "835",
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
R2OO_HNd_2H
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 {Cs,O} 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
""",
  kf = Arrhenius(A=(4.79E+07,A_UNITS,"+-",0.0),
                 n=(1.46,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "837",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
R2OO_2H_2H
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
""",
  kf = Arrhenius(A=(1.56E+07,A_UNITS,"+-",0.0),
                 n=(1.69,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "836",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
R2OO_NdNd_2H
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {1,S}
8 H 0 {2,S}
9 H 0 {2,S}
""",
  kf = Arrhenius(A=(5.06E+08,A_UNITS,"+-",0.0),
                 n=(1.19,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "838",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
R2OO_2H_HNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 H 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(9.79E+08,A_UNITS,"+-",0.0),
                 n=(1.17,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.1,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "839",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
R2OO_HNd_HNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 {Cs,O} 0 {1,S}
8 H 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(1.65E+09,A_UNITS,"+-",0.0),
                 n=(1.01,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "840",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
R2OO_NdNd_HNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {1,S}
8 H 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(6.48E+10,A_UNITS,"+-",0.0),
                 n=(0.57,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "841",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
R2OO_2H_NdNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 H 0 {1,S}
8 {Cs,O} 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(7.48E+09,A_UNITS,"+-",0.0),
                 n=(1.08,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(29.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "842",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
R2OO_HNd_NdNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 H 0 {1,S}
7 {Cs,O} 0 {1,S}
8 {Cs,O} 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(8.11E+14,A_UNITS,"+-",0.0),
                 n=(-0.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(30.4,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "843",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
rate(
  group1 = 
"""
R2OO_NdNd_NdNd
1 *1 C 0 {2,S}, {5,S}, {6,S}, {7,S}
2 *2 C 0 {1,S}, {3,S}, {8,S}, {9,S}
3 *3 O 0 {2,S}, {4,S}
4 *4 O 1 {3,S}
5 *5 H 0 {1,S}
6 {Cs,O} 0 {1,S}
7 {Cs,O} 0 {1,S}
8 {Cs,O} 0 {2,S}
9 {Cs,O} 0 {2,S}
""",
  kf = Arrhenius(A=(3.10E+19,A_UNITS,"+-",0.0),
                 n=(-1.78,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(31.7,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "844",
  short_comment = "Sumathy\'s CBS-QB3 calculations. Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d\') level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


