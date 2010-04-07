# encoding: utf-8
header = """
intra_OH_migration

RnOOH -> HORO.


Reverse name: none

(1) BREAK_BOND		{*2,S,*3}
(2) FORM_BOND		{*1,S,*3}
(3) GAIN_RADICAL	{*2,1}
(4) LOSE_RADICAL 	{*1,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "intra_OH_migration"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f33: intra hydroxyl migration
// Originally from rate_library_4.txt, Cath, 03/07/28

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f33_intra_hydroxyl_migration
//			RnOOH		Y_rad_out			Temp.		A			n	a		E0		DA		Dn		Da		DE0		Rank	Comments
	
"""

# Set some units for all the rates in this file

A_UNITS = "cm^3/mol/s"
E_UNITS = "kcal/mol"

# And these are the rates...


# Number 1
rate(
  group1 = 
"""
RnOOH
Union {ROOH, R2OOH, R3OOH, R4OOH}
""",
  group2 = 
"""
Y_rad_out
1 *1 {Cd,Cs,Sid,Sis} 1
""",
  kf = Arrhenius(A=(1E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "826",
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
R2OOH_S
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 *2 O 0 {2,S}, {4,S}
4 *3 O 0 {3,S}, {5,S}
5 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(3.39E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "827",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 3
rate(
  group1 = 
"""
R2OOH_S
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 *2 O 0 {2,S}, {4,S}
4 *3 O 0 {3,S}, {5,S}
5 H 0 {4,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.69E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(24.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "828",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
rate(
  group1 = 
"""
R3OOH_SS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 *2 O 0 {3,S}, {5,S}
5 *3 O 0 {4,S}, {6,S}
6 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(4.47E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(27.8,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "829",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 5
rate(
  group1 = 
"""
R3OOH_SS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 *2 O 0 {3,S}, {5,S}
5 *3 O 0 {4,S}, {6,S}
6 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(2.88E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "830",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 6
rate(
  group1 = 
"""
R3OOH_SS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 *2 O 0 {3,S}, {5,S}
5 *3 O 0 {4,S}, {6,S}
6 H 0 {5,S}
""",
  group2 = 
"""
C_rad_out_Cs2
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(4.79E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(26.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "831",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
rate(
  group1 = 
"""
R4OOH_SSS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 {Cd,Cs} 0 {3,S}, {5,S}
5 *2 O 0 {4,S}, {6,S}
6 *3 O 0 {5,S}, {7,S}
7 H 0 {6,S}
""",
  group2 = 
"""
C_rad_out_2H
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 H 0 {1,S}
""",
  kf = Arrhenius(A=(1.12E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(16.3,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "832",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 8
rate(
  group1 = 
"""
R4OOH_SSS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 {Cd,Cs} 0 {3,S}, {5,S}
5 *2 O 0 {4,S}, {6,S}
6 *3 O 0 {5,S}, {7,S}
7 H 0 {6,S}
""",
  group2 = 
"""
C_rad_out_H/NonDeC
1 *1 C 1 {2,S}, {3,S}
2 H 0 {1,S}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(1.00E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.6,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "833",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 9
rate(
  group1 = 
"""
R4OOH_SSS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 {Cd,Cs} 0 {3,S}, {5,S}
5 *2 O 0 {4,S}, {6,S}
6 *3 O 0 {5,S}, {7,S}
7 H 0 {6,S}
""",
  group2 = 
"""
C_rad_out_Cs2
1 *1 C 1 {2,S}, {3,S}
2 Cs 0 {1,S}
3 Cs 0 {1,S}
""",
  kf = Arrhenius(A=(8.71E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.9,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "834",
  short_comment = "CBS-QB3 calculations (Catherina Wijaya). Treatment of hindered rotor included; hindered rotor PES are done at B3LYP/6-31g(d) level.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


