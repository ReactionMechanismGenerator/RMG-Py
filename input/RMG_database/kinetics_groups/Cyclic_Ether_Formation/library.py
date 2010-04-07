# encoding: utf-8
header = """
Cyclic_Ether_Formation

RnOOH -> RO + OH


Reverse name: OH+CyclicEther_Form_Alkyl-hydroperoxyl

(1) BREAK_BOND		{*2,S,*3}
(2) FORM_BOND		{*1,S,*2}
(3) GAIN_RADICAL	{*3,1}
(4) LOSE_RADICAL 	{*1,1}


Generated on 7th April 2010 at 17:08
"""

reaction_family_name = "Cyclic_Ether_Formation"

# These lines were in the RMG library file but were not translated into anything useful:
unread_lines= """
// rate library for f31: cyclic ether formation reaction
// original from rate_library_4.txt, Cath, 03/07/28

// jing, define key word for format of the rate: either Arrhenius or Arrhenius_EP
Arrhenius_EP

//f31_intermolecular_HA

// Catherina Wijaya Thesis, pg 159.

//No.		RnOOH		Y_rad_intra			Temp.		A			n		a		E0		A		Dn		Da		DE0		Rank	Comments

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
Union {R2OOH, R3OOH, R4OOH, R5OOH}
""",
  group2 = 
"""
Y_rad_intra
1 *1 R 1
""",
  kf = Arrhenius(A=(1E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(10,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 0,
  old_id = "812",
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
C_pri_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 *4 C 0 {1,S}
""",
  kf = Arrhenius(A=(3.98E+12,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(1.3,None,"+-",0.3),
                 E0=(37.0,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "813",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
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
C_sec_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 *4 C 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(1.38E+12,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(1.3,None,"+-",0.3),
                 E0=(37.0,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "814",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 4
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
C_ter_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.09E+12,A_UNITS,"*/",1.2),
                 n=(0,None,"+-",0.0),
                 alpha=(1.3,None,"+-",0.3),
                 E0=(37.0,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "815",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
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
C_pri_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 *4 C 0 {1,S}
""",
  kf = Arrhenius(A=(4.47E+11,A_UNITS,"*/",1.74),
                 n=(0,None,"+-",0.0),
                 alpha=(1.0,None,"+-",0.1),
                 E0=(38.2,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "816",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
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
C_sec_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 *4 C 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(2.04E+11,A_UNITS,"*/",1.74),
                 n=(0,None,"+-",0.0),
                 alpha=(1.0,None,"+-",0.1),
                 E0=(38.2,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "817",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 7
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
C_ter_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.31E+11,A_UNITS,"*/",1.74),
                 n=(0,None,"+-",0.0),
                 alpha=(1.0,None,"+-",0.1),
                 E0=(38.2,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "818",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
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
C_pri_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 *4 C 0 {1,S}
""",
  kf = Arrhenius(A=(5.13E+10,A_UNITS,"*/",1.41),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(14.8,E_UNITS,"+-",2.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "819",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
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
C_sec_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 *4 C 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(3.63E+10,A_UNITS,"*/",1.41),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(13.0,E_UNITS,"+-",2.5)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "820",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 10
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
C_ter_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 {R!H} 0 {1,S}
4 {R!H} 0 {1,S}
""",
  kf = Arrhenius(A=(2.57E+10,A_UNITS,"*/",1.41),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(11.5,E_UNITS,"+-",3.0)
                 ),
  temperature_range = (300,1500),
  rank = 2,
  old_id = "821",
  short_comment = "CBS-QB3 and BH&HLYP calculations (Catherina Wijaya & Sumathy Raman). Including treatment of hindered rotor.",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 11
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
Cs_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(6.00E+11,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(22.000,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "822",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 12
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
Cs_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(7.50E+10,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(15.250,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "823",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 13
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
Cs_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(9.38E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(7.000,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "824",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)

# Number 14
rate(
  group1 = 
"""
R5OOH_SSSS
1 *1 {Cd,Cs} 1 {2,S}
2 *4 {Cd,Cs} 0 {1,S}, {3,S}
3 {Cd,Cs} 0 {2,S}, {4,S}
4 {Cd,Cs} 0 {3,S}, {5,S}
5 {Cd,Cs} 0 {4,S}, {6,S}
6 *2 O 0 {5,S}, {7,S}
7 *3 O 0 {6,S}, {8,S}
8 H 0 {7,S}
""",
  group2 = 
"""
Cs_rad_intra
1 *1 C 1 {2,S}, {3,S}, {4,S}
2 *4 C 0 {1,S}
3 R 0 {1,S}
4 R 0 {1,S}
""",
  kf = Arrhenius(A=(1.17E+09,A_UNITS,"+-",0.0),
                 n=(0,None,"+-",0.0),
                 alpha=(0,None,"+-",0.0),
                 E0=(1.800,E_UNITS,"+-",0.0)
                 ),
  temperature_range = (300,1500),
  rank = 5,
  old_id = "825",
  short_comment = "Curran\'s [8] estimation in reaction type 19, QOOH = cyclic ether + OH",
  long_comment = 
"""
""",
   history = [("2010-04-07","Generated from current RMG library.","rwest@mit.edu")]
)


